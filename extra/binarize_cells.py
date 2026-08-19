#!/usr/bin/env python3
"""
binarize_cells.py
=================

Adds explicit positive/negative binary calls to a PIPEX cell_data.csv,
using three independent thresholding methods applied directly to the
cell-level mean intensity distribution of each marker.

For every marker M in the input, three new columns are written:

    M_otsu_pos      (0/1)   2-class Otsu on cell means
    M_triangle_pos  (0/1)   Triangle threshold on cell means
    M_gmm_pos       (0/1)   2-component GMM crossover on cell means

This is the first stage of the PIPEX v2.0 extra/ folder pipeline. The
second stage (`pipex_marker_annotator.py`) adds a QC framework on top of
these binary calls. The third stage (`pipex_score_comparison.py`)
benchmarks PIPEX's internal scores against the consensus.

Why this script exists
----------------------

PIPEX v2.0 already writes three continuous per-cell scores at
segmentation time (`_otsu3`, `_triangle_score`, `_gmm_prob`). Those
scores are computed on the whole-image pixel intensity distribution,
which separates signal pixels from background pixels rather than
positive cells from negative cells. As a result, thresholding those
columns with a naive "score > 0" cutoff overcalls positivity for any
marker where most cells have some signal above background.

binarize_cells.py fits the same three thresholding methods to the
distribution of cell-level mean intensities instead. On this scale the
thresholds sit between the negative and positive cell populations, so
the resulting binary columns can be used directly for cell typing,
gating, or as filters for downstream analyses.

Usage
-----

By default the script reads `cell_data.csv` from the current working
directory. Run it once (with your PIPEX virtual environment activated)
from any folder that already contains a PIPEX cell_data.csv:

    python binarize_cells.py

To point it at a file elsewhere, or to restrict the panel to a subset
of markers:

    python binarize_cells.py --input path/to/cell_data.csv
    python binarize_cells.py --markers CD8 CD31 Ki67

For a full list of options, including transform choice and output
folder:

    python binarize_cells.py --help

Output
------

Two files are written next to the input CSV (or into --output-dir if
given):

    <stem>_binarized.csv     All original columns plus the three new
                             binary columns per marker.
    <stem>_binarize_qc.csv   Per-marker record of the thresholds each
                             method landed on and the resulting
                             positivity rates.

Dependencies
------------

    Required: numpy, pandas, scipy, scikit-learn, scikit-image
    Optional: anndata (only if you pass an .h5ad file as input)

Author: Mariya Mardamshina, Lundberg lab, Stanford University.
Part of the PIPΣX v2.0 extra/ folder.
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from skimage.filters import threshold_otsu, threshold_triangle
from sklearn.mixture import GaussianMixture

# Quiet sklearn / skimage warnings during per-marker fitting.
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#
# These constants set the default behaviour of the script. All of them
# can be overridden on the command line, so the typical workflow is to
# edit this block once for your project and then run the script with no
# arguments.

# Path to the PIPEX cell_data.csv you want to process. The default is
# a relative path, so if you run the script from a folder that contains
# a cell_data.csv it will pick it up automatically. Otherwise, either
# edit this constant or pass --input on the command line.
DEFAULT_INPUT_PATH = Path("cell_data.csv")

# Intensity transform applied before thresholding. Choices:
#   "arcsinh"  -- recommended for most fluorescence-scale spatial
#                 proteomics data. Cofactor is the 20th percentile of
#                 nonzero values (marker-adaptive).
#   "log1p"    -- classical log transform, good for counts-like data.
#   "none"     -- pass raw intensities straight to the thresholders.
DEFAULT_TRANSFORM = "arcsinh"

# Optional restricted marker panel. When empty, the script auto-detects
# every marker in the input CSV that has PIPEX helper columns. When
# populated, only these markers are processed and everything else in the
# CSV is left untouched. Useful for excluding technical channels (DAPI,
# Histone, BetaActin, and other segmentation-reference stains) from the
# analytical panel while keeping them in the CSV for QC purposes.
#
# Example: 46 biological markers, excluding DAPI/Histone/BetaActin
#   MARKERS_TO_ANALYZE = [
#       "ARL13B", "ATM", "Bcl2", "CD107a", "CD11c", "CD14", "CD163",
#       "CD20", "CD31", "CD34", "CD39", "CD3e", "CD40", "CD44", "CD45",
#       "CD4", "CD56", "CD68", "CD8", "Caveolin", "CollagenIV",
#       "Ecadherin", "ER", "EpCam", "GP100", "HIF1A", "HLAA", "HLADR",
#       "ICOS", "IDO1", "KRT5", "KRT8", "Ki67", "PCNA", "PCNT", "PD1",
#       "PDL1", "PanCK", "Podoplanin", "SMA", "SOX2", "TOX", "VISTA",
#       "Vimentin", "bCatenin1", "iNOS",
#   ]
MARKERS_TO_ANALYZE: list[str] = []


# ---------------------------------------------------------------------------
# Internal constants (rarely need editing)
# ---------------------------------------------------------------------------

# Standard PIPEX morphology columns to never treat as markers.
PIPEX_MORPH_COLS = {
    "cell_id", "size", "x", "y", "solidity", "eccentricity", "memref",
    "tile", "tile_id", "sample", "sample_id", "donor", "donor_id",
}

# Helper-column suffixes PIPEX writes per marker.
PIPEX_HELPER_SUFFIXES = (
    "_local_90", "_ratio_pixels", "_otsu3", "_triangle_score", "_gmm_prob",
)

# Clustering / annotation columns commonly produced by PIPEX downstream steps.
CLUSTERING_PREFIXES = ("leiden", "kmeans")

# Random seed for reproducible GMM initialisation.
RANDOM_STATE = 42

# Minimum number of distinct positive values required to fit thresholds.
MIN_NONZERO = 30


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(verbose: bool = False) -> logging.Logger:
    log = logging.getLogger("pipex_binarize")
    log.handlers.clear()
    h = logging.StreamHandler(sys.stdout)
    h.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    log.addHandler(h)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    return log


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------

def load_input(path: Path, log: logging.Logger) -> pd.DataFrame:
    """Load CSV or h5ad into a DataFrame with cells as rows."""
    suffix = path.suffix.lower()
    if suffix == ".csv":
        log.info(f"Reading CSV: {path}")
        df = pd.read_csv(path)
    elif suffix in (".h5ad", ".h5"):
        try:
            import anndata as ad
        except ImportError as e:
            raise ImportError(
                "Reading .h5ad requires the `anndata` package. "
                "Install with: pip install anndata"
            ) from e
        log.info(f"Reading AnnData: {path}")
        a = ad.read_h5ad(path)
        X = np.asarray(a.X.todense()) if hasattr(a.X, "todense") else np.asarray(a.X)
        df_x = pd.DataFrame(X, columns=list(a.var_names), index=a.obs_names)
        df = pd.concat([a.obs.reset_index(drop=True),
                        df_x.reset_index(drop=True)], axis=1)
    else:
        raise ValueError(f"Unsupported input extension: {suffix}")
    log.info(f"Loaded {df.shape[0]:,} cells x {df.shape[1]:,} columns")
    return df


# ---------------------------------------------------------------------------
# Marker detection
# ---------------------------------------------------------------------------

def detect_markers(df: pd.DataFrame,
                   user_markers: Optional[list[str]],
                   log: logging.Logger) -> list[str]:
    """Identify marker columns in a PIPEX dataframe.

    A column M is treated as a marker if it is numeric, not a known
    morphology / clustering / helper column, and either the user
    supplied an explicit list, or `M_local_90` also exists (PIPEX v2.0
    fingerprint).
    """
    cols = df.columns.tolist()
    if user_markers:
        missing = [m for m in user_markers if m not in cols]
        if missing:
            raise ValueError(f"Markers not found in input: {missing}")
        log.info(f"Using user-supplied marker list ({len(user_markers)} markers)")
        return list(user_markers)

    markers: list[str] = []
    for c in cols:
        if c in PIPEX_MORPH_COLS:
            continue
        if any(c.startswith(p) for p in CLUSTERING_PREFIXES):
            continue
        if any(c.endswith(s) for s in PIPEX_HELPER_SUFFIXES):
            continue
        if not pd.api.types.is_numeric_dtype(df[c]):
            continue
        # PIPEX v2.0 fingerprint: helper column for the marker exists.
        if f"{c}_local_90" in cols:
            markers.append(c)

    if not markers:
        # Fallback: no PIPEX helper columns present (e.g. h5ad input).
        # Treat all numeric, non-morph, non-clustering columns as markers.
        log.warning("No PIPEX helper columns found; falling back to "
                    "all non-morphology numeric columns as markers.")
        for c in cols:
            if c in PIPEX_MORPH_COLS:
                continue
            if any(c.startswith(p) for p in CLUSTERING_PREFIXES):
                continue
            if any(c.endswith(s) for s in PIPEX_HELPER_SUFFIXES):
                continue
            if pd.api.types.is_numeric_dtype(df[c]):
                markers.append(c)

    log.info(f"Auto-detected {len(markers)} marker columns")
    log.debug(f"Markers: {markers}")
    return markers


# ---------------------------------------------------------------------------
# Intensity transformation
# ---------------------------------------------------------------------------

def transform_values(x: np.ndarray, kind: str) -> np.ndarray:
    """Optional intensity transform. Same convention as the annotator:
    arcsinh cofactor = 20th percentile of nonzero values."""
    if kind == "none":
        return x
    if kind == "log1p":
        return np.log1p(np.clip(x, 0, None))
    if kind == "arcsinh":
        nz = x[x > 0]
        if nz.size == 0:
            return np.zeros_like(x)
        cofactor = max(np.percentile(nz, 20), 1e-3)
        return np.arcsinh(x / cofactor)
    raise ValueError(f"Unknown transform: {kind}")


# ---------------------------------------------------------------------------
# Per-method thresholding on cell means
# ---------------------------------------------------------------------------

def fit_otsu_2class(x: np.ndarray) -> tuple[float, np.ndarray]:
    """Two-class Otsu threshold on cell means. Returns (threshold, positive-mask).

    This is the standard positive-vs-negative Otsu, DIFFERENT from PIPEX's
    core `_otsu3` which uses the lower of a 3-class multi-Otsu split on
    the whole-image pixel distribution.
    """
    tau = float(threshold_otsu(x))
    pos = (x > tau).astype(np.int8)
    return tau, pos


def fit_triangle_2class(x: np.ndarray) -> tuple[float, np.ndarray]:
    """Triangle threshold on cell means. Returns (threshold, positive-mask).

    Distinct from PIPEX's core `_triangle_score` because it is computed
    on the cell-mean distribution rather than on the whole-image pixel
    distribution. On markers with many background pixels, the pixel-level
    Triangle threshold sits near the background peak and calls most cells
    positive; the cell-level Triangle threshold sits closer to the
    biological boundary between the negative and positive cell populations.
    """
    tau = float(threshold_triangle(x))
    pos = (x > tau).astype(np.int8)
    return tau, pos


def fit_gmm_2comp(x: np.ndarray) -> tuple[float, np.ndarray]:
    """Two-component GMM on cell means. Returns (crossover-threshold, positive-mask).

    The threshold is the posterior crossover point BETWEEN the two
    component means, not the raw predict_proba at 0.5 over the whole
    intensity range. This avoids a known GMM pathology on zero-inflated
    markers, where cells with intensity = 0 lie in the joint tail of
    both components and can be erroneously assigned to whichever
    component has the wider tail. Constraining the boundary to live
    between the two means guarantees biologically sensible calls.

    Distinct from PIPEX's core `_gmm_prob` which fits the GMM on
    whole-image pixels and returns the raw posterior; the cell-level
    fit is more robust for downstream binary calling because the cell
    means already integrate over within-cell pixel heterogeneity.
    """
    xx = x.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=RANDOM_STATE,
                          n_init=3, max_iter=200).fit(xx)
    means = gmm.means_.flatten()
    high = int(np.argmax(means))
    low = 1 - high
    mu_low, mu_high = float(means[low]), float(means[high])

    if mu_high > mu_low:
        # Find the posterior crossover (p_high = 0.5) in the valley.
        grid = np.linspace(mu_low, mu_high, 2001).reshape(-1, 1)
        p_grid_high = gmm.predict_proba(grid)[:, high]
        sign_changes = np.where(np.diff(np.sign(p_grid_high - 0.5)))[0]
        if sign_changes.size:
            tau = float(grid[sign_changes[0], 0])
        else:
            tau = (mu_low + mu_high) / 2.0
    else:
        # Degenerate: both means equal → no separation possible.
        tau = float(np.median(x))

    pos = (x > tau).astype(np.int8)
    return tau, pos


# ---------------------------------------------------------------------------
# Per-marker pipeline
# ---------------------------------------------------------------------------

def binarize_marker(df: pd.DataFrame,
                    marker: str,
                    transform: str,
                    log: logging.Logger) -> tuple[pd.DataFrame, dict]:
    """Run the three thresholding methods on one marker.
    Returns (new_columns_dataframe, qc_row)."""

    x_raw = df[marker].to_numpy(dtype=float, copy=True)
    finite = np.isfinite(x_raw)
    x_raw = np.where(finite, x_raw, 0.0)  # treat NaN as 0 (PIPEX uses zeros)

    n_nonzero = int((x_raw > 0).sum())
    if n_nonzero < MIN_NONZERO or len(np.unique(x_raw)) < 5:
        log.warning(f"  {marker}: too few non-zero values "
                    f"({n_nonzero}); marking all cells negative.")
        new = pd.DataFrame({
            f"{marker}_otsu_pos": np.zeros(len(df), dtype=np.int8),
            f"{marker}_triangle_pos": np.zeros(len(df), dtype=np.int8),
            f"{marker}_gmm_pos": np.zeros(len(df), dtype=np.int8),
        })
        qc = dict(marker=marker, n_cells=len(df), n_nonzero=n_nonzero,
                  transform=transform,
                  threshold_otsu=np.nan, threshold_triangle=np.nan,
                  threshold_gmm=np.nan,
                  pct_pos_otsu=0.0, pct_pos_triangle=0.0, pct_pos_gmm=0.0,
                  usable=False, reason="too few non-zero values")
        return new, qc

    x_used = transform_values(x_raw, transform)

    tau_o, pos_o = fit_otsu_2class(x_used)
    tau_t, pos_t = fit_triangle_2class(x_used)
    tau_g, pos_g = fit_gmm_2comp(x_used)

    new = pd.DataFrame({
        f"{marker}_otsu_pos": pos_o,
        f"{marker}_triangle_pos": pos_t,
        f"{marker}_gmm_pos": pos_g,
    })

    pct_pos_o = float(100.0 * pos_o.mean())
    pct_pos_t = float(100.0 * pos_t.mean())
    pct_pos_g = float(100.0 * pos_g.mean())

    qc = dict(
        marker=marker,
        n_cells=len(df),
        n_nonzero=n_nonzero,
        transform=transform,
        threshold_otsu=tau_o,
        threshold_triangle=tau_t,
        threshold_gmm=tau_g,
        pct_pos_otsu=pct_pos_o,
        pct_pos_triangle=pct_pos_t,
        pct_pos_gmm=pct_pos_g,
        usable=True,
        reason="",
    )

    log.info(
        f"  {marker:>14s}  Otsu+={pct_pos_o:5.1f}%  "
        f"Tri+={pct_pos_t:5.1f}%  GMM+={pct_pos_g:5.1f}%  "
        f"tau_o={tau_o:.3f} tau_t={tau_t:.3f} tau_g={tau_g:.3f}"
    )

    return new, qc


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute explicit positive/negative binary calls per "
                    "cell for each marker in a PIPEX v2.0 cell_data.csv, "
                    "using 2-class Otsu, Triangle, and 2-component GMM "
                    "thresholding on cell-level mean intensities.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", "-i", default=None, type=Path,
                   help=f"Path to PIPEX cell_data.csv (or .h5ad). "
                        f"Default: {DEFAULT_INPUT_PATH} (relative to current directory).")
    p.add_argument("--output-dir", "-o", default=None, type=Path,
                   help="Directory to write outputs. Default: input's folder.")
    p.add_argument("--markers", "-m", default=None, nargs="+",
                   help="Optional explicit list of marker columns. "
                        "If omitted, markers are auto-detected from PIPEX "
                        "helper-column fingerprints (M_local_90 must exist).")
    p.add_argument("--transform", "-t",
                   choices=["none", "log1p", "arcsinh"],
                   default=DEFAULT_TRANSFORM,
                   help=f"Intensity transform applied before thresholding "
                        f"(default: {DEFAULT_TRANSFORM}).")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Verbose logging.")
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    log = setup_logger(args.verbose)

    # Resolve input path from the default if not given on CLI.
    if args.input is None:
        in_path: Path = DEFAULT_INPUT_PATH.expanduser().resolve()
        log.info(f"No --input given; using default: {in_path}")
    else:
        in_path = args.input.expanduser().resolve()

    if not in_path.is_file():
        log.error(f"Input not found: {in_path}")
        log.error("Either run this script from a folder that contains "
                  "cell_data.csv, edit DEFAULT_INPUT_PATH at the top of "
                  "the script, or pass --input explicitly.")
        return 2

    out_dir: Path = (args.output_dir or in_path.parent).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = in_path.stem

    df = load_input(in_path, log)

    # Marker selection precedence: explicit CLI --markers wins; else
    # project MARKERS_TO_ANALYZE; else auto-detect all.
    if args.markers is not None:
        marker_list = args.markers
        log.info(f"Using CLI --markers list ({len(marker_list)} markers)")
    elif MARKERS_TO_ANALYZE:
        marker_list = list(MARKERS_TO_ANALYZE)
        log.info(f"Using project MARKERS_TO_ANALYZE list ({len(marker_list)} markers)")
    else:
        marker_list = None  # detect_markers will auto-detect
    markers = detect_markers(df, marker_list, log)
    if not markers:
        log.error("No markers found. Use --markers to supply an explicit list.")
        return 2

    log.info(f"Binarizing {len(markers)} markers (transform={args.transform})...")
    bin_frames: list[pd.DataFrame] = []
    qc_rows: list[dict] = []
    for marker in markers:
        new_cols, qc = binarize_marker(df, marker, args.transform, log)
        bin_frames.append(new_cols)
        qc_rows.append(qc)

    bin_df = pd.concat(bin_frames, axis=1)
    out_df = pd.concat([df.reset_index(drop=True),
                        bin_df.reset_index(drop=True)], axis=1)
    out_csv = out_dir / f"{stem}_binarized.csv"
    out_df.to_csv(out_csv, index=False)
    log.info(f"Wrote binarized table: {out_csv}  ({out_df.shape[0]:,} x {out_df.shape[1]:,})")

    qc_df = pd.DataFrame(qc_rows)
    qc_csv = out_dir / f"{stem}_binarize_qc.csv"
    qc_df.to_csv(qc_csv, index=False)
    log.info(f"Wrote per-marker QC report: {qc_csv}")

    n_unusable = int((~qc_df["usable"]).sum())
    if n_unusable:
        bad = qc_df.loc[~qc_df["usable"], ["marker", "reason"]].to_dict("records")
        log.warning(f"{n_unusable} marker(s) had too few non-zero values "
                    f"and were marked all-negative; see {qc_csv.name}:")
        for row in bad:
            log.warning(f"  - {row['marker']}: {row['reason']}")

    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
