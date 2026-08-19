#!/usr/bin/env python3
"""
pipex_marker_annotator.py
=========================

Adds a per-cell confidence layer, cross-method consensus voting, and
per-marker reliability flagging on top of the binary calls produced by
`binarize_cells.py`. This script does not recompute the binary calls
itself. It reads them from a prior binarize_cells.py run, along with
the thresholds recorded in that run's QC file, and adds a QC layer on
top.

This is the second stage of the PIPEX v2.0 extra/ folder pipeline. The
first stage (`binarize_cells.py`) writes the binary calls. The third
stage (`pipex_score_comparison.py`) benchmarks PIPEX's internal scores
against the consensus that this script produces.

What this script adds
---------------------

For every marker M in the input, the following per-cell columns are
appended:

    M_otsu_conf        [0.5, 1.0]   distance-from-threshold confidence
    M_triangle_conf    [0.5, 1.0]   distance-from-threshold confidence
    M_gmm_conf         [0.5, 1.0]   distance-from-threshold confidence
    M_consensus_pos    (0/1)        majority vote across the 3 methods
    M_n_methods_agree  (2 or 3)     how many methods agree with consensus
    M_unanimous        (0/1)        set if all 3 methods agreed
    M_combined_conf    [0.0, 1.0]   mean confidence x agreement fraction

A separate per-marker QC report records, for each marker, Sarle's
bimodality coefficient, the BIC preference for a 2-component over
1-component GMM fit, whether the marker was flagged as reliable, and a
free-text reason if not.

A note on confidence scoring
----------------------------

Confidence for all three methods is computed uniformly as
sigmoid(|x_transformed - threshold| / MAD), which puts the three methods
on a common [0.5, 1.0] scale. This makes their combined confidence
score meaningful. The GMM's shape-awareness is instead captured at the
marker level, through the bimodality coefficient and BIC preference in
the QC report, which flag markers where the GMM fit should not be
trusted overall.

Reliability flagging
--------------------

A marker is flagged as unreliable if any of the following holds:

    1. The cell-level distribution is essentially unimodal
       (Sarle's BC < 0.45 AND 2-component GMM BIC preference < 10).
    2. The consensus call is near-uniform (< 0.5% or > 99.5% positive
       cells), meaning there is no biological signal to threshold.
    3. One of the three methods disagrees with the median of the other
       two by more than 30 percentage points.

Each failure mode points to a specific pathology, so the reason string
tells you what went wrong rather than just that something did.

Usage
-----

By default the script reads `cell_data_binarized.csv` from the current
working directory and auto-locates the QC file
(`cell_data_binarize_qc.csv`) next to it. Run it once from the folder
that contains the binarize_cells.py output:

    python pipex_marker_annotator.py

To point at specific files, restrict to a subset of markers, or disable
the diagnostic plots:

    python pipex_marker_annotator.py --input path/to/binarized.csv
    python pipex_marker_annotator.py --markers CD8 CD31 Ki67
    python pipex_marker_annotator.py --no-plots

For a full list of options:

    python pipex_marker_annotator.py --help

Output
------

Two files are written next to the input CSV (or into --output-dir if
given), plus one folder of per-marker diagnostic plots when --plots is
enabled:

    <stem>_annotated.csv        Input columns plus the annotation
                                columns listed above.
    <stem>_reliability.csv      Per-marker reliability report with
                                bimodality metrics, positivity rates,
                                agreement stats, and reliable / reason.
    <stem>_diagnostics/         One PNG per marker showing the raw
                                and transformed histograms with the
                                three thresholds overlaid.

Dependencies
------------

    Required: numpy, pandas, scipy, scikit-learn
    Optional: matplotlib (for --plots), anndata (for .h5ad input)

Author: Mariya Mardamshina, Lundberg lab, Stanford University.
Part of the PIPΣX v2.0 extra/ folder of standalone post-processing scripts.
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
from scipy.special import expit  # numerically stable sigmoid
from sklearn.mixture import GaussianMixture

# Quiet sklearn warnings during per-marker QC.
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

# Path to the binarize_cells.py output CSV you want to annotate. The
# default is a relative path, so if you run the script from a folder
# that contains cell_data_binarized.csv it will pick it up. Otherwise,
# edit this constant or pass --input on the command line.
DEFAULT_INPUT_PATH = Path("cell_data_binarized.csv")

# Filename convention for the QC report produced by binarize_cells.py.
# The script will auto-locate this file next to the input CSV. If your
# QC file lives elsewhere or has a different name, pass --qc-report.
DEFAULT_QC_NAME = "cell_data_binarize_qc.csv"

# Write a per-marker diagnostic PNG (raw + transformed histogram with
# the three thresholds overlaid) into <output-dir>/<stem>_diagnostics/.
# Set to False to skip; useful when re-running on many datasets and
# only the tabular outputs matter.
DEFAULT_PLOTS = True

# Optional restricted marker panel. When empty, every marker in the QC
# report gets annotated. When populated, only these markers are
# processed. Note: if you already set MARKERS_TO_ANALYZE upstream in
# binarize_cells.py, the restriction propagates automatically (only
# those markers appear in the QC report), so this list can stay empty.
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

# Suffixes written by binarize_cells.py (the upstream stage).
BINARIZE_POS_SUFFIXES = ("_otsu_pos", "_triangle_pos", "_gmm_pos")

# Clustering / annotation columns commonly produced by PIPEX downstream steps.
CLUSTERING_PREFIXES = ("leiden", "kmeans")

# Random seed for reproducible GMM initialisation (for BIC ΔBIC computation).
RANDOM_STATE = 42

# Minimum number of distinct positive values required for QC computation.
MIN_NONZERO = 30


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(verbose: bool = False) -> logging.Logger:
    log = logging.getLogger("pipex_annot")
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


def find_qc_report(input_path: Path,
                   explicit_qc_path: Optional[Path],
                   log: logging.Logger) -> Path:
    """Locate the binarize_cells.py QC report.

    Explicit path (--qc-report) wins. Otherwise we look for a file with
    the same folder as the input and a name matching the pattern
    <stem>_binarize_qc.csv, where <stem> is the input stem with the
    "_binarized" suffix removed (e.g. cell_data_binarized.csv ->
    cell_data_binarize_qc.csv).
    """
    if explicit_qc_path is not None:
        p = explicit_qc_path.expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"QC report not found: {p}")
        log.info(f"Using QC report: {p}")
        return p

    stem = input_path.stem
    if stem.endswith("_binarized"):
        base_stem = stem[: -len("_binarized")]
    else:
        base_stem = stem
    candidate = input_path.parent / f"{base_stem}_binarize_qc.csv"
    if candidate.is_file():
        log.info(f"Auto-located QC report: {candidate}")
        return candidate

    raise FileNotFoundError(
        f"Could not find QC report next to {input_path}. "
        f"Expected: {candidate}. Provide it explicitly with --qc-report."
    )


def load_qc_report(path: Path, log: logging.Logger) -> pd.DataFrame:
    """Load the binarize_cells.py QC report and index it by marker."""
    qc = pd.read_csv(path)
    required = {"marker", "transform", "threshold_otsu",
                "threshold_triangle", "threshold_gmm"}
    missing = required - set(qc.columns)
    if missing:
        raise ValueError(
            f"QC report {path} is missing required columns: {sorted(missing)}. "
            f"Was it produced by binarize_cells.py?"
        )
    log.info(f"Loaded QC report for {len(qc)} markers")
    return qc.set_index("marker")


# ---------------------------------------------------------------------------
# Marker detection
# ---------------------------------------------------------------------------

def detect_markers(df: pd.DataFrame,
                   qc_index: pd.Index,
                   user_markers: Optional[list[str]],
                   log: logging.Logger) -> list[str]:
    """Identify marker columns to annotate.

    A column M is treated as a marker if:
      - all three of `M_otsu_pos`, `M_triangle_pos`, `M_gmm_pos` are
        present (binarize_cells.py fingerprint), AND
      - M itself is present and numeric (raw mean intensity), AND
      - M appears in the QC report.
    """
    cols = set(df.columns)
    qc_markers = set(qc_index)

    if user_markers:
        missing = [m for m in user_markers if m not in cols]
        if missing:
            raise ValueError(f"Markers not found in input: {missing}")
        candidates = list(user_markers)
    else:
        # Auto-detect from binarize_cells.py fingerprint.
        candidates = []
        for c in df.columns:
            if c in PIPEX_MORPH_COLS:
                continue
            if any(c.startswith(p) for p in CLUSTERING_PREFIXES):
                continue
            if any(c.endswith(s) for s in PIPEX_HELPER_SUFFIXES):
                continue
            if any(c.endswith(s) for s in BINARIZE_POS_SUFFIXES):
                continue
            if not pd.api.types.is_numeric_dtype(df[c]):
                continue
            # Fingerprint: all three binarize_cells.py columns exist.
            if all(f"{c}{sfx}" in cols for sfx in BINARIZE_POS_SUFFIXES):
                candidates.append(c)

    markers = []
    skipped_no_qc = []
    skipped_no_pos = []
    for m in candidates:
        if m not in qc_markers:
            skipped_no_qc.append(m)
            continue
        if not all(f"{m}{sfx}" in cols for sfx in BINARIZE_POS_SUFFIXES):
            skipped_no_pos.append(m)
            continue
        markers.append(m)

    if skipped_no_qc:
        log.warning(f"{len(skipped_no_qc)} marker(s) skipped: not in QC "
                    f"report ({skipped_no_qc[:5]}{'...' if len(skipped_no_qc) > 5 else ''})")
    if skipped_no_pos:
        log.warning(f"{len(skipped_no_pos)} marker(s) skipped: missing "
                    f"_pos columns from binarize_cells.py "
                    f"({skipped_no_pos[:5]}{'...' if len(skipped_no_pos) > 5 else ''})")

    log.info(f"Annotating {len(markers)} markers")
    log.debug(f"Markers: {markers}")
    return markers


# ---------------------------------------------------------------------------
# Transform & robust scale (must match binarize_cells.py conventions)
# ---------------------------------------------------------------------------

def transform_values(x: np.ndarray, kind: str) -> np.ndarray:
    """Reproduce the transform used by binarize_cells.py."""
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


def robust_scale(x: np.ndarray) -> float:
    """MAD-based scale (in input units), with a fallback to std."""
    mad = stats.median_abs_deviation(x, scale="normal")  # ~ std for Gaussian
    if not np.isfinite(mad) or mad <= 0:
        s = float(np.std(x))
        return s if s > 0 else 1.0
    return float(mad)


# ---------------------------------------------------------------------------
# Bimodality QC (per marker, computed on cell means)
# ---------------------------------------------------------------------------

def sarle_bimodality(x: np.ndarray) -> float:
    """Sarle's bimodality coefficient. BC > 5/9 ~ 0.555 suggests bimodality."""
    n = x.size
    if n < 8:
        return np.nan
    g = stats.skew(x, bias=False)
    k = stats.kurtosis(x, fisher=True, bias=False)  # excess kurtosis
    denom = k + (3.0 * (n - 1) ** 2) / ((n - 2) * (n - 3))
    if denom <= 0:
        return np.nan
    return (g ** 2 + 1.0) / denom


def gmm_bic_delta(x: np.ndarray) -> float:
    """BIC of 1-component minus BIC of 2-component GMM. Positive = bimodal."""
    xx = x.reshape(-1, 1)
    try:
        bic1 = GaussianMixture(n_components=1, random_state=RANDOM_STATE).fit(xx).bic(xx)
        bic2 = GaussianMixture(n_components=2, random_state=RANDOM_STATE,
                               n_init=2).fit(xx).bic(xx)
        return float(bic1 - bic2)
    except Exception:
        return np.nan


# ---------------------------------------------------------------------------
# Confidence & consensus
# ---------------------------------------------------------------------------

def distance_confidence(x_used: np.ndarray, tau: float, scale: float) -> np.ndarray:
    """Sigmoid distance-from-threshold confidence in [0.5, 1.0).

    Confidence approaches 0.5 at the threshold and approaches 1.0 as
    the cell moves away in either direction.
    """
    return expit(np.abs(x_used - tau) / max(scale, 1e-9))


def consensus_calls(otsu_pos: np.ndarray,
                    tri_pos: np.ndarray,
                    gmm_pos: np.ndarray,
                    otsu_conf: np.ndarray,
                    tri_conf: np.ndarray,
                    gmm_conf: np.ndarray
                    ) -> dict[str, np.ndarray]:
    """Compute majority-vote consensus + agreement + combined confidence."""
    stack = np.stack([otsu_pos, tri_pos, gmm_pos], axis=1)
    votes = stack.sum(axis=1)
    consensus = (votes >= 2).astype(np.int8)
    n_agree = np.where(consensus == 1, votes, 3 - votes).astype(np.int8)  # 2 or 3
    unanimous = (n_agree == 3).astype(np.int8)

    mean_conf = (otsu_conf + tri_conf + gmm_conf) / 3.0
    combined = mean_conf * (n_agree / 3.0)

    return {
        "consensus_pos": consensus,
        "n_methods_agree": n_agree,
        "unanimous": unanimous,
        "combined_conf": combined,
    }


# ---------------------------------------------------------------------------
# Optional plotting
# ---------------------------------------------------------------------------

def plot_marker(marker: str,
                x_raw: np.ndarray,
                x_used: np.ndarray,
                tau_otsu: float,
                tau_tri: float,
                tau_gmm: float,
                qc_row: dict,
                out_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))

    axes[0].hist(x_raw, bins=80, color="0.7", edgecolor="white")
    axes[0].set_title(f"{marker} - raw mean intensity")
    axes[0].set_xlabel("intensity")
    axes[0].set_ylabel("# cells")

    axes[1].hist(x_used, bins=80, color="0.7", edgecolor="white")
    if np.isfinite(tau_otsu):
        axes[1].axvline(tau_otsu, color="C0", lw=1.5, label=f"Otsu={tau_otsu:.2f}")
    if np.isfinite(tau_tri):
        axes[1].axvline(tau_tri, color="C1", lw=1.5, label=f"Triangle={tau_tri:.2f}")
    if np.isfinite(tau_gmm):
        axes[1].axvline(tau_gmm, color="C3", lw=1.5, label=f"GMM={tau_gmm:.2f}")
    axes[1].set_title(
        f"thresholds (BC={qc_row.get('bimodality_coef', np.nan):.2f}, "
        f"deltaBIC={qc_row.get('gmm_bic_delta', np.nan):.0f})"
    )
    axes[1].set_xlabel("intensity (post-transform)")
    axes[1].legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_dir / f"{marker}.png", dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main per-marker pipeline
# ---------------------------------------------------------------------------

def annotate_marker(df: pd.DataFrame,
                    marker: str,
                    qc: pd.DataFrame,
                    plots_dir: Optional[Path],
                    log: logging.Logger) -> tuple[pd.DataFrame, dict]:
    """Add confidence + consensus + reliability for one marker.
    Returns (new_columns_dataframe, reliability_qc_row)."""

    # Pull the raw intensity and the three binary calls.
    x_raw = df[marker].to_numpy(dtype=float, copy=True)
    finite = np.isfinite(x_raw)
    x_raw = np.where(finite, x_raw, 0.0)

    pos_o = df[f"{marker}_otsu_pos"].to_numpy(dtype=np.int8)
    pos_t = df[f"{marker}_triangle_pos"].to_numpy(dtype=np.int8)
    pos_g = df[f"{marker}_gmm_pos"].to_numpy(dtype=np.int8)

    # Pull the thresholds and transform from the QC report.
    qc_row = qc.loc[marker]
    transform = str(qc_row["transform"])
    tau_o = float(qc_row["threshold_otsu"])
    tau_t = float(qc_row["threshold_triangle"])
    tau_g = float(qc_row["threshold_gmm"])

    # Reproduce the transformed intensity to compute distance confidence.
    x_used = transform_values(x_raw, transform)
    scale = robust_scale(x_used)

    # Distance-based confidence for each method. If a threshold is NaN
    # (binarize_cells.py skipped this method for this marker), the
    # confidence is 0.5 (uninformative) rather than NaN, so consensus
    # arithmetic stays well-defined.
    conf_o = distance_confidence(x_used, tau_o, scale) if np.isfinite(tau_o) else np.full(len(df), 0.5)
    conf_t = distance_confidence(x_used, tau_t, scale) if np.isfinite(tau_t) else np.full(len(df), 0.5)
    conf_g = distance_confidence(x_used, tau_g, scale) if np.isfinite(tau_g) else np.full(len(df), 0.5)

    cons = consensus_calls(pos_o, pos_t, pos_g, conf_o, conf_t, conf_g)

    new = pd.DataFrame({
        f"{marker}_otsu_conf": np.round(conf_o, 4),
        f"{marker}_triangle_conf": np.round(conf_t, 4),
        f"{marker}_gmm_conf": np.round(conf_g, 4),
        f"{marker}_consensus_pos": cons["consensus_pos"],
        f"{marker}_n_methods_agree": cons["n_methods_agree"],
        f"{marker}_unanimous": cons["unanimous"],
        f"{marker}_combined_conf": np.round(cons["combined_conf"], 4),
    })

    # Per-marker QC.
    n_nonzero = int((x_raw > 0).sum())
    bc = sarle_bimodality(x_used) if n_nonzero >= MIN_NONZERO else np.nan
    bic_delta = gmm_bic_delta(x_used) if n_nonzero >= MIN_NONZERO else np.nan

    pct_pos_o = float(100.0 * pos_o.mean())
    pct_pos_t = float(100.0 * pos_t.mean())
    pct_pos_g = float(100.0 * pos_g.mean())
    pct_pos_cons = float(100.0 * cons["consensus_pos"].mean())
    pct_unan = float(100.0 * cons["unanimous"].mean())
    mean_combined = float(cons["combined_conf"].mean())

    # Reliability heuristic (same three failure modes as before).
    pcts = np.array([pct_pos_o, pct_pos_t, pct_pos_g])
    method_names = np.array(["otsu", "triangle", "gmm"])
    deltas = np.array([abs(pcts[i] - np.median(np.delete(pcts, i)))
                       for i in range(3)])
    outlier_idx = int(np.argmax(deltas))
    outlier_delta = float(deltas[outlier_idx])
    method_outlier = method_names[outlier_idx] if outlier_delta > 30 else ""

    unreliable_unimodal = (
        np.isfinite(bc) and bc < 0.45
    ) and (
        np.isfinite(bic_delta) and bic_delta < 10
    )
    unreliable_extreme = (pct_pos_cons < 0.5) or (pct_pos_cons > 99.5)
    unreliable_outlier = bool(method_outlier)
    reliable = not (unreliable_unimodal or unreliable_extreme or unreliable_outlier)
    reasons = []
    if unreliable_unimodal:
        reasons.append("unimodal distribution")
    if unreliable_extreme:
        reasons.append("near-uniform call (<0.5% or >99.5% positive)")
    if unreliable_outlier:
        reasons.append(f"{method_outlier} disagrees by {outlier_delta:.0f}pp")
    reason = "; ".join(reasons)

    qc_out = dict(
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
        pct_pos_consensus=pct_pos_cons,
        pct_unanimous=pct_unan,
        mean_combined_conf=mean_combined,
        bimodality_coef=bc,
        gmm_bic_delta=bic_delta,
        method_outlier=method_outlier,
        outlier_delta_pp=outlier_delta,
        reliable=reliable,
        reason=reason,
    )

    if plots_dir is not None:
        plot_marker(marker, x_raw, x_used, tau_o, tau_t, tau_g, qc_out, plots_dir)

    log.info(
        f"  {marker:>14s}  Otsu+={pct_pos_o:5.1f}%  "
        f"Tri+={pct_pos_t:5.1f}%  GMM+={pct_pos_g:5.1f}%  "
        f"Cons+={pct_pos_cons:5.1f}%  unanim={pct_unan:5.1f}%  "
        f"BC={bc:4.2f}  reliable={reliable}"
    )

    return new, qc_out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Add per-cell confidence, cross-method consensus, "
                    "and per-marker reliability flagging on top of the "
                    "binary calls produced by binarize_cells.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", "-i", default=None, type=Path,
                   help=f"Path to the binarize_cells.py output "
                        f"(typically cell_data_binarized.csv or .h5ad). "
                        f"Default: {DEFAULT_INPUT_PATH} (relative to current directory).")
    p.add_argument("--qc-report", "-q", default=None, type=Path,
                   help=f"Path to the binarize_cells.py QC report. "
                        f"If omitted, we auto-locate {DEFAULT_QC_NAME} "
                        f"next to --input.")
    p.add_argument("--output-dir", "-o", default=None, type=Path,
                   help="Directory to write outputs. Default: input's folder.")
    p.add_argument("--markers", "-m", default=None, nargs="+",
                   help="Optional explicit list of marker columns. "
                        "If omitted, markers are auto-detected from the "
                        "binarize_cells.py fingerprint (_otsu_pos, "
                        "_triangle_pos, _gmm_pos all present).")
    p.add_argument("--plots", action="store_true", default=DEFAULT_PLOTS,
                   help=f"Write a per-marker diagnostic histogram PNG "
                        f"to <output-dir>/<stem>_diagnostics/ "
                        f"(default: {DEFAULT_PLOTS}). Use --no-plots to disable.")
    p.add_argument("--no-plots", action="store_false", dest="plots",
                   help="Disable per-marker diagnostic plots.")
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
        log.error("Run binarize_cells.py first to produce "
                  "cell_data_binarized.csv, or run this script from a "
                  "folder that already contains one, or pass --input "
                  "explicitly.")
        return 2

    qc_path = find_qc_report(in_path, args.qc_report, log)

    out_dir: Path = (args.output_dir or in_path.parent).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = in_path.stem

    df = load_input(in_path, log)
    qc = load_qc_report(qc_path, log)
    # Marker selection precedence: explicit CLI --markers wins; else
    # project MARKERS_TO_ANALYZE; else auto-detect all from QC report.
    if args.markers is not None:
        marker_list = args.markers
        log.info(f"Using CLI --markers list ({len(marker_list)} markers)")
    elif MARKERS_TO_ANALYZE:
        marker_list = list(MARKERS_TO_ANALYZE)
        log.info(f"Using project MARKERS_TO_ANALYZE list ({len(marker_list)} markers)")
    else:
        marker_list = None  # detect_markers will auto-detect
    markers = detect_markers(df, qc.index, marker_list, log)
    if not markers:
        log.error(
            "No markers found. The input should be produced by "
            "binarize_cells.py and contain _otsu_pos, _triangle_pos, "
            "and _gmm_pos columns for each marker, matching the QC report."
        )
        return 2

    plots_dir: Optional[Path] = None
    if args.plots:
        plots_dir = out_dir / f"{stem}_diagnostics"
        plots_dir.mkdir(parents=True, exist_ok=True)
        log.info(f"Writing diagnostic plots to: {plots_dir}")

    log.info(f"Annotating {len(markers)} markers...")
    annot_frames: list[pd.DataFrame] = []
    qc_rows: list[dict] = []
    for marker in markers:
        new_cols, qc_row = annotate_marker(df, marker, qc, plots_dir, log)
        annot_frames.append(new_cols)
        qc_rows.append(qc_row)

    annot = pd.concat(annot_frames, axis=1)
    out_df = pd.concat([df.reset_index(drop=True),
                        annot.reset_index(drop=True)], axis=1)
    out_csv = out_dir / f"{stem}_annotated.csv"
    out_df.to_csv(out_csv, index=False)
    log.info(f"Wrote annotated table: {out_csv}  ({out_df.shape[0]:,} x {out_df.shape[1]:,})")

    qc_df = pd.DataFrame(qc_rows)
    qc_csv = out_dir / f"{stem}_reliability.csv"
    qc_df.to_csv(qc_csv, index=False)
    log.info(f"Wrote per-marker reliability report: {qc_csv}")

    n_unreliable = int((~qc_df["reliable"]).sum())
    if n_unreliable:
        bad = qc_df.loc[~qc_df["reliable"], ["marker", "reason"]].to_dict("records")
        log.warning(f"{n_unreliable} marker(s) flagged as unreliable; "
                    f"see {qc_csv.name}:")
        for row in bad:
            log.warning(f"  - {row['marker']}: {row['reason']}")

    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
