#!/usr/bin/env python3
"""
pipex_score_comparison.py
=========================

Benchmarks the three continuous per-marker scores that PIPEX writes at
segmentation time (`_otsu3`, `_triangle_score`, `_gmm_prob`) against the
cell-level consensus produced by the upstream stages of the extra/
folder pipeline.

This is the third stage of the PIPEX v2.0 extra/ folder pipeline. Stage
1 (`binarize_cells.py`) writes cell-level binary calls; stage 2
(`pipex_marker_annotator.py`) adds confidence and consensus; this
script takes the annotator's output and asks how well PIPEX's internal
scores agree with that consensus.

Two views of agreement
----------------------

For every (marker, method) pair the script reports two complementary
numbers, because they answer different questions.

  1. ROC AUC of the PIPEX score against the annotator's consensus call.
     Cutoff-free. AUC = 1.0 means the two methods rank cells identically
     regardless of where any threshold is drawn. AUC = 0.5 means the
     PIPEX score contains no information about the consensus. This is
     the most defensible summary when the correct PIPEX cutoff is
     unknown.

  2. Cutoff-dependent metrics: confusion matrix, agreement, Cohen's
     kappa, McNemar p, and direction-of-bias for a specific cutoff.
     The cutoff can come from three sources:
       --auto-cutoff (default)         Per-marker per-method optimum via
                                       Youden's J on the ROC curve.
       --no-auto-cutoff                Falls back to the fixed defaults
                                       below.
       --otsu3-cutoff / --triangle-    User-supplied fixed cutoffs (only
       cutoff / --gmm-prob-cutoff      used when --no-auto-cutoff is set).

Why the auto-cutoff can fall back
---------------------------------

Youden's J is meaningless when the score contains no information about
the consensus. In that case the ROC curve sits on the diagonal, `argmax`
of a flat J returns an endpoint, and the "optimal cutoff" is either
below every observation (calling everything positive) or above every
observation (calling everything negative). The script detects this and
falls back to the fixed default cutoff, recording the reason in a
`degenerate_reason` column. Four reasons are possible:

    "annotator all one class"     Consensus is uniform; ROC undefined.
    "AUC not above chance"        Score has no rank agreement.
    "cutoff outside score range"  Youden pick sits below or above all
                                  observed scores.
    "roc_curve failed"            sklearn raised on the ROC fit.

Usage
-----

By default the script reads `cell_data_binarized_annotated.csv` from
the current working directory. Run it once from the folder that
contains the annotator's output:

    python pipex_score_comparison.py

Other common invocations:

    python pipex_score_comparison.py --input path/to/annotated.csv
    python pipex_score_comparison.py --no-auto-cutoff
    python pipex_score_comparison.py --markers CD8 CD31 Ki67

For a full list of options:

    python pipex_score_comparison.py --help

Output
------

Two files are written next to the input CSV (or into --output-dir if
given), plus a folder of per-marker diagnostic plots when --plots is
enabled:

    <stem>_pipex_comparison.csv          One row per (marker, method),
                                         with ROC AUC, applied cutoff,
                                         positive rates, agreement,
                                         Cohen's kappa, McNemar p,
                                         direction of bias, and any
                                         degenerate_reason.
    <stem>_pipex_comparison_summary.txt  Human-readable summary,
                                         per method, listing every
                                         degenerate marker with reason.
    <stem>_comparison_diagnostics/       One PNG per marker showing
                                         PIPEX score against raw
                                         intensity, coloured by
                                         agreement status.

Dependencies
------------

    Required: numpy, pandas, scipy, scikit-learn
    Optional: matplotlib (only when --plots is enabled)

Author: Mariya Mardamshina, Lundberg lab, Stanford University.
Part of the PIPΣX v2.0 extra/ folder.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import cohen_kappa_score, confusion_matrix, roc_auc_score, roc_curve


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#
# These constants set the default behaviour of the script. All of them
# can be overridden on the command line, so the typical workflow is to
# edit this block once for your project and then run the script with no
# arguments.

# Path to the pipex_marker_annotator.py output CSV you want to
# benchmark. The default is a relative path, so if you run the script
# from a folder that contains cell_data_binarized_annotated.csv it will
# pick it up. Otherwise, edit this constant or pass --input on the
# command line.
DEFAULT_INPUT_PATH = Path("cell_data_binarized_annotated.csv")

# Use per-marker per-method Youden's J to select the cutoff for the
# cutoff-dependent metrics (see the module docstring for what this
# means). When set to False, the script uses the fixed cutoffs set via
# --otsu3-cutoff, --triangle-cutoff, --gmm-prob-cutoff.
DEFAULT_AUTO_CUTOFF = True

# Write per-marker diagnostic PNGs (PIPEX score against raw intensity,
# coloured by agreement status) into
# <output-dir>/<stem>_comparison_diagnostics/. Set to False to skip.
DEFAULT_PLOTS = True

# Optional restricted marker panel. When empty, every marker present in
# the annotated CSV is benchmarked. When populated, only these markers
# are processed. Note: if you already set MARKERS_TO_ANALYZE upstream
# in binarize_cells.py or pipex_marker_annotator.py, the restriction
# propagates automatically, so this list can stay empty.
MARKERS_TO_ANALYZE: list[str] = []


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(verbose: bool = False) -> logging.Logger:
    log = logging.getLogger("pipex_compare")
    log.handlers.clear()
    h = logging.StreamHandler(sys.stdout)
    h.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    log.addHandler(h)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    return log


# ---------------------------------------------------------------------------
# Method definition
# ---------------------------------------------------------------------------

# (method_name, pipex_column_suffix, annotator_pos_suffix, default_cutoff,
#  comparator: ">" or ">=").
METHODS = [
    ("otsu",     "_otsu3",          "_otsu_pos",     0.0, ">"),
    ("triangle", "_triangle_score", "_triangle_pos", 0.0, ">"),
    ("gmm",      "_gmm_prob",       "_gmm_pos",      0.5, ">="),
]


# ---------------------------------------------------------------------------
# Marker detection from an annotated dataframe
# ---------------------------------------------------------------------------

def detect_markers(df: pd.DataFrame,
                   user_markers: Optional[list[str]],
                   log: logging.Logger) -> list[str]:
    """A marker is a column M such that BOTH PIPEX's score columns AND
    the annotator's *_pos columns exist.

    If user_markers is provided, restrict to that list (intersected
    with what actually has both column families). Missing / mismatched
    names are logged and skipped rather than raising.
    """
    cols = set(df.columns)

    def has_all_families(col: str) -> bool:
        if not pd.api.types.is_numeric_dtype(df[col]):
            return False
        for _, suffix, ann_suffix, _, _ in METHODS:
            if (col + suffix) not in cols or (col + ann_suffix) not in cols:
                return False
        return True

    if user_markers:
        markers, missing = [], []
        for m in user_markers:
            if m in cols and has_all_families(m):
                markers.append(m)
            else:
                missing.append(m)
        if missing:
            log.warning(f"{len(missing)} requested marker(s) skipped "
                        f"(missing PIPEX or annotator columns): "
                        f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
        log.info(f"Using {len(markers)} of {len(user_markers)} requested markers.")
        return markers

    candidates = [c for c in df.columns if has_all_families(c)]
    log.info(f"Detected {len(candidates)} markers with both PIPEX and "
             f"annotator outputs.")
    return candidates


# ---------------------------------------------------------------------------
# Per-marker, per-method comparison
# ---------------------------------------------------------------------------

# Below this AUC, Youden's J is essentially flat along the ROC diagonal
# and picking argmax returns a meaningless endpoint. Treat as degenerate
# and refuse to return an auto-cutoff.
MIN_AUC_FOR_AUTOCUTOFF = 0.55


def auc_and_optimal_cutoff(score: np.ndarray,
                           ann_call: np.ndarray
                           ) -> tuple[float, float, str]:
    """ROC AUC of score predicting ann_call, plus Youden's-J optimal
    cutoff. Returns (auc, optimal_cutoff, degenerate_reason).

    degenerate_reason is an empty string when the auto-cutoff is
    trustworthy. It is populated with a short label when the cutoff
    should NOT be used:
      - "annotator all one class": annotator called 0 or 100% positive,
        so ROC is undefined.
      - "AUC not above chance":    score has no discriminative power
        against the annotator; Youden's J flattens along the diagonal.
      - "cutoff outside score range": Youden's J picked a threshold
        outside the observed score range, which never separates cells.
      - "roc_curve failed":        sklearn raised on the ROC fit.

    In all degenerate cases the returned cutoff is NaN.
    """
    y = ann_call.astype(int)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan"), float("nan"), "annotator all one class"
    try:
        auc = float(roc_auc_score(y, score))
    except Exception:
        return float("nan"), float("nan"), "roc_curve failed"

    if not np.isfinite(auc) or auc < MIN_AUC_FOR_AUTOCUTOFF:
        return auc, float("nan"), "AUC not above chance"

    try:
        fpr, tpr, thr = roc_curve(y, score)
        # Youden's J = TPR - FPR; pick the threshold that maximises it.
        # roc_curve returns thresholds with thr[0] = max+1 (sentinel);
        # skip it so we don't pick the sentinel by accident.
        j = tpr - fpr
        best_idx = int(np.argmax(j[1:])) + 1 if len(j) > 1 else 0
        cutoff = float(thr[best_idx])
    except Exception:
        return auc, float("nan"), "roc_curve failed"

    # Reject cutoffs that sit outside the observed score range — these
    # never separate cells regardless of what Youden's J says.
    s_min, s_max = float(np.min(score)), float(np.max(score))
    if not np.isfinite(cutoff) or cutoff <= s_min or cutoff >= s_max:
        return auc, float("nan"), "cutoff outside score range"

    return auc, cutoff, ""


def compare_one(pipex_call: np.ndarray,
                ann_call: np.ndarray) -> dict:
    """Compute confusion + agreement metrics for two binary call vectors."""
    n = len(pipex_call)
    # Confusion matrix with both labels enumerated (handles all-0 / all-1).
    cm = confusion_matrix(pipex_call, ann_call, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()  # PIPEX rows, annotator cols
    # Treating PIPEX as the "row" and annotator as "column":
    #   tn = both negative
    #   fp = PIPEX negative, annotator positive  ("annotator-only positives")
    #   fn = PIPEX positive, annotator negative  ("PIPEX-only positives")
    #   tp = both positive
    pipex_pos = float((tp + fn) / n) if n else float("nan")
    ann_pos = float((tp + fp) / n) if n else float("nan")
    agreement = float((tp + tn) / n) if n else float("nan")

    # Cohen's kappa, with safe handling for the degenerate cases where
    # one or both vectors is constant.
    try:
        kappa = float(cohen_kappa_score(pipex_call, ann_call))
        if not np.isfinite(kappa):
            kappa = np.nan
    except Exception:
        kappa = np.nan

    # McNemar exact test on the discordant pairs (fp, fn). Only the
    # discordant pairs matter; the test asks whether
    # P(annotator+ | PIPEX-) == P(annotator- | PIPEX+).
    n_disc = int(fp + fn)
    if n_disc == 0:
        mcnemar_p = 1.0
    else:
        # Two-sided exact binomial test: under H0 the smaller of (fp,
        # fn) is Binomial(n_disc, 0.5).
        k = int(min(fp, fn))
        mcnemar_p = float(2.0 * stats.binom.cdf(k, n_disc, 0.5))
        mcnemar_p = min(mcnemar_p, 1.0)

    # Direction-of-bias label.
    if mcnemar_p > 0.05 or n_disc < 10:
        direction = "no asymmetric bias"
    elif fp > fn:
        direction = "annotator more inclusive"
    else:
        direction = "PIPEX more inclusive"

    return dict(
        n=int(n),
        pipex_pos_rate=pipex_pos,
        annotator_pos_rate=ann_pos,
        both_negative=int(tn),
        both_positive=int(tp),
        pipex_only_positive=int(fn),
        annotator_only_positive=int(fp),
        agreement=agreement,
        cohen_kappa=kappa,
        mcnemar_p=mcnemar_p,
        direction=direction,
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_marker(marker: str,
                df: pd.DataFrame,
                cutoffs: dict[str, float],
                comparators: dict[str, str],
                results: dict[str, dict],
                out_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    x_raw = df[marker].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharex=True)
    for ax, (name, suffix, ann_suffix, _, _) in zip(axes, METHODS):
        score = df[marker + suffix].to_numpy(dtype=float)
        cutoff = cutoffs[name]
        comp = comparators[name]
        pipex_pos = (score > cutoff) if comp == ">" else (score >= cutoff)
        ann_pos = df[marker + ann_suffix].to_numpy(dtype=int).astype(bool)

        # Four cells categories.
        both_neg = (~pipex_pos) & (~ann_pos)
        both_pos = pipex_pos & ann_pos
        pipex_only = pipex_pos & (~ann_pos)
        ann_only = (~pipex_pos) & ann_pos

        ax.scatter(x_raw[both_neg], score[both_neg], s=4, alpha=0.35,
                   color="0.7", label=f"both - ({both_neg.sum()})")
        ax.scatter(x_raw[both_pos], score[both_pos], s=4, alpha=0.6,
                   color="C2", label=f"both + ({both_pos.sum()})")
        ax.scatter(x_raw[pipex_only], score[pipex_only], s=8, alpha=0.7,
                   color="C3", label=f"PIPEX-only ({pipex_only.sum()})")
        ax.scatter(x_raw[ann_only], score[ann_only], s=8, alpha=0.7,
                   color="C0", label=f"annotator-only ({ann_only.sum()})")

        ax.axhline(cutoff, color="k", lw=0.8, ls="--")
        r = results[name]
        ax.set_title(
            f"{name}: agree {r['agreement']:.1%}, κ={r['cohen_kappa']:.2f}\n"
            f"{r['direction']}",
            fontsize=9
        )
        ax.set_xlabel("raw intensity")
        ax.set_ylabel(f"{name} score")
        ax.legend(fontsize=7, frameon=False, loc="best")

    fig.suptitle(f"{marker}: PIPEX score vs annotator binary call",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / f"{marker}.png", dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare PIPEX v2.0 internal scores against "
                    "pipex_marker_annotator.py binary calls.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", "-i", default=None, type=Path,
                   help=f"Annotated CSV from pipex_marker_annotator.py "
                        f"(must contain both PIPEX score columns and "
                        f"annotator *_pos columns). "
                        f"Default: {DEFAULT_INPUT_PATH} (relative to current directory).")
    p.add_argument("--output-dir", "-o", default=None, type=Path,
                   help="Directory to write outputs. Default: input's folder.")
    p.add_argument("--otsu3-cutoff", type=float, default=0.0,
                   help="Threshold on PIPEX's _otsu3 score above which a "
                        "cell is positive (default: > 0).")
    p.add_argument("--triangle-cutoff", type=float, default=0.0,
                   help="Threshold on PIPEX's _triangle_score above which "
                        "a cell is positive (default: > 0).")
    p.add_argument("--gmm-prob-cutoff", type=float, default=0.5,
                   help="Threshold on PIPEX's _gmm_prob at or above which "
                        "a cell is positive (default: >= 0.5).")
    p.add_argument("--auto-cutoff", action="store_true", default=DEFAULT_AUTO_CUTOFF,
                   help=f"Per (marker, method), pick the cutoff that "
                        f"maximises Youden's J on the ROC curve "
                        f"(annotator call as the reference). Overrides "
                        f"the per-method cutoff flags (default: {DEFAULT_AUTO_CUTOFF}). "
                        f"Use --no-auto-cutoff to disable.")
    p.add_argument("--no-auto-cutoff", action="store_false", dest="auto_cutoff",
                   help="Use the fixed --*-cutoff values instead of Youden's J.")
    p.add_argument("--plots", action="store_true", default=DEFAULT_PLOTS,
                   help=f"Write a per-marker comparison plot (default: {DEFAULT_PLOTS}). "
                        f"Use --no-plots to disable.")
    p.add_argument("--no-plots", action="store_false", dest="plots",
                   help="Disable per-marker comparison plots.")
    p.add_argument("--markers", "-m", default=None, nargs="+",
                   help="Optional explicit list of marker columns to "
                        "benchmark. If omitted, uses the project "
                        "MARKERS_TO_ANALYZE list (or all detected markers "
                        "if that is empty).")
    p.add_argument("--verbose", "-v", action="store_true")
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
        log.error("Run binarize_cells.py then pipex_marker_annotator.py "
                  "first to produce cell_data_binarized_annotated.csv, "
                  "or run this script from a folder that already contains "
                  "one, or pass --input explicitly.")
        return 2
    out_dir: Path = (args.output_dir or in_path.parent).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = in_path.stem.replace("_annotated", "")

    log.info(f"Reading annotated table: {in_path}")
    df = pd.read_csv(in_path)
    log.info(f"Loaded {df.shape[0]:,} cells x {df.shape[1]:,} columns")

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
        log.error("No markers found with both PIPEX and annotator columns. "
                  "Did you run pipex_marker_annotator.py first?")
        return 2

    cutoffs = {
        "otsu":     args.otsu3_cutoff,
        "triangle": args.triangle_cutoff,
        "gmm":      args.gmm_prob_cutoff,
    }
    comparators = {name: comp for name, _, _, _, comp in METHODS}

    if args.auto_cutoff:
        log.info("Auto-cutoff mode: per-marker, per-method optimal "
                 "cutoff from Youden's J on the ROC.")
    else:
        log.info(f"Fixed PIPEX cutoffs: otsu3 > {cutoffs['otsu']}, "
                 f"triangle_score > {cutoffs['triangle']}, "
                 f"gmm_prob >= {cutoffs['gmm']}")

    plots_dir: Optional[Path] = None
    if args.plots:
        plots_dir = out_dir / f"{stem}_comparison_diagnostics"
        plots_dir.mkdir(parents=True, exist_ok=True)
        log.info(f"Writing per-marker plots to: {plots_dir}")

    rows: list[dict] = []
    n_fallback = 0
    for marker in markers:
        per_method_results: dict[str, dict] = {}
        for name, suffix, ann_suffix, _, comp in METHODS:
            score = df[marker + suffix].to_numpy(dtype=float)
            ann = df[marker + ann_suffix].to_numpy(dtype=int)

            # Cutoff-free: ROC AUC + Youden's-J optimal cutoff (with
            # degenerate-case detection).
            auc, opt_cut, degenerate = auc_and_optimal_cutoff(score, ann)

            # Pick the cutoff for the binarised comparison.
            if args.auto_cutoff and np.isfinite(opt_cut):
                cut = opt_cut
                cutoff_source = "auto"
                # Auto-cutoff is a threshold the score must exceed.
                pipex_call = (score > cut)
            elif args.auto_cutoff:
                # Auto mode requested but Youden's J was degenerate:
                # fall back to the fixed default for this method so we
                # still produce a computable comparison row.
                cut = cutoffs[name]
                cutoff_source = f"fallback:{degenerate}"
                pipex_call = (score > cut) if comp == ">" else (score >= cut)
                n_fallback += 1
            else:
                cut = cutoffs[name]
                cutoff_source = "fixed"
                pipex_call = (score > cut) if comp == ">" else (score >= cut)

            r = compare_one(pipex_call.astype(int), ann.astype(int))
            r["marker"] = marker
            r["method"] = name
            r["pipex_cutoff"] = float(cut)
            r["pipex_cutoff_source"] = cutoff_source
            r["roc_auc"] = auc
            r["optimal_cutoff_youden"] = opt_cut
            r["degenerate_reason"] = degenerate
            rows.append(r)
            per_method_results[name] = r
            log.info(
                f"  {marker:>14s} | {name:>8s}: "
                f"AUC={auc:5.3f}  "
                f"PIPEX+={r['pipex_pos_rate']:5.1%}  "
                f"Annot+={r['annotator_pos_rate']:5.1%}  "
                f"agree={r['agreement']:5.1%}  "
                f"κ={r['cohen_kappa']:5.2f}  "
                f"cut={cut:8.5f}({cutoff_source})"
            )
        if plots_dir is not None:
            # Use the cutoffs we actually applied.
            applied_cutoffs = {n: per_method_results[n]["pipex_cutoff"]
                               for n in per_method_results}
            plot_marker(marker, df, applied_cutoffs, comparators,
                        per_method_results, plots_dir)

    if n_fallback:
        log.warning(
            f"{n_fallback} (marker, method) pairs fell back from auto-cutoff "
            f"to the fixed default because Youden's J was degenerate "
            f"(see the `degenerate_reason` column in the output CSV)."
        )

    out = pd.DataFrame(rows)
    cols_order = ["marker", "method",
                  "roc_auc", "optimal_cutoff_youden",
                  "pipex_cutoff", "pipex_cutoff_source", "degenerate_reason", "n",
                  "pipex_pos_rate", "annotator_pos_rate",
                  "both_negative", "both_positive",
                  "pipex_only_positive", "annotator_only_positive",
                  "agreement", "cohen_kappa", "mcnemar_p", "direction"]
    out = out[cols_order]
    out_csv = out_dir / f"{stem}_pipex_comparison.csv"
    out.to_csv(out_csv, index=False)
    log.info(f"Wrote per-marker comparison: {out_csv}")

    # Human-readable summary. Force UTF-8 so this works on Windows
    # regardless of the system code page.
    summary_path = out_dir / f"{stem}_pipex_comparison_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("PIPEX vs annotator comparison summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Cells:   {df.shape[0]:,}\n")
        f.write(f"Markers: {len(markers)}\n")
        if args.auto_cutoff:
            f.write("Cutoffs: per-marker, per-method (Youden's J).\n\n")
        else:
            f.write(f"Cutoffs: otsu3 > {cutoffs['otsu']}, "
                    f"triangle_score > {cutoffs['triangle']}, "
                    f"gmm_prob >= {cutoffs['gmm']} (fixed).\n\n")
        for name, _, _, _, _ in METHODS:
            sub = out[out.method == name]
            mean_auc = sub.roc_auc.mean(skipna=True)
            mean_agree = sub.agreement.mean()
            mean_kappa = sub.cohen_kappa.mean(skipna=True)
            n_degenerate = int(sub["degenerate_reason"].astype(bool).sum())
            n_informative = len(sub) - n_degenerate
            f.write(f"{name.upper()}\n")
            f.write(f"  mean ROC AUC (cutoff-free):       {mean_auc:.3f}\n")
            f.write(f"  mean agreement (at chosen cutoff): {mean_agree:.1%}\n")
            f.write(f"  mean Cohen's kappa:               {mean_kappa:.3f}\n")
            f.write(f"  informative markers (auto-cutoff valid): {n_informative} / {len(sub)}\n")
            if n_degenerate:
                f.write(f"  degenerate markers (fallback to fixed cutoff): {n_degenerate}\n")
                degen = sub[sub["degenerate_reason"].astype(bool)]
                for _, r in degen.iterrows():
                    f.write(f"    {r['marker']:>14s}  AUC={r['roc_auc']:.3f}  "
                            f"reason: {r['degenerate_reason']}\n")
            f.write(f"  bias: PIPEX more inclusive in "
                    f"{int((sub.direction == 'PIPEX more inclusive').sum())} markers, "
                    f"annotator more inclusive in "
                    f"{int((sub.direction == 'annotator more inclusive').sum())} markers, "
                    f"no bias in "
                    f"{int((sub.direction == 'no asymmetric bias').sum())} markers\n")
            worst = sub.nsmallest(5, "roc_auc")[
                ["marker", "roc_auc", "agreement", "cohen_kappa", "direction"]]
            f.write(f"  5 worst markers by ROC AUC:\n")
            for _, r in worst.iterrows():
                f.write(f"    {r.marker:>14s}  AUC={r.roc_auc:.3f}  "
                        f"agree={r.agreement:.1%}  kappa={r.cohen_kappa:.2f}  "
                        f"{r.direction}\n")
            f.write("\n")
    log.info(f"Wrote summary: {summary_path}")

    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
