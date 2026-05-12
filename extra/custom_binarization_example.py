"""
custom_binarization_example.py
================================
Applies marker +/- binarization to the continuous scores already computed
by the PIPEX segmentation step and stored in cell_data.csv.

No image processing is required, this works entirely from the CSV.
For each detected marker a new column {marker}_custom_binarization is added
with values '{marker}+' or '{marker}-'.

Optionally generates per-marker threshold diagnostic plots (PROMINENT style):
intensity histogram with GMM, Triangle, and Otsu thresholds overlaid, the
fitted GMM mixture curve, the GMM separation score, and the cascade method
that was actually used for the call.

Binarization cascade (mirrors the logic in new_binarizing_methods.py):
  1. _gmm_prob >= GMM_PROB_THRESHOLD   (if column present and not all-NaN)
  2. _triangle_score >= 0              (fallback: signed distance above threshold)
  3. _otsu3 >= 0                       (final fallback: signed distance above threshold)

Overcall protection: if on a re-run the positive fraction increases more than
OVERCALL_FOLD times relative to the previous labels AND exceeds OVERCALL_PCT%
of all cells, the existing labels are kept and a warning is printed.
On a first run there are no previous labels, so this check is a no-op.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

# ── Configuration ──────────────────────────────────────────────────────────────

CELL_DATA_PATH = os.path.join("data", "analysis", "cell_data.csv")
DIAG_DIR       = os.path.join("data", "analysis", "threshold_diagnostics")
MAKE_PLOTS     = True

# Minimum GMM posterior probability to call a cell positive.
# 0.5 = equal probability of belonging to either population.
GMM_PROB_THRESHOLD = 0.5

# Overcall protection: both conditions must be exceeded to trigger.
OVERCALL_FOLD = 5.0   # fold increase vs. previous run labels
OVERCALL_PCT  = 15.0  # percentage of all cells labelled positive

# Optional sample column. If present in the CSV, diagnostic plots are
# faceted by sample (one subplot per sample, in a grid).
SAMPLE_COL = "sample_id"

# ── Score cascade definition ───────────────────────────────────────────────────

# Each entry: (column suffix, decision threshold).
# The cascade tries them in order and uses the first one with valid data.
# For _triangle_score and _otsu3 a threshold of 0 means "above the reference
# threshold computed at segmentation time", which is the natural cut point.
SCORE_CASCADE = [
    ("_gmm_prob",       GMM_PROB_THRESHOLD),
    ("_triangle_score", 0.0),
    ("_otsu3",          0.0),
]

# ── Helpers ────────────────────────────────────────────────────────────────────

def detect_markers(df):
    """Return markers that have at least one recognised score column in the CSV."""
    score_suffixes = {sfx for sfx, _ in SCORE_CASCADE} | {"_ratio_pixels", "_local_90"}
    candidates = set()
    for col in df.columns:
        for sfx in score_suffixes:
            if col.endswith(sfx):
                candidates.add(col[:-len(sfx)])
    return sorted(m for m in candidates if m in df.columns)


def check_overcall(n_pos_new, n_pos_old, n_total):
    """Return (is_overcall, fold) compared to previous label count."""
    if n_pos_old <= 0:
        return False, 1.0
    fold = n_pos_new / max(n_pos_old, 1)
    pct  = n_pos_new / n_total * 100
    return (fold > OVERCALL_FOLD and pct > OVERCALL_PCT), fold


def binarize_marker(marker, df):
    """
    Walk the score cascade and return (pos_mask, method_name) for one marker.
    pos_mask is a boolean Series aligned to df.index.
    Returns (None, None) if no usable score column exists.
    """
    for suffix, threshold in SCORE_CASCADE:
        col = marker + suffix
        if col not in df.columns:
            continue
        valid = df[col].notna()
        if not valid.any():
            continue
        pos_mask = pd.Series(False, index=df.index)
        pos_mask[valid] = df.loc[valid, col] >= threshold
        return pos_mask, suffix.lstrip("_")

    return None, None


# ── Diagnostic plotting ────────────────────────────────────────────────────────

def recover_pipex_thresholds(df_sub, marker):
    """
    Recover the per-marker thresholds PIPEX used, from the stored score columns.

      triangle threshold = intensity - _triangle_score   (signed distance)
      otsu     threshold = intensity - _otsu3            (signed distance)
      gmm      threshold = lowest intensity where _gmm_prob >= 0.5

    Returns a dict with any subset of {'gmm', 'triangle', 'otsu'}.
    """
    out = {}
    if marker not in df_sub.columns:
        return out

    tri_col = f"{marker}_triangle_score"
    if tri_col in df_sub.columns:
        diffs = (df_sub[marker] - df_sub[tri_col]).dropna()
        if len(diffs):
            out["triangle"] = float(diffs.median())

    otsu_col = f"{marker}_otsu3"
    if otsu_col in df_sub.columns:
        diffs = (df_sub[marker] - df_sub[otsu_col]).dropna()
        if len(diffs):
            out["otsu"] = float(diffs.median())

    gmm_col = f"{marker}_gmm_prob"
    if gmm_col in df_sub.columns:
        sub = df_sub[[marker, gmm_col]].dropna().sort_values(marker)
        if len(sub):
            pos = sub[sub[gmm_col] >= GMM_PROB_THRESHOLD]
            if len(pos):
                out["gmm"] = float(pos[marker].iloc[0])

    return out


def fit_gmm_overlay(values, n_components=2):
    """
    Fit a 2-component GMM for the visual overlay only (not used for calls).
    Returns (xs, total_pdf, separation) or (None, None, None) on failure.
    Separation = |mu1 - mu0| / sqrt(sd0^2 + sd1^2).
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 50 or np.ptp(v) == 0:
        return None, None, None
    try:
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(v.reshape(-1, 1))
        order   = np.argsort(gmm.means_.flatten())
        means   = gmm.means_.flatten()[order]
        sds     = np.sqrt(gmm.covariances_.flatten())[order]
        weights = gmm.weights_[order]
        xs = np.linspace(v.min(), np.percentile(v, 99.5), 500)
        pdf = sum(
            w * (1.0 / (s * np.sqrt(2 * np.pi))) *
            np.exp(-0.5 * ((xs - m) / s) ** 2)
            for w, m, s in zip(weights, means, sds)
        )
        sep = abs(means[1] - means[0]) / np.sqrt(sds[0] ** 2 + sds[1] ** 2)
        return xs, pdf, float(sep)
    except Exception:
        return None, None, None


def _plot_one_panel(ax, values, marker, thresholds, method_used, panel_title):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        ax.set_axis_off()
        return

    xmax = np.percentile(v, 99.5)
    plot_v = v[v <= xmax] if xmax > v.min() else v

    ax.hist(plot_v, bins=80, density=True, alpha=0.55,
            color="steelblue", edgecolor="none")

    xs, pdf, sep = fit_gmm_overlay(v)
    if xs is not None:
        ax.plot(xs, pdf, color="salmon", linewidth=1.2)

    styles = {
        "gmm":      dict(color="red",    linestyle="-",  linewidth=2.2),
        "triangle": dict(color="orange", linestyle="--", linewidth=1.6),
        "otsu":     dict(color="purple", linestyle=":",  linewidth=1.8),
    }
    label_short = {"gmm": "gmm", "triangle": "tri", "otsu": "otsu"}
    for name, t in thresholds.items():
        ax.axvline(t, label=f"{label_short[name]}: {t:.2f}", **styles[name])

    sep_str = f", sep={sep:.2f}" if sep is not None else ""
    ax.set_title(f"{panel_title}\n({method_used}{sep_str})", fontsize=10)
    ax.set_xlabel(marker)
    ax.set_ylabel("density")
    ax.legend(fontsize=8, frameon=True)


def make_marker_diagnostic(df, marker, method_used, out_dir):
    """Save a per-marker PNG. If SAMPLE_COL is present, faceted by sample."""
    os.makedirs(out_dir, exist_ok=True)

    if SAMPLE_COL in df.columns:
        samples = sorted(df[SAMPLE_COL].dropna().unique().tolist())
    else:
        samples = [None]

    n = len(samples)
    ncols = min(4, n) if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)
    axes = axes.flatten()

    for i, sid in enumerate(samples):
        sub = df if sid is None else df[df[SAMPLE_COL] == sid]
        thr = recover_pipex_thresholds(sub, marker)
        panel_title = marker if sid is None else str(sid)
        _plot_one_panel(axes[i], sub[marker], marker, thr,
                        method_used, panel_title)

    for j in range(len(samples), len(axes)):
        axes[j].set_axis_off()

    fig.suptitle(f"{marker}, threshold diagnostics",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    out_path = os.path.join(out_dir, f"{marker}_thresholds.png")
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if not os.path.exists(CELL_DATA_PATH):
        raise FileNotFoundError(f"cell_data.csv not found at: {CELL_DATA_PATH}")

    print(f"Reading {CELL_DATA_PATH} ...", flush=True)
    df = pd.read_csv(CELL_DATA_PATH)
    print(f"  {len(df):,} cells loaded", flush=True)

    markers = detect_markers(df)
    if not markers:
        print("No markers with score columns found, nothing to do.")
        return
    print(f"  {len(markers)} markers detected: {', '.join(markers)}\n", flush=True)

    methods_used = {}

    for marker in markers:
        out_col = f"{marker}_custom_binarization"
        pos_label = f"{marker}+"
        neg_label = f"{marker}-"

        # Snapshot previous labels for overcall comparison on re-runs
        n_pos_old = 0
        if out_col in df.columns:
            n_pos_old = int(df[out_col].astype(str).str.endswith("+").sum())

        pos_mask, method = binarize_marker(marker, df)
        if pos_mask is None:
            print(f"  [{marker}] no usable score column, skipped")
            continue

        n_pos   = int(pos_mask.sum())
        n_total = len(df)
        pct_pos = n_pos / n_total * 100

        # Overcall check (only meaningful when previous labels exist)
        if n_pos_old > 0:
            is_overcall, fold = check_overcall(n_pos, n_pos_old, n_total)
            if is_overcall:
                print(f"  [{marker}] WARNING: overcall detected "
                      f"({fold:.1f}x vs previous run, {pct_pos:.1f}% positive), "
                      f"previous labels kept")
                methods_used[marker] = "kept-previous"
                continue

        df[out_col] = pos_mask.map({True: pos_label, False: neg_label})
        methods_used[marker] = method
        print(f"  [{marker}] {n_pos:,} positive ({pct_pos:.1f}%), method: {method}")

    df.to_csv(CELL_DATA_PATH, index=False)
    print(f"\nSaved: {CELL_DATA_PATH}", flush=True)

    # ── Diagnostic plots ───────────────────────────────────────────────────────

    if MAKE_PLOTS:
        print(f"\nWriting threshold diagnostics to {DIAG_DIR} ...", flush=True)
        for marker in markers:
            method = methods_used.get(marker)
            if method is None or method == "kept-previous":
                continue
            if marker not in df.columns:
                print(f"  [{marker}] raw intensity column missing, plot skipped")
                continue
            try:
                out = make_marker_diagnostic(df, marker, method, DIAG_DIR)
                print(f"  [{marker}] {os.path.basename(out)}")
            except Exception as e:
                print(f"  [{marker}] plot failed: {e}")


if __name__ == "__main__":
    main()
