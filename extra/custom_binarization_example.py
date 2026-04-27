"""
custom_binarization_example.py
================================
Applies marker +/- binarization to the continuous scores already computed
by the PIPEX segmentation step and stored in cell_data.csv.

No image processing is required — this works entirely from the CSV.
For each detected marker a new column {marker}_custom_binarization is added
with values '{marker}+' or '{marker}-'.

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

# ── Configuration ──────────────────────────────────────────────────────────────

CELL_DATA_PATH = os.path.join("data", "analysis", "cell_data.csv")

# Minimum GMM posterior probability to call a cell positive.
# 0.5 = equal probability of belonging to either population.
GMM_PROB_THRESHOLD = 0.5

# Overcall protection: both conditions must be exceeded to trigger.
OVERCALL_FOLD = 5.0   # fold increase vs. previous run labels
OVERCALL_PCT  = 15.0  # percentage of all cells labelled positive

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


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if not os.path.exists(CELL_DATA_PATH):
        raise FileNotFoundError(f"cell_data.csv not found at: {CELL_DATA_PATH}")

    print(f"Reading {CELL_DATA_PATH} ...", flush=True)
    df = pd.read_csv(CELL_DATA_PATH)
    print(f"  {len(df):,} cells loaded", flush=True)

    markers = detect_markers(df)
    if not markers:
        print("No markers with score columns found — nothing to do.")
        return
    print(f"  {len(markers)} markers detected: {', '.join(markers)}\n", flush=True)

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
            print(f"  [{marker}] no usable score column — skipped")
            continue

        n_pos   = int(pos_mask.sum())
        n_total = len(df)
        pct_pos = n_pos / n_total * 100

        # Overcall check (only meaningful when previous labels exist)
        if n_pos_old > 0:
            is_overcall, fold = check_overcall(n_pos, n_pos_old, n_total)
            if is_overcall:
                print(f"  [{marker}] WARNING: overcall detected "
                      f"({fold:.1f}x vs previous run, {pct_pos:.1f}% positive) "
                      f"— previous labels kept")
                continue

        df[out_col] = pos_mask.map({True: pos_label, False: neg_label})
        print(f"  [{marker}] {n_pos:,} positive ({pct_pos:.1f}%) — method: {method}")

    df.to_csv(CELL_DATA_PATH, index=False)
    print(f"\nSaved: {CELL_DATA_PATH}", flush=True)


if __name__ == "__main__":
    main()