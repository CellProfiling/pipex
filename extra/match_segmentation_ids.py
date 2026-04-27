"""
match_segmentation_ids.py
=========================
Matches cells between two PIPEX cell_data.csv files based on spatial proximity
(x, y centroid coordinates) and writes a matched_cell_id column into the target file.

Typical use case: you ran PIPEX twice on the same image — once with whole-cell
segmentation (source) and once with nuclear segmentation (target) — and you want
to know which whole-cell ID corresponds to each nucleus so you can merge marker
data across both runs.

Usage
-----
Edit the three variables in the "Configuration" block below, then run:

    python match_segmentation_ids.py

Output
------
The target CSV is overwritten in-place with a new column `matched_cell_id`
inserted immediately after `cell_id`.  A value of 0 means no source cell was
found within `minimum_distance_allowed` pixels of that target centroid.

Dependencies
------------
    pip install pandas scipy
"""

import pandas as pd
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Configuration — edit these three lines
# ---------------------------------------------------------------------------
SOURCE_FILE              = "../data/analysis/cell_data.csv"       # whole-cell segmentation
TARGET_FILE              = "../data/analysis/cell_data_nuclear.csv"  # nuclear segmentation
MINIMUM_DISTANCE_ALLOWED = 10  # pixels; matches farther than this are rejected
# ---------------------------------------------------------------------------


def match_segmentation_ids(source_file, target_file, minimum_distance_allowed):
    print(f"Reading source file : {source_file}")
    source = pd.read_csv(source_file)

    print(f"Reading target file : {target_file}")
    target = pd.read_csv(target_file)

    # Build a KD-tree from the source (x, y) centroids for fast nearest-neighbour lookup.
    source_tree = cKDTree(source[["x", "y"]].values)

    # Query: for every target centroid find the closest source centroid.
    distances, indices = source_tree.query(target[["x", "y"]].values, k=1)

    # Map to source cell_id; set to 0 when the nearest neighbour is too far away.
    matched_ids = source["cell_id"].iloc[indices].values
    matched_ids[distances > minimum_distance_allowed] = 0

    # Insert matched_cell_id right after cell_id (drop first if re-running).
    if "matched_cell_id" in target.columns:
        target.drop(columns=["matched_cell_id"], inplace=True)
    insert_position = target.columns.get_loc("cell_id") + 1
    target.insert(insert_position, "matched_cell_id", matched_ids)

    # Report match statistics before saving.
    n_total   = len(target)
    n_matched = (target["matched_cell_id"] != 0).sum()
    print(f"Matched {n_matched:,} / {n_total:,} cells "
          f"({100 * n_matched / n_total:.1f}%) within {minimum_distance_allowed} px")

    print(f"Writing result to   : {target_file}")
    target.to_csv(target_file, index=False)
    print("Done.")


if __name__ == "__main__":
    match_segmentation_ids(SOURCE_FILE, TARGET_FILE, MINIMUM_DISTANCE_ALLOWED)
