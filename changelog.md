Changelog
=========

### Version 2.0

#### Segmentation

- **`{marker}_triangle_score` continuous score**
  Signed distance of each cell's normalized mean intensity from the global Triangle threshold, computed per marker. Follows the same scale and sign convention as the existing `_otsu3` score and is stored as a new column in `cell_data.csv`. Particularly useful for markers with strongly right-skewed distributions, where the Triangle method tends to place a higher threshold than Otsu; comparing the two scores can flag markers that are genuinely difficult to threshold.

- **`{marker}_gmm_prob` continuous score**
  Posterior probability (0–1) that a cell belongs to the positive population, derived from a two-component Gaussian Mixture Model fit to the normalized marker image. Unlike the distance-based scores, this captures the shape and width of both populations, so a cell sitting at the boundary gets ~0.5 regardless of how spread out the distributions are. Stored as a new column in `cell_data.csv`; set to `NaN` when the GMM fit is not trusted (see `-gmm_min_separation`).

- **`-gmm_min_separation` parameter**
  Controls the minimum separation between the two GMM components (in combined standard deviation units) required to trust the fit. When the marker distribution is not clearly bimodal and the separation falls below this value, `_gmm_prob` is set to `NaN` for all cells of that marker rather than reporting a potentially misleading probability. Default is 0.5; increase it to demand a cleaner bimodal shape before trusting the score.

#### Preprocessing

- **`-tophat_radius` parameter**
  Applies a morphological top-hat background subtraction to the marker image before segmentation, using a circular structuring element of the specified radius. This removes slowly-varying background signal such as autofluorescence and uneven staining while preserving cell-scale structures. Set to 0 (default) to disable; a good starting point is roughly 1.5× the expected cell radius in pixels.

#### Analysis

- **`-leiden_res` parameter**
  Exposes the Leiden clustering resolution as a user-configurable parameter. Higher values produce more, finer-grained clusters; lower values produce fewer, broader ones. The typical range for spatial proteomics panels of 10–40 markers is 0.1–2.0, with 0.3–0.5 recommended as a starting point. Default is 0.5.

- **`-k_estimation` parameter**
  When kmeans clustering is enabled, runs 20 successive fits (k = 1 to 19) and produces three diagnostic plots — distortion (elbow method), inertia, and silhouette score per k — to help choose the number of clusters before committing to a final run. Higher silhouette scores indicate better-separated clusters.

- **Multiple binarized suffixes via `-use_bin`**
  Previously, `-use_bin` accepted a single column suffix and substituted it for all marker intensity inputs to clustering. It now accepts a comma-separated list of suffixes, applying each one to its corresponding marker. This allows mixing score types across markers in the same analysis run.

- **Multiple parallel cluster refinements**
  The cluster refinement step (`cell_types.csv`) now supports multiple independent refinements in a single run by assigning distinct `ref_id` values to groups of rules. Each group produces its own output column and JSON report, runs in parallel against the same clustering results, and does not filter the output of previous groups. A typical use is a first strict pass for well-defined populations and a second looser pass for ambiguous clusters.

- **Neighborhood analysis**
  New full-featured spatial analysis activated via `-neigh_cluster_id`. For each k value in `-neigh_k_values` (default: 1, 5, 10), PIPEX computes the cell type composition of the k nearest neighbors of every cell, producing a heatmap and stacked bar chart per k. Independently, DBSCAN clustering is run on each cell type's spatial coordinates to classify it as scattered sparsely, clustered sparsely, or clustered densely — with the density threshold controlled by `-neigh_density_threshold`. Results are saved to `cell_data.csv` and a dedicated spatial distribution CSV. See Annex 5 in the README for full details.

#### Filter masks

- **LMD export**
  New output mode for `generate_filtered_masks.py` that produces an XML cutting file compatible with Leica's Laser Microdissection software. Four parameters control the output geometry: `-shape_dilation` expands each cell outline by a given number of pixels, `-convolution_smoothing` controls contour smoothness, `-path_optimization` selects the cutting path order strategy (none, Hilbert, or greedy), and `-distance_heuristic` merges nearby shapes into a single cutting group to reduce stage movements.

#### TissUUmaps export

- **`-launch` parameter for `generate_tissuumaps.py`**
  When set to `yes`, automatically starts a local HTTP server serving the `TissUUmaps_webexport` folder after export completes and opens the result in the default browser. The server runs on the first available port starting at 8080 and keeps running until the process is interrupted with `Ctrl+C`. Requires `include_html=yes` to have generated the webexport folder first.

#### Extra scripts

- **`extra/` folder**
  New folder with standalone Python scripts for post-processing tasks outside the core pipeline. Scripts work directly from `cell_data.csv` — no images required — and are designed as ready-to-adapt examples.

- **`extra/custom_binarization_example.py`**
  Re-applies marker +/- binarization to an existing `cell_data.csv` using the continuous scores already computed by segmentation (`_gmm_prob`, `_triangle_score`, `_otsu3`). Useful to re-tune thresholds or regenerate binarization columns without re-running the full pipeline.

- **`extra/match_segmentation_ids.py`**
  Matches cells between two `cell_data.csv` files from different segmentation runs on the same image (e.g. whole-cell vs. nuclear) by nearest-neighbour spatial proximity. Writes a `matched_cell_id` column into the target file so the two tables can be joined downstream.

#### Documentation

- **Annex 3: Continuous marker scores**
  New annex in the README documenting all five per-cell continuous scores (`_local_90`, `_ratio_pixels`, `_otsu3`, `_triangle_score`, `_gmm_prob`), including their scale, interpretation, and guidance on which score to use for different marker types and downstream tasks.

- **Annex 5: Neighborhood cell type analysis**
  New annex in the README (renumbered from the former Annex 4 slot) providing a full explanation of the neighbor composition analysis and the spatial distribution classification methodology, including the DBSCAN approach and the density threshold definition.