# Installing PIPEX 2.0

This guide covers a clean installation of PIPEX v2.0 (released April 2026) on Windows with PyCharm. Linux and macOS notes are included where they differ. PIPEX v2.0 must be installed in a separate virtual environment from any earlier version. The two installations cannot share dependencies.

## Requirements

- Windows 10/11, Linux, or macOS
- Python 3.12 (specified by the developer). Earlier versions are untested with v2.0; later versions (3.13, 3.14) lack TensorFlow wheels and will fail at install time.
- 8 GB RAM for images up to 6k resolution, 16 GB for images up to 12k
- ~6 GB free disk space for the environment and dependencies
- GPU optional. A CUDA-capable GPU accelerates Stardist segmentation but is not required.

## Step 1: Confirm Python 3.12 is available

```bash
py -3.12 --version       # Windows
python3.12 --version     # Linux / macOS
```

If 3.12 is not installed, download from https://www.python.org/downloads/ and install. On Windows, do not tick "Add Python to PATH" if another Python version already controls PATH. Use the `py` launcher to invoke 3.12 explicitly. Existing Python installations on the system are not affected.

## Step 2: Download the v2.0 source

Two options:

**Option A, release archive**
1. Open https://github.com/CellProfiling/pipex/releases
2. Under "PIPEX 2.0" (tag `2.0`), expand "Assets"
3. Download the `Source code (zip)` and extract to a permanent location, for example `C:\Users\<username>\PycharmProjects\PIPEX_v2\`

**Option B, git clone**
```bash
git clone --branch 2.0 https://github.com/CellProfiling/pipex.git PIPEX_v2
```

Use a folder name that makes the version explicit. Do not place v1.0 and v2.0 source in the same parent directory without clear naming.

## Step 3: Create an isolated PyCharm project

1. Open PyCharm. File > Open and select the `PIPEX_v2` folder, or File > New Project pointing to the same location.
2. Configure the interpreter: Settings > Project > Python Interpreter > Add Interpreter > Add Local Interpreter > Virtualenv Environment > New.
3. Set Base interpreter to Python 3.12, typically `C:\Python312\python.exe` on Windows or `/usr/bin/python3.12` on Linux.
4. Leave "Inherit global site-packages" unchecked. This is critical to prevent the v2.0 environment from seeing v1.0 packages.
5. Click OK. PyCharm creates a `.venv` folder inside the project.

## Step 4: Install system packages (Linux only)

```bash
sudo apt-get install python3.12-dev libvips python3-opencv gcc python3-tk
```

On macOS, use Homebrew:

```bash
brew install vips opencv
```

Skip this step on Windows.

## Step 5: Activate the venv and install Python dependencies

Open the PyCharm Terminal. Confirm the prompt starts with `(.venv)`. If not:

```bash
.venv\Scripts\activate.bat        # Windows
source .venv/bin/activate         # Linux / macOS
```

Verify the active interpreter is 3.12:

```bash
python --version
```

Then:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Installation typically takes 5 to 15 minutes.

## Step 6: Verify the installation

```bash
python -c "import stardist, tensorflow, skimage, scanpy, anndata; print('PIPEX 2.0 dependencies OK')"
```

Launch the GUI to confirm a working Tk environment:

```bash
python pipexGUI.py        # Windows
sh pipex.sh               # Linux / macOS
```

## Step 7: Run a minimal test

1. Place a small set of TIFF images named after their markers (for example `DAPI.tif`, `CDH1.tif`) in a working folder such as `D:\PIPEX_test\data\`.
2. Edit `pipex_batch_list.txt` in the project root with a minimal segmentation command using v2.0 features:

```
segmentation.py -data=D:/PIPEX_test/data -nuclei_marker=DAPI -nuclei_diameter=20 -nuclei_expansion=10 -measure_markers=DAPI,CDH1 -gmm_min_separation=0.5
```

3. From the project root, with the venv active:

```bash
python pipex.py
```

Output appears in `D:\PIPEX_test\data\analysis\`. A successful v2.0 run produces `cell_data.csv` containing the new continuous score columns: `_local_90`, `_ratio_pixels`, `_otsu3`, `_triangle_score`, and `_gmm_prob` for each marker.

## What is new in v2.0 (relevant to installation and first run)

- **New per-cell continuous scores**: `_triangle_score` and `_gmm_prob` columns are added to `cell_data.csv` for every marker. See Annex 3 of the README for interpretation.
- **New segmentation parameter `-gmm_min_separation`**: controls when the GMM-based score is trusted. Default 0.5.
- **New preprocessing parameter `-tophat_radius`**: morphological top-hat background subtraction. Set to 0 to disable.
- **New analysis parameter `-leiden_res`**: exposes Leiden clustering resolution. Default 0.5; range 0.3 to 0.5 recommended for spatial proteomics panels of 10 to 40 markers.
- **New analysis parameter `-k_estimation`**: produces elbow, inertia, and silhouette plots over k = 1 to 19 to guide kmeans cluster number selection.
- **Multiple suffixes accepted in `-use_bin`**: comma-separated list, allowing different score types per marker in the same clustering run.
- **Multiple parallel cluster refinements**: distinct `ref_id` values in `cell_types.csv` produce independent refinement outputs in one pass.
- **Neighborhood analysis**: activated via `-neigh_cluster_id`, produces neighbor composition heatmaps and stacked bar charts at configurable k values, plus DBSCAN-based spatial distribution classification per cell type.
- **LMD export from `generate_filtered_masks.py`**: produces XML cutting files compatible with Leica Laser Microdissection software, with shape dilation, contour smoothing, and cutting path optimization.
- **TissUUmaps `-launch` parameter**: starts a local HTTP server and opens the export in the default browser automatically.
- **`extra/` folder**: standalone post-processing scripts (`custom_binarization_example.py`, `match_segmentation_ids.py`) that operate on `cell_data.csv` without requiring the full pipeline.

If migrating from v1.0, do not assume `cell_data.csv` schemas are interchangeable. v2.0 adds columns; downstream scripts written for v1.0 may need updating.

## Common issues

- **`Could not find a version that satisfies the requirement tensorflow`**: Python version mismatch. Confirm `python --version` inside the active venv returns 3.12.x.
- **Stardist installation fails on Windows**: ensure Microsoft C++ Build Tools are installed. Download from https://visualstudio.microsoft.com/visual-cpp-build-tools/.
- **`ModuleNotFoundError: No module named 'stardist'`**: the virtual environment is not active. Reactivate before running any script.
- **OpenCV import errors on Linux**: install `python3-opencv` via apt.
- **`_gmm_prob` column is all NaN for a marker**: the GMM fit was not trusted because the two components were too close. Lower `-gmm_min_separation` (for example to 0.3) if the marker has weak but real bimodality, or accept the NaN if the marker distribution is genuinely unimodal.
- **Out-of-memory errors during segmentation on large images**: on Linux, prepend `swap 8` to `pipex_batch_list.txt`. On Windows or macOS, use `max_res 8000` to downscale internally.

## Coexistence with PIPEX 1.0

Each version of PIPEX must run from its own folder with its own `.venv`. The two cannot share an environment.

- Use separate PyCharm projects, one per version. File > Recent Projects to switch.
- Never activate two venvs in the same terminal session. The `(.venv)` prefix in the prompt indicates which environment is active. Use `where python` (Windows) or `which python` (Linux/macOS) to confirm.
- Use separate data folders per analysis. v2.0 adds columns to `cell_data.csv`; mixing outputs across versions will break downstream code that expects a fixed schema.
- Pin dependency versions for reproducibility:

```bash
pip freeze > requirements_v2_locked.txt
```

Commit this file alongside any analysis project that depends on a specific v2.0 install state.

## Reference

- PIPEX repository: https://github.com/CellProfiling/pipex
- v2.0 release: https://github.com/CellProfiling/pipex/releases/tag/2.0
- Full documentation: see `README.md` and `PIPEX.pdf` in the project root.
- Annex 3 (continuous marker scores) and Annex 5 (neighborhood analysis) in the README cover the v2.0 additions in detail.
