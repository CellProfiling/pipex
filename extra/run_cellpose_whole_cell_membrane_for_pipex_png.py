# run_cellpose_wholecell_for_pipex.py  (PNG input version)
# ============================================================================
# CellPose WHOLE-CELL segmentation for PIPEX, using DAPI (nucleus) +
# a COMPOSITE membrane channel built from multiple markers.
# INPUT IMAGES: .png (also supports .tif/.tiff).
# Mask output stays .tif (int32) for PIPEX; overlays are .png.
#
#     conda activate cellpose
#     python run_cellpose_wholecell_for_pipex.py
#
# PIPEX Segmentation GUI:
#     Custom segmentation      = cellpose_mask.tif   (int32, cv2-readable)
#     Custom segmentation type = Full        (whole-cell)
#     Nuclei expansion         = 1
#     Measure markers          = (paste string from Step 5)
# ============================================================================

import os
import glob
import numpy as np
import tifffile
from skimage.transform import resize
from skimage import exposure
from skimage.io import imread as sk_imread   # reads PNG and TIFF

from cellpose import models
try:
    import cellpose
    CELLPOSE_VERSION = cellpose.version
except Exception:
    CELLPOSE_VERSION = "unknown"


# ============================================================================
# CONFIG
# ============================================================================
DATA_DIR = os.path.dirname(os.path.abspath(__file__))

# Input image extension(s). PNG first, but TIFF also allowed.
INPUT_EXTENSIONS = (".png", ".tif", ".tiff")
INPUT_EXT = ".png"          # extension used to build marker filenames

NUCLEAR_FILE = "DAPI" + INPUT_EXT
IMAGE_MPP    = 0.3125

# --- Composite membrane: markers that outline cell boundaries ---------------
# Names WITHOUT extension. Only files that exist are used.
MEMBRANE_MARKERS = [
    "b-Catenin1", "E-cadherin", "Pan-Cytokeratin", "CD45", "Caveolin"
]
MEMBRANE_COMBINE = "max"     # "max" (recommended) or "sum"

# --- CellPose parameters ----------------------------------------------------
CELL_DIAMETER      = 25
MODEL_TYPE         = "cyto3"
FLOW_THRESHOLD     = 0.4
CELLPROB_THRESHOLD = 0.0

# Anti-tiling / anti-stripe
TILE_OVERLAP = 0.5
AUGMENT      = True
RESAMPLE     = True

USE_GPU  = False
USE_CLAHE = False

# --- Boundary overlay (QC image) --------------------------------------------
OVERLAY_COLOR      = (1, 0, 0)   # RED
OVERLAY_LINE_THICK = 1
OVERLAY_ON_MEMBRANE = True

# --- PIPEX / output ---------------------------------------------------------
COMPARTMENT_LABEL = "Full"
PIPEX_MAX_RESOLUTION = int(os.environ.get("PIPEX_MAX_RESOLUTION", 32768))
PRE_DOWNSCALE = False

# Mask outputs stay TIFF (int32) for PIPEX; overlays are PNG.
OUTPUT_NPY = os.path.join(DATA_DIR, "cellpose_mask.npy")
OUTPUT_TIF = os.path.join(DATA_DIR, "cellpose_mask.tif")
COMPOSITE_TIF = os.path.join(DATA_DIR, "composite_membrane_QC.tif")
OUTPUT_OVERLAY = os.path.join(DATA_DIR, "cellpose_segmentation_boundaries.png")
OUTPUT_OVERLAY_MEM = os.path.join(DATA_DIR, "cellpose_boundaries_on_membrane.png")

# --- Step 5 options ---------------------------------------------------------
EXCLUDE_NUCLEAR_INPUT  = False
EXCLUDE_MEMBRANE_INPUT = False
MEMBRANE_FILE = "b-Catenin1" + INPUT_EXT
ALWAYS_EXCLUDE = {"cellpose_mask", "mesmer_mask", "mesmer_mask_nuclear",
                  "composite_membrane_QC",
                  "segmentation_mask", "segmentation_binary_mask"}
MASK_PREFIXES = ("cellpose_mask", "mesmer_mask", "composite_membrane_QC")


# ============================================================================
# Image reading helper (handles PNG + TIFF, forces 2D grayscale)
# ============================================================================
def read_image_2d(path):
    """Read PNG/TIFF as a 2D grayscale array. Drops alpha / averages RGB."""
    img = sk_imread(path)
    img = np.asarray(img)
    if img.ndim == 3:
        # RGB or RGBA PNG -> take luminance (or first channel if identical)
        if img.shape[2] == 4:      # RGBA
            img = img[..., :3]
        # if all channels equal (grayscale saved as RGB) just take one
        if img.shape[2] == 3 and np.allclose(img[..., 0], img[..., 1]) \
                and np.allclose(img[..., 1], img[..., 2]):
            img = img[..., 0]
        else:
            img = img.mean(axis=2)  # convert RGB to grayscale
    return img


def image_shape_2d(path):
    """Return 2D (H, W) shape of a PNG/TIFF, ignoring channels."""
    img = read_image_2d(path)
    return img.shape


# ============================================================================
# Helpers
# ============================================================================
def normalize_channel(img, low=1, high=99, use_clahe=False):
    img = img.astype("float32")
    lo, hi = np.percentile(img, [low, high])
    if hi > lo:
        img = np.clip((img - lo) / (hi - lo), 0.0, 1.0)
    if use_clahe:
        img = exposure.equalize_adapthist(img, clip_limit=0.01)
    return img.astype("float32")


def save_boundary_overlay(bg_img, mask, out_path,
                          boundary_color=(1, 0, 0), line_thickness=1,
                          low=1, high=99.5):
    from skimage.segmentation import find_boundaries
    from skimage.io import imsave

    d = bg_img.astype("float32")
    lo, hi = np.percentile(d, [low, high])
    if hi > lo:
        d = np.clip((d - lo) / (hi - lo), 0, 1)
    rgb = np.stack([d, d, d], axis=-1)

    boundaries = find_boundaries(mask, mode="outer")
    if line_thickness > 1:
        from skimage.morphology import binary_dilation, disk
        boundaries = binary_dilation(boundaries, disk(line_thickness - 1))
    rgb[boundaries] = boundary_color
    imsave(out_path, (rgb * 255).astype(np.uint8))
    print("Saved boundary overlay:", out_path)


def build_composite_membrane(data_dir, marker_names, ref_shape, combine="max"):
    ref_H, ref_W = ref_shape
    used, missing, stack = [], [], []
    for m in marker_names:
        path = os.path.join(data_dir, m + INPUT_EXT)
        if not os.path.exists(path):
            missing.append(m); continue
        img = read_image_2d(path)
        if img.ndim != 2 or img.shape != (ref_H, ref_W):
            missing.append(m + " (shape)"); continue
        stack.append(normalize_channel(img))
        used.append(m)
    if not stack:
        raise SystemExit(">>> ERROR: no valid membrane markers found for composite.")
    arr = np.stack(stack, axis=0)
    composite = arr.max(axis=0) if combine == "max" else arr.sum(axis=0)
    composite = normalize_channel(composite)
    print(f">>> Composite membrane from {len(used)} markers ({combine}): {used}")
    if missing:
        print(f">>> Skipped (missing/shape): {missing}")
    return composite


def verify_marker_shapes(data_dir, reference_file, extensions=INPUT_EXTENSIONS):
    ref_path = os.path.join(data_dir, reference_file)
    if not os.path.exists(ref_path):
        print(f">>> ERROR: reference file not found: {ref_path}")
        return False, None
    ref_shape = image_shape_2d(ref_path)
    if len(ref_shape) != 2:
        print(f">>> ERROR: reference {reference_file} is not 2D (shape={ref_shape}).")
        return False, None
    ref_H, ref_W = ref_shape
    print(f">>> Reference: {reference_file}  shape=({ref_H}, {ref_W})")
    print("-" * 60)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(data_dir, f"*{ext}")))
    files = sorted(set(f for f in files
                       if not os.path.basename(f).startswith(MASK_PREFIXES)))
    all_ok = True
    for f in files:
        name = os.path.basename(f)
        try:
            shape = image_shape_2d(f)
        except Exception as e:
            print(f"  [READ ERROR]    {name}: {e}"); all_ok = False; continue
        if len(shape) != 2:
            print(f"  [NOT 2D]        {name}: shape={shape}"); all_ok = False
        elif shape != (ref_H, ref_W):
            print(f"  [SIZE MISMATCH] {name}: shape={shape}"); all_ok = False
        else:
            print(f"  [OK]            {name}: shape={shape}")
    print("-" * 60)
    print(">>> All images 2D and matching." if all_ok else ">>> PROBLEMS FOUND.")
    return all_ok, (ref_H, ref_W)


def compute_pipex_scale_factor(h, w, max_res):
    if h <= max_res and w <= max_res:
        return 1
    factor, i = 0, 2
    while factor == 0:
        factor = i if max(h, w) / i <= max_res else 0
        if factor == 0:
            i *= 2
    return factor


def pipex_downscale(np_img, factor):
    if factor <= 1:
        return np_img
    return resize(np_img, (np_img.shape[0]//factor, np_img.shape[1]//factor),
                  order=0, preserve_range=True, anti_aliasing=False)


def generate_measure_markers(data_dir, reference_shape, extensions=INPUT_EXTENSIONS):
    ref_H, ref_W = reference_shape
    exclude = set(ALWAYS_EXCLUDE)
    if EXCLUDE_NUCLEAR_INPUT:  exclude.add(os.path.splitext(NUCLEAR_FILE)[0])
    if EXCLUDE_MEMBRANE_INPUT: exclude.add(os.path.splitext(MEMBRANE_FILE)[0])
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(data_dir, f"*{ext}")))
    files = sorted(set(f for f in files
                       if not os.path.basename(f).startswith(MASK_PREFIXES)))
    markers, skipped = [], []
    for f in files:
        name = os.path.splitext(os.path.basename(f))[0]
        if name in exclude:
            skipped.append((name, "excluded")); continue
        try:
            shape = image_shape_2d(f)
        except Exception:
            skipped.append((name, "unreadable")); continue
        if len(shape) != 2 or shape != (ref_H, ref_W):
            skipped.append((name, f"shape {shape}")); continue
        markers.append(name)
    return markers, skipped


def run_cellpose_2ch(nuc, mem):
    major = 0
    try:
        major = int(str(CELLPOSE_VERSION).split(".")[0])
    except Exception:
        pass
    print(f">>> CellPose version: {CELLPOSE_VERSION} (major={major})")
    print(f">>> diameter={CELL_DIAMETER}, tile_overlap={TILE_OVERLAP}, "
          f"augment={AUGMENT}, resample={RESAMPLE}")

    img = np.stack([mem, nuc], axis=-1)   # chan1=membrane, chan2=nucleus
    anti = dict(tile_overlap=TILE_OVERLAP, augment=AUGMENT, resample=RESAMPLE)

    if major >= 4:
        model = models.CellposeModel(gpu=USE_GPU)
        try:
            out = model.eval(img, diameter=CELL_DIAMETER,
                             flow_threshold=FLOW_THRESHOLD,
                             cellprob_threshold=CELLPROB_THRESHOLD, **anti)
        except TypeError:
            out = model.eval(img, diameter=CELL_DIAMETER,
                             flow_threshold=FLOW_THRESHOLD,
                             cellprob_threshold=CELLPROB_THRESHOLD)
        masks = out[0]
    else:
        try:
            model = models.Cellpose(gpu=USE_GPU, model_type=MODEL_TYPE)
            try:
                out = model.eval(img, diameter=CELL_DIAMETER, channels=[1, 2],
                                 flow_threshold=FLOW_THRESHOLD,
                                 cellprob_threshold=CELLPROB_THRESHOLD, **anti)
            except TypeError:
                out = model.eval(img, diameter=CELL_DIAMETER, channels=[1, 2],
                                 flow_threshold=FLOW_THRESHOLD,
                                 cellprob_threshold=CELLPROB_THRESHOLD)
            masks = out[0]
        except Exception as e:
            print(f">>> fallback CellposeModel ({e})")
            model = models.CellposeModel(gpu=USE_GPU, model_type=MODEL_TYPE)
            out = model.eval(img, diameter=CELL_DIAMETER, channels=[1, 2],
                             flow_threshold=FLOW_THRESHOLD,
                             cellprob_threshold=CELLPROB_THRESHOLD)
            masks = out[0]
    return np.asarray(masks).astype(np.int32)


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70); print("STEP 1: Verifying image shapes"); print("=" * 70)
    ok, ref_shape = verify_marker_shapes(DATA_DIR, NUCLEAR_FILE)
    if not ok:
        raise SystemExit(">>> Aborting: fix image shape problems first.")
    H, W = ref_shape

    print("\n" + "=" * 70); print("STEP 1b: Building composite membrane"); print("=" * 70)
    dapi_raw = read_image_2d(os.path.join(DATA_DIR, NUCLEAR_FILE))
    nuc = normalize_channel(dapi_raw, use_clahe=USE_CLAHE)
    mem = build_composite_membrane(DATA_DIR, MEMBRANE_MARKERS, (H, W),
                                   combine=MEMBRANE_COMBINE)
    tifffile.imwrite(COMPOSITE_TIF, (mem * 65535).astype(np.uint16))
    print(">>> Saved composite membrane QC image:", COMPOSITE_TIF)

    print("\n" + "=" * 70); print("STEP 2: Downscale guard"); print("=" * 70)
    sf = compute_pipex_scale_factor(H, W, PIPEX_MAX_RESOLUTION)
    if sf > 1 and PRE_DOWNSCALE:
        nuc = pipex_downscale(nuc, sf).astype("float32")
        mem = pipex_downscale(mem, sf).astype("float32")
        dapi_raw = pipex_downscale(dapi_raw, sf)
        print(f">>> Pre-downscaled by {sf} -> {nuc.shape}")
    elif sf > 1:
        print(f">>> Large image (factor {sf}); PIPEX will downscale mask.")
    else:
        print(f">>> Image {H}x{W} within limit. No downscaling.")

    print("\n" + "=" * 70); print("STEP 3: Running CellPose WHOLE-CELL"); print("=" * 70)
    print(">>> augment=True -> slower on CPU; please wait...")
    mask = run_cellpose_2ch(nuc, mem)
    print(f">>> Cells segmented: {int(mask.max())}")

    print("\n" + "=" * 70); print("STEP 4: Saving mask outputs"); print("=" * 70)
    np.save(OUTPUT_NPY, mask.astype(np.int32)); print("Saved:", OUTPUT_NPY)
    tifffile.imwrite(OUTPUT_TIF, mask.astype(np.int32))
    print("Saved (int32 tif):", OUTPUT_TIF)

    print("\n" + "=" * 70)
    print("STEP 4b: Saving DAPI + segmentation boundary overlay (QC, RED)")
    print("=" * 70)
    save_boundary_overlay(dapi_raw, mask, OUTPUT_OVERLAY,
                          boundary_color=OVERLAY_COLOR,
                          line_thickness=OVERLAY_LINE_THICK)
    if OVERLAY_ON_MEMBRANE:
        save_boundary_overlay((mem * 65535).astype(np.uint16), mask,
                              OUTPUT_OVERLAY_MEM,
                              boundary_color=OVERLAY_COLOR,
                              line_thickness=OVERLAY_LINE_THICK)

    print("\n" + "=" * 70); print("STEP 5: Measure markers string"); print("=" * 70)
    markers, skipped = generate_measure_markers(DATA_DIR, (H, W))
    if skipped:
        print("Skipped:", [s[0] for s in skipped])
    ms = ",".join(markers)
    print(f">>> {len(markers)} markers detected.\n"); print(ms, "\n")
    with open(os.path.join(DATA_DIR, "measure_markers.txt"), "w", encoding="utf-8") as f:
        f.write(ms)
    print("Saved: measure_markers.txt")

    print("\n" + "=" * 70)
    print("DONE. PIPEX Segmentation GUI:")
    print("  Custom segmentation      =", os.path.basename(OUTPUT_TIF), "(int32 tif)")
    print("  Custom segmentation type =", COMPARTMENT_LABEL)
    print("  Nuclei expansion         = 1")
    print("  Measure markers          = (paste from measure_markers.txt)")
    print("  QC overlay (DAPI)        =", os.path.basename(OUTPUT_OVERLAY))
    if OVERLAY_ON_MEMBRANE:
        print("  QC overlay (membrane)    =", os.path.basename(OUTPUT_OVERLAY_MEM))
    print("=" * 70)


if __name__ == "__main__":
    main()