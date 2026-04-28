import sys
import os
import argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_USE_LEGACY_KERAS'] = '1'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
import datetime

import numpy as np
import pandas as pd
import tensorflow as tf
from pipex_utils import iter_marker_images, log, sanitize_marker_list, validate_marker_files

from stardist.models import StarDist2D
from stardist.plot import render_label

import cv2

from skimage.io import imsave, imread
from skimage.filters import threshold_multiotsu, threshold_triangle
from sklearn.mixture import GaussianMixture
from skimage.measure import regionprops
from skimage.segmentation import watershed, mark_boundaries, expand_labels, relabel_sequential
from skimage.feature import peak_local_max
from skimage.transform import resize


pipex_max_resolution = 32768
if "PIPEX_MAX_RESOLUTION" in os.environ:
    pipex_max_resolution = int(os.environ.get('PIPEX_MAX_RESOLUTION'))
pipex_scale_factor = 0
data_folder = os.environ.get('PIPEX_DATA')
stardist_tile_threshold = 4096
watershed_tile_threshold = 2048
watershed_tile_size = 2048

nuclei_marker = ""
nuclei_diameter = 0
nuclei_expansion = 0
nuclei_definition = 0
nuclei_closeness = 0
nuclei_area_limit = 0
membrane_marker = ""
membrane_diameter = 0
membrane_compactness = 0.9
membrane_keep = "no"
adjust_images = 0
custom_segmentation = ""
custom_segmentation_type = "full"
measure_markers = ""
gmm_min_separation = 0.5


def downscale_images(np_img):
    if len(np_img) > pipex_max_resolution or len(np_img[0]) > pipex_max_resolution:
        global pipex_scale_factor
        if pipex_scale_factor == 0:
            i = 2
            while pipex_scale_factor == 0:
                if max(len(np_img), len(np_img[0])) / i <= pipex_max_resolution:
                    pipex_scale_factor = i
                else:
                    i = i * 2
            global nuclei_diameter
            nuclei_diameter = nuclei_diameter / pipex_scale_factor
            global nuclei_expansion
            nuclei_expansion = nuclei_expansion / pipex_scale_factor
            global membrane_diameter
            membrane_diameter = membrane_diameter / pipex_scale_factor
        return resize(np_img, (len(np_img) // pipex_scale_factor, len(np_img[0]) // pipex_scale_factor), order=0, preserve_range=True, anti_aliasing=False).astype('uint16')
    return np_img


def upscale_results(df):
    if pipex_scale_factor > 0:
        for fname in ["segmentation_mask.tif", "segmentation_binary_mask.tif"]:
            path = os.path.join(data_folder, "analysis", fname)
            image = cv2.imread(path, cv2.IMREAD_UNCHANGED)
            image = cv2.resize(image, (image.shape[1] * pipex_scale_factor, image.shape[0] * pipex_scale_factor), interpolation=cv2.INTER_NEAREST)
            cv2.imwrite(path, image)

        labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
        labels = labels.repeat(pipex_scale_factor, axis=0).repeat(pipex_scale_factor, axis=1)
        np.save(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'), labels)

        path = os.path.join(data_folder, "analysis", "segmentation_mask_show.jpg")
        image = cv2.imread(path, cv2.IMREAD_UNCHANGED)
        image = cv2.resize(image, (image.shape[1] * pipex_scale_factor, image.shape[0] * pipex_scale_factor), interpolation=cv2.INTER_NEAREST)
        cv2.imwrite(path, image)

        pipex_scale_factor_n2 = pow(pipex_scale_factor, 2)
        df['x'] = df['x'] * pipex_scale_factor
        df['y'] = df['y'] * pipex_scale_factor
        df['size'] = df['size'].apply(lambda x: int(x * pipex_scale_factor_n2))


def cell_segmentation(nuclei_img_orig, membrane_img_orig, custom_img_orig):
    if custom_segmentation == "" or custom_segmentation_type != "nuc":
        img_min = np.amin(nuclei_img_orig)
        img_range = np.amax(nuclei_img_orig) - img_min
        if img_range == 0:
            print(">>> WARNING: nuclei image has uniform intensity, cannot segment", flush=True)
            return np.zeros_like(nuclei_img_orig, dtype=np.int32), set()
        nuclei_img = (nuclei_img_orig - img_min) / img_range

        try:
            gpus = tf.config.list_physical_devices('GPU')
            if gpus:
                try:
                    for gpu in gpus:
                        tf.config.experimental.set_memory_growth(gpu, True)
                    print(f">>> GPU(s) detected: {len(gpus)} device(s) available", flush=True)
                except RuntimeError as e:
                    print(f">>> GPU configuration error: {e}", flush=True)
            else:
                print(">>> No GPU detected, using CPU", flush=True)
        except ImportError:
            print(">>> TensorFlow not configured for GPU detection", flush=True)

        #run stardist over nuclei image
        model = StarDist2D.from_pretrained('2D_versatile_fluo')
        sd_labels = None
        #for big images (>stardist_tile_threshold), run predict_instances_big method using 2048 square tiles
        if max(len(nuclei_img), len(nuclei_img[0])) > stardist_tile_threshold:
            sd_labels, _ = model.predict_instances_big(nuclei_img,axes='YX',block_size=1024,min_overlap=128,prob_thresh=(nuclei_definition if nuclei_definition > 0 else None), nms_thresh=(nuclei_closeness if nuclei_closeness > 0 else None))
        else:
            sd_labels, _ = model.predict_instances(nuclei_img,axes='YX',prob_thresh=(nuclei_definition if nuclei_definition > 0 else None), nms_thresh=(nuclei_closeness if nuclei_closeness > 0 else None))
        log("Stardist prediction done")

        if nuclei_area_limit > 0:
            detections = regionprops(sd_labels)
            for curr_detection in detections:
                if curr_detection.area > nuclei_area_limit:
                    sd_labels[sd_labels == curr_detection.label] = 0

        imsave(os.path.join(data_folder, "analysis", "quality_control", "stardist_result.jpg"), (render_label(sd_labels, img=None)[..., :3] * 255).astype(np.uint8))
        log("Stardist base result image saved")

        #if nuclei_expansion parameter is required, expand labelled regions (avoiding overlap) to specified size
        if nuclei_expansion > 0:
            sd_labels_expanded = expand_labels(sd_labels, distance=nuclei_expansion)
            imsave(os.path.join(data_folder, "analysis", "quality_control", "stardist_result_expanded.jpg"), np.uint8(mark_boundaries(nuclei_img_orig, sd_labels_expanded) * 255))
            log("Stardist expanded result image saved")
        else:
            sd_labels_expanded = sd_labels.copy()
    else:
        sd_labels_expanded = custom_img_orig

    affected_by_membrane = set()
    if membrane_diameter > 0 or custom_segmentation_type == "mem":
        if custom_segmentation == "" or custom_segmentation_type != "mem":
            #if membrane marker is provided, run custom watershed segmentation
            mem_min = np.amin(membrane_img_orig)
            mem_range = np.amax(membrane_img_orig) - mem_min
            if mem_range == 0:
                print(">>> WARNING: membrane image has uniform intensity, skipping membrane segmentation", flush=True)
                membrane_img = np.zeros_like(membrane_img_orig, dtype=np.float64)
            else:
                membrane_img = (membrane_img_orig - mem_min) / mem_range

            membrane_keep_index = -1
            membrane_intensity_mean = threshold_multiotsu(membrane_img, 5)[0]
            tiles = []
            if len(membrane_img) > watershed_tile_threshold or len(membrane_img[0]) > watershed_tile_threshold:
                num_rows = int(len(membrane_img) / watershed_tile_size)
                if (len(membrane_img) % watershed_tile_size != 0):
                    num_rows = num_rows + 1
                num_columns = int(len(membrane_img[0]) / watershed_tile_size)
                if (len(membrane_img[0]) % watershed_tile_size != 0):
                    num_columns = num_columns + 1
                for row in range(num_rows):
                    tiles.append([])
                    for column in range(num_columns):
                        tiles[row].append(membrane_img[(row * watershed_tile_size):((row + 1) * watershed_tile_size), (column * watershed_tile_size):((column + 1) * watershed_tile_size)])
            else:
                tiles.append([membrane_img])

            for tile_row in range(len(tiles)):
                for tile_column in range(len(tiles[tile_row])):
                    tile = tiles[tile_row][tile_column]
                    tile_x = tile_row * watershed_tile_size
                    tile_y = tile_column * watershed_tile_size
                    tile_desc = str(tile_row) + "_" + str(tile_column)
                    tile_orig = membrane_img[(tile_row * watershed_tile_size):((tile_row + 1) * watershed_tile_size), (tile_column * watershed_tile_size):((tile_column + 1) * watershed_tile_size)]
                    # LoG preprocessing: Gaussian smooth then Laplacian for edge-enhanced height map
                    sigma = max(1.0, float(membrane_diameter) / 6.0)
                    tile_f32 = tile.astype(np.float32)
                    tile_smooth = cv2.GaussianBlur(tile_f32, (0, 0), sigma)
                    tile_lap = cv2.Laplacian(tile_smooth, cv2.CV_32F)
                    tile_lap_neg = np.clip(-tile_lap, 0.0, None).astype(np.float64)
                    lap_max = tile_lap_neg.max()
                    if lap_max > 0:
                        tile_lap_neg /= lap_max
                    height_map = tile_smooth.astype(np.float64) * 0.5 + tile_lap_neg * 0.5

                    # Pass 1: nucleus-seeded watershed — nuclei flood outward, bright membrane acts as barrier
                    tile_nuc_slice = sd_labels[tile_x:tile_x + len(tile), tile_y:tile_y + len(tile[0])]
                    watershed_mask = (tile >= membrane_intensity_mean) | (tile_nuc_slice > 0)
                    ws_labels = watershed(height_map, markers=tile_nuc_slice, mask=watershed_mask, compactness=membrane_compactness)
                    log("Watershed of tile " + tile_desc + " done")

                    imsave(os.path.join(data_folder, "analysis", "quality_control", "wathershed_tile_" + tile_desc + "_result.jpg"), np.uint8(mark_boundaries(tile_orig, ws_labels) * 255))
                    log("Watershed of tile " + tile_desc + " result image saved")

                    # Apply membrane-guided cell boundaries to sd_labels_expanded
                    tile_exp = sd_labels_expanded[tile_x:tile_x + len(ws_labels), tile_y:tile_y + len(ws_labels[0])]
                    affected_by_membrane.update(np.unique(ws_labels[ws_labels > 0]).tolist())
                    tile_exp[ws_labels > 0] = ws_labels[ws_labels > 0]
                    tile_exp[(ws_labels == 0) & (tile_nuc_slice == 0)] = 0

                    # Pass 2: detect membrane-only cells in unclaimed high-signal areas
                    if membrane_keep == 'yes':
                        unclaimed_mask = (ws_labels == 0) & (tile >= membrane_intensity_mean)
                        if unclaimed_mask.any():
                            max_coords = peak_local_max(np.where(unclaimed_mask, tile, 0.0), min_distance=max(1, int(membrane_diameter)), labels=unclaimed_mask.astype(np.int32))
                            if len(max_coords) > 0:
                                seed_labels = np.zeros_like(ws_labels, dtype=np.int32)
                                for i, (r, c) in enumerate(max_coords):
                                    seed_labels[r, c] = membrane_keep_index - i
                                membrane_keep_index -= len(max_coords)
                                ws_keep = watershed(height_map, markers=seed_labels, mask=unclaimed_mask, compactness=membrane_compactness)
                                tile_exp[ws_keep != 0] = ws_keep[ws_keep != 0]
                                log("Membrane-only watershed of tile " + tile_desc + " done")

                    imsave(os.path.join(data_folder, "analysis", "quality_control", "wathershed_tile_" + tile_desc + "_result_merged_by_nuclei.jpg"), np.uint8(mark_boundaries(tile_orig, tile_exp) * 255))
                    log("Watershed of tile " + tile_desc + " final result image saved")

        else:
            ws_labels = custom_img_orig
            # merge resulting segments so they don't cut nuclei (not expanded)
            ws_regions = {}
            unique_mem_labels = np.unique(ws_labels)
            unique_mem_labels = unique_mem_labels[unique_mem_labels != 0]
            for mem_label in unique_mem_labels:
                mask = ws_labels == mem_label
                nuc_vals = sd_labels[mask]
                ws_regions[mem_label] = set(nuc_vals[nuc_vals != 0].tolist())

            # merge resulting segments that contain same nuclei and/or nothing
            value_to_first = {}
            for region in ws_regions:
                key = frozenset(ws_regions[region])
                if key not in value_to_first:
                    value_to_first[key] = region
            ws_regions_merged = {region: value_to_first[frozenset(ws_regions[region])] for region in ws_regions}

            if ws_regions_merged:
                lookup = np.arange(max(ws_regions_merged) + 1, dtype=ws_labels.dtype)
                for label, merged in ws_regions_merged.items():
                    lookup[label] = merged
                nz = ws_labels != 0
                ws_labels[nz] = lookup[ws_labels[nz]]

            imsave(os.path.join(data_folder, "analysis", "quality_control", "wathershed_result_merged_by_nuclei.jpg"),
                   np.uint8(mark_boundaries(nuclei_img_orig, ws_labels) * 255))
            print(">>> Watershed preliminary nuclei filter result image saved =",
                  datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

            # cut expanded nuclei that collide with watershed segments
            in_expansion = (sd_labels != sd_labels_expanded) | (sd_labels_expanded == 0)
            candidates = in_expansion & (ws_labels != 0)
            for mem_label in np.unique(ws_labels[candidates]):
                pixel_mask = candidates & (ws_labels == mem_label)
                region = ws_regions[mem_label]
                exp_vals = sd_labels_expanded[pixel_mask]
                if len(region) == 0:
                    affected_by_membrane.update(exp_vals[exp_vals > 0].tolist())
                    sd_labels_expanded[pixel_mask] = 0
                else:
                    not_in_region = np.array([e == 0 or e not in region for e in exp_vals])
                    if not_in_region.any():
                        if membrane_keep == 'yes':
                            membrane_only_label = sum(n for n in region if n < 0)
                            new_vals = exp_vals.copy()
                            new_vals[not_in_region] = membrane_only_label
                            sd_labels_expanded[pixel_mask] = new_vals
                        else:
                            affected_by_membrane.update(exp_vals[not_in_region & (exp_vals > 0)].tolist())
                            new_vals = exp_vals.copy()
                            new_vals[not_in_region] = 0
                            sd_labels_expanded[pixel_mask] = new_vals

        top_positive_label = np.max(sd_labels_expanded)
        negative_label_mask = sd_labels_expanded < 0
        affected_by_membrane.update(np.unique(sd_labels_expanded[negative_label_mask]).tolist())
        sd_labels_expanded[negative_label_mask] = top_positive_label - sd_labels_expanded[negative_label_mask]


        #find rare disjointed segmented cells and using their associated convex hull instead
        segment_properties = regionprops(sd_labels_expanded)
        for segment in segment_properties:
            contours, hierarchy = cv2.findContours(segment.image.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if len(contours) > 1:
                log("Found disjointed segment " + str(segment.label) + ", using convex hull instead")
                bbox = segment.bbox
                filling_image = segment.image_convex
                if segment.solidity > 0.0:
                    rows, cols = np.where(filling_image > 0)
                    sd_labels_expanded[rows + bbox[0], cols + bbox[1]] = segment.label

    del sd_labels
    if custom_segmentation != "" and custom_segmentation_type == "full":
        sd_labels_expanded = custom_img_orig
    else:
        sd_labels_expanded = relabel_sequential(sd_labels_expanded)[0]

    np.save(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'), sd_labels_expanded)
    log("Final joined segmentation result numpy binary data saved")

    imsave(os.path.join(data_folder, "analysis", "segmentation_mask_show.jpg"), np.uint8(mark_boundaries(nuclei_img_orig, sd_labels_expanded) * 255))
    log("Final joined segmentation result image over nuclei saved")

    black_canvas = np.zeros((len(nuclei_img_orig), len(nuclei_img_orig[0]), 3), dtype=np.uint8)
    boundaries_rgb = np.uint8(mark_boundaries(black_canvas, sd_labels_expanded, color=(0, 1, 0)) * 255)
    boundaries_rgba = np.dstack([boundaries_rgb, np.where(boundaries_rgb.any(axis=2), 255, 0).astype(np.uint8)])
    imsave(os.path.join(data_folder, "analysis", "segmentation_boundaries.png"), boundaries_rgba)
    log("Final segmentation boundaries overlay saved")

    sdLabels_expanded_binary = np.copy(sd_labels_expanded)
    sdLabels_expanded_binary[sdLabels_expanded_binary > 0] = 1
    imsave(os.path.join(data_folder, "analysis", "segmentation_binary_mask.tif"), np.uint8(sdLabels_expanded_binary * 255))
    del sdLabels_expanded_binary

    if np.amax(sd_labels_expanded) <= 255:
        imsave(os.path.join(data_folder, "analysis", "segmentation_mask.tif"), sd_labels_expanded.astype(np.uint8))
    elif np.amax(sd_labels_expanded) <= 65535:
        imsave(os.path.join(data_folder, "analysis", "segmentation_mask.tif"), sd_labels_expanded.astype(np.uint16))
    else:
        imsave(os.path.join(data_folder, "analysis", "segmentation_mask.tif"), sd_labels_expanded.astype(np.uint32))
    log("Final joined segmentation result image saved")

    return sd_labels_expanded, affected_by_membrane


#Function to calculate the marker intensities for each cell
def marker_calculation(marker, marker_img, cellLabels, data_table):
    #applying segmentation mask over the marker image
    marker_img_min, marker_img_max = np.percentile(marker_img, (1, 99.5))
    marker_img_range = marker_img_max - marker_img_min
    if marker_img_range == 0:
        marker_img_range = 1
    marker_img_norm = np.clip((marker_img - marker_img_min) / marker_img_range, 0.0, 1.0)
    try:
        c_otsu = threshold_multiotsu(marker_img_norm, 3)
        cell_binarized_threshold = c_otsu[0]
    except Exception:
        cell_binarized_threshold = float(np.mean(marker_img_norm))
        print(f">>> WARNING: threshold_multiotsu failed for {marker}, using mean {cell_binarized_threshold:.4f} as fallback", flush=True)
    log("Marker " + marker + " binarize threshold " + str(cell_binarized_threshold))

    try:
        tri_threshold = float(threshold_triangle(marker_img_norm))
    except Exception:
        tri_threshold = None

    gmm_model = None
    gmm_pos_idx = 1
    try:
        vals = marker_img_norm.flatten()
        vals = vals[np.isfinite(vals)]
        if len(vals) >= 10 and np.std(vals) > 0:
            gmm = GaussianMixture(n_components=2, random_state=42, max_iter=300)
            gmm.fit(vals.reshape(-1, 1))
            means = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_.flatten())
            order = np.argsort(means)
            sep = (means[order[1]] - means[order[0]]) / (stds[order[0]] + stds[order[1]])
            if sep >= gmm_min_separation:
                gmm_model = gmm
                gmm_pos_idx = int(order[1])
    except Exception:
        pass

    marker_properties = regionprops(cellLabels, marker_img)
    #obtaining mean intensity per cell
    for cell in marker_properties:
        data_table[cell.label][marker] = cell.intensity_mean
        cell_image = cell.image_intensity
        cell_image = cell_image[(cell_image != 0) & (~np.isnan(cell_image))]
        data_table[cell.label][marker + '_local_90'] = np.quantile(cell_image, 0.9) if len(cell_image) > 0 else 0
        cell_image_norm = np.maximum((cell_image - marker_img_min) / marker_img_range, 0.0)
        data_table[cell.label][marker + '_ratio_pixels'] = np.count_nonzero(cell_image_norm >= cell_binarized_threshold) / cell.area
        cell_mean_norm = np.maximum((cell.intensity_mean - marker_img_min) / marker_img_range, 0.0)
        data_table[cell.label][marker + '_otsu3'] = cell_mean_norm - cell_binarized_threshold
        data_table[cell.label][marker + '_triangle_score'] = (cell_mean_norm - tri_threshold) if tri_threshold is not None else np.nan
        data_table[cell.label][marker + '_gmm_prob'] = (
            float(gmm_model.predict_proba([[cell_mean_norm]])[0][gmm_pos_idx])
            if gmm_model is not None else np.nan)

    log("Marker " + marker + " calculated")


#Function to handle the command line parameters passed
def options(argv):
    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='segmentation.py')
    parser.add_argument('--data', default=os.environ.get('PIPEX_DATA'),
        help='path to images folder : example -> -data=/lab/projectX/images')
    parser.add_argument('--nuclei_marker', default='',
        help='name before . in image file : example -> -nuclei_marker=DAPI1')
    parser.add_argument('--nuclei_diameter', type=int, default=0,
        help='number of pixels : example -> -nuclei_diameter=20')
    parser.add_argument('--nuclei_expansion', type=int, default=0,
        help='number of pixels, can be 0 : example -> -nuclei_expansion=20')
    parser.add_argument('--nuclei_definition', type=float, default=0,
        help='gradation between 0.001 and 0.999 : example -> -nuclei_definition=0.1')
    parser.add_argument('--nuclei_closeness', type=float, default=0,
        help='gradation between 0.001 and 0.999 : example -> -nuclei_closeness=0.6')
    parser.add_argument('--nuclei_area_limit', type=float, default=0,
        help='number of pixels : example -> -nuclei_area_limit=3200')
    parser.add_argument('--membrane_marker', default='',
        help='name before . in image file : example -> -membrane_marker=CDH1')
    parser.add_argument('--membrane_diameter', type=int, default=0,
        help='number of pixels : example -> -membrane_diameter=25')
    parser.add_argument('--membrane_compactness', type=float, default=0.9,
        help='"squareness" of the membrane, gradation between 0.001 and 0.999 : example -> -membrane_compactness=0.5')
    parser.add_argument('--membrane_keep', choices=['yes', 'no'], default='no',
        help='keep segmented membranes without nuclei : example -> -membrane_keep=no')
    parser.add_argument('--custom_segmentation', default='',
        help='file path to a pre-made custom segmentation : example -> -custom_segmentation=/data/custom_seg.npy')
    parser.add_argument('--custom_segmentation_type', choices=['full', 'nuc', 'mem'], default='full',
        help='type of the custom segmentation : example -> -custom_segmentation_type=full')
    parser.add_argument('--measure_markers', type=lambda s: [x.strip() for x in s.split(',')], default=[],
        help='list of marker names to measure : example -> -measure_markers=AMY2A,SST,GORASP2')
    parser.add_argument('--gmm_min_separation', type=float, default=0.5,
        help='minimum separation between GMM components (in combined std units) to trust the fit and compute gmm_prob : example -> -gmm_min_separation=0.5')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)

if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    nuclei_marker = args.nuclei_marker
    nuclei_diameter = args.nuclei_diameter
    nuclei_expansion = args.nuclei_expansion
    nuclei_definition = args.nuclei_definition
    nuclei_closeness = args.nuclei_closeness
    nuclei_area_limit = args.nuclei_area_limit
    membrane_marker = args.membrane_marker
    membrane_diameter = args.membrane_diameter
    membrane_compactness = args.membrane_compactness
    membrane_keep = args.membrane_keep
    custom_segmentation = args.custom_segmentation
    custom_segmentation_type = args.custom_segmentation_type
    measure_markers = sanitize_marker_list(args.measure_markers)
    gmm_min_separation = args.gmm_min_separation

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
    with open(os.path.join(data_folder, 'log_settings_segmentation.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time segmentation = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))

    log("Start time segmentation")

    if measure_markers:
        validate_marker_files(data_folder, measure_markers)

    os.makedirs(os.path.join(data_folder, 'analysis', 'quality_control'), exist_ok=True)

    #finding nuclei and membrane image filenames by marker name
    nuclei_img = None
    membrane_img = None
    custom_img = None
    if custom_segmentation == "" or custom_segmentation_type != "full":
        loaded = {}
        def collect_images(marker, raw_img):
            if marker == nuclei_marker:
                loaded['nuclei'] = downscale_images(raw_img)
            elif marker == membrane_marker:
                loaded['membrane'] = downscale_images(raw_img)

        markers_to_load = [m for m in [nuclei_marker, membrane_marker] if m]
        iter_marker_images(data_folder, markers_to_load, collect_images)
        nuclei_img = loaded.get('nuclei')
        membrane_img = loaded.get('membrane')
    else:
        try_image = False
        try:
            custom_img = downscale_images(np.load(os.path.join(data_folder, custom_segmentation)))
        except Exception as e:
            try_image = True
            print('>>> ', e, flush=True)
        if try_image:
            try:
                custom_img = downscale_images(cv2.imread(os.path.join(data_folder, custom_segmentation), cv2.IMREAD_UNCHANGED))
            except Exception as e:
                print('>>> Could not custom segmentation file ' + custom_segmentation, flush=True)
                print('>>> ', e, flush=True)

    #performing segmentation
    if custom_segmentation == "" or custom_segmentation_type != "full":
        if nuclei_img is None:
            print(">>> ERROR: nuclei marker '" + nuclei_marker + "' not found in data folder, cannot segment", flush=True)
            sys.exit(1)
        if membrane_marker != "" and membrane_img is None:
            print(">>> WARNING: membrane marker '" + membrane_marker + "' not found in data folder, continuing without membrane refinement", flush=True)
    elif custom_img is None:
        print(">>> ERROR: custom segmentation file '" + custom_segmentation + "' could not be loaded, cannot segment", flush=True)
        sys.exit(1)
    cellLabels, membraneAffected = cell_segmentation(nuclei_img, membrane_img, custom_img)

    del nuclei_img
    del membrane_img
    del custom_img
    #creating data table with segmented cell information
    cellProperties = regionprops(cellLabels)
    data_table = {}
    for cell in cellProperties:
        data_cell = {}
        data_cell['cell_id'] = cell.label
        data_cell['size'] = cell.area
        data_cell['x'] = int(cell.centroid[1])
        data_cell['y'] = int(cell.centroid[0])
        data_cell['solidity'] = cell.solidity
        data_cell['eccentricity'] = cell.eccentricity
        data_cell['memref'] = int(cell.label in membraneAffected)
        data_table[cell.label] = data_cell

    #calculating marker intensities per cell
    def handle_marker(marker, raw_img):
        marker_calculation(marker, downscale_images(raw_img), cellLabels, data_table)

    iter_marker_images(data_folder, measure_markers, handle_marker)


    #dumpming data_table in cell_data.csv file
    df = pd.DataFrame.from_dict(data_table, orient='index')
    upscale_results(df)
    binarized_marker_columns = []
    for marker in measure_markers:
        binarized_marker_columns.append(marker + "_local_90")
        binarized_marker_columns.append(marker + "_ratio_pixels")
        binarized_marker_columns.append(marker + "_otsu3")
        binarized_marker_columns.append(marker + "_triangle_score")
        binarized_marker_columns.append(marker + "_gmm_prob")
    measure_markers.extend(binarized_marker_columns)
    measure_markers.insert(0, 'memref')
    measure_markers.insert(0, 'eccentricity')
    measure_markers.insert(0, 'solidity')
    measure_markers.insert(0, 'y')
    measure_markers.insert(0, 'x')
    measure_markers.insert(0, 'size')
    measure_markers.insert(0, 'cell_id')
    df = df.reindex(measure_markers, axis=1).fillna(0)
    df.to_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'), index=False)

    log("End time segmentation")
