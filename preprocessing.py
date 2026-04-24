import sys
import os
import argparse
import datetime
import numpy as np
from pipex_utils import iter_marker_images, log, sanitize_marker_list, validate_marker_files

from skimage.filters import threshold_multiotsu
from skimage.io import imsave
from skimage.transform import resize

import cv2


pipex_max_resolution = 30000
if "PIPEX_MAX_RESOLUTION" in os.environ:
    pipex_max_resolution = int(os.environ.get('PIPEX_MAX_RESOLUTION'))
pipex_scale_factor = 0
data_folder = './data'
preprocess_markers = []

thres_min = 0.0
thres_max = 1.0
otsu_threshold_levels = -1
bin_min = 0
bin_max = 2

exposure = 1.0

tile_size = 0
light_gradient = 0
flatten_spots = 'no'
balance_tiles = 'no'
stitch_size = 0
tophat_radius = 0


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
        return resize(np_img, (len(np_img) // pipex_scale_factor, len(np_img[0]) // pipex_scale_factor), order=0, preserve_range=True, anti_aliasing=False).astype('uint16')
    return np_img


def upscale_results(marker):
    if pipex_scale_factor > 0:
        path = os.path.join(data_folder, "preprocessed", marker + ".tif")
        image = cv2.imread(path, cv2.IMREAD_UNCHANGED)
        image = cv2.resize(image, (image.shape[1] * pipex_scale_factor, image.shape[0] * pipex_scale_factor), interpolation=cv2.INTER_NEAREST)
        cv2.imwrite(path, image)


def rescale_tile_intensity(tile, mean_in, mean_factor, dev_in, dev_factor, bins):
    mask = tile != 0
    result = np.clip(tile + mean_in * (mean_factor - 1), 0.0, 1.0)
    result_dev = (tile - mean_in) * (dev_factor - 1)
    result = np.where(tile >= mean_in,
                      np.clip(result + result_dev, 0.0, 1.0),
                      np.clip(-result_dev, 0.0, 1.0))
    above_max = result > bins[bin_max]
    result = np.where(above_max,
                      bins[bin_max] + np.power(np.maximum(result - bins[bin_max], 0) * 100, 0.6) / 100,
                      result)
    below_min = result < bins[bin_min]
    result = np.where(below_min,
                      bins[bin_min] - np.power(np.maximum(bins[bin_min] - result, 0) * 100, 0.6) / 100,
                      result)
    return np.where(mask, result, tile)


def apply_tile_compensation(f_name, np_img, bins, tile_data, local_tile_size, local_stitch_size):
    num_rows = int(len(np_img) / local_tile_size)
    num_columns = int(len(np_img[0]) / local_tile_size)

    for row in range(num_rows):
        for column in range(num_columns):
            if tile_data['tiles'][row][column]['samples'] > 0:
                tile = np_img[(row * local_tile_size):((row + 1) * local_tile_size), (column * local_tile_size):((column + 1) * local_tile_size)]

                mean_factor = tile_data['chosen']['mean'] / tile_data['tiles'][row][column]['mean']
                dev_factor = tile_data['chosen']['deviation'] / max(0.001, tile_data['tiles'][row][column]['deviation'])

                values = tile[np.nonzero(tile)]
                curr_mean = np.mean(values)
                curr_dev = np.std(values)

                np_img[(row * local_tile_size):((row + 1) * local_tile_size), (column * local_tile_size):((column + 1) * local_tile_size)] = rescale_tile_intensity(tile, curr_mean, mean_factor, curr_dev, dev_factor, bins)

    if tile_size == local_tile_size:
        log("Balanced tiles for image " + f_name)

    if local_stitch_size > 0:
        for row in range(num_rows):
            for column in range(num_columns):
                if row < num_rows - 1:
                    for gradient in range(3):
                        streg_min_y = int((row + 1) * local_tile_size - local_stitch_size / pow(2, (gradient + 1)))
                        streg_max_y = int((row + 1) * local_tile_size + local_stitch_size / pow(2, (gradient + 1)))
                        streg_min_x = int(column * local_tile_size)
                        streg_max_x = int((column + 1) * local_tile_size)
                        stitch = np_img[streg_min_y:streg_max_y,streg_min_x:streg_max_x]
                        kernel_size = 1 + gradient * 2
                        np_img[streg_min_y:streg_max_y, streg_min_x:streg_max_x] = cv2.GaussianBlur(stitch, (kernel_size, kernel_size), 0)
                if (column < num_columns - 1):
                    for gradient in range(3):
                        streg_min_y = int(row * local_tile_size - (local_stitch_size / pow(2, (gradient + 1)) if row > 0 else 0))
                        streg_max_y = int((row + 1) * local_tile_size + (local_stitch_size / pow(2, (gradient + 1)) if row < num_rows - 1 else 0))
                        streg_min_x = int((column + 1) * local_tile_size - local_stitch_size / pow(2, (gradient + 1)))
                        streg_max_x = int((column + 1) * local_tile_size + local_stitch_size / pow(2, (gradient + 1)))
                        stitch = np_img[streg_min_y:streg_max_y, streg_min_x:streg_max_x]
                        kernel_size = 1 + gradient * 2
                        np_img[streg_min_y:streg_max_y, streg_min_x:streg_max_x] = cv2.GaussianBlur(stitch, (kernel_size, kernel_size), 0)

        if tile_size == local_tile_size:
            log("Smoothed stitched lines for image " + f_name)


def generate_tile_compensation_data(np_img, bins, local_tile_size):
    tile_data = {}
    tile_data['tiles'] = []

    tiles = []
    num_rows = int(len(np_img) / local_tile_size)
    num_columns = int(len(np_img[0]) / local_tile_size)
    for row in range(num_rows):
        tiles.append([])
        tile_data['tiles'].append([])
        for column in range(num_columns):
            tiles[row].append(np_img[(row * local_tile_size):((row + 1) * local_tile_size), (column * local_tile_size):((column + 1) * local_tile_size)])
            tile_data['tiles'][row].append([])

    for tile_row in range(len(tiles)):
        for tile_column in range(len(tiles[tile_row])):
            tile = tiles[tile_row][tile_column]

            bin_count = np.histogram(np.ravel(tile), bins)[0]

            values = tile[np.nonzero(tile)]
            values = values[values >= bins[bin_min]]
            values = values[values <= bins[bin_max]]

            curr_subtile_data = {}
            curr_subtile_data['row'] = tile_row
            curr_subtile_data['column'] = tile_column
            curr_subtile_data['samples'] = bin_count[bin_max]
            if (curr_subtile_data['samples'] == 0):
                curr_subtile_data['mean'] = 0
                curr_subtile_data['deviation'] = 0
            else:
                curr_subtile_data['mean'] = np.mean(values)
                curr_subtile_data['deviation'] = np.std(values)

            tile_data['tiles'][tile_row][tile_column] = curr_subtile_data

    valid_tiles = [curr_data for row in tile_data['tiles'] for curr_data in row if curr_data['samples'] > 0]
    if valid_tiles:
        median_mean = np.median([t['mean'] for t in valid_tiles])
        tile_data['chosen'] = min(valid_tiles, key=lambda t: abs(t['mean'] - median_mean))

    if tile_size == local_tile_size:
        log("Calculated reference tiles balance")

    return tile_data

def apply_tile_gradient_compensation(f_name, np_img, bins, gradient_data):
    num_rows = int(len(np_img) / tile_size)
    num_columns = int(len(np_img[0]) / tile_size)
    kernel_size = int(tile_size / light_gradient)

    for row in range(num_rows):
        for column in range(num_columns):
            tile_ratio = gradient_data[row][column]['ratio_gradient']
            if tile_ratio > 0:
                for row_kernel in range(light_gradient):
                    for column_kernel in range(light_gradient):
                        kernel_ratio = gradient_data[row][column]['kernels'][row_kernel][column_kernel]['ratio_gradient']
                        if kernel_ratio > 0 and kernel_ratio < tile_ratio:
                            tile_kernel_sy = row * tile_size + row_kernel * kernel_size
                            tile_kernel_sx = column * tile_size + column_kernel * kernel_size

                            final_ratio = 1 + (bins[bin_max] * (tile_ratio / kernel_ratio)) / bins[bin_max]

                            np_img[tile_kernel_sy:tile_kernel_sy + kernel_size, tile_kernel_sx:tile_kernel_sx + kernel_size] = np.clip(
                                np_img[tile_kernel_sy:tile_kernel_sy + kernel_size, tile_kernel_sx:tile_kernel_sx + kernel_size] * final_ratio, 0, 1.0)

                tile = np_img[(row * tile_size):((row + 1) * tile_size), (column * tile_size):((column + 1) * tile_size)]
                subtile_data = generate_tile_compensation_data(tile, bins, kernel_size)
                apply_tile_compensation('', tile, bins, subtile_data, kernel_size, kernel_size / 10)

    imsave(os.path.join(data_folder, "preprocessed", os.path.splitext(f_name)[0] + "_gradient.jpg"), np.uint8(np_img * 255))

    log("Applied gradient fix for image " + f_name)


def generate_tile_gradient_data(np_img, bins, tile_size):
    num_rows = int(len(np_img) / tile_size)
    num_columns = int(len(np_img[0]) / tile_size)
    kernel_size = int(tile_size / light_gradient)
    gradient_data = []
    for row in range(num_rows):
        gradient_data.append([])
        for column in range(num_columns):
            gradient_data[row].append([])
            tile_gradient = {}
            tile_gradient_data = []
            max_gradient = 0
            ratio_gradient = 0
            max_gradient_kernel = [0, 0]
            for row_kernel in range(light_gradient):
                tile_gradient_data.append([])
                for column_kernel in range(light_gradient):
                    tile_kernel_sy = row * tile_size + row_kernel * kernel_size
                    tile_kernel_sx = column * tile_size + column_kernel * kernel_size
                    tile_kernel_ey = kernel_size + ((tile_size % light_gradient) if row_kernel == light_gradient - 1 else 0)
                    tile_kernel_ex = kernel_size + ((tile_size % light_gradient) if column_kernel == light_gradient - 1 else 0)
                    tile_kernel = np_img[tile_kernel_sy:(tile_kernel_sy + tile_kernel_ey), tile_kernel_sx:(tile_kernel_sx + tile_kernel_ex)]
                    bin_count = np.histogram(np.ravel(tile_kernel), bins)[0]
                    samples = bin_count[bin_min] + bin_count[bin_max]
                    top_samples = bin_count[bin_max]
                    middle_samples = bin_count[bin_min]
                    tile_gradient_info = {}
                    tile_gradient_info['bins'] = bin_count
                    if top_samples > 0:
                        tile_gradient_info['ratio_gradient'] = max(0.01, top_samples / middle_samples) if middle_samples > 0 else 1
                    elif middle_samples > 0:
                        tile_gradient_info['ratio_gradient'] = 0.01
                    else:
                        samples = 1
                        tile_gradient_info['ratio_gradient'] = -1
                    tile_gradient_data[row_kernel].append(tile_gradient_info)
                    if max_gradient == 0 or tile_gradient_info['ratio_gradient'] > 0 and max_gradient < tile_gradient_info['ratio_gradient'] * np.log(samples):
                        max_gradient = tile_gradient_info['ratio_gradient'] * np.log(samples)
                        ratio_gradient = tile_gradient_info['ratio_gradient']
                        max_gradient_kernel = [row_kernel, column_kernel]

            tile_gradient['kernels'] = tile_gradient_data
            tile_gradient['max_gradient'] = max_gradient
            tile_gradient['ratio_gradient'] = ratio_gradient
            tile_gradient['max_gradient_kernel'] = max_gradient_kernel
            gradient_data[row][column] = tile_gradient

    log("Calculated reference tiles gradient")

    return gradient_data


def apply_thresholds(marker, np_img, threshold_min, threshold_max):
    if threshold_min > 0:
        np_thres = np_img.copy()
        np_thres[np_thres >= threshold_min] = 0
        imsave(os.path.join(data_folder, "preprocessed", marker + "_threshold_bottom.jpg"), np.uint8(np_thres * 255))
        del np_thres
        np_img[np_img < threshold_min] = 0
    if threshold_max < 1:
        np_thres = np_img.copy()
        np_thres[np_img <= threshold_max] = 0
        imsave(os.path.join(data_folder, "preprocessed", marker + "_threshold_top.jpg"), np.uint8(np_thres * 255))
        del np_thres
        np_img[np_img > threshold_max] = 0
    log("Applied min and max thresholds")


def preprocess_image(marker, marker_img):
    log("Preprocessing " + marker + " start")
    p_low, p_high = np.percentile(marker_img, (1, 99.5))
    if p_high == p_low:
        print(">>> WARNING: marker", marker, "has uniform intensity, skipping", flush=True)
        return
    np_img = np.clip((marker_img - p_low) / (p_high - p_low), 0.0, 1.0)

    apply_thresholds(marker, np_img, thres_min, thres_max)

    if tophat_radius > 0:
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (tophat_radius * 2 + 1, tophat_radius * 2 + 1))
        np_img = np.clip(cv2.morphologyEx(np_img.astype(np.float32), cv2.MORPH_TOPHAT, kernel).astype(np.float64), 0.0, 1.0)
        log("Applied top-hat background subtraction")

    global otsu_threshold_levels
    if otsu_threshold_levels >=0:
        np_img_denoised = cv2.GaussianBlur(np_img.astype(np.float32), (5, 5), 0).astype(np.float64)
        final_c_otsu = []
        if otsu_threshold_levels == 0:
            c_otsu = threshold_multiotsu(np_img_denoised, 3)
            final_c_otsu = c_otsu.copy()
            c_otsu = np.insert(c_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=c_otsu)
            regions = np.reshape(np.array([c_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_3.jpg"), np.uint8(regions * 255))
    
            c_otsu = threshold_multiotsu(np_img_denoised, 4)
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [2])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_4_1-2.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_4_1-3.jpg"),
                   np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_4_2-3.jpg"), np.uint8(regions * 255))
    
            c_otsu = threshold_multiotsu(np_img_denoised, 5)
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [2, 3])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_1-2.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1, 3])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_1-3.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1, 2])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_1-4.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0, 3])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_2-3.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0, 2])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_2-4.jpg"), np.uint8(regions * 255))
    
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0, 1])
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu_5_3-4.jpg"), np.uint8(regions * 255))
        else:
            final_c_otsu = threshold_multiotsu(np_img_denoised, otsu_threshold_levels)
            temp_otsu = final_c_otsu.copy()
            delete_intervals = []
            if bin_min > 1:
                delete_intervals = list(range(0, bin_min - 1))
            if bin_max < otsu_threshold_levels - 1:
                delete_intervals.extend(range(bin_max, otsu_threshold_levels - 1))
            temp_otsu = np.delete(temp_otsu, delete_intervals)
            temp_otsu = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=temp_otsu)
            regions = np.reshape(np.array([temp_otsu[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0])))
            imsave(os.path.join(data_folder, "preprocessed", marker + "_otsu.jpg"), np.uint8(regions * 255))
    
        np_img[np_img < final_c_otsu[bin_min]] = 0
        if flatten_spots == "yes":
            np_img[np_img > final_c_otsu[bin_max]] = final_c_otsu[bin_max]
        log("Otsu thresholded result image calculated")

        if balance_tiles == "yes":
            final_c_otsu = np.insert(final_c_otsu, 0, 0.0)
            final_c_otsu = np.append(final_c_otsu, 1.0)
            if light_gradient > 1:
                gradient_data = generate_tile_gradient_data(np_img, final_c_otsu, tile_size)
                apply_tile_gradient_compensation(marker, np_img, final_c_otsu, gradient_data)
            np.nan_to_num(np_img, copy=False)
            tile_data = generate_tile_compensation_data(np_img, final_c_otsu, tile_size)
            apply_tile_compensation(marker, np_img, final_c_otsu, tile_data, tile_size, stitch_size)

    if exposure != 1.0:
        np_img = np.clip(np_img * exposure, 0.0, 1.0)
        log("Exposure calculated")

    imsave(os.path.join(data_folder, "preprocessed", marker + ".tif"), np.uint16(np_img * 65535))
    log("Preprocessed result image " + marker + ".tif saved")


#Function to handle the command line parameters passed
def options(argv):
    def _parse_otsu(s):
        if ':' in s:
            parts = s.split(':')
            return int(parts[0]), int(parts[1]), int(parts[2])
        n = int(s)
        return n, 0, n - 1

    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='preprocessing.py')
    parser.add_argument('--data', default='./data',
        help="path to images folder : example -> -data=/lab/projectX/images")
    parser.add_argument('--preprocess_markers', type=lambda s: [x.strip() for x in s.split(',')], default=[],
        help='list of markers to preprocess : example -> -preprocess_markers=DAPI,CTNNB1,AMY2A,SST')
    parser.add_argument('--threshold_min', type=lambda s: float(s) / 100, default=0.0,
        help='percentage of intensity : example -> -threshold_min=1')
    parser.add_argument('--threshold_max', type=lambda s: float(s) / 100, default=1.0,
        help='percentage of intensity : example -> -threshold_max=99')
    parser.add_argument('--otsu_threshold_levels', type=_parse_otsu, default=None, metavar='LEVELS[:MIN:MAX]',
        help='otsu classes, e.g. 3 or with bin filtering 5:1:2 : example -> -otsu_threshold_levels=3')
    parser.add_argument('--flatten_spots', choices=['yes', 'no'], default='no',
        help='flatten spots : example -> -flatten_spots=no')
    parser.add_argument('--tile_size', type=int, default=0,
        help='number of pixels : example -> -tile_size=1844')
    parser.add_argument('--light_gradient', type=lambda s: 2 ** int(s), default=0,
        help='light gradient correction level : example -> -light_gradient=3')
    parser.add_argument('--balance_tiles', choices=['yes', 'no'], default='no',
        help='balance tiles : example -> -balance_tiles=yes')
    parser.add_argument('--stitch_size', type=int, default=0,
        help='number of pixels : example -> -stitch_size=20')
    parser.add_argument('--tophat_radius', type=int, default=0,
        help='tophat filter radius in pixels : example -> -tophat_radius=10')
    parser.add_argument('--exposure', type=lambda s: int(s) / 100, default=1.0,
        help='percentage of base intensity : example -> -exposure=300')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)


if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    preprocess_markers = sanitize_marker_list(args.preprocess_markers)
    thres_min = args.threshold_min
    thres_max = args.threshold_max
    if args.otsu_threshold_levels is not None:
        otsu_threshold_levels, bin_min, bin_max = args.otsu_threshold_levels
    flatten_spots = args.flatten_spots
    tile_size = args.tile_size
    light_gradient = args.light_gradient
    balance_tiles = args.balance_tiles
    stitch_size = args.stitch_size
    tophat_radius = args.tophat_radius
    exposure = args.exposure

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
    with open(os.path.join(data_folder, 'log_settings_preprocessing.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time preprocessing = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))

    log("Start time preprocessing")

    validate_marker_files(data_folder, preprocess_markers)

    try:
        os.mkdir(os.path.join(data_folder, 'preprocessed'))
    except OSError as error:
        print('>>> preprocessed folder already exists, overwriting results', flush=True)

    def handle_preprocess(marker, raw_img):
        preprocess_image(marker, downscale_images(raw_img))
        upscale_results(marker)

    iter_marker_images(data_folder, preprocess_markers, handle_preprocess)

    log("End time preprocessing")
