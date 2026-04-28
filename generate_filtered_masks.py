import sys
import os
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from skimage.io import imsave
from skimage.segmentation import relabel_sequential, mark_boundaries
from lmd.lib import SegmentationLoader
import cv2
from pipex_utils import log


data_folder = os.environ.get('PIPEX_DATA')
field = ""
values = ""
tile_size = 0
tile_overlap = 0
tile_overlap_percentage = 0
tile_relabel = 'no'
extend_tile = 'no'
lmd = "no"
shape_dilation = 0
convolution_smoothing = 15
path_optimization = "none"
distance_heuristic = 300


def closest_cell(df, x, y):
    dists = np.hypot(df['x'].to_numpy() - x, df['y'].to_numpy() - y)
    idx = np.argmin(dists)
    row = df.iloc[idx]
    return row['cell_id'], [row['y'], row['x']]



#Function to handle the command line parameters passed
def options(argv):
    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='generate_filtered_masks.py')
    parser.add_argument('--data', default=os.environ.get('PIPEX_DATA'),
        help='path to images folder : example -> -data=/lab/projectX/images')
    parser.add_argument('--field', default='',
        help='column in cell_data.csv to filter cells by : example -> -field=leiden_id')
    parser.add_argument('--values', default='',
        help='comma-separated values in the selected column to filter by : example -> -values=3,6,7')
    parser.add_argument('--lmd', choices=['yes', 'no'], default='no',
        help='create LMD XML file : example -> -lmd=no')
    parser.add_argument('--shape_dilation', type=int, default=0,
        help='dilation of shapes in pixels : example -> -shape_dilation=0')
    parser.add_argument('--convolution_smoothing', type=int, default=15,
        help='number of datapoints for smoothing shapes : example -> -convolution_smoothing=15')
    parser.add_argument('--path_optimization', choices=['none', 'hilbert', 'greedy'], default='none',
        help='optimization of cutting path between shapes : example -> -path_optimization=none')
    parser.add_argument('--distance_heuristic', type=int, default=300,
        help='nearest neighbour heuristic distance for merging shapes : example -> -distance_heuristic=300')
    parser.add_argument('--tile_size', type=int, default=0,
        help='number of pixels per square tile : example -> -tile_size=2048')
    parser.add_argument('--tile_overlap', type=int, default=0,
        help='pixels of surrounding overlap per tile : example -> -tile_overlap=128')
    parser.add_argument('--tile_percentage_overlap', type=int, default=0, dest='tile_overlap_percentage',
        help="tile size's percentage of surrounding overlap : example -> -tile_percentage_overlap=10")
    parser.add_argument('--tile_relabel', choices=['yes', 'no'], default='no',
        help='relabel sequentially the tile segments : example -> -tile_relabel=yes')
    parser.add_argument('--extend_tile', choices=['yes', 'no'], default='no',
        help='extend border tiles : example -> -extend_tile=no')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)


if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    field = args.field
    values = args.values
    lmd = args.lmd
    shape_dilation = args.shape_dilation
    convolution_smoothing = args.convolution_smoothing
    path_optimization = args.path_optimization
    distance_heuristic = args.distance_heuristic
    tile_size = args.tile_size
    tile_overlap = args.tile_overlap
    tile_overlap_percentage = args.tile_overlap_percentage
    tile_relabel = args.tile_relabel
    extend_tile = args.extend_tile
    df_filtered = None
    value_array = []

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
    with open(os.path.join(data_folder, 'log_settings_filter.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time filter = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))

    log("Start time generate_filtered_masks")

    #Load segmentation data in numpy array format
    labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
    df = pd.read_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'))

    if field != "" and values != "":
        value_array = values.split(",")

        df_filtered = df[df[field].astype(str).isin(value_array)]
        label_list = df_filtered['cell_id'].tolist()

        labels = np.where(np.isin(labels, label_list), labels, 0)

        np.save(os.path.join(data_folder, 'analysis', 'segmentation_data_filtered.npy'), labels)
        log("Filtered segmentation result numpy binary data saved")

        labels_binary = np.copy(labels)
        labels_binary[labels_binary > 0] = 1
        imsave(os.path.join(data_folder, 'analysis', 'segmentation_filtered_binary_mask.tif'), np.uint8(labels_binary * 255))
        log("Filtered segmentation result binary image saved")
        del labels_binary

        black_canvas = np.zeros((len(labels), len(labels[0]), 3), dtype=np.uint8)
        boundaries_rgb = np.uint8(mark_boundaries(black_canvas, labels, color=(0, 1, 0)) * 255)
        boundaries_rgba = np.dstack([boundaries_rgb, np.where(boundaries_rgb.any(axis=2), 255, 0).astype(np.uint8)])
        imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_boundaries.png"), boundaries_rgba)
        log("Final segmentation boundaries overlay saved")

        if np.amax(labels) <= 255:
            imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_mask.tif"), labels.astype(np.uint8))
        elif np.amax(labels) <= 65535:
            imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_mask.tif"), labels.astype(np.uint16))
        else:
            imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_mask.tif"), labels.astype(np.uint32))
        log("Filtered segmentation result image saved")

    if lmd == "yes":
        _1, top_left_point = closest_cell(df, df['x'].quantile(.25), df['y'].quantile(.25))
        _2, top_right_point = closest_cell(df, df['x'].quantile(.75), df['y'].quantile(.25))
        _3, bottom_right_point = closest_cell(df, df['x'].quantile(.75), df['y'].quantile(.75))
        log("LMD calibration cell/points: cell_id " + str(int(_1)) + " - " + str(top_left_point) + ", cell_id " + str(int(_2)) + " - " + str(top_right_point) + ", cell_id " + str(int(_3)) + " - " + str(bottom_right_point) + " chosen")

        cells = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
        cells = np.where(np.isin(cells, [_1, _2, _3]), cells, 0)
        cells[cells > 0] = 1

        overlay_img = np.zeros((len(cells), len(cells[0]), 4), np.uint8)
        overlay_img[:,:,1] = np.where(cells > 0, 255, 0)
        overlay_img[:,:,2] = np.where(cells > 0, 255, 0)
        overlay_img[:,:,3] = np.where(cells > 0, 32, 0)
        overlay_img = cv2.line(overlay_img, (int(top_left_point[1]), int(max(0, top_left_point[0] - 10))), (int(top_left_point[1]), int(min(len(cells), top_left_point[0] + 10))), (255, 0, 0, 255), 1)
        overlay_img = cv2.line(overlay_img, (int(max(0, top_left_point[1] - 10)), int(top_left_point[0])), (int(min(len(cells[0]), top_left_point[1] + 10)), int(top_left_point[0])), (255, 0, 0, 255), 1)
        overlay_img = cv2.line(overlay_img, (int(top_right_point[1]), int(max(0, top_right_point[0] - 10))), (int(top_right_point[1]), int(min(len(cells), top_right_point[0] + 10))), (255, 0, 0, 255), 1)
        overlay_img = cv2.line(overlay_img, (int(max(0, top_right_point[1] - 10)), int(top_right_point[0])), (int(min(len(cells[0]), top_right_point[1] + 10)), int(top_right_point[0])), (255, 0, 0, 255), 1)
        overlay_img = cv2.line(overlay_img, (int(bottom_right_point[1]), int(max(0, bottom_right_point[0] - 10))), (int(bottom_right_point[1]), int(min(len(cells), bottom_right_point[0] + 10))), (255, 0, 0, 255), 1)
        overlay_img = cv2.line(overlay_img, (int(max(0, bottom_right_point[1] - 10)), int(bottom_right_point[0])), (int(min(len(cells[0]), bottom_right_point[1] + 10)), int(bottom_right_point[0])), (255, 0, 0, 255), 1)
        cv2.imwrite(os.path.join(data_folder, "analysis", "lmd_calibration_points_overlay.png"), overlay_img)
        log("LMD calibration points transparent image saved")

        calibration_points = np.array([top_left_point, top_right_point, bottom_right_point])
        loader_config = { 'orientation_transform': np.array([[0, -1], [1, 0]]),
                         'shape_dilation': shape_dilation,
                         'convolution_smoothing': convolution_smoothing,
                         'path_optimization': path_optimization,
                         'distance_heuristic': distance_heuristic }
        sl = SegmentationLoader(config=loader_config)

        if field != "" and values != "":
            unique_labels = set(np.unique(labels[labels > 0]).tolist())
            for i in value_array:
                df_group = df_filtered[df_filtered[field].astype(str) == i]
                label_group = list(filter(unique_labels.__contains__, df_group['cell_id'].tolist()))
                cell_sets = [{"classes": np.unique(label_group), "well": "A" + str(i)}]
                shape_collection = sl(labels, cell_sets, calibration_points)
                shape_collection.plot(fig_size=(10, 10), save_name=os.path.join(data_folder, "analysis", "lmd_shapes_plot_" + str(i) + ".jpg"))
                log("LMD shapes plot for group " + str(i) + " saved")
                log("LMD shapes statistics for group " + str(i) + " block start")
                print(shape_collection.stats())
                log("LMD shapes statistics for group " + str(i) + " block end")
                shape_collection.save(os.path.join(data_folder, "analysis", "lmd_shapes_" + str(i) + ".xml"))
                log("LMD XML file for group " + str(i) + " saved")
        else:
            cell_sets = [{"classes": np.unique(labels[labels > 0]), "well": "1"}]
            shape_collection = sl(labels, cell_sets, calibration_points)
            shape_collection.plot(fig_size=(10, 10), save_name=os.path.join(data_folder, "analysis", "lmd_shapes_plot.jpg"))
            log("LMD shapes plot saved")
            log("LMD shapes statistics block start")
            print(shape_collection.stats())
            log("LMD shapes statistics block end")
            shape_collection.save(os.path.join(data_folder, "analysis", "lmd_shapes.xml"))
            log("LMD XML file saved")




    if tile_size > 0:
        if tile_overlap == 0 and tile_overlap_percentage > 0:
            tile_overlap = int(tile_size * tile_overlap_percentage / 100)

        num_rows = int(len(labels) / tile_size)
        if len(labels) % tile_size != 0 and extend_tile == 'no':
            num_rows = num_rows + 1
        num_columns = int(len(labels[0]) / tile_size)
        if len(labels[0]) % tile_size != 0 and extend_tile == 'no':
            num_columns = num_columns + 1
        for row in range(num_rows):
            for column in range(num_columns):
                min_y = (row * tile_size)
                max_y = ((row + 1) * tile_size)
                min_x = (column * tile_size)
                max_x = ((column + 1) * tile_size)
                if tile_overlap > 0:
                    if row > 0:
                        min_y = min_y - tile_overlap
                    if column > 0:
                        min_x = min_x - tile_overlap
                    if row < num_rows - 1:
                        max_y = max_y + tile_overlap
                    if column < num_columns - 1:
                        max_x = max_x + tile_overlap
                if extend_tile == 'yes':
                    if row == num_rows - 1:
                        max_y = len(labels)
                    if column == num_columns - 1:
                        max_x = len(labels[0])

                tile = labels[min_y:max_y, min_x:max_x]
                if tile_relabel == "yes":
                    tile = relabel_sequential(tile)[0]
                tile_desc = str(row) + '_' + str(column)
                np.save(os.path.join(data_folder, 'analysis', 'segmentation_data_filtered_tile_' + tile_desc + '.npy'), tile)
                log("Filtered segmentation tile " + tile_desc + " result numpy binary data saved")

                tileBinary = np.copy(tile)
                tileBinary[tileBinary > 0] = 1
                imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_tile_" + tile_desc + "_binary_mask.tif"), np.uint8(tileBinary * 255))
                log("Filtered segmentation tile " + tile_desc + " result binary image saved")
                del tileBinary

                tile_canvas = np.zeros((len(tile), len(tile[0]), 3), dtype=np.uint8)
                tile_boundaries_rgb = np.uint8(mark_boundaries(tile_canvas, tile, color=(0, 1, 0)) * 255)
                tile_boundaries_rgba = np.dstack([tile_boundaries_rgb, np.where(tile_boundaries_rgb.any(axis=2), 255, 0).astype(np.uint8)])
                imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_tile_" + tile_desc + "_boundaries.png"), tile_boundaries_rgba)
                log("Final segmentation " + tile_desc + " boundaries overlay saved")

                if np.amax(tile) <= 255:
                    imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_tile_" + tile_desc + "_mask.tif"), tile.astype(np.uint8))
                elif np.amax(tile) <= 65535:
                    imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_tile_" + tile_desc + "_mask.tif"), tile.astype(np.uint16))
                else:
                    imsave(os.path.join(data_folder, "analysis", "segmentation_filtered_tile_" + tile_desc + "_mask.tif"), tile.astype(np.uint32))
                log("Filtered segmentation tile " + tile_desc + " result image saved")


    log("End time generate_filtered_masks")
