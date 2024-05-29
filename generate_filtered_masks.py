import math
import sys
import os
import datetime
import pandas as pd
import numpy as np
from skimage.io import imsave
from skimage.segmentation import relabel_sequential
from lmd.lib import SegmentationLoader, Collection

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
    min_dist = 1000000
    cell_id = None
    cell_id_coord = None
    for index, row in df.iterrows():
        curr_dist = math.dist([x, y], [row['x'], row['y']])
        if curr_dist < min_dist:
            min_dist = curr_dist
            cell_id = row['cell_id']
            cell_id_coord = [row['y'], row['x']]

    return cell_id, cell_id_coord



#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print('generate_filtered_masks.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-field=<optional, name of the column in cell_data.csv to filter the cells by> : example -> -field=leiden_id\n\t-values=<optional, values, comma-separated, present in the selected colum of cell_data.csv to filter the cells by> : example -> -values=3,6,7\n\t-lmd=<yes or no to create LMD XML file> : example -> -lmd=no\n\t-shape_dilation=<optional, dilation of shapes in pixels> : example -> -shape_dilation=0\n\t-convolution_smoothing=<optional, number of datapoints used for smoothing the shapes> : example -> -convolution_smoothing=15\n\t-path_optimization=<optional, "none"/"hilbert"/"greedy" optimization of cutting path between shapes> : example -> -path_optimization=none\n\t-distance_heuristic=<optional, nearest neighbour heuristic distance for merging shapes> : example -> -distance_heuristic=300\n\t-tile_size=<optional, number of pixels of each square tile segmented> : example -> -tile_size=2048\n\t-tile_overlap=<optional, number of pixels of surrounding overlap of each square tile segmented> : example -> -tile_overlap=128\n\t-tile_percentage_overlap=<optional, tile size\'s percentage of surrounding overlap of each square tile segmented> : example -> -tile_percentage_overlap=10\n\t-tile_relabel=<optional, yes or no to relabel sequentially the tile segments> : example -> -tile_relabel=yes\n\t-extend_tile=<yes or no to have bigger border tiles> : example -> -extend_tile=no', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print('generate_filtered_masks.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-field=<optional, name of the column in cell_data.csv to filter the cells by> : example -> -field=leiden_id\n\t-values=<optional, values, comma-separated, present in the selected colum of cell_data.csv to filter the cells by> : example -> -values=3,6,7\n\t-lmd=<yes or no to create LMD XML file> : example -> -lmd=no\n\t-shape_dilation=<optional, dilation of shapes in pixels> : example -> -shape_dilation=0\n\t-convolution_smoothing=<optional, number of datapoints used for smoothing the shapes> : example -> -convolution_smoothing=15\n\t-path_optimization=<optional, "none"/"hilbert"/"greedy" optimization of cutting path between shapes> : example -> -path_optimization=none\n\t-distance_heuristic=<optional, nearest neighbour heuristic distance for merging shapes> : example -> -distance_heuristic=300\n\t-tile_size=<optional, number of pixels of each square tile segmented> : example -> -tile_size=2048\n\t-tile_overlap=<optional, number of pixels of surrounding overlap of each square tile segmented> : example -> -tile_overlap=128\n\t-tile_percentage_overlap=<optional, tile size\'s percentage of surrounding overlap of each square tile segmented> : example -> -tile_percentage_overlap=10\n\t-tile_relabel=<optional, yes or no to relabel sequentially the tile segments> : example -> -tile_relabel=yes\n\t-extend_tile=<yes or no to have bigger border tiles> : example -> -extend_tile=no', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-field='):
                global field
                field = arg[7:]
            elif arg.startswith('-values='):
                global values
                values = arg[8:]
            elif arg.startswith('-lmd='):
                global lmd
                lmd = arg[5:]
            elif arg.startswith('-shape_dilation='):
                global shape_dilation
                shape_dilation = int(arg[16:])
            elif arg.startswith('-convolution_smoothing='):
                global convolution_smoothing
                convolution_smoothing = int(arg[23:])
            elif arg.startswith('-path_optimization='):
                global path_optimization
                path_optimization = arg[19:]
            elif arg.startswith('-distance_heuristic='):
                global distance_heuristic
                distance_heuristic = int(arg[20:])
            elif arg.startswith('-tile_size='):
                global tile_size
                tile_size = int(arg[11:])
            elif arg.startswith('-tile_overlap='):
                global tile_overlap
                tile_overlap = int(arg[14:])
            elif arg.startswith('-tile_percentage_overlap='):
                global tile_overlap_percentage
                tile_overlap_percentage = int(arg[25:])
            elif arg.startswith('-tile_relabel='):
                global tile_relabel
                tile_relabel = arg[14:]
            elif arg.startswith('-extend_tile='):
                global extend_tile
                extend_tile = arg[13:]


if __name__ =='__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()
    with open(data_folder + '/log_settings_filter.txt', 'w+', encoding='utf-8') as f:
        f.write(">>> Start time filter = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))
        f.close()

    print(">>> Start time generate_filtered_masks =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    #Load segmentation data in numpy array format
    labels = np.load(data_folder + '/analysis/segmentation_data.npy')
    df = pd.read_csv(data_folder + '/analysis/cell_data.csv')

    if (field != "" and values != ""):
        value_array = values.split(",")
        value_array = [int(i) for i in value_array]

        df_filtered = df[df[field].isin(value_array)]
        label_list = df_filtered['cell_id'].tolist()

        ix = np.in1d(labels.ravel(), label_list).reshape(labels.shape)
        labels = np.where(ix, labels, 0)

        np.save(data_folder + '/analysis/segmentation_data_filtered.npy', labels)
        print(">>> Filtered segmentation result numpy binary data saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

        labelsBinary = np.copy(labels)
        labelsBinary[labelsBinary > 0] = 1
        imsave(data_folder + "/analysis/segmentation_filtered_binary_mask.tif", np.uint16(labelsBinary * 65535))
        print(">>> Filtered segmentation result binary image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        del labelsBinary

        imsave(data_folder + "/analysis/segmentation_filtered_mask.tif", np.uint16(labels * 65535))
        print(">>> Filtered segmentation result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    if (lmd == "yes"):
        _1, top_left_point = closest_cell(df, df['x'].quantile(.25), df['y'].quantile(.25))
        _2, top_right_point = closest_cell(df, df['x'].quantile(.75), df['y'].quantile(.25))
        _3, bottom_right_point = closest_cell(df, df['x'].quantile(.75), df['y'].quantile(.75))
        print(">>> LMD calibration cell/points",_1,"-",top_left_point,",",_2,"-",top_right_point,",",_3,"-",bottom_right_point,"chosen =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        calibration_points = np.array([top_left_point, top_right_point, bottom_right_point])
        loader_config = { 'orientation_transform': np.array([[0, -1], [1, 0]]),
                         'shape_dilation': shape_dilation,
                         'convolution_smoothing': convolution_smoothing,
                         'path_optimization': path_optimization,
                         'distance_heuristic': distance_heuristic }
        sl = SegmentationLoader(config=loader_config)

        if (field != "" and values != ""):
            cell_sets = []
            value_array = values.split(",")
            value_array = [int(i) for i in value_array]
            for i in value_array:
                df_group = df_filtered[df_filtered[field] == i]
                label_group = list(filter(np.unique(labels).__contains__, df_group['cell_id'].tolist()))
                cell_sets = [{"classes": np.unique(label_group), "well": "A" + str(i)}]
                shape_collection = sl(labels, cell_sets, calibration_points)
                shape_collection.plot(fig_size=(10, 10), save_name=data_folder + "/analysis/lmd_shapes_plot_" + str(i) + ".jpg")
                print(">>> LMD shapes plot for group",str(i),"saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                print(">>> LMD shapes statistics for group",str(i)," block start=", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                print(shape_collection.stats())
                print(">>> LMD shapes statistics for group",str(i)," block end=", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                shape_collection.save(data_folder + "/analysis/lmd_shapes_" + str(i) + ".xml")
                print(">>> LMD XML file for group",str(i)," saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        else:
            cell_sets = [{"classes": np.unique(labels), "well": "1"}]
            shape_collection = sl(labels, cell_sets, calibration_points)
            shape_collection.plot(fig_size=(10, 10), save_name=data_folder + "/analysis/lmd_shapes_plot.jpg")
            print(">>> LMD shapes plot saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            print(">>> LMD shapes statistics block start=", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            shape_collection.stats()
            print(">>> LMD shapes statistics block end=", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            shape_collection.save(data_folder + "/analysis/lmd_shapes.xml")
            print(">>> LMD XML file saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)




    if (tile_size > 0):
        if (tile_overlap == 0 and tile_overlap_percentage > 0):
            tile_overlap = int(tile_size * tile_overlap_percentage / 100)

        num_rows = int(len(labels) / tile_size)
        if (len(labels) % tile_size != 0 and extend_tile == 'no'):
            num_rows = num_rows + 1
        num_columns = int(len(labels[0]) / tile_size)
        if (len(labels[0]) % tile_size != 0 and extend_tile == 'no'):
            num_columns = num_columns + 1
        for row in range(num_rows):
            for column in range(num_columns):
                min_y = (row * tile_size)
                max_y = ((row + 1) * tile_size)
                min_x = (column * tile_size)
                max_x = ((column + 1) * tile_size)
                if (tile_overlap > 0):
                    if (row > 0):
                        min_y = min_y - tile_overlap
                    if (column > 0):
                        min_x = min_x - tile_overlap
                    if (row < num_rows - 1):
                        max_y = max_y + tile_overlap
                    if (column < num_columns - 1):
                        max_x = max_x + tile_overlap
                if (extend_tile == 'yes'):
                    if (row == num_rows - 1):
                        max_y = len(labels)
                    if (column == num_columns - 1):
                        max_x = len(labels[0])

                tile = labels[min_y:max_y, min_x:max_x]
                if tile_relabel == "yes":
                    tile = relabel_sequential(tile)[0]
                tile_desc = str(row) + '_' + str(column)
                np.save(data_folder + '/analysis/segmentation_data_filtered_tile_' + tile_desc + '.npy', tile)
                print(">>> Filtered segmentation tile ",tile_desc," result numpy binary data saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

                tileBinary = np.copy(tile)
                tileBinary[tileBinary > 0] = 1
                imsave(data_folder + "/analysis/segmentation_filtered_tile_" + tile_desc + "_binary_mask.tif", np.uint16(tileBinary * 65535))
                print(">>> Filtered segmentation tile ",tile_desc," result binary image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                del tileBinary

                imsave(data_folder + "/analysis/segmentation_filtered_tile_" + tile_desc + "_mask.tif", np.uint16(tile * 65535))
                print(">>> Filtered segmentation tile ",tile_desc," result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)


    print(">>> End time generate_filtered_masks =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
