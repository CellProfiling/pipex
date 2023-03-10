import sys
import os
import datetime
import pandas as pd
import numpy as np
from skimage.io import imsave
from skimage.segmentation import relabel_sequential


data_folder = os.environ.get('PIPEX_DATA')
field = ""
values = ""
tile_size = 0
tile_overlap = 0
tile_overlap_percentage = 0
tile_relabel = 'no'
extend_tile = 'no'


#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print('generate_filtered_masks.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-field=<optional, name of the column in cell_data.csv to filter the cells by> : example -> -field=leiden_id\n\t-values=<optional, values, comma-separated, present in the selected colum of cell_data.csv to filter the cells by> : example -> -values=3,6,7\n\t-tile_size=<optional, number of pixels of each square tile segmented> : example -> -tile_size=2048\n\t-tile_overlap=<optional, number of pixels of surrounding overlap of each square tile segmented> : example -> -tile_overlap=128\n\t-tile_percentage_overlap=<optional, tile size\'s percentage of surrounding overlap of each square tile segmented> : example -> -tile_percentage_overlap=10\n\t-tile_relabel=<optional, yes or no to relabel sequentially the tile segments> : example -> -tile_relabel=yes\n\t-extend_tile=<yes or no to have bigger border tiles> : example -> -extend_tile=no', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print('generate_filtered_masks.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-field=<optional, name of the column in cell_data.csv to filter the cells by> : example -> -field=leiden_id\n\t-values=<optional, values, comma-separated, present in the selected colum of cell_data.csv to filter the cells by> : example -> -values=3,6,7\n\t-tile_size=<optional, number of pixels of each square tile segmented> : example -> -tile_size=2048\n\t-tile_overlap=<optional, number of pixels of surrounding overlap of each square tile segmented> : example -> -tile_overlap=128\n\t-tile_percentage_overlap=<optional, tile size\'s percentage of surrounding overlap of each square tile segmented> : example -> -tile_percentage_overlap=10\n\t-tile_relabel=<optional, yes or no to relabel sequentially the tile segments> : example -> -tile_relabel=yes\n\t-extend_tile=<yes or no to have bigger border tiles> : example -> -extend_tile=no', flush=True)
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
                tile_relabel = int(arg[14:])
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
                                
    print(">>> Start time generate_filtered_masks =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    
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
        print(">>> Filtered segmentation result numpy binary data saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

        labelsBinary = np.copy(labels)
        labelsBinary[labelsBinary > 0] = 1
        imsave(data_folder + "/analysis/segmentation_filtered_binary_mask.tif", np.uint16(labelsBinary * 65535))
        print(">>> Filtered segmentation result binary image saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
        del labelsBinary

        imsave(data_folder + "/analysis/segmentation_filtered_mask.tif", np.uint16(labels * 65535))
        print(">>> Filtered segmentation result image saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)  
        
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
                print(">>> Filtered segmentation tile ",tile_desc," result numpy binary data saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

                tileBinary = np.copy(tile)
                tileBinary[tileBinary > 0] = 1
                imsave(data_folder + "/analysis/segmentation_filtered_tile_" + tile_desc + "_binary_mask.tif", np.uint16(tileBinary * 65535))
                print(">>> Filtered segmentation tile ",tile_desc," result binary image saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
                del tileBinary

                imsave(data_folder + "/analysis/segmentation_filtered_tile_" + tile_desc + "_mask.tif", np.uint16(tile * 65535))
                print(">>> Filtered segmentation tile ",tile_desc," result image saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)  
                
           
    print(">>> End time generate_filtered_masks =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

