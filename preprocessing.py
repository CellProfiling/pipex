import sys
import os
import datetime
import argparse
import fnmatch
import PIL
import PIL.Image
import numpy as np
import math
from tifffile import TiffFile
from xml.etree import ElementTree

from skimage.filters import threshold_multiotsu
from skimage.io import imsave, imread

import cv2

import scipy.stats as st


PIL.Image.MAX_IMAGE_PIXELS = 10000000000
data_folder = './data'

thres_min = 0.0
thres_max = 1.0

tile_size = 0
bright_levels = 0
light_gradient = 0
flatten_spots = 'no'
balance_tiles = 'no'
stitch_size = 0

exposure = 1.0
heat_map='no'

bin_min = 1
bin_max = 2
       
       
def rescale_tile_intensity(x, mean_in, mean_factor, dev_in, dev_factor, bins):
    result = min(max(0, x + (mean_in * mean_factor) - mean_in), 1.0)
    result_dev = ((x - mean_in) * dev_factor) - (x - mean_in)
    result = min(max(0, result + result_dev if x >= mean_in else - result_dev), 1.0) 
    if (result > bins[bin_max]):
        result = bins[bin_max] + math.pow((result - bins[bin_max]) * 100, 0.6) / 100
    elif (result < bins[bin_min]):
        result = bins[bin_min] - math.pow((bins[bin_min] - result) * 100, 0.6) / 100

    return result
    
       
def apply_tile_compensation(f_name, np_img, bins, tile_data, local_tile_size, local_stitch_size):
    num_rows = int(len(np_img) / local_tile_size)
    num_columns = int(len(np_img[0]) / local_tile_size)
            
    for row in range(num_rows):
        for column in range(num_columns):
            if (tile_data['tiles'][row][column]['samples'] > 0):         
                tile = np_img[(row * local_tile_size):((row + 1) * local_tile_size), (column * local_tile_size):((column + 1) * local_tile_size)]
   
                mean_factor = tile_data['chosen']['mean'] / tile_data['tiles'][row][column]['mean']
                dev_factor = tile_data['chosen']['deviation'] / max(0.001, tile_data['tiles'][row][column]['deviation']) 
           
                values = tile[np.nonzero(tile)]
                curr_mean = np.mean(values)
                curr_dev = np.std(values)
            
                np_img[(row * local_tile_size):((row + 1) * local_tile_size), (column * local_tile_size):((column + 1) * local_tile_size)] = np.reshape(np.array([rescale_tile_intensity(x, curr_mean, mean_factor, curr_dev, dev_factor, bins) for x in np.ravel(tile)]), (local_tile_size, local_tile_size)) 
                    
    if (tile_size == local_tile_size):                
        print(">>> Balanced tiles for image",f_name,"=",datetime.datetime.now().strftime("%H:%M:%S"),flush=True)
            
    if (local_stitch_size > 0):
        for row in range(num_rows):
            for column in range(num_columns):
                if (row < num_rows - 1):
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
                        
        if (tile_size == local_tile_size):                
            print(">>> Smoothed stitched lines for image",f_name,"=", datetime.datetime.now().strftime("%H:%M:%S"),flush=True)
                               
    
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
            tile_x = tile_row * local_tile_size
            tile_y = tile_column * local_tile_size
	    
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

    for row in tile_data['tiles']:
        for curr_data in row:
            if curr_data['samples'] == 0:
                continue
            if not 'chosen' in tile_data:
                tile_data['chosen'] = curr_data
            else:
                if curr_data['samples'] > tile_data['chosen']['samples']:
                    tile_data['chosen'] = curr_data
                    
    if (tile_size == local_tile_size):                
        print(">>> Calculated reference tiles balance =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    
    return tile_data   
   
def rescale_gradient_intensity(intensity, ratio):
    return min(max(0, intensity * ratio), 1.0)             
    
           
def apply_tile_gradient_compensation(f_name, np_img, bins, gradient_data):
    num_rows = int(len(np_img) / tile_size)
    num_columns = int(len(np_img[0]) / tile_size)
    kernel_size = int(tile_size / light_gradient)
                              
    for row in range(num_rows):
        for column in range(num_columns):
            tile_ratio = gradient_data[row][column]['ratio_gradient']
            if (tile_ratio > 0): 
                for row_kernel in range(light_gradient):
                    for column_kernel in range(light_gradient):
                        kernel_ratio = gradient_data[row][column]['kernels'][row_kernel][column_kernel]['ratio_gradient']
                        if (kernel_ratio > 0 and kernel_ratio < tile_ratio):                
                            tile_kernel_sy = row * tile_size + row_kernel * kernel_size
                            tile_kernel_sx = column * tile_size + column_kernel * kernel_size
                            tile_kernel_ey = kernel_size + (tile_size % light_gradient) if row_kernel == light_gradient - 1 else 0
                            tile_kernel_ex = kernel_size + (tile_size % light_gradient) if column_kernel == light_gradient - 1 else 0
                            tile_kernel = np_img[tile_kernel_sy:(tile_kernel_sy + tile_kernel_ey), tile_kernel_sx:(tile_kernel_sx + tile_kernel_ex)]
                        
                            final_ratio = 1 + (bins[bin_max] * (tile_ratio / kernel_ratio)) / bins[bin_max]
                        
                            for i1 in range(kernel_size):
                                for i2 in range(kernel_size):
                                    np_img[tile_kernel_sy + i1][tile_kernel_sx + i2] = rescale_gradient_intensity(np_img[tile_kernel_sy + i1][tile_kernel_sx + i2], final_ratio)
                
                tile = np_img[(row * tile_size):((row + 1) * tile_size), (column * tile_size):((column + 1) * tile_size)]
                subtile_data = generate_tile_compensation_data(tile, bins, kernel_size)   
                apply_tile_compensation('', tile, bins, subtile_data, kernel_size, kernel_size / 10)

    imsave(data_folder + "/preprocessed/" + os.path.splitext(file)[0] + "_gradient.jpg", np.uint8(np_img * 255))
    
    print(">>> Applied gradient fix for image",f_name,"=", datetime.datetime.now().strftime("%H:%M:%S"),flush=True)
                  
                
def generate_tile_gradient_data(np_img, bins, tile_size):    
    num_rows = int(len(np_img) / tile_size)
    num_columns = int(len(np_img[0]) / tile_size)
       
    kernel_size = int(tile_size / light_gradient)
    gradient_data = []     
    max_gradient = 0
    max_gradient_kernel = []
    ratio_gradient = 0
    for row in range(num_rows):
        gradient_data.append([])
        for column in range(num_columns):
            gradient_data[row].append([])  
            tile = np_img[(row * tile_size):((row + 1) * tile_size), (column * tile_size):((column + 1) * tile_size)] 
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
                    tile_kernel_ey = kernel_size + (tile_size % light_gradient) if row_kernel == light_gradient - 1 else 0
                    tile_kernel_ex = kernel_size + (tile_size % light_gradient) if column_kernel == light_gradient - 1 else 0
                    tile_kernel = np_img[tile_kernel_sy:(tile_kernel_sy + tile_kernel_ey), tile_kernel_sx:(tile_kernel_sx + tile_kernel_ex)]
                    bin_count = np.histogram(np.ravel(tile_kernel), bins)[0]
                    samples = 0
                    samples = bin_count[bin_min] + bin_count[bin_max]
                    top_samples = bin_count[bin_max]
                    middle_samples = bin_count[bin_min]
                    tile_gradient_info = {}
                    tile_gradient_info['bins'] = bin_count  
                    if (top_samples > 0): 
                        tile_gradient_info['ratio_gradient'] = max(0.01, top_samples / middle_samples) if middle_samples > 0 else 1       
                    elif (middle_samples > 0):
                        tile_gradient_info['ratio_gradient'] = 0.01 
                    else:                              
                        samples = 1                        
                        tile_gradient_info['ratio_gradient'] = -1
                    tile_gradient_data[row_kernel].append(tile_gradient_info)  
                    if (max_gradient == 0 or tile_gradient_info['ratio_gradient'] > 0 and max_gradient < tile_gradient_info['ratio_gradient'] * math.log(samples)):
                        max_gradient = tile_gradient_info['ratio_gradient'] * math.log(samples)
                        ratio_gradient = tile_gradient_info['ratio_gradient']
                        max_gradient_kernel = [row_kernel, column_kernel]
                 
            tile_gradient['kernels'] = tile_gradient_data  
            tile_gradient['max_gradient'] = max_gradient   
            tile_gradient['ratio_gradient'] = ratio_gradient   
            tile_gradient['max_gradient_kernel'] = max_gradient_kernel     
            gradient_data[row][column] = tile_gradient
    
    print(">>> Calculated reference tiles gradient =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    
    return gradient_data


def apply_thresholds(f_name, np_img, threshold_min, threshold_max): 
    if (threshold_min > 0 or threshold_max < 100):
        np_img[np_img < threshold_min] = 0
        np_img[np_img > threshold_max] = 0
        
    if (f_name != ''):                
       print(">>> Applied min and max thresholds =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)


#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print ('preprocessing.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to \'./data\'> : example -> -data=/lab/projectX/images\n\t-threshold_min=<number, percentage of intensity> : example -> -threshold_min=1\n\t-threshold_max=<number, percentage of intensity> : example -> -threshold_max=99\n\t-tile_size=<number of pixels> : example -> -tile_size=1844\n\t-bright_levels=<number> : example -> -bright_levels=3\n\t-flatten_spots=<yes or no> : example -> -flatten_spots=no\n\t-light_gradient=<number> : example -> -light_gradient=3\n\t-balance_tiles=<yes or no> : example -> -balance_tiles=yes\n\t-stitch_size=<number of pixels> : example -> -stitch_size=20\n\t-exposure=<number, percentage of the base intensity> : example -> -exposure=300\n\t-heat_map=<yes or no> : example -> -heat_map=yes', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print ('preprocessing.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to \'./data\'> : example -> -data=/lab/projectX/images\n\t-threshold_min=<number, percentage of intensity> : example -> -threshold_min=1\n\t-threshold_max=<number, percentage of intensity> : example -> -threshold_max=99\n\t-tile_size=<number of pixels> : example -> -tile_size=1844\n\t-bright_levels=<number> : example -> -bright_levels=3\n\t-flatten_spots=<yes or no> : example -> -flatten_spots=no\n\t-light_gradient=<number> : example -> -light_gradient=3\n\t-balance_tiles=<yes or no> : example -> -balance_tiles=yes\n\t-stitch_size=<number of pixels> : example -> -stitch_size=20\n\t-exposure=<number, percentage of the base intensity> : example -> -exposure=300\n\t-heat_map=<yes or no> : example -> -heat_map=yes', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-threshold_min='):
                global thres_min
                thres_min = float(arg[15:]) / 100
            elif arg.startswith('-threshold_max='):
                global thres_max
                thres_max = float(arg[15:]) / 100
            elif arg.startswith('-tile_size='):
                global tile_size
                tile_size = int(arg[11:])
            elif arg.startswith('-bright_levels='):
                global bright_levels
                global bin_min
                global bin_max
                par_bl = arg[15:]
                if (":" in par_bl):
                    bright_levels = int(par_bl.split(":")[0])
                    bin_min = int(par_bl.split(":")[1])
                    bin_max = int(par_bl.split(":")[2])
                else: 
                    bright_levels = int(par_bl)
            elif arg.startswith('-flatten_spots='):
                global flatten_spots
                flatten_spots = arg[15:]
            elif arg.startswith('-light_gradient='):
                global light_gradient
                light_gradient = 2 ** int(arg[16:])
            elif arg.startswith('-balance_tiles='):
                global balance_tiles
                balance_tiles = arg[15:]
            elif arg.startswith('-stitch_size='):
                global stitch_size
                stitch_size = int(arg[13:])
            elif arg.startswith('-exposure='):
                global exposure
                exposure = int(arg[10:]) / 100
            elif arg.startswith('-heat_map='):
                global heat_map
                heat_map = arg[10:]
           
           
def preprocess_image(f_name, numpy_img):           
    #normalizing images
    np_img_min = np.amin(numpy_img)
    np_img_max = np.amax(numpy_img)
    np_img = (numpy_img - np_img_min) / (np_img_max - np_img_min)   
        
    apply_thresholds(f_name, np_img, thres_min, thres_max)
                
    if (balance_tiles == "yes"):
        global bright_levels
        if (flatten_spots == 'yes'): 
            c_otsu = threshold_multiotsu(np_img, 3 if bright_levels == 0 else bright_levels)
            bins = np.insert(c_otsu, 0, 0.0)
            bins = np.append(bins, 1.0)
            spot_threshold = bins[bin_max]
            np_img = np.reshape(np.array([x if x <= spot_threshold else spot_threshold + math.pow((x - spot_threshold) * 100, 0.6) / 100 for x in np.ravel(np_img)]), (len(np_img), len(np_img[0])))
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_flattened.jpg", np.uint8(np_img * 255))
                
        if bright_levels == 0:
            c_otsu = threshold_multiotsu(np_img, 3)
            bins = np.insert(c_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_3.jpg", np.uint8(regions * 255))
        
            c_otsu = threshold_multiotsu(np_img, 4)
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [2])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_4_1-2.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_4_1-3.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_4_2-3.jpg", np.uint8(regions * 255))
        
            c_otsu = threshold_multiotsu(np_img, 5)
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [2,3])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_1-2.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1,3])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_1-3.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [1,2])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_1-4.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0,3])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_2-3.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0,2])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_2-4.jpg", np.uint8(regions * 255))
        
            temp_otsu = c_otsu.copy()
            temp_otsu = np.delete(temp_otsu, [0,1])
            bins = np.insert(temp_otsu, 0, 0.0)
            regions = np.digitize(np_img, bins=bins)
            regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
            imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu_5_3-4.jpg", np.uint8(regions * 255))
          
            bright_levels = 3
                  
        c_otsu = threshold_multiotsu(np_img, bright_levels)
        bins = np.insert(c_otsu, 0, 0.0)
        bins = np.append(bins, 1.0)
        regions = np.digitize(np_img, bins=bins)
        regions = np.reshape(np.array([bins[x - 1] for x in np.ravel(regions)]), (len(np_img), len(np_img[0]))) 
        imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + "_otsu.jpg", np.uint8(regions * 255))
        
        if (light_gradient > 1):
            gradient_data = generate_tile_gradient_data(np_img, bins, tile_size)

        tile_data = generate_tile_compensation_data(np_img, bins, tile_size)

        if (light_gradient > 1):
            apply_tile_gradient_compensation(f_name, np_img, bins, gradient_data)
                
        apply_tile_compensation(f_name, np_img, bins, tile_data, tile_size, stitch_size)
        
    if (exposure != 1):
        np_img = np.reshape(np.array([max(0, min(1, x * exposure)) for x in np.ravel(np_img)]), (len(np_img), len(np_img[0])))
    imsave(data_folder + "/preprocessed/" + f_name, np.uint16(np_img * 65535))
    print(">>> Preprocessed result image",f_name,"saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
            
    if (heat_map == 'yes'):    
        imsave(data_folder + "/preprocessed/" + os.path.splitext(f_name)[0] + '_heatmap.jpg', cv2.applyColorMap(np.uint8(np.reshape(np.array([1 - x for x in np.ravel(np_img)]), (len(np_img), len(np_img[0]))) * 255), cv2.COLORMAP_PARULA))
        print(">>> Preprocessed heatmap for image ",f_name,"saved =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
                       

if __name__ =='__main__':
    options(sys.argv[1:])
    
    pidfile_filename = './RUNNING'    
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        
    print(">>> Start time preprocessing =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    try: 
        os.mkdir(data_folder + '/preprocessed') 
    except OSError as error: 
        print('>>> preprocessed folder already exists, overwriting results', flush=True)  
    
    tile_data = {}
    gradient_data = []
    for file in sorted(os.listdir(data_folder)):
        file_path = data_folder + '/' + file
        if fnmatch.fnmatch(file, '*Empty.*') or fnmatch.fnmatch(file, '*Blank.*') or os.path.isdir(file_path):
            continue

        next_try = False
        try:
            with TiffFile(file_path) as tif:
                if len(tif.series[0].pages) == 1:
                    preprocess_image(file, next(iter(tif.series[0].pages)).asarray())
                else:
                    for page in tif.series[0].pages:                        
                        biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                        preprocess_image(os.path.splitext(file)[0] + '_' + biomarker + '.tif', page.asarray())
        except:
            print('>>> checking type of ' + file_path + ', not QPTIFF', flush=True)  
            next_try = True

        if next_try:
            next_try = False
            try:
                np_img = imread(file_path)
                preprocess_image(file, np_img)
            except:
                next_try = True

        if next_try:
            try:
                np_img = np.array(PIL.Image.open(file_path))
                preprocess_image(file, np_img)
            except:
                print('>>> Could not read image ' + file_path, flush=True)  

                      
    print(">>> End time preprocessing =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
        
