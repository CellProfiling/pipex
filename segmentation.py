import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import datetime
import fnmatch
import PIL
import numpy as np
import math
import pandas as pd
from tifffile import TiffFile
from xml.etree import ElementTree

from stardist.models import StarDist2D
from stardist.plot import render_label

import cv2

from skimage.io import imsave, imread
from skimage.exposure import is_low_contrast, equalize_adapthist
from skimage.filters import threshold_multiotsu
from skimage.measure import regionprops
from skimage.segmentation import watershed, mark_boundaries, expand_labels, relabel_sequential
from skimage.transform import resize


PIL.Image.MAX_IMAGE_PIXELS = 10000000000
pipex_max_resolution = 30000
pipex_scale_factor = 0
pipex_img_bit_depth = 255
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
measure_markers = ""
bin_percentage = 100.0


def downscale_images(np_img):
    if (len(np_img) > pipex_max_resolution or len(np_img[0]) > pipex_max_resolution):
        global pipex_scale_factor
        if (pipex_scale_factor == 0):
            i = 2
            while (pipex_scale_factor == 0):
                if (max(len(np_img), len(np_img[0])) / i <= pipex_max_resolution):
                    pipex_scale_factor = i
                else:
                    i = i * 2
            global nuclei_diameter
            nuclei_diameter = nuclei_diameter / pipex_scale_factor
            global nuclei_expansion
            nuclei_expansion = nuclei_expansion / pipex_scale_factor
            global membrane_diameter
            membrane_diameter = membrane_diameter / pipex_scale_factor
        return resize(np_img, (len(np_img) / pipex_scale_factor, len(np_img[0]) / pipex_scale_factor), order=0, preserve_range=True, anti_aliasing=False).astype('uint16')
    return np_img


def upscale_results(df):
    if (pipex_scale_factor > 0):
        image = PIL.Image.open(data_folder + "/analysis/segmentation_binary_mask.tif")
        image = image.resize((image.size[0] * pipex_scale_factor, image.size[1] * pipex_scale_factor))
        image.save(data_folder + "/analysis/segmentation_binary_mask.tif")

        image = PIL.Image.open(data_folder + "/analysis/segmentation_mask_show.jpg")
        image = image.resize((image.size[0] * pipex_scale_factor, image.size[1] * pipex_scale_factor))
        image.save(data_folder + "/analysis/segmentation_mask_show.jpg")

        labels = np.load(data_folder + '/analysis/segmentation_data.npy')
        labels = labels.repeat(pipex_scale_factor, axis=0).repeat(pipex_scale_factor, axis=1)
        np.save(data_folder + '/analysis/segmentation_data.npy', labels)

        pipex_scale_factor_n2 = pow(pipex_scale_factor, 2)
        df['x'] = df['x'] * pipex_scale_factor
        df['y'] = df['y'] * pipex_scale_factor
        df['size'] = df['size'].apply(lambda x: int(x * pipex_scale_factor_n2))


def cell_segmentation(nuclei_img_orig, membrane_img_orig):
    #normalizing images
    nuclei_img = (nuclei_img_orig - np.amin(nuclei_img_orig)) / (np.amax(nuclei_img_orig) - np.amin(nuclei_img_orig))
    if membrane_diameter > 0:
        membrane_img = (membrane_img_orig - np.amin(membrane_img_orig)) / (np.amax(membrane_img_orig) - np.amin(membrane_img_orig))

    #if adjust_images parameter is required, check and enhance contrast and gamma levels of nuclei and membrane images
    if adjust_images == 'yes':
        if is_low_contrast(nuclei_img, fraction_threshold=0.05):
            print(">>> Nuclei image has low contrast, equalizing =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            nuclei_img = equalize_adapthist(nuclei_img, clip_limit=0.03)
            im = PIL.Image.fromarray((nuclei_img * 255).astype(np.uint8))
            im = im.convert('RGB')
            im.save(data_folder + "/analysis/quality_control/nuclei_equalized.jpg")
        nuclei_hist = np.histogram(nuclei_img, bins=3)[0]
        gamma = (nuclei_hist[1] + nuclei_hist[2]) / (nuclei_hist[0] / 3)
        if gamma < 0.01:
            print(">>> Nuclei image has low gamma, adjusting =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            nuclei_img = nuclei_img * np.log(0.01 / gamma)
            nuclei_img = np.clip(nuclei_img, 0, 1)
            im = PIL.Image.fromarray((nuclei_img * 255).astype(np.uint8))
            im = im.convert('RGB')
            im.save(data_folder + "/analysis/quality_control/nuclei_adjusted.jpg")

        if membrane_diameter > 0:
            if is_low_contrast(membrane_img, fraction_threshold=0.075):
                print(">>> Membrane image has low contrast, equalizing =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                membrane_img = equalize_adapthist(membrane_img, clip_limit=0.03)
                im = PIL.Image.fromarray((membrane_img * 255).astype(np.uint8))
                im = im.convert('RGB')
                im.save(data_folder + "/analysis/quality_control/membrane_equalized.jpg")
            membrane_hist = np.histogram(membrane_img, bins=3)[0]
            gamma = (membrane_hist[1] + membrane_hist[2]) / (membrane_hist[0] / 3)
            if gamma < 0.01:
                print(">>> Membrane image has low gamma, adjusting =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                membrane_img = membrane_img * np.log(0.01 / gamma)
                np.clip(membrane_img, 0, 1)
                im = PIL.Image.fromarray((membrane_img * 255).astype(np.uint8))
                im = im.convert('RGB')
                im.save(data_folder + "/analysis/quality_control/membrane_adjusted.jpg")

    #run stardist over nuclei image
    model = StarDist2D.from_pretrained('2D_versatile_fluo')
    sdLabels = None
    _ = None
    #for big images (>stardist_tile_threshold), run predict_instances_big method using 2048 square tiles
    if max(len(nuclei_img), len(nuclei_img[0])) > stardist_tile_threshold:
        sdLabels, _ = model.predict_instances_big(nuclei_img,axes='YX',block_size=2048,min_overlap=128,prob_thresh=(nuclei_definition if nuclei_definition > 0 else None), nms_thresh=(nuclei_closeness if nuclei_closeness > 0 else None))
    else:
        sdLabels, _ = model.predict_instances(nuclei_img,axes='YX',prob_thresh=(nuclei_definition if nuclei_definition > 0 else None), nms_thresh=(nuclei_closeness if nuclei_closeness > 0 else None))
    print(">>> Stardist prediction done =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    if nuclei_area_limit > 0:
        detections = regionprops(sdLabels)
        for curr_detection in detections:
            if curr_detection.area > nuclei_area_limit:
                sdLabels[sdLabels == curr_detection.label] = 0

    im = PIL.Image.fromarray((render_label(sdLabels, img=None) * 255).astype(np.uint8))
    im = im.convert('RGB')
    im.save(data_folder + "/analysis/quality_control/stardist_result.jpg")
    print(">>> Stardist base result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    #if nuclei_expansion parameter is required, expand labelled regions (avoiding overlap) to specified size
    if nuclei_expansion >= 0:
        sdLabelsExpanded = expand_labels(sdLabels, distance=nuclei_expansion)
        imsave(data_folder + "/analysis/quality_control/stardist_result_expanded.jpg", np.uint8(mark_boundaries(nuclei_img_orig, sdLabelsExpanded) * 255))
        print(">>> Stardist expanded result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
    else:
        sdLabelsExpanded = sdLabels


    #if membrane marker is provided, run custom watershed segmentation
    affected_by_membrane = set()
    if membrane_diameter > 0:
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
                #run a basic watershed with segments approximatelly dimensioned by membrane_diameter and high compactness
                num_markers = (len(tile) / membrane_diameter) * (len(tile[0]) / membrane_diameter)
                wsLabels = watershed(tile * 255, markers=num_markers, compactness=membrane_compactness)
                print(">>> Watershed of tile ",tile_desc," done =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

                imsave(data_folder + "/analysis/quality_control/wathershed_tile_" + tile_desc + "_result.jpg", np.uint8(mark_boundaries(tile_orig, wsLabels) * 255))
                print(">>> Watershed of tile ",tile_desc," result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

                membrane_keep_intensity = 0
                membrane_properties = {}
                membrane_region_list = regionprops(wsLabels, tile_orig)
                num_membrane_regions = 0
                for curr_membrane in membrane_region_list:
                    if curr_membrane.intensity_mean > membrane_intensity_mean and curr_membrane.area > membrane_diameter * membrane_diameter / 10:
                        if membrane_keep == "yes":
                            membrane_properties[curr_membrane.label] = curr_membrane.intensity_mean
                            num_membrane_regions = num_membrane_regions + 1
                    else:
                        wsLabels[wsLabels == curr_membrane.label] = 0
                if membrane_keep == "yes" and num_membrane_regions > 0:
                    membrane_keep_intensity = membrane_intensity_mean

                #merge resulting segments so they don't cut nuclei (not expanded)
                wsRegions = {}
                for row in range(len(wsLabels)):
                    for column in range(len(wsLabels[row])):
                        memLabel = wsLabels[row][column]
                        if memLabel == 0:
                            continue
                        if not memLabel in wsRegions:
                            wsRegions[memLabel] = set()
                        nucLabel = sdLabels[row + tile_x][column + tile_y]
                        if nucLabel != 0:
                            if not nucLabel in wsRegions[memLabel]:
                                wsRegions[memLabel].add(nucLabel)
                        elif membrane_keep == 'yes' and memLabel in membrane_properties:
                            if membrane_properties[memLabel] >= membrane_keep_intensity:
                                wsRegions[memLabel].add(membrane_keep_index)
                                membrane_keep_index = membrane_keep_index - 1
                            del membrane_properties[memLabel]

                #merge resulting segments that contain same nuclei and/or nothing
                wsRegionsMerged = {}
                for region in wsRegions:
                    regionValue = wsRegions[region]
                    found = False
                    for region2 in wsRegions:
                        if wsRegions[region2] == regionValue:
                            found = True
                            wsRegionsMerged[region] = region2
                            break
                    if not found:
                        wsRegionsMerged[region] = region

                for row in range(len(wsLabels)):
                    for column in range(len(wsLabels[row])):
                        if wsLabels[row][column] == 0:
                            continue
                        wsLabels[row][column] = wsRegionsMerged[wsLabels[row][column]]

                imsave(data_folder + "/analysis/quality_control/wathershed_tile_" + tile_desc +"_result_merged_by_nuclei.jpg", np.uint8(mark_boundaries(tile_orig, wsLabels) * 255))
                print(">>> Watershed of tile ",tile_desc," preliminary nuclei filter result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

                #cut expanded nuclei that collide with watershed segments
                for row in range(len(wsLabels)):
                    for column in range(len(wsLabels[row])):
                        nucLabel = sdLabels[row + tile_x][column + tile_y]
                        expLabel = sdLabelsExpanded[row + tile_x][column + tile_y]
                        if nucLabel != expLabel or expLabel == 0:
                            memLabel = wsLabels[row][column]
                            if memLabel == 0:
                                continue
                            if len(wsRegions[memLabel]) == 0:
                                sdLabelsExpanded[row + tile_x][column + tile_y] = 0
                                if expLabel > 0:
                                    affected_by_membrane.add(expLabel)
                            elif expLabel == 0 or not expLabel in wsRegions[memLabel]:
                                if membrane_keep == 'yes':
                                    membrane_only_label = 0 + sum(number for number in wsRegions[memLabel] if number < 0)
                                    sdLabelsExpanded[row + tile_x][column + tile_y] = membrane_only_label
                                else:
                                    sdLabelsExpanded[row + tile_x][column + tile_y] = 0
                                    if expLabel > 0:
                                        affected_by_membrane.add(expLabel)

        top_positive_label = np.max(sdLabelsExpanded)
        negative_label_mask = sdLabelsExpanded < 0
        sdLabelsExpanded[negative_label_mask] = top_positive_label - sdLabelsExpanded[negative_label_mask]
        affected_by_membrane.update(list(np.unique(sdLabelsExpanded[sdLabelsExpanded < 0])))


        #find rare disjointed segmented cells and using their associated convex hull instead
        segmentProperties = regionprops(sdLabelsExpanded)
        for segment in segmentProperties:
            contours, hierarchy = cv2.findContours(segment.image.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if len(contours) > 1:
                print(">>> Found disjointed segment " + str(segment.label) + ", using convex hull instead =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                bbox = segment.bbox
                filling_image = segment.image_convex
                if segment.solidity > 0.0:
                    for row in range(len(filling_image)):
                        for column in range(len(filling_image[row])):
                           filled_pixel = filling_image[row][column]
                           if (filled_pixel > 0):
                               sdLabelsExpanded[row + bbox[0]][column + bbox[1]] = segment.label

    else:
        del sdLabels

    sdLabelsExpanded = relabel_sequential(sdLabelsExpanded)[0]

    np.save(data_folder + '/analysis/segmentation_data.npy', sdLabelsExpanded)
    print(">>> Final joined segmentation result numpy binary data saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    imsave(data_folder + "/analysis/segmentation_mask_show.jpg", np.uint8(mark_boundaries(nuclei_img_orig, sdLabelsExpanded) * 255))
    print(">>> Final joined segmentation result image over nuclei saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    sdLabelsExpandedBinary = np.copy(sdLabelsExpanded)
    sdLabelsExpandedBinary[sdLabelsExpandedBinary > 0] = 1
    imsave(data_folder + "/analysis/segmentation_binary_mask.tif", np.uint8(sdLabelsExpandedBinary * 255))
    del sdLabelsExpandedBinary

    if (np.amax(sdLabelsExpanded) <= 255):
        imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint8(sdLabelsExpanded * 255))
    elif (np.amax(sdLabelsExpanded) <= 65535):
        imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint16(sdLabelsExpanded * 65535))
    else:
        imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint32(sdLabelsExpanded * 4294967296))
    print(">>> Final joined segmentation result image saved =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    return sdLabelsExpanded, affected_by_membrane


#Function to calculate the marker intensities for each cell
def marker_calculation(marker, marker_img, cellLabels, data_table):
    #applying segmentation mask over the marker image
    marker_img_otsu = (marker_img - np.amin(marker_img)) / (np.amax(marker_img) - np.amin(marker_img))
    marker_img_otsu = np.uint8(marker_img_otsu * 255)
    c_otsu = threshold_multiotsu(marker_img_otsu, 3)
    cell_binarized_threshold = ((max(0.0, c_otsu[0] - 1) / 255.0) * pipex_img_bit_depth) * bin_percentage / 100.0
    print(">>> Marker " + marker + " binarize threshold " + str(cell_binarized_threshold) + " =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
    markerProperties = regionprops(cellLabels, marker_img)
    #obtaining mean intensity per cell
    for cell in markerProperties:
        data_table[cell.label][marker] = cell.intensity_mean
        cell_image = cell.image_intensity
        cell_image = cell_image[(cell_image != 0) & (~np.isnan(cell_image))]
        data_table[cell.label][marker + '_local_90'] = np.quantile(cell_image, 0.9) if len(cell_image) > 0 else 0
        data_table[cell.label][marker + '_ratio_pixels'] = np.count_nonzero(cell_image >= cell_binarized_threshold) / cell.area
        data_table[cell.label][marker + '_bin'] = 1 if cell.intensity_mean > cell_binarized_threshold else 0

    print(">>> Marker " + marker + " calculated =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)


#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print('segmentation.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-nuclei_marker=<name before . in image file> : example, from image filename "reg001_cyc001_ch001_DAPI1.tif"-> -nuclei_marker=DAPI1\n\t-nuclei_diameter=<number of pixels> : example -> -nuclei_diameter=20\n\t-nuclei_expansion=<number of pixels, can be 0> : example -> -nuclei_expansion=20\n\t-nuclei_definition=<optional, gradation between 0.001 and 0.999> : example -> -nuclei_definition=0.1\n\t-nuclei_closeness=<optional, gradation between 0.001 and 0.999> : example -> -nuclei_closeness=0.6\n\t-nuclei_area_limit=<optional, number of pixels> : example -> -nuclei_area_limit=3200\n\t-membrane_marker=<optional, name before . in image file> : example, from image filename "reg001_cyc008_ch003_CDH1.tif" -> -membrane_marker=CDH1\n\t-membrane_diameter=<optional, number of pixels> : example -> -membrane_diameter=25\n\t-membrane_compactness=<optional, \"squareness\" of the membrane, gradation between 0.001 and 0.999> : example -> -membrane_compactness=0.5\n\t-membrane_keep=<yes or no to keep segmented membranes without nuclei> : example -> -membrane_keep=no\n\t-adjust_images=<yes or no to enhance poor images> : example -> -adjust_images=yes\n\t-custom_segmentation=<optional, file path to a pre-made custom segmentation> : example -> -custom_segmentation=/data/custom_seg.npy\n\t-measure_markers=<list of markers names before . in image files> : example -> measure_markers=AMY2A,SST,GORASP2\n\t-bin_percentage=<optional, multiplying factor for binarize threshold> : example -> -bin_percentage=10', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print ('segmentation.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-nuclei_marker=<name before . in image file> : example, from image filename "reg001_cyc001_ch001_DAPI1.tif"-> -nuclei_marker=DAPI1\n\t-nuclei_diameter=<number of pixels> : example -> -nuclei_diameter=20\n\t-nuclei_expansion=<number of pixels, can be 0> : example -> -nuclei_expansion=20\n\t-nuclei_definition=<optional, gradation between 0.001 and 0.999> : example -> -nuclei_definition=0.1\n\t-nuclei_closeness=<optional, gradation between 0.001 and 0.999> : example -> -nuclei_closeness=0.6\n\t-nuclei_area_limit=<optional, number of pixels> : example -> -nuclei_area_limit=3200\n\t-membrane_marker=<optional, name before . in image file> : example, from image filename "reg001_cyc008_ch003_CDH1.tif" -> -membrane_marker=CDH1\n\t-membrane_diameter=<optional, number of pixels> : example -> -membrane_diameter=25\n\t-membrane_compactness=<optional, \"squareness\" of the membrane, gradation between 0.001 and 0.999> : example -> -membrane_compactness=0.5\n\t-membrane_keep=<yes or no to keep segmented membranes without nuclei> : example -> -membrane_keep=no\n\t-adjust_images=<yes or no to enhance poor images> : example -> -adjust_images=yes\n\t-custom_segmentation=<optional, file path to a pre-made custom segmentation> : example -> -custom_segmentation=/data/custom_seg.npy\n\t-measure_markers=<list of markers names before . in image files> : example -> measure_markers=AMY2A,SST,GORASP2\n\t-bin_percentage=<optional, multiplying factor for binarize threshold> : example -> -bin_percentage=10', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-nuclei_marker='):
                global nuclei_marker
                nuclei_marker = arg[15:]
            elif arg.startswith('-nuclei_diameter='):
                global nuclei_diameter
                nuclei_diameter = int(arg[17:])
            elif arg.startswith('-nuclei_expansion='):
                global nuclei_expansion
                nuclei_expansion = int(arg[18:])
            elif arg.startswith('-nuclei_definition='):
                global nuclei_definition
                nuclei_definition = float(arg[19:])
            elif arg.startswith('-nuclei_closeness='):
                global nuclei_closeness
                nuclei_closeness = float(arg[18:])
            elif arg.startswith('-nuclei_area_limit='):
                global nuclei_area_limit
                nuclei_area_limit = float(arg[19:])
            elif arg.startswith('-membrane_marker='):
                global membrane_marker
                membrane_marker = arg[17:]
            elif arg.startswith('-membrane_diameter='):
                global membrane_diameter
                membrane_diameter = int(arg[19:])
            elif arg.startswith('-membrane_compactness='):
                global membrane_compactness
                membrane_compactness = float(arg[22:])
            elif arg.startswith('-membrane_keep='):
                global membrane_keep
                membrane_keep = arg[15:]
            elif arg.startswith('-adjust_images='):
                global adjust_images
                adjust_images = arg[15:]
            elif arg.startswith('-custom_segmentation='):
                global custom_segmentation
                custom_segmentation = arg[21:]
            elif arg.startswith('-measure_markers='):
                global measure_markers
                measure_markers = arg[17:].split(",")
            elif arg.startswith('-bin_percentage='):
                global bin_percentage
                bin_percentage = float(arg[16:])


if __name__ =='__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()
    with open(data_folder + '/log_settings_segmentation.txt', 'w+', encoding='utf-8') as f:
        f.write(">>> Start time segmentation = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))
        f.close()

    print(">>> Start time segmentation =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    try:
        os.mkdir(data_folder + '/analysis')
    except OSError as error:
        print('>>> analysis folder already exists, overwriting results', flush=True)
    try:
        os.mkdir(data_folder + '/analysis/quality_control')
    except OSError as error:
        print('>>> quality_control folder already exists, overwriting results', flush=True)

    #finding nuclei and membrane image filenames by marker name
    nuclei_img = None
    membrane_img = None
    if custom_segmentation == "":
        for file in os.listdir(data_folder):
            file_path = data_folder + '/' + file
            if os.path.isdir(file_path):
                continue

            next_try = False
            try:
                with TiffFile(file_path) as tif:
                    if len(tif.series[0].pages) == 1:
                        if fnmatch.fnmatch(file, '*' + nuclei_marker + '.*'):
                            nuclei_img = downscale_images(next(iter(tif.series[0].pages)).asarray())
                        if membrane_marker != "" and fnmatch.fnmatch(file, '*' + membrane_marker + '.*'):
                            membrane_img = downscale_images(next(iter(tif.series[0].pages)).asarray())
                    else:
                        for page in tif.series[0].pages:
                            biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                            if biomarker == nuclei_marker:
                                nuclei_img = downscale_images(page.asarray())
                            if biomarker == membrane_marker:
                                membrane_img = downscale_images(page.asarray())
            except Exception as e:
                print('>>> checking type of ' + file_path + ', not QPTIFF', flush=True)
                print('>>> ', e, flush=True)
                next_try = True

            if next_try:
                next_try = False
                try:
                    if fnmatch.fnmatch(file, '*' + nuclei_marker + '.*'):
                        curr_image = imread(file_path)
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        nuclei_img = downscale_images(curr_image)
                    if membrane_marker != "" and fnmatch.fnmatch(file, '*' + membrane_marker + '.*'):
                        curr_image = imread(file_path)
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        curr_image = curr_image[:, :, 0]
                        membrane_img = downscale_images(curr_image)
                except Exception as e:
                    next_try = True
                    print('>>> ', e, flush=True)

            if next_try:
                try:
                    if fnmatch.fnmatch(file, '*' + nuclei_marker + '.*'):
                        curr_image = np.array(PIL.Image.open(file_path))
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        nuclei_img = downscale_images(curr_image)
                    if membrane_marker != "" and fnmatch.fnmatch(file, '*' + membrane_marker + '.*'):
                        curr_image = np.array(PIL.Image.open(file_path))
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        membrane_img = downscale_images(curr_image)
                except Exception as e:
                    print('>>> Could not read image ' + file_path, flush=True)
                    print('>>> ', e, flush=True)

        #performing segmentation
        cellLabels, membraneAffected = cell_segmentation(nuclei_img, membrane_img)
    else:
        try_image = False
        label_data = None
        membraneAffected = set()
        try:
            cellLabels = downscale_images(np.load(data_folder + '/' + custom_segmentation))
            np.save(data_folder + '/analysis/segmentation_data.npy', cellLabels)
            print(">>> Custom segmentation result numpy binary data saved =",
                  datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

            cellLabelsBinary = np.copy(cellLabels)
            cellLabelsBinary[cellLabelsBinary > 0] = 1
            imsave(data_folder + "/analysis/segmentation_binary_mask.tif", np.uint8(cellLabelsBinary * 255))
            del cellLabelsBinary

            if (np.amax(cellLabels) <= 255):
                imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint8(cellLabels * 255))
            elif (np.amax(cellLabels) <= 65535):
                imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint16(cellLabels * 65535))
            else:
                imsave(data_folder + "/analysis/segmentation_mask.tif", np.uint32(cellLabels * 4294967296))
            print(">>> Custom segmentation binary result image saved =",
                  datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        except Exception as e:
            try_image = True
            print('>>> ', e, flush=True)

        if try_image:
            try:
                cellLabels = downscale_images(np.array(PIL.Image.open(data_folder + '/' + custom_segmentation)))
            except Exception as e:
                print('>>> Could not custom segmentation file ' + custom_segmentation, flush=True)
                print('>>> ', e, flush=True)

    del nuclei_img
    del membrane_img
    #creating data table with segmented cell information
    cellProperties = regionprops(cellLabels)
    data_table = {}
    for cell in cellProperties:
        data_cell = {}
        data_cell['cell_id'] = cell.label
        data_cell['size'] = cell.area
        data_cell['x'] = int(cell.centroid[1])
        data_cell['y'] = int(cell.centroid[0])
        data_cell['memref'] = int(cell.label in membraneAffected)
        data_table[cell.label] = data_cell

    #calculating marker intensities per cell
    for file in os.listdir(data_folder):
        file_path = data_folder + '/' + file
        if os.path.isdir(file_path):
            continue

        next_try = False
        try:
            with TiffFile(file_path) as tif:
                if len(tif.series[0].pages) == 1:
                    for marker in measure_markers:
                        if marker + '.' in file:
                            marker_calculation(marker, downscale_images(next(iter(tif.series[0].pages)).asarray()), cellLabels, data_table)
                            break
                else:
                    for page in tif.series[0].pages:
                        biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                        if biomarker in measure_markers:
                            marker_calculation(biomarker, downscale_images(page.asarray()), cellLabels, data_table)
        except Exception as e:
            print('>>> checking type of ' + file_path + ', not QPTIFF', flush=True)
            print('>>> ', e, flush=True)
            next_try = True

        if next_try:
            next_try = False
            try:
                for marker in measure_markers:
                    if marker + '.' in file:
                        curr_image = imread(file_path)
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        marker_calculation(marker, downscale_images(curr_image), cellLabels, data_table)
                        break
            except Exception as e:
                next_try = True
                print('>>> ', e, flush=True)

        if next_try:
            try:
                for marker in measure_markers:
                    if marker + '.' in file:
                        curr_image = np.array(PIL.Image.open(file_path))
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        marker_calculation(marker, downscale_images(curr_image), cellLabels, data_table)
                        break
            except Exception as e:
                print('>>> Could not read image ' + file_path, flush=True)
                print('>>> ', e, flush=True)


    #dumpming data_table in cell_data.csv file
    df = pd.DataFrame.from_dict(data_table, orient='index')
    upscale_results(df)
    binarized_marker_columns = []
    for marker in measure_markers:
        binarized_marker_columns.append(marker + "_local_90")
        binarized_marker_columns.append(marker + "_ratio_pixels")
        binarized_marker_columns.append(marker + "_bin")
    measure_markers.extend(binarized_marker_columns)
    measure_markers.insert(0, 'memref')
    measure_markers.insert(0, 'y')
    measure_markers.insert(0, 'x')
    measure_markers.insert(0, 'size')
    measure_markers.insert(0, 'cell_id')
    df = df.reindex(measure_markers, axis=1)
    df.to_csv(data_folder + '/analysis/cell_data.csv', index=False)

    print(">>> End time segmentation =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
