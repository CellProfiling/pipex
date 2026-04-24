import sys
import os
import argparse
import datetime
import pandas as pd
import numpy as np
import diplib as dip
import geojson
from pipex_utils import log, sanitize_marker_list, validate_marker_columns


data_folder = os.environ.get('PIPEX_DATA')
included_markers = []
cluster_id = ""
cluster_color = ""


#Function to handle the command line parameters passed
def options(argv):
    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='generate_geojson.py')
    parser.add_argument('--data', default=os.environ.get('PIPEX_DATA'),
        help='path to images folder : example -> -data=/lab/projectX/images')
    parser.add_argument('--included_markers', type=lambda s: s.split(','), default=[],
        help='list of markers to include : example -> -included_markers=AMY2A,SST,GORASP2')
    parser.add_argument('--cluster_id', default='',
        help='column for cluster id information : example -> -cluster_id=kmeans')
    parser.add_argument('--cluster_color', default='',
        help='column for cluster color information : example -> -cluster_color=kmeans_color')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)

if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    included_markers = sanitize_marker_list(args.included_markers)
    cluster_id = args.cluster_id
    cluster_color = args.cluster_color

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
    with open(os.path.join(data_folder, 'log_settings_geojson.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time geojson = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))

    log("Start time generate_geojson")

    #Load segmentation data in numpy array format
    labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
    df = pd.read_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'))

    #Getting the list of marker names
    markers = list(df.columns.values)
    markers = markers[(df.columns.get_loc("y") + 1):]
    #saveguard if analysis.py has been executed before and cluster_id + cluster_color already exists
    if 'cluster_id' in markers:
        markers = markers[:-(len(df.columns) - df.columns.get_loc("cluster_id"))]
    elif 'leiden' in markers:
        markers = markers[:-(len(df.columns) - df.columns.get_loc("leiden"))]
    elif 'kmeans' in markers:
        markers = markers[:-(len(df.columns) - df.columns.get_loc("kmeans"))]

    # If a specific list of markers is informed, we use it
    if len(included_markers) > 0:
        markers = included_markers

    validate_marker_columns(df, markers)

    #calculate segmentation polygons via fast chaincodes from diplib
    chaincodes = dip.GetImageChainCodes(labels.astype('uint32'))
    borders = {}
    for chaincode in chaincodes:
        borders[chaincode.objectID] = np.array(chaincode.Polygon()).tolist()
    log("Labelled regions to approximate polygons conversion finished")

    df = df.set_index('cell_id')

    #generating geojson data to import in qupath
    GEOdata = []
    for label in borders:
        if label not in df.index:
            continue
        cell = df.loc[label]
        final_coords = borders[label]
        final_coords.append(final_coords[0])
        cell_data = {}
        cell_data["id"] = "cell" + str(label)
        cell_data["type"] = "Feature"
        cell_data["geometry"] = {
            "type" : "Polygon",
            "coordinates" : [final_coords],
        }
        cell_data["properties"] = {
            "name" : "cell" + str(label),
            "object_type" : "detection",
            "isLocked" : "false",
            "collectionIndex": 0,
        }

        cell_data["properties"]["measurements"] = []
        for marker in markers:
            cell_data["properties"]["measurements"].append({
                "name" : marker,
                "value" : float(cell[marker])
                })

        #if cluster_id parameter is selected, add cluster_id and cluster_color
        if cluster_id != '' and cluster_id in df.columns:
            cell_data["properties"]["classification"] = {
                "name": str(cell[cluster_id]),
                "colorRGB": str(cell[cluster_color] if cluster_color != '' else '')
                }
        GEOdata.append(cell_data)

    #dump GEOdata variable to json file
    with open(os.path.join(data_folder, 'analysis', 'cell_segmentation_geo.json'), 'w') as outfile:
        geojson.dump(GEOdata, outfile)

    log("End time generate_geojson")
