import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import datetime
import pandas as pd
import numpy as np
import diplib as dip
import skimage as sk
import geojson


data_folder = os.environ.get('PIPEX_DATA')
included_markers = []
cluster_id = ""
cluster_color = ""


#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print('generate_geojson.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-included_markers=<optional, list of present specific markers to include> : example -> -included_markers=AMY2A,SST,GORASP2\n\t-cluster_id=<optional, name of the column to add as cluster id information from cell_data.csv> : example -> -cluster_id=kmeans\n\t-cluster_color=<optional, name of the column to add as cluster information color from cell_data.csv> : example -> -cluster_color=kmeans_color', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print('generate_geojson.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-included_markers=<optional, list of present specific markers to include> : example -> -included_markers=AMY2A,SST,GORASP2\n\t-cluster_id=<optional, name of the column to add as cluster id information from cell_data.csv> : example -> -cluster_id=kmeans\n\t-cluster_color=<optional, name of the column to add as cluster information color from cell_data.csv> : example -> -cluster_color=kmeans_color', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-included_markers='):
                global included_markers
                included_markers = arg[18:].split(",")
            elif arg.startswith('-cluster_id='):
                global cluster_id
                cluster_id = arg[12:]
            elif arg.startswith('-cluster_color='):
                global cluster_color
                cluster_color = arg[15:]


if __name__ =='__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()

    print(">>> Start time generate_geojson =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #Load segmentation data in numpy array format
    labels = np.load(data_folder + '/analysis/segmentation_data.npy')
    df = pd.read_csv(data_folder + '/analysis/cell_data.csv')

    markers = []
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

    #calculate segmentation polygons via fast chaincodes from diplib
    chaincodes = dip.GetImageChainCodes(labels.astype('uint32'))
    borders = {}
    for chaincode in chaincodes:
        borders[chaincode.objectID] = np.array(chaincode.Polygon()).tolist()
    print(">>> Labelled regions to approximate polygons conversion finished =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #generating geojson data to import in qupath
    GEOdata = []
    for label in borders:
        if (label not in df['cell_id'].values):
            continue
        cell_row = (df['cell_id'] == label)
        final_coords = borders[label]
        final_coords.append(final_coords[0])
        cell_data = {}
        cell_data["id"] = "cell" + str(label)
        cell_data["type"] = "Feature"
        cell_data["geometry"] = {
            "type" : "Polygon",
            "coordinates" : [final_coords]}
        cell_data["properties"] = {
            "name" : "cell" + str(label),
            "object_type" : "detection",
            "isLocked" : "false",
            }

        cell_data["properties"]["measurements"] = []
        for marker in markers:
            cell_data["properties"]["measurements"].append({
                "name" : marker,
                "value" : str(df[cell_row][marker].values[0])
                })

        #if cluster_id parameter is selected, add cluster_id and cluster_color
        if cluster_id != '' and len(df[cell_row][cluster_id]) > 0:
            cell_data["properties"]["classification"] = {
                "name": str(df[cell_row][cluster_id].values[0]),
                "colorRGB": str(df[cell_row][cluster_color].values[0] if cluster_color != '' else '')
                }
        GEOdata.append(cell_data)

    #dump GEOdata variable to json file
    with open(data_folder + '/analysis/cell_segmentation_geo.json', 'w') as outfile:
        geojson.dump(GEOdata, outfile)

    print(">>> End time generate_geojson =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
