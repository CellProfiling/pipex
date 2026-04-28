import os
os.environ['PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION'] = 'python'
import scanpy as sc
import scipy
import fnmatch
import json
import copy
import datetime
import sys
import argparse
import threading
import webbrowser
import http.server
import socketserver
from skimage.measure import approximate_polygon
import numpy as np
from pipex_utils import log

data_folder = os.environ.get('PIPEX_DATA')
include_marker_images = "no"
include_geojson = "no"
compress_geojson = "no"
include_html = "no"
launch = "no"

def find_marker_file(folder, marker):
    """Return the filename of the first image in folder whose name ends with marker."""
    for fname in sorted(os.listdir(folder)):
        if not os.path.isdir(os.path.join(folder, fname)) and fnmatch.fnmatch(fname, f'*{marker}.*'):
            return fname
    return None


def exporting_tissuumaps ():
    # Check if the required files are present
    global include_geojson
    if include_marker_images == "no" and include_geojson == "yes":
        print(">>> Impossible to display geojson without a background image", flush=True)
        include_geojson = "no"
    if include_geojson == "yes":
        if not os.path.exists(os.path.join(data_folder, 'analysis', 'cell_segmentation_geo.json')):
            print(">>> Impossible to display geojson without a cell segmentation file", flush=True)
            include_geojson = "no"
    
    adata = sc.read_h5ad(os.path.join(data_folder, 'analysis/downstream/anndata.h5ad'))
    log("AnnData loaded")

    # Make sure that the X matrix is in the compressed sparse column (CSC) format (required by TissUUmaps)
    adata.X = scipy.sparse.csc_matrix(adata.X)

    # Add image layers and cell segmentation geoJSON file to the AnnData object:
    if include_marker_images == "yes":
        markers = adata.var_names
    elif include_marker_images == "no":
        markers = []
    else:
        markers = include_marker_images.split(",")
    if include_geojson == "yes" and include_marker_images != "no":
        if compress_geojson == "yes":
            import geobuf
            with open(os.path.join(data_folder, 'analysis', 'cell_segmentation_geo.json'), 'r') as geojson_file:
                geojson_data = json.load(geojson_file)
                for feature in geojson_data:
                    if "measurements" in feature["properties"]:
                        del feature["properties"]["measurements"]
                geojson_data_approx = []
                for feature in geojson_data:
                    feature_approx = copy.deepcopy(feature)
                    for i in range(len(feature["geometry"]["coordinates"])):
                        feature_approx["geometry"]["coordinates"][i] = approximate_polygon(np.array(feature["geometry"]["coordinates"][i]), tolerance=0.75).tolist()
                    geojson_data_approx.append(feature_approx)
                geojson_data_approx = {
                    "type":"FeatureCollection",
                    "features":geojson_data_approx
                }
                pbf = geobuf.encode(geojson_data_approx, 3)
                data = geobuf.geobuf_pb2.Data()
                data.ParseFromString(pbf)
                data.precision = 3
                pbf = data.SerializeToString()
                cell_segmentation_path = "../cell_segmentation_geo.pbf"
                with open(os.path.join(data_folder, 'analysis', 'cell_segmentation_geo.pbf'), 'wb') as f:
                    f.write(pbf)
            log("GeoJSON compressed to pbf")
        else:
            cell_segmentation_path = "../cell_segmentation_geo.json"
        regionFiles = [
            {
                "name": "Cell segmentation",
                "path": cell_segmentation_path,
                "title": "Cell segmentation",
                "settings":[
                    {
                        "module": "regionUtils",
                        "function": "_regionStrokeWidth",
                        "value": "0.5"
                    },
                    {
                        "module": "regionUtils",
                        "function": "_regionStrokeAdaptOnZoom",
                        "value": True
                    },
                    {
                        "module": "glUtils",
                        "function": "_regionShowOnTop",
                        "value": False
                    }
                ]
            }
        ]
    else:
        regionFiles = []
    layers = []
    for marker in markers:
        fname = find_marker_file(data_folder, marker)
        if fname is None:
            print(f">>> Warning: no image found for marker {marker}, skipping layer", flush=True)
            continue
        layers.append({"name": marker, "tileSource": f"../../{fname}.dzi"})

    adata.uns["tmap"] = json.dumps({
        "layers": layers,
        "regionFiles": regionFiles,
        "plugins": ["Feature_Space","InteractionQC","Spot_Inspector"],
        "settings": [
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["Spot_Inspector",
                [
                    {"name": "_layer_format", "value":"{layout-row6}"},
                    {"name": "_cmap", "value":"undefined"},
                ],False]
            },
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["InteractionQC",[],False]
            },
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["Feature_Space",[],False]
            }
        ],
    })
    adata.write_h5ad(os.path.join(data_folder, 'analysis', 'downstream', 'anndata_TissUUmaps.h5ad'))
    log("AnnData TissUUmaps file saved")

    if include_html == "yes":
        import tissuumaps

        state = tissuumaps.read_h5ad.h5ad_to_tmap("", os.path.join(data_folder, 'analysis', 'downstream', 'anndata_TissUUmaps.h5ad'))
        tissuumaps.views.exportToStatic(
            json.dumps(state), 
            os.path.join(data_folder, 'analysis', 'downstream', 'TissUUmaps_webexport'),
            os.path.join(data_folder, 'analysis', 'downstream')
        )
        log("HTML static export saved")
        from urllib.request import urlretrieve
        for plugin in ["Feature_Space","InteractionQC","Spot_Inspector"]:
            url = f"https://tissuumaps.github.io/TissUUmaps/plugins/latest/{plugin}.js"
            filename = os.path.join(data_folder, 'analysis', 'downstream', 'TissUUmaps_webexport', 'plugins', f"{plugin}.js")
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            urlretrieve(url, filename)
        log("TissUUmaps plugins downloaded")

#Function to handle the command line parameters passed
def options(argv):
    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='generate_tissuumaps.py')
    parser.add_argument('--data', default=os.environ.get('PIPEX_DATA'),
        help='path to images folder : example -> -data=/lab/projectX/images')
    parser.add_argument('--include_marker_images', default='no',
        help='yes/no or comma-separated marker list for image layers : example -> -include_marker_images=DAPI,SST,GORASP2')
    parser.add_argument('--include_geojson', choices=['yes', 'no'], default='no',
        help='include cell segmentation as regions : example -> -include_geojson=yes')
    parser.add_argument('--compress_geojson', choices=['yes', 'no'], default='no',
        help='compress geojson regions into pbf : example -> -compress_geojson=yes')
    parser.add_argument('--include_html', choices=['yes', 'no'], default='no',
        help='export html page for web sharing : example -> -include_html=yes')
    parser.add_argument('--launch', choices=['yes', 'no'], default='no',
        help='launch local web server and open browser after export : example -> -launch=yes')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)

if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    include_marker_images = args.include_marker_images
    include_geojson = args.include_geojson
    compress_geojson = args.compress_geojson
    include_html = args.include_html
    launch = args.launch

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
    with open(os.path.join(data_folder, 'log_settings_tissuumaps.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time tissuumaps = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))

    log("Start time exporting tissuumaps")

    exporting_tissuumaps()

    if launch == "yes":
        webexport_path = os.path.join(data_folder, 'analysis', 'downstream', 'TissUUmaps_webexport')
        if not os.path.isdir(webexport_path):
            print(">>> WARNING: TissUUmaps webexport directory not found, cannot launch", flush=True)
        else:
            port = 8080
            while port < 8200:
                try:
                    handler = http.server.SimpleHTTPRequestHandler
                    handler.log_message = lambda *args: None
                    httpd = socketserver.TCPServer(("", port), handler)
                    break
                except OSError:
                    port += 1
            else:
                print(">>> WARNING: could not find a free port to launch TissUUmaps", flush=True)
                httpd = None
            if httpd:
                os.chdir(webexport_path)
                thread = threading.Thread(target=httpd.serve_forever, daemon=True)
                thread.start()
                url = f"http://localhost:{port}"
                print(f">>> TissUUmaps running at {url} — press Ctrl+C to stop", flush=True)
                webbrowser.open(url)
                try:
                    thread.join()
                except KeyboardInterrupt:
                    httpd.shutdown()
                    print(">>> TissUUmaps server stopped", flush=True)

    log("End time exporting tissuumaps")