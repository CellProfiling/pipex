import os
import sys
import fnmatch
import datetime
import cv2
from tifffile import TiffFile
from xml.etree import ElementTree


def log(msg, prefix=">>>"):
    print(prefix + " " + msg + " = " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)


def sanitize_marker_list(markers):
    """Strip whitespace, drop empty strings and duplicates, preserving order."""
    seen = []
    for m in markers:
        m = m.strip()
        if m and m not in seen:
            seen.append(m)
    return seen


def validate_marker_files(data_folder, markers):
    """Check each marker has at least one matching file in data_folder.
    Logs missing markers and exits if any are not found."""
    available = set()
    for file in os.listdir(data_folder):
        file_path = os.path.join(data_folder, file)
        if os.path.isdir(file_path):
            continue
        try:
            with TiffFile(file_path) as tif:
                if len(tif.series[0].pages) == 1:
                    for m in markers:
                        if fnmatch.fnmatch(file, '*' + m + '.*'):
                            available.add(m)
                else:
                    for page in tif.series[0].pages:
                        try:
                            biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                            if biomarker in markers:
                                available.add(biomarker)
                        except Exception:
                            pass
        except Exception:
            for m in markers:
                if fnmatch.fnmatch(file, '*' + m + '.*'):
                    available.add(m)
    missing = [m for m in markers if m not in available]
    if missing:
        log("ERROR: the following markers were not found as image files in " + data_folder + ": " + ", ".join(missing))
        sys.exit(1)


def validate_marker_columns(df, markers, source="cell_data.csv"):
    """Check each marker exists as a column in df.
    Logs missing markers and exits if any are not found."""
    missing = [m for m in markers if m not in df.columns]
    if missing:
        log("ERROR: the following markers were not found as columns in " + source + ": " + ", ".join(missing))
        sys.exit(1)


def iter_marker_images(data_folder, marker_names, callback):
    """
    Scan data_folder for images matching any name in marker_names and call
    callback(marker_name, raw_image_array) for each match.

    Handles three formats:
      - Single-channel TIFF: matched by filename glob '*<marker>.*'
      - Akoya qptiff multi-page: matched by the Biomarker XML tag on each page
      - Any other readable image (cv2 fallback): matched by filename glob

    For single-page files the loop breaks after the first matching marker so
    each file is only loaded once.
    """
    marker_set = set(marker_names)
    for file in os.listdir(data_folder):
        file_path = os.path.join(data_folder, file)
        if os.path.isdir(file_path):
            continue

        next_try = False
        try:
            with TiffFile(file_path) as tif:
                if len(tif.series[0].pages) == 1:
                    for marker in marker_names:
                        if fnmatch.fnmatch(file, '*' + marker + '.*'):
                            callback(marker, next(iter(tif.series[0].pages)).asarray())
                            break
                else:
                    for page in tif.series[0].pages:
                        biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                        if biomarker in marker_set:
                            callback(biomarker, page.asarray())
        except Exception as e:
            print('>>> checking type of ' + file_path + ', not QPTIFF', flush=True)
            print('>>> ', e, flush=True)
            next_try = True

        if next_try:
            try:
                for marker in marker_names:
                    if fnmatch.fnmatch(file, '*' + marker + '.*'):
                        curr_image = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
                        if len(curr_image.shape) > 2:
                            curr_image = curr_image[:, :, 0]
                        callback(marker, curr_image)
                        break
            except Exception as e:
                print('>>> Could not read image ' + file_path, flush=True)
                print('>>> ', e, flush=True)