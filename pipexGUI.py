import os
import sys
import datetime
import psutil
import FreeSimpleGUI as sg
import pipex


data_folder = os.path.join(os.getcwd(), "data")
if "PIPEX_DATA" in os.environ:
    data_folder = os.path.abspath(os.environ.get('PIPEX_DATA'))

pidfile_filename = './RUNNING'
if "PIPEX_WORK" in os.environ:
    pidfile_filename = './work/RUNNING'

try:
    with open(pidfile_filename,'r') as f:
        lines = f.readlines()
        if psutil.pid_exists(int(lines[0])):
            print(">>> Another PIPEX process seems to be running, exiting =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            sys.exit()
except IOError:
    pass

info_icon = b'iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAAC8HpUWHRSYXcgcHJvZmlsZSB0eXBlIGV4aWYAAHja7ZddkuwmDIXfWUWWYEkIieXwW5UdZPk5GLdzu3vuQydPqWooAwb5IPQBPRPGX3/O8AcSZU8hqnnKKR1IMcfMBQ0/dipnSUc8yzPFawjvT/3hHmB0CWrZr54u+0c/3QK7KmjpL0LeroH6PJCvGdhfhHhXsjxa7X4J5UtIeA/QJVD2so6U3X5dQh277o+V+H7CKqI/u/32boheV8wjzENIDpQsvh2Q9XCQgobsEoYkEW2YomRJlxgC8lOc7pTh0RwXinejJyp3i37uD6+0Il8m8hLkdNc/9gfSlwG55+Gn/eNXi5/7TbZUOF6iv545u89zzVhFiQmhTteiHks5W7CrmGJN7QF66TA8Cgk7c0Z27OqGrdCPdlTkRpkYDCZF6lRo0jjrRg0uRh6BDQ3mxnJ2uhhnbrL5IdNkkyxdHBQbsAt6+faFzmnz0cI5m2PmTjBlghitffFpDp9+MOc6CkSH37GCX8wr2HBjkVslzECE5hVUPQP8yK9pcRUQ1BXldUQyAlu3RFX65yaQE7TAUFHvM0jWLwGECFMrnCEBAVAjUUp0GLMRIZAOQAWus0SuIECq3OEkR5EENs5ranxidJqyMroD+nGZgYRKEgObLAWwYlTsH4uOPVRUNKpqUlPXrCVJiklTSpbWpVhMLAZTS2bmlq24eHT15Obu2UvmLLg0Nads2XPOpWDOAuWCrwsMSqlcpcaqoaZq1WuupWH7tNi0pWbNW26lc5eO+6Onbt177mXQwFYacehIw4aPPMrEVpsSZpw607TpM89yU7uwvuUPqNFFjU9Sy9Buaug1e0jQuk50MQMwDpFA3BYCbGhezA6nGHmRW8yOzDgVynBSF7NOixgIxkGskx7sAm+ii9x/4hYsPnHjf0suLHQfknvn9hO1vn6G2klsn8IV1ENw+jA+vLCX9WP3VoffDXxaf4W+Ql+hr9BX6Cv0FfofCU388bD+C/wbJaGnjctq5JEAAA0aaVRYdFhNTDpjb20uYWRvYmUueG1wAAAAAAA8P3hwYWNrZXQgYmVnaW49Iu+7vyIgaWQ9Ilc1TTBNcENlaGlIenJlU3pOVGN6a2M5ZCI/Pgo8eDp4bXBtZXRhIHhtbG5zOng9ImFkb2JlOm5zOm1ldGEvIiB4OnhtcHRrPSJYTVAgQ29yZSA0LjQuMC1FeGl2MiI+CiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogIDxyZGY6RGVzY3JpcHRpb24gcmRmOmFib3V0PSIiCiAgICB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIKICAgIHhtbG5zOnN0RXZ0PSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvc1R5cGUvUmVzb3VyY2VFdmVudCMiCiAgICB4bWxuczpkYz0iaHR0cDovL3B1cmwub3JnL2RjL2VsZW1lbnRzLzEuMS8iCiAgICB4bWxuczpHSU1QPSJodHRwOi8vd3d3LmdpbXAub3JnL3htcC8iCiAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyIKICAgIHhtbG5zOnhtcD0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wLyIKICAgeG1wTU06RG9jdW1lbnRJRD0iZ2ltcDpkb2NpZDpnaW1wOjZiYTJkMjY3LTViNmQtNGQ1MS1hMzI3LWE3YzY5MTdjZTQ4NCIKICAgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDo4ZDQ3YzcxYy04ZGU0LTRkYmQtODQ5MS00NDg0OGNjNGQwMDQiCiAgIHhtcE1NOk9yaWdpbmFsRG9jdW1lbnRJRD0ieG1wLmRpZDozYzEwOWNjMS1mODI2LTRmNjMtYWY2YS00ZGUyYzViNjBkMDgiCiAgIGRjOkZvcm1hdD0iaW1hZ2UvcG5nIgogICBHSU1QOkFQST0iMi4wIgogICBHSU1QOlBsYXRmb3JtPSJMaW51eCIKICAgR0lNUDpUaW1lU3RhbXA9IjE2NDYzMTMwMjYwNTkwNjUiCiAgIEdJTVA6VmVyc2lvbj0iMi4xMC4yOCIKICAgdGlmZjpPcmllbnRhdGlvbj0iMSIKICAgeG1wOkNyZWF0b3JUb29sPSJHSU1QIDIuMTAiPgogICA8eG1wTU06SGlzdG9yeT4KICAgIDxyZGY6QmFnPgogICAgIDxyZGY6bGkKICAgICAgc3RFdnQ6YWN0aW9uPSJzYXZlZCIKICAgICAgc3RFdnQ6Y2hhbmdlZD0iLyIKICAgICAgc3RFdnQ6aW5zdGFuY2VJRD0ieG1wLmlpZDo5ZjBmMGVlOS1lNjhjLTRhZTMtYTY4Mi1lOTU2YmRkNjRmNTUiCiAgICAgIHN0RXZ0OnNvZnR3YXJlQWdlbnQ9IkdpbXAgMi4xMCAoTGludXgpIgogICAgICBzdEV2dDp3aGVuPSIyMDIyLTAzLTAzVDE0OjEwOjI2KzAxOjAwIi8+CiAgICA8L3JkZjpCYWc+CiAgIDwveG1wTU06SGlzdG9yeT4KICA8L3JkZjpEZXNjcmlwdGlvbj4KIDwvcmRmOlJERj4KPC94OnhtcG1ldGE+CiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAKPD94cGFja2V0IGVuZD0idyI/Pr5MA1AAAAGFaUNDUElDQyBwcm9maWxlAAAokX2RPUjDQBiG36aKP7Q42EHEIUJ1siAq4qhVKEKFUCu06mBy6R80aUhSXBwF14KDP4tVBxdnXR1cBUHwB8TRyUnRRUr8Lim0iPGO4x7e+96Xu+8AoV5mmtUxDmi6baYScTGTXRW7XhGkGUYPhmVmGXOSlITv+LpHgO93MZ7lX/fnCKs5iwEBkXiWGaZNvEE8vWkbnPeJI6woq8TnxGMmXZD4keuKx2+cCy4LPDNiplPzxBFisdDGShuzoqkRTxFHVU2nfCHjscp5i7NWrrLmPfkLQzl9ZZnrtIaQwCKWIEGEgipKKMNGjHadFAspOo/7+Addv0QuhVwlMHIsoAINsusH/4PfvbXykxNeUigOdL44zscI0LULNGqO833sOI0TIPgMXOktf6UOzHySXmtp0SOgbxu4uG5pyh5wuQMMPBmyKbtSkJaQzwPvZ/RNWaD/Fuhd8/rWPMfpA5CmXiVvgINDYLRA2es+7+5u79u/Nc3+/QAHZnJ84H+0RgAAAAZiS0dEAAAAAAAA+UO7fwAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB+YDAw0KGlaJBaIAAAOPSURBVHja7ZvPS1VBFMc/6cvwB0Wgm1pFi9ASDPsBFSYEglA86PU3SOuE+gd0r0vXutKFJLxd9GOh8HBhlLqxyB6oC9uI5MJ82aJ5dBlm3tzmztwfvPuFi8zz3jPfOXfmnDNnzoUcOXLkyNG8OBNTP33AMDAA9APdwBWgVfy/BmwD+8A6sAZ8ADayrNxBYArYBU4tr11gGriVpYGXgEqEQeuuCvAszQO/C6x6GLh8rYq+UoMOMU0bkT4EXgOvgBHgOlAIyCiI30aAl+LeQ4PMadF3orgqDJWOZBkoAm0Wss+KZ8sN5G8KDolgCDjQEJsHeh321Stkqvo6EFxixSPgSEFmS7g7XxgWfcj9HglOsb151eDngM4Y+u8EZjVKGIpjzaum/XgCS/CFZjl4swkdGoM3ZiFnEqiK56uibWPRxxR8Nnx5h2kHb74DWNEYs2VL4uMaF+k8yJE7mbWQM2nw7ZOW/OYUspwGS6sKa29j8KoGBVQjGMYtRcToLLaXidq6ujChbhQXKcsquVBARRHk2MI0A75H5Dqv2EBF3tLKJKNEeCYbMOEgYpRlDrq0/GUHrlTnBVYcuS957zAVRZiczCg6iicmHMUBKhQVSRXrNJa8pW0j/WhTbKX7dDe3GKxqEG+B4wwo4Bh4ZxhLKAXclNrLGUrLyVwHdDcWGgi5IbU/WpKZEX9bgXbgEtDF34z0qfjtAtADnAPuOAhi1qR2v40CeqT2jiWZsQRmgMy1W3djo3OBE/7l7espqhMLMs+B38B54DHw0HC/ixlQAH4F2jXDyw4VtrrCJ0MwdNtRP6H4tyQwPb+myVo2UkDtP+xF2lAwjCWUAral9rUMKUDm+s1GAftS+3KGFCBz/WGjgHVDYJRmyFw/2yhADibuZUgB910Ecb42Q4ue3aCzzdAmsBdodwGjGXj7o4JrHXtiLFZxwEIKwtqooXeUFJ7zlJjvJeA8JQZuk6K+FeA8KQpu0+I+FTCMp7Q4uDsYAXhjUMADC5leD0ZAfTQ2F/LZVuGaLgJPgJ8GBSyKBEa72IKHgfejMbA/HJ0hWlGUCbEcjoL98bhPBcR6PA7pKpBQvXmvBRJ1pKFERrXmYymRqUNXJPWFJiiSCs6Epi2TC9qEJAslN0iwUDLoHcKUyi4RvlR2iYyUysrBUlMWS8t4ir9y+RIZQv2DiR2ifzAx6Itkmj+Zed8ok5MjR44cOXJExx+XlqahG0Iq8wAAAABJRU5ErkJggg=='

tooltip_1 = '<name before . in image file>: \nexample, from image filename \"reg001_cyc001_ch001_DAPI1.tif\"-> DAPI1'
tooltip_2 = '<number of pixels>: \nexample -> 20'
tooltip_3 = '<number of pixels, can be 0>: \nexample -> 20'
tooltip_4 = '<optional, name before . in image file>: \nexample, from image filename "reg001_cyc008_ch003_CDH1.tif" -> CDH1'
tooltip_5 = '<optional, number of pixels>: \nexample -> 25'
tooltip_6 = '<optional, "squareness" of the membrane, gradation between 0.001 and 0.999>: \nexample -> 0.9'
tooltip_7 = '<yes or no to enhance poor images>: \nexample -> yes'
tooltip_8 = '<list of markers names before . in image files>: \nexample -> DAPI1,CDH1,AMY2A,SST,GORASP2'
tooltip_9 = '<optional, one-side approximate resolution>: \nexample -> 1000'
tooltip_10 = '<optional, name of the column to add as cluster id information from cell_data.csv>: \nexample -> kmeans'
tooltip_11 = '<name of the column in cell_data.csv to filter the cells by>: \nexample -> cluster_id'
tooltip_12 = '<values, comma-separated, present in the selected colum of cell_data.csv to filter the cells by>: \nexample -> 3,6,7'
tooltip_13 = '<number of pixels>: \nexample -> 2048'
tooltip_14 = '<optional, number of pixels to be added around>: \nexample -> 128'
tooltip_15 = '<optional, percentage of tile size to be added around>: \nexample -> 10'
tooltip_16 = '<yes or no to have bigger border tiles>: \nexample -> no'
tooltip_17 = '<intensities below this percentage will be deleted>: \nexample -> 1'
tooltip_18 = '<intensities above this percentage will be deleted>: \nexample -> 99'
tooltip_19 = '<percentage of the image base intensity to apply>: \nexample -> 150'
tooltip_20 = '<number of pixels of a tile>: \nexample -> 1844'
tooltip_21 = '<number of pixels to use as smooth region for the tile cut>: \nexample -> 20'
tooltip_22 = '<yes or no to generate heat maps per image>: \nexample -> yes'
tooltip_23 = '<optional, "strictness" of stardist detections, gradation between 0.001 and 0.999>: \nexample -> 0.5'
tooltip_24 = '<number, main levels of intensity in the image [normally 3-5] and their 2 selected focus>: \nexample -> 4:1:3'
tooltip_25 = '<number, factor of complexity of light issues [1 to 4 should be enough]>: \nexample -> 3'
tooltip_26 = '<yes or no to reduce artifacts or foldings to a main intensity level>: \nexample -> yes'
tooltip_27 = '<optional, list of present specific markers to analyze>: \nexample -> AMY2A,SST,GORASP2'
tooltip_28 = '<percentage of biggest cells to remove>: \nexample -> 5'
tooltip_29 = '<yes or no to perform leiden clustering>: \nexample -> yes'
tooltip_30 = '<yes or no to perform kmeans clustering>: \nexample -> yes'
tooltip_31 = '<yes or no to show elbow analysis for kmeans>: \nexample -> yes'
tooltip_32 = '<force k number of cluster in kmeans>: \nexample -> 10'
tooltip_33 = '<yes or no to keep segmented membranes without nuclei as cells>: \nexample -> no'
tooltip_34 = '<percentage of smallest cells to remove>: \nexample -> 5'
tooltip_35 = '<optional, yes or no to apply custom Cell Profiling lab\'s biomarkers filtering>: \nexample -> yes'
tooltip_36 = '<optional, yes or no to apply log1p normalization to the markers>: \nexample -> yes'
tooltip_37 = '<optional, yes or no to apply quantile normalization to the markers>: \nexample -> yes'
tooltip_38 = '<optional, name of the column in cell_data.csv to perform batch correction by>: \nexample -> batch_id'
tooltip_39 = '<optional, yes or no to apply 0 to 1 re-scale normalization>: \nexample -> yes'
tooltip_40 = '<optional, yes or no to use binarized columns for clustering>: \nexample -> no'
tooltip_41 = '<optional, list of present specific markers to inlcude>: \nexample -> AMY2A,SST,GORASP2'
tooltip_42 = '<optional, name of the column to add as cluster color information from cell_data.csv>: \nexample -> kmeans_color'
tooltip_43 = '<optional, "strictness" of stardist detections proximity, gradation between 0.001 and 0.999>: \nexample -> 0.5'
tooltip_44 = '<optional, maximum allowed pixel area for initial Stardist detections>: \nexample -> 1600'
tooltip_45 = '<optional, yes or no to relabel sequentially the tile segments>: \nexample -> yes'
tooltip_46 = '<optional, file path to a pre-made custom segmentation>`: \nexample -> -custom_segmentation=/data/custom_seg.npy'
tooltip_47 = '<optional, yes or no to refine the cluster results through cell_types.csv data>: \nexample -> yes'
tooltip_48 = '<optional, yes or no or list of present specific markers to display as image layers> : example -> DAPI,SST,GORASP2'
tooltip_49 = '<optional, yes or no to include cell segmentation as regions> : example -> yes'
tooltip_50 = '<optional, yes or no to compress geojson regions into pbf> : example -> yes'
tooltip_51 = '<optional, yes or no to export html page for sharing the TissUUmaps project on the web> : example -> yes'
tooltip_52 = '<optional, dilation of shapes in pixels>: \nexample -> 0'
tooltip_53 = '<optional, number of datapoints used for smoothing the shapes>: \nexample -> 15'
tooltip_54 = '<optional, "none"/"hilbert"/"greedy" optimization of cutting path between shapes> : example -> none'
tooltip_55 = '<optional, nearest neighbour heuristic distance for merging shapes>: \nexample -> 300'



sg.theme('LightBrown10')

column = [[sg.Text('PIPEX data folder:', font='any 12'), sg.In(default_text=data_folder, size=(50,1), key='-DATA_FOLDER-'), sg.FolderBrowse(initial_folder=data_folder)],
          [sg.Text('Choose the sequence of commands you want PIPEX to perform:', font='any 12')],
          [sg.Text('_'*85)],
          [sg.Checkbox('Preprocessing', font='any 12 bold', key='-PREPROCESS-', enable_events=True)],
          [sg.Text('  - Min. threshold:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20, disabled=True, key='-PREPROCESS_THRMIN-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_17)],
          [sg.Text('  - Max. threshold:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='100',s=20, disabled=True, key='-PREPROCESS_THRMAX-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_18)],
          [sg.Text('  - Balance tiles', pad=((20,0), (0,0))), sg.Checkbox('',key='-PREPROCESS_TILFIX-', disabled=True, enable_events=True)],
          [sg.Text('  - Tile size:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='1844',s=20,disabled=True, key='-PREPROCESS_TILSIZ-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_20)],
          [sg.Text('  - Bright levels:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='3',s=20,disabled=True, key='-PREPROCESS_TILOTS-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_24)],
          [sg.Text('  - Flatten spots', pad=((20,0), (0,0))), sg.Checkbox('',key='-PREPROCESS_TILFLA-', disabled=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_26)],
          [sg.Text('  - Light gradient:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='3',s=20,disabled=True, key='-PREPROCESS_TILKER-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_25)],
          [sg.Text('  - Stitch size:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='40',s=20,disabled=True, key='-PREPROCESS_TILSTI-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_21)],
          [sg.Text('  - Exposure:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='100',s=20, disabled=True, key='-PREPROCESS_EXPOSU-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_19)],
          [sg.Text('  - Generate heat maps', pad=((20,0), (0,0))), sg.Checkbox('',key='-PREPROCESS_HEAMAP-', disabled=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_22)],
          [sg.Text('_'*85)],
          [sg.Checkbox('Cell segmentation', font='any 12 bold', key='-SEGMENTATION-', enable_events=True)],
          [sg.Text('  - NUCLEI marker:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='DAPI1',s=20,disabled=True, key='-SEGMENTATION_NUCMARK-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_1)],
          [sg.Text('  - NUCLEI diameter:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='20',s=20,disabled=True, key='-SEGMENTATION_NUCDIAM-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_2)],
          [sg.Text('  - NUCLEI expansion:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='20',s=20,disabled=True, key='-SEGMENTATION_NUCEXPA-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_3)],
          [sg.Text('  - NUCLEI definition:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-SEGMENTATION_NUCDEFI-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_23)],
          [sg.Text('  - NUCLEI closeness:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-SEGMENTATION_NUCCLOS-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_43)],
          [sg.Text('  - NUCLEI area limit:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-SEGMENTATION_NUCARLI-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_44)],
          [sg.Text('  - Use membrane marker', pad=((20,0), (0,0))), sg.Checkbox('',key='-SEGMENTATION_MEMUSE-', disabled=True, enable_events=True)],
          [sg.Text('  - MEMBRANE marker:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='HSP90B1',s=20,disabled=True, key='-SEGMENTATION_MEMMARK-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_4)],
          [sg.Text('  - MEMBRANE diameter:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='25',s=20,disabled=True, key='-SEGMENTATION_MEMDIAM-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_5)],
          [sg.Text('  - MEMBRANE compactness:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0.9',s=20,disabled=True, key='-SEGMENTATION_MEMCOMP-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_6)],
          [sg.Text('  - Keep membrane without nuclei:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',disabled=True, key='-SEGMENTATION_MEMKEEP-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_33)],
          [sg.Text('  - Adjust images:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',disabled=True, key='-SEGMENTATION_ADJUST-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_7)],
          [sg.Text('  - Custom segmentation:', s=35, pad=((20,0), (0,0))), sg.In(default_text="", size=(30,1), key='-SEGMENTATION_CUSSEG-'), sg.FileBrowse(initial_folder=data_folder), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_46)],
          [sg.Text('  - Measure markers, comma-separated:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='GORASP2,AMY2A',s=40,disabled=True, key='-SEGMENTATION_MEASURE-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_8)],
          [sg.Text('_'*85)],
          [sg.Checkbox('Downstream analysis', font='any 12 bold', key='-ANALYSIS-', enable_events=True)],
          [sg.Text(' NOTE: requires previous \'Segmentation\' results')],
          [sg.Text('  - Image size:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='1000',s=20,disabled=True, key='-ANALYSIS_SIZE-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_9)],
          [sg.Text('  - Analysis markers, comma-separated:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='',s=40,disabled=True, key='-ANALYSIS_MARKER-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_27)],
          [sg.Text('  - Use binarized markers', pad=((20,0), (0,0))), sg.Checkbox('',key='-ANALYSIS_USEBIN-', disabled=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_40)],
          [sg.Text('  - Cell size top crop:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-ANALYSIS_TOPTHR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_28)],
          [sg.Text('  - Cell size bottom crop:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-ANALYSIS_BOTTHR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_34)],
          [sg.Text('  - Custom Cell Profiling filtering:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',default=False,disabled=True, key='-ANALYSIS_CUSFIL-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_35)],
          [sg.Text('  - log1p normalization:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',default=False,disabled=True, key='-ANALYSIS_LOGNOR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_37)],
          [sg.Text('  - Standard normalization:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',default=True,disabled=True, key='-ANALYSIS_STDNOR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_39)],
          [sg.Text('  - Batch correction by column:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='',s=20,disabled=True, key='-ANALYSIS_BATCOR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_38)],
          [sg.Text('  - Quantile normalization:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',disabled=True, key='-ANALYSIS_QUANOR-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_36)],
          [sg.Text('  - Perform leiden cluster', pad=((20,0), (0,0))), sg.Checkbox('',key='-ANALYSIS_LEIDEN-', disabled=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_29)],
          [sg.Text('  - Perform Kmeans cluster', pad=((20,0), (0,0))), sg.Checkbox('',key='-ANALYSIS_KMEANS-', disabled=True, enable_events=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_30)],
          [sg.Text('  - K clusters:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='10',s=20,disabled=True, key='-ANALYSIS_KCLUST-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_32)],
          [sg.Text('  - Calculate elbow method', pad=((20,0), (0,0))), sg.Checkbox('',key='-ANALYSIS_ELBOW-', disabled=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_31)],
          [sg.Text('  - Refine clusters', pad=((20,0), (0,0))), sg.Checkbox('',key='-ANALYSIS_REFINE-', disabled=True, enable_events=True), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_47)],
          [sg.Text('_'*85)],
          [sg.Checkbox('QuPath GeoJSON', font='any 12 bold', key='-QUPATH-', enable_events=True)],
          [sg.Text(' NOTE: requires previous \'Segmentation\' results', pad=((20,0), (0,0)))],
          [sg.Text(' NOTE: requires a previous \'Downstream analysis\' results if you want clustering data', pad=((20,0), (0,0)))],
          [sg.Text('  - Included markers, comma-separated:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='',s=40,disabled=True, key='-QUPATH_MARKER-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_41)],
          [sg.Text('  - Cluster id column name:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='',s=20,disabled=True, key='-QUPATH_CLUID-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_10)],
          [sg.Text('  - Cluster color column name:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='',s=20,disabled=True, key='-QUPATH_CLUCOL-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_42)],
          [sg.Text('_'*85)],
          [sg.Checkbox('Filter segmentation', font='any 12 bold', key='-FILTERED-', enable_events=True)],
          [sg.Text(' NOTE: requires previous \'Segmentation\' results', pad=((20,0), (0,0)))],
          [sg.Text(' NOTE: requires a previous \'Downstream analysis\' results if you want cluster filtering', pad=((20,0), (0,0)))],
          [sg.Text('  - Perform cluster filtering', pad=((20,0), (0,0))), sg.Checkbox('',key='-FILTERED_CLUFIL-', disabled=True, enable_events=True)],
          [sg.Text('  - Cluster column name:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='cluster_id',s=20,disabled=True, key='-FILTERED_FIELD-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_11)],
          [sg.Text('  - Cluster id(s), comma-separated:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='3,4,8',s=40,disabled=True, key='-FILTERED_VALUE-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_12)],
          [sg.Text('  - LMD export', pad=((20,0), (0,0))), sg.Checkbox('',key='-FILTERED_LMD-', disabled=True, enable_events=True)],
          [sg.Text('  - Shape dilation:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-FILTERED_LMDDIL-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_52)],
          [sg.Text('  - Shape smoothing:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='15',s=20,disabled=True, key='-FILTERED_LMDSMO-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_53)],
          [sg.Text('  - Path optimization:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='none',s=20,disabled=True, key='-FILTERED_LMDPTH-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_54)],
          [sg.Text('  - Distance heuristic:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='300',s=20,disabled=True, key='-FILTERED_LMDDIS-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_55)],
          [sg.Text('  - Perform tiling', pad=((20,0), (0,0))), sg.Checkbox('',key='-FILTERED_TILING-', disabled=True, enable_events=True)],
          [sg.Text('  - Tile size:',s=35, pad=((20,0), (0,0))), sg.Input(default_text='2048',s=20,disabled=True, key='-FILTERED_TILSIZ-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_13)],
          [sg.Text('  - Tile overlap, in pixels:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-FILTERED_TILOVE-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_14)],
          [sg.Text('  - Tile overlap, in percentage:',s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='0',s=20,disabled=True, key='-FILTERED_TILPER-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_15)],
          [sg.Text('  - Tile relabel:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',disabled=True, key='-FILTERED_TILLAB-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_45)],
          [sg.Text('  - Extend tiles:',s=35, pad=((20,0), (0,0))), sg.Checkbox('',disabled=True, key='-FILTERED_TILEXT-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_16)],
          [sg.Text('_'*85)],
          [sg.Checkbox('Export to TissUUmaps', font='any 12 bold', key='-TISSUUMAPS-', enable_events=True)],
          [sg.Text(' NOTE: requires a previous \'Downstream analysis\' results', pad=((20,0), (0,0)))],
          [sg.Text(' NOTE: requires previous \'QuPath GeoJSON\' results if you want to include regions', pad=((20,0), (0,0)))],
          [sg.Text('  - Include marker images, comma-separated:', s=(35,1), pad=((20,0), (0,0))), sg.Input(default_text='',s=40,disabled=True, key='-TISSUUMAPS_MARKER-'), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_48)],
          [sg.Text('  - Include regions', pad=((20,0), (0,0))), sg.Checkbox('',key='-TISSUUMAPS_REGION-', disabled=True, enable_events=True, default=False), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_49)],
          [sg.Text('  - Compress regions', pad=((20,0), (0,0))), sg.Checkbox('',key='-TISSUUMAPS_COMPRESS_REGIONS-', disabled=True, enable_events=True, default=False), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_50)],
          [sg.Text('  - Include html', pad=((20,0), (0,0))), sg.Checkbox('',key='-TISSUUMAPS_HTML-', disabled=True, enable_events=True, default=False), sg.Image(data=info_icon,subsample=3,tooltip=tooltip_51)]]

layout = [[sg.Column(column, scrollable=True,  vertical_scroll_only=True, size=(800,700))],
          [sg.Text('_'*85)],
          [sg.Button('Run'), sg.Button('Batch mode'), sg.Button('Cancel')]]

window = sg.Window('PIPEX GUI', layout)

cancel = False
batch_mode = False

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Cancel':
        cancel = True
        break
    if event == 'Batch mode':
        batch_mode = True
        break
    if event == 'Run':
        break
    if event == '-PREPROCESS-':
        window['-PREPROCESS_THRMIN-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_THRMAX-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_EXPOSU-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILFIX-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILSIZ-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILOTS-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILFLA-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILKER-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_TILSTI-'].update(disabled=(not values['-PREPROCESS-']))
        window['-PREPROCESS_HEAMAP-'].update(disabled=(not values['-PREPROCESS-']))
        if (values['-PREPROCESS-']):
            window['-PREPROCESS_TILSIZ-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
            window['-PREPROCESS_TILOTS-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
            window['-PREPROCESS_TILFLA-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
            window['-PREPROCESS_TILKER-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
            window['-PREPROCESS_TILSTI-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
    if event == '-PREPROCESS_TILFIX-':
        window['-PREPROCESS_TILSIZ-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
        window['-PREPROCESS_TILOTS-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
        window['-PREPROCESS_TILFLA-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
        window['-PREPROCESS_TILKER-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
        window['-PREPROCESS_TILSTI-'].update(disabled=(not values['-PREPROCESS_TILFIX-']))
    if event == '-SEGMENTATION-':
        window['-SEGMENTATION_NUCMARK-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_NUCDIAM-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_NUCEXPA-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_NUCDEFI-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_NUCCLOS-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_NUCARLI-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEMUSE-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEMMARK-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEMDIAM-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEMCOMP-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEMKEEP-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_ADJUST-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_CUSSEG-'].update(disabled=(not values['-SEGMENTATION-']))
        window['-SEGMENTATION_MEASURE-'].update(disabled=(not values['-SEGMENTATION-']))
        if (values['-SEGMENTATION-']):
            window['-SEGMENTATION_MEMMARK-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
            window['-SEGMENTATION_MEMDIAM-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
            window['-SEGMENTATION_MEMCOMP-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
            window['-SEGMENTATION_MEMKEEP-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
    if event == '-SEGMENTATION_MEMUSE-':
        window['-SEGMENTATION_MEMMARK-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
        window['-SEGMENTATION_MEMDIAM-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
        window['-SEGMENTATION_MEMCOMP-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
        window['-SEGMENTATION_MEMKEEP-'].update(disabled=(not values['-SEGMENTATION_MEMUSE-']))
    if event == '-ANALYSIS-':
        window['-ANALYSIS_SIZE-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_MARKER-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_USEBIN-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_TOPTHR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_BOTTHR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_CUSFIL-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_LOGNOR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_STDNOR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_QUANOR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_BATCOR-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_LEIDEN-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_KMEANS-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_ELBOW-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_KCLUST-'].update(disabled=(not values['-ANALYSIS-']))
        window['-ANALYSIS_REFINE-'].update(disabled=(not values['-ANALYSIS-']))
        if (values['-ANALYSIS-']):
            window['-ANALYSIS_ELBOW-'].update(disabled=(not values['-ANALYSIS_KMEANS-']))
            window['-ANALYSIS_KCLUST-'].update(disabled=(not values['-ANALYSIS_KMEANS-']))
    if event == '-ANALYSIS_KMEANS-':
            window['-ANALYSIS_ELBOW-'].update(disabled=(not values['-ANALYSIS_KMEANS-']))
            window['-ANALYSIS_KCLUST-'].update(disabled=(not values['-ANALYSIS_KMEANS-']))
    if event == '-QUPATH-':
        window['-QUPATH_MARKER-'].update(disabled=(not values['-QUPATH-']))
        window['-QUPATH_CLUID-'].update(disabled=(not values['-QUPATH-']))
        window['-QUPATH_CLUCOL-'].update(disabled=(not values['-QUPATH-']))
    if event == '-FILTERED-':
        window['-FILTERED_CLUFIL-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_FIELD-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_LMD-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_LMDDIL-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_LMDSMO-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_LMDPTH-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_LMDDIS-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_VALUE-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILING-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILSIZ-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILOVE-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILPER-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILLAB-'].update(disabled=(not values['-FILTERED-']))
        window['-FILTERED_TILEXT-'].update(disabled=(not values['-FILTERED-']))
        if (values['-FILTERED-']):
            window['-FILTERED_FIELD-'].update(disabled=(not values['-FILTERED_CLUFIL-']))
            window['-FILTERED_VALUE-'].update(disabled=(not values['-FILTERED_CLUFIL-']))
            window['-FILTERED_LMDDIL-'].update(disabled=(not values['-FILTERED_LMD-']))
            window['-FILTERED_LMDSMO-'].update(disabled=(not values['-FILTERED_LMD-']))
            window['-FILTERED_LMDPTH-'].update(disabled=(not values['-FILTERED_LMD-']))
            window['-FILTERED_LMDDIS-'].update(disabled=(not values['-FILTERED_LMD-']))
            window['-FILTERED_TILSIZ-'].update(disabled=(not values['-FILTERED_TILING-']))
            window['-FILTERED_TILOVE-'].update(disabled=(not values['-FILTERED_TILING-']))
            window['-FILTERED_TILPER-'].update(disabled=(not values['-FILTERED_TILING-']))
            window['-FILTERED_TILLAB-'].update(disabled=(not values['-FILTERED_TILING-']))
            window['-FILTERED_TILEXT-'].update(disabled=(not values['-FILTERED_TILING-']))
    if event == '-FILTERED_CLUFIL-':
        window['-FILTERED_FIELD-'].update(disabled=(not values['-FILTERED_CLUFIL-']))
        window['-FILTERED_VALUE-'].update(disabled=(not values['-FILTERED_CLUFIL-']))
    if event == '-FILTERED_LMD-':
        window['-FILTERED_LMDDIL-'].update(disabled=(not values['-FILTERED_LMD-']))
        window['-FILTERED_LMDSMO-'].update(disabled=(not values['-FILTERED_LMD-']))
        window['-FILTERED_LMDPTH-'].update(disabled=(not values['-FILTERED_LMD-']))
        window['-FILTERED_LMDDIS-'].update(disabled=(not values['-FILTERED_LMD-']))
    if event == '-FILTERED_TILING-':
        window['-FILTERED_TILSIZ-'].update(disabled=(not values['-FILTERED_TILING-']))
        window['-FILTERED_TILOVE-'].update(disabled=(not values['-FILTERED_TILING-']))
        window['-FILTERED_TILPER-'].update(disabled=(not values['-FILTERED_TILING-']))
        window['-FILTERED_TILLAB-'].update(disabled=(not values['-FILTERED_TILING-']))
        window['-FILTERED_TILEXT-'].update(disabled=(not values['-FILTERED_TILING-']))
    if event == '-TISSUUMAPS-':
        window['-TISSUUMAPS_MARKER-'].update(disabled=(not values['-TISSUUMAPS-']))
        window['-TISSUUMAPS_REGION-'].update(disabled=(not values['-TISSUUMAPS-']))
        window['-TISSUUMAPS_COMPRESS_REGIONS-'].update(disabled=(not values['-TISSUUMAPS-']))
        window['-TISSUUMAPS_HTML-'].update(disabled=(not values['-TISSUUMAPS-']))

window.close()

if cancel:
    sys.exit()

if batch_mode:
    if "PIPEX_WORK" in os.environ:
        os.system("sudo docker-compose up pipex")
    else:
        pipex.batch_processor()
    sys.exit()

batch_filename = './pipex_batch_list.txt'
batch_data = values['-DATA_FOLDER-']
custom_segmentation_file = values['-SEGMENTATION_CUSSEG-']
if "PIPEX_WORK" in os.environ:
    batch_filename = os.environ['PIPEX_WORK'] + '/pipex_batch_list.txt'
    batch_data = batch_data.replace(os.environ['PIPEX_WORK'], './work')
    custom_segmentation_file = custom_segmentation_file.replace(os.environ['PIPEX_WORK'], './work')

batch_list = ''
if values['-PREPROCESS-']:
    batch_list = (batch_list + '\n' +
        'preprocessing.py -data=' + batch_data +
        ' -threshold_min=' + values['-PREPROCESS_THRMIN-'] +
        ' -threshold_max=' + values['-PREPROCESS_THRMAX-'] +
        ' -exposure=' + values['-PREPROCESS_EXPOSU-'] +
        ' -heat_map=' + ('yes' if values['-PREPROCESS_HEAMAP-'] else 'no'))
    if (values['-PREPROCESS_TILFIX-']):
        batch_list = (batch_list +
            ' -tile_size=' + values['-PREPROCESS_TILSIZ-'] +
            ' -bright_levels=' + values['-PREPROCESS_TILOTS-'] +
            ' -flatten_spots=' + ('yes' if values['-PREPROCESS_TILFLA-'] else 'no') +
            ' -light_gradient=' + values['-PREPROCESS_TILKER-'] +
            ' -balance_tiles=yes' +
            ' -stitch_size=' + values['-PREPROCESS_TILSTI-'])

if values['-SEGMENTATION-']:
    batch_list = (batch_list + '\n' +
        'segmentation.py -data=' + batch_data +
        ' -nuclei_marker=' + values['-SEGMENTATION_NUCMARK-'] +
        ' -nuclei_diameter=' + values['-SEGMENTATION_NUCDIAM-'] +
        ' -nuclei_expansion=' + values['-SEGMENTATION_NUCEXPA-'] +
        ' -nuclei_definition=' + values['-SEGMENTATION_NUCDEFI-'] +
        ' -nuclei_closeness=' + values['-SEGMENTATION_NUCCLOS-'] +
        ' -nuclei_area_limit=' + values['-SEGMENTATION_NUCARLI-'] +
        ' -adjust_images=' + ('yes' if values['-SEGMENTATION_ADJUST-'] else 'no'))
    if (values['-SEGMENTATION_MEMUSE-']):
        batch_list = (batch_list +
            ' -membrane_marker=' + values['-SEGMENTATION_MEMMARK-'] +
            ' -membrane_diameter=' + values['-SEGMENTATION_MEMDIAM-'] +
            ' -membrane_compactness=' + values['-SEGMENTATION_MEMCOMP-'] +
            ' -membrane_keep=' + ('yes' if values['-SEGMENTATION_MEMKEEP-'] else 'no'))
    if (values['-SEGMENTATION_CUSSEG-'] != ''):
        batch_list = (batch_list +
            ' -custom_segmentation=' + custom_segmentation_file)
    if (values['-SEGMENTATION_MEASURE-'] != ''):
        batch_list = (batch_list +
            ' -measure_markers=' + values['-SEGMENTATION_MEASURE-'])

if values['-ANALYSIS-']:
    batch_list = (batch_list + '\n' +
        'analysis.py -data=' + batch_data +
        ' -image_size=' + values['-ANALYSIS_SIZE-'] +
        ' -cellsize_max=' + values['-ANALYSIS_TOPTHR-'] +
        ' -cellsize_min=' + values['-ANALYSIS_BOTTHR-'] +
        ' -custom_filter=' + ('yes' if values['-ANALYSIS_CUSFIL-'] else 'no') +
        ' -log_norm=' + ('yes' if values['-ANALYSIS_LOGNOR-'] else 'no') +
        ' -std_norm=' + ('yes' if values['-ANALYSIS_STDNOR-'] else 'no') +
        ' -quantile_norm=' + ('yes' if values['-ANALYSIS_QUANOR-'] else 'no') +
        ' -batch_corr=' + values['-ANALYSIS_BATCOR-'] +
        ' -use_bin=' + ('yes' if values['-ANALYSIS_USEBIN-'] else 'no'))
    if (values['-ANALYSIS_LEIDEN-']):
        batch_list = (batch_list +
            ' -leiden=' + ('yes' if values['-ANALYSIS_LEIDEN-'] else 'no'))
    if (values['-ANALYSIS_KMEANS-']):
        batch_list = (batch_list +
            ' -kmeans=' + ('yes' if values['-ANALYSIS_KMEANS-'] else 'no') +
            ' -k_clusters=' + values['-ANALYSIS_KCLUST-'] +
            ' -elbow=' + ('yes' if values['-ANALYSIS_ELBOW-'] else 'no'))
    if (values['-ANALYSIS_REFINE-']):
        batch_list = (batch_list +
            ' -refine_clusters=' + ('yes' if values['-ANALYSIS_REFINE-'] else 'no'))
    if (values['-ANALYSIS_MARKER-'] != ''):
        batch_list = (batch_list +
            ' -analysis_markers=' + values['-ANALYSIS_MARKER-'])

if values['-QUPATH-']:
    batch_list = (batch_list + '\n' +
        'generate_geojson.py -data=' + batch_data)
    if (values['-QUPATH_MARKER-'] != ''):
        batch_list = (batch_list +
                      ' -included_markers=' + values['-QUPATH_MARKER-'])
    if (values['-QUPATH_CLUID-'] != ''):
        batch_list = (batch_list +
                      ' -cluster_id=' + values['-QUPATH_CLUID-'])
    if (values['-QUPATH_CLUCOL-'] != ''):
        batch_list = (batch_list +
                      ' -cluster_color=' + values['-QUPATH_CLUCOL-'])

if values['-FILTERED-']:
    batch_list = (batch_list + '\n' +
        'generate_filtered_masks.py -data=' + batch_data)
    if (values['-FILTERED_CLUFIL-']):
        batch_list = (batch_list +
            ' -field=' + values['-FILTERED_FIELD-'] +
            ' -values=' + values['-FILTERED_VALUE-'])
    if (values['-FILTERED_LMD-']):
        batch_list = (batch_list +
            ' -lmd=yes' +
            ' -shape_dilation=' + values['-FILTERED_LMDDIL-'] +
            ' -convolution_smoothing=' + values['-FILTERED_LMDSMO-'] +
            ' -path_optimization=' + values['-FILTERED_LMDPTH-'] +
            ' -distance_heuristic=' + values['-FILTERED_LMDDIS-'])
    if (values['-FILTERED_TILING-']):
        batch_list = (batch_list +
            ' -tile_size=' + values['-FILTERED_TILSIZ-'] +
            ' -tile_overlap=' + values['-FILTERED_TILOVE-'] +
            ' -tile_percentage_overlap=' + values['-FILTERED_TILPER-'] +
            ' -tile_relabel=' + ('yes' if values['-FILTERED_TILLAB-'] else 'no') +
            ' -extend_tile=' + ('yes' if values['-FILTERED_TILEXT-'] else 'no'))

if values['-TISSUUMAPS-']:
    batch_list = (batch_list + '\n' +
        'generate_tissuumaps.py -data=' + batch_data +
        ' -include_marker_images=' + values['-TISSUUMAPS_MARKER-'] +
        ' -include_geojson=' + ('yes' if values['-TISSUUMAPS_REGION-'] else 'no') +
        ' -compress_geojson=' + ('yes' if values['-TISSUUMAPS_COMPRESS_REGIONS-'] else 'no') +
        ' -include_html=' + ('yes' if values['-TISSUUMAPS_HTML-'] else 'no'))

if (batch_list != ''):
    batch_list = '#Auto-generated by PIPEX GUI' + batch_list + '\n'
    with open(batch_filename,'w') as f:
        f.writelines(batch_list)

    if "PIPEX_WORK" in os.environ:
        os.system("sudo docker-compose up pipex")
    else:
        pipex.batch_processor()
