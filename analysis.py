import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import datetime
import numpy as np
import PIL
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import qnorm
from combat.pycombat import pycombat
import json
import statistics


data_folder = os.environ.get('PIPEX_DATA')
image_size = 1000
analysis_markers = []
cellsize_max = 0
cellsize_min = 0
custom_filter = "no"
std_norm = "no"
log_norm = "no"
quantile_norm = "no"
batch_corr = ""
use_bin = "no"
leiden = "no"
kmeans = "no"
elbow = "no"
k_clusters = 10
refine_clusters = "no"
neigh_cluster_id = ""

max_samples = 200000


#Function to perform all data filtering, normalization and derived calculations
def data_calculations():
    #Reading the cell segmentation csv file
    df_norm = pd.read_csv(data_folder + '/analysis/cell_data.csv')

    markers = []
    #Getting the list of marker names
    markers = list(df_norm.columns.values)
    markers = markers[(df_norm.columns.get_loc("y") + 1):]
    #saveguard if analysis.py has been executed before and cluster_id + cluster_color already exists
    if any("_bin_thres" in s for s in markers):
        markers = markers[:-(len(df_norm.columns) - df_norm.columns.get_loc(list(filter(lambda x: '_bin_thres' in x, markers))[0]))]
    else:
        if 'leiden' in markers:
            markers = markers[:-(len(df_norm.columns) - df_norm.columns.get_loc("leiden"))]
        elif 'kmeans' in markers:
            markers = markers[:-(len(df_norm.columns) - df_norm.columns.get_loc("kmeans"))]

    #If a specific list of markers is informed, we use it
    if len(analysis_markers) > 0:
        markers = analysis_markers

    print(">>> List of markers to analyze gathered ",markers," =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    for marker in markers:
        df_norm[marker] = pd.to_numeric(df_norm[marker]).fillna(0)

    #We filter biggest and smallest cells
    filter_set = set()
    cellsize_df = df_norm['size']

    if cellsize_min > 0 or cellsize_max > 0:
        size_quantiles = cellsize_df.quantile([0.0000000001 + cellsize_min, 0.9999999999 - cellsize_max])

        filter_set.update(df_norm[df_norm['size'] > size_quantiles[0.9999999999 - cellsize_max]].index.values.tolist())
        filter_set.update(df_norm[df_norm['size'] < size_quantiles[0.0000000001 + cellsize_min]].index.values.tolist())

    if custom_filter == 'yes':
        if 'DAPI' in markers:
            filter_set.update(df_norm[df_norm['DAPI'] > df_norm['DAPI'].quantile(.99)].index.values.tolist())
        if 'CDH1' in markers:
            filter_set.update(df_norm[df_norm['CDH1'] > df_norm['CDH1'].quantile(.99)].index.values.tolist())
        if 'CTNNB1' in markers:
            filter_set.update(df_norm[df_norm['CTNNB1'] > df_norm['CTNNB1'].quantile(.99)].index.values.tolist())

    df_norm.drop(index=filter_set, axis=0, inplace=True)

    if cellsize_min > 0 or cellsize_max > 0 or custom_filter == 'yes':
        print(">>> Cells filtered =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #We normalize all the markers through min-max
    for marker in markers:
        if log_norm == 'yes':
            df_norm[marker] = np.log1p(df_norm[marker])
        if std_norm == 'yes':
            marker_min = df_norm[marker].min()
            marker_max = df_norm[marker].max()
            df_norm[marker] = df_norm[marker].apply(lambda x: (x - marker_min) / (marker_max - marker_min))

    if log_norm == 'yes' or std_norm == 'yes':
        print(">>> Markers normalized =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #Alternative normalization via z-score
    #for marker in markers:
    #    df_norm[marker] = pd.to_numeric(df[marker])
    #    if log_norm == 'yes':
    #        df_norm[marker] = np.log1p(df[marker])
    #    marker_mean = df_norm[marker].mean()
    #    marker_std = df_norm[marker].std()
    #    df_norm[marker] = df_norm[marker].apply(lambda x: (x - marker_mean) / marker_std)

    #ComBat batch correction
    if batch_corr != '':
        batch = []
        for batch_id in df_norm[batch_corr].unique():
            df_batch = df_norm[(df_norm[batch_corr] == batch_id)]
            batch.extend([batch_id for _ in range(len(df_batch))])

        df_norm[markers] = pycombat(df_norm[markers].transpose(), batch).transpose()

        print(">>> ComBat batch correction performed =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #Quantile normalization
    if quantile_norm == 'yes':
        df_norm[markers] = qnorm.quantile_normalize(df_norm[markers].transpose()).transpose()

        print(">>> Quantile normalization performed =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    df_norm.to_csv(data_folder + '/analysis/downstream/cell_data_norm.csv', index=False)

    #We calculate and plot the correlations between all markers
    df_corr = df_norm.copy()
    df_corr = df_corr[markers].corr()

    plt.figure()
    fig1, ax1 = plt.subplots(figsize=(7,5))
    sns_heatmap = sns.heatmap(df_corr, annot=True, annot_kws={"fontsize":5}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = False, linewidths=.1, cbar=False, ax=ax1)
    plt.savefig(data_folder + '/analysis/downstream/correlation_heatmap.jpg')
    plt.clf()
    plt.close()

    #We calculate and plot the dendogram of the correlations clustermap
    plt.figure()
    sns_clustermap = sns.clustermap(df_corr, figsize=(7,5))
    plt.savefig(data_folder + '/analysis/downstream/correlation_dendogram.jpg')
    plt.clf()
    plt.close()

    #We mark the zeroes as NaN to calculate mean, std, quartiles, etc... only over the cells that have the marker expressed
    #for marker in markers:
    #    df_norm.loc[marker] = np.NaN

    columns = ['marker', 'num_cells', 'percent_cells_50+', 'mean', 'median', 'std', 'q10', 'q25', 'q50', 'q75', 'q90'] + markers

    #We create an extra pandas dataframe with global information about each marker. This will be a new output in PIPEX analysis step as a csv file called 'cell_data_markers.csv'
    df_ext = pd.DataFrame(columns=columns)

    for marker in markers:
        marker_df = df_norm[marker]
        quantiles = marker_df.quantile([.10, .25, .50, .75, .90])

        x_tiles = int(df_norm['x'].max() / 1844 + 1 if df_norm['x'].max() % 1844 != 0 else 0)
        y_tiles = int(df_norm['y'].max() / 1844 + 1 if df_norm['y'].max() % 1844 != 0 else 0)
        qif = np.zeros((y_tiles, x_tiles))
        for i in range(len(qif)):
            for j in range(len(qif[0])):
                tile_min_x = j * 1844
                tile_min_y = i * 1844
                tile_max_x = (j + 1) * 1844
                tile_max_y = (i + 1) * 1844
                marker_tile = df_norm.loc[(df_norm['x'] > tile_min_x) & (df_norm['y'] > tile_min_y) & (df_norm['x'] <= tile_max_x) & (df_norm['y'] <= tile_max_y)]

                qif[i][j] = marker_tile[marker].sum() / (1844 ** 2)
        row = {'marker': marker,
               'num_cells': len(marker_df),
               'percent_cells_50+': (len(marker_df[(marker_df >= 0.5)]) / len(marker_df)) * 100,
               'mean': marker_df.mean(),
               'median': marker_df.median(),
               'std': marker_df.std(),
               'q10': quantiles[.10],
               'q25': quantiles[.25],
               'q50': quantiles[.50],
               'q75': quantiles[.75],
               'q90': quantiles[.90]}

        for i in range(len(qif)):
            for j in range(len(qif[0])):
                row['QIF_' + str(i) + '_' + str(j)] = qif[i][j]

        df_ext = pd.concat([df_ext, pd.DataFrame([row])], ignore_index=True)

    for marker in markers:
        df_ext[marker] = df_corr[marker].values

    #We generate a boxplot with each marker. Note: this shows calculations with cells expressing the marker and ignoring the other that don't
    plt.figure()
    fig2, ax2 = plt.subplots(figsize=(7,5))
    sns_boxplot = sns.boxplot(x = "variable", y = "value", data = pd.melt(df_norm[markers]), showfliers = False, color="skyblue")
    sns_boxplot.set_xticklabels(sns_boxplot.get_xticklabels(), rotation = 45)
    sns_boxplot.set(xlabel=None)
    sns_boxplot.set(ylabel=None)
    plt.savefig(data_folder + '/analysis/downstream/markers_boxplot.jpg')
    plt.clf()
    plt.close()

    df_ext.to_csv(data_folder + '/analysis/downstream/cell_data_markers.csv', index=False)
    print(">>> Markers information calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    del df_corr
    del df_ext

    fill_surface_html_template(markers, df_norm)
    print(">>> Markers surface plot generated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    return df_norm, markers


#Function pass data to the surface html template
def fill_surface_html_template(markers, df_norm):
    html_template = open("markers_surface.html", "r")
    html_content = html_template.read()
    markers_formatted = "["
    for marker in markers:
        markers_formatted = markers_formatted + '\"' + marker + '\",'
    markers_formatted = markers_formatted[:-1] + "]"
    html_content = html_content.replace("$$$MARKERS$$$", markers_formatted)
    minDim = min(min(df_norm['x'].values), min(df_norm['y'].values))
    maxDim = max(max(df_norm['x'].values), max(df_norm['y'].values))
    factor = int(maxDim / 100) + 1 if maxDim % 100 > 0 else 0
    ticktext = [*range(minDim, maxDim, factor * 20)]
    ticktext_formatted = "["
    for tick in ticktext:
        ticktext_formatted = ticktext_formatted + str(tick) + ','
    ticktext_formatted = ticktext_formatted[:-1] + "]"
    html_content = html_content.replace("$$$TICKS$$$", ticktext_formatted)
    df_surface = df_norm
    if std_norm != 'yes':
        df_surface = df_norm.copy()
        for marker in markers:
            marker_min = df_surface[marker].min()
            marker_max = df_surface[marker].max()
            df_surface[marker] = df_surface[marker].apply(lambda x: (x - marker_min) / (marker_max - marker_min))
    z = []
    for marker in markers:
        z.append([])
    for x in range(minDim, maxDim, factor):
        for m in range(len(markers)):
            z_row = []
            for y in range(minDim, maxDim, factor):
                df_cell = df_surface[(df_surface['x'] >= x) & (df_surface['y'] >= y) & (df_surface['x'] < x + factor) & (df_surface['y'] < y + factor)]
                if df_cell.empty:
                    z_row.append(0)
                else:
                    z_row.append(df_cell[markers[m]].max())
            z[m].append(z_row)
    z_formatted = "["
    for z_marker in z:
        z_formatted = z_formatted + '['
        for z_row in z_marker:
            z_formatted = z_formatted + '['
            for z_elem in z_row:
                z_formatted = z_formatted + str(z_elem) + ','
            z_formatted = z_formatted[:-1] + "],"
        z_formatted = z_formatted[:-1] + "],"
    z_formatted = z_formatted[:-1] + "]"
    html_content = html_content.replace("$$$DATA$$$", z_formatted)
    f = open(data_folder + '/analysis/downstream/markers_surface.html', 'w')
    f.write(html_content)
    f.close()


#Function to convert regular RGB to packed integer
def generate_leiden_color(leiden_id, leiden_color_list):
    try:
        rgb = PIL.ImageColor.getcolor(leiden_color_list[int(leiden_id)], "RGB")
        return (rgb[0] << 16) + (rgb[1] << 8) + rgb[2]
    except Exception as e:
        return ''


def check_cell_type_threshold(curr_marker, curr_rule, curr_score, high_threshold, low_threshold):
    if curr_rule == 'high':
        if curr_score >= high_threshold:
            return 100
        elif curr_score > low_threshold:
            return 50
    elif curr_rule == 'low':
        if curr_score <= low_threshold:
            return 100
        elif curr_score < high_threshold:
            return 50
    else:
        if curr_score > low_threshold and curr_score < high_threshold:
            return 100
        else:
            return 50
    return 0


def check_cell_type(row, cluster_id, clustering_merge_data):
    high_threshold = clustering_merge_data['scores'][cluster_id]['q75']
    low_threshold = clustering_merge_data['scores'][cluster_id]['q25']
    final_score = 0
    num_rules = 0
    if not pd.isnull(row['marker1']):
        curr_marker = row['marker1']
        curr_rule = row['rule1']
        num_rules = num_rules + 1
        if curr_marker in clustering_merge_data['scores'][cluster_id]['markers']:
            curr_score = clustering_merge_data['scores'][cluster_id]['markers'][curr_marker]
            final_score = check_cell_type_threshold(curr_marker, curr_rule, curr_score, high_threshold, low_threshold)
    if not pd.isnull(row['marker2']):
        curr_marker = row['marker2']
        curr_rule = row['rule2']
        num_rules = num_rules + 1
        if curr_marker in clustering_merge_data['scores'][cluster_id]['markers']:
            curr_score = clustering_merge_data['scores'][cluster_id]['markers'][curr_marker]
            final_score = final_score + check_cell_type_threshold(curr_marker, curr_rule, curr_score, high_threshold, low_threshold)
    if not pd.isnull(row['marker3']):
        curr_marker = row['marker3']
        curr_rule = row['rule3']
        num_rules = num_rules + 1
        if curr_marker in clustering_merge_data['scores'][cluster_id]['markers']:
            curr_score = clustering_merge_data['scores'][cluster_id]['markers'][curr_marker]
            final_score = final_score + check_cell_type_threshold(curr_marker, curr_rule, curr_score, high_threshold, low_threshold)
    if final_score > 0:
        return final_score / num_rules

    return None


def calculate_cluster_info(adata, cluster_type):
    plt.figure()
    sc.pl.umap(adata, color=[cluster_type], show=False, save='_' + cluster_type)
    plt.clf()
    plt.close()

    plt.figure()
    sc.pl.spatial(adata, color=cluster_type, spot_size=20, show=False, save='_spatial_' + cluster_type)
    plt.clf()
    plt.close()

    try:
        sq.gr.spatial_neighbors(adata, coord_type="generic")
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_type)
        plt.figure()
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_type, method="single", show=False,
                               save='nhood_enrichment_' + cluster_type + '.jpg')
        plt.clf()
        plt.close()
    except Exception as e:
        print(e)
        print('>>> Neighborhood calculations failed for cluster ' + cluster_type, flush=True)

    try:
        sq.gr.interaction_matrix(adata, cluster_key=cluster_type)
        plt.figure()
        sq.pl.interaction_matrix(adata, cluster_key=cluster_type, show=False, save='interaction_matrix_' + cluster_type + '.jpg')
        plt.clf()
        plt.close()
    except Exception as e:
        print(e)
        print('>>> Interaction matrix analysis failed for cluster ' + cluster_type, flush=True)

    try:
        sc.tl.rank_genes_groups(adata, cluster_type, method='t-test')
        sc.settings.set_figure_params(format='jpg', figsize=(image_size / 200, image_size / 200))
        plt.figure()
        sc.pl.rank_genes_groups(adata, n_genes=len(markers), sharey=False, show=False, save='')
        plt.clf()
        plt.close()
    except Exception as e:
        print(e)
        print('>>> Rank genes groups analysis failed for cluster ' + cluster_type, flush=True)

    sc.settings.set_figure_params(format='jpg', figsize=(image_size / 100, image_size / 100))
    plt.figure()
    sc.pl.heatmap(adata, markers, groupby=cluster_type, swap_axes=True, cmap='viridis', dendrogram=False, show=False,
                  save='_' + cluster_type)
    plt.clf()
    plt.close()


def refine_clustering(adata, cluster_type):
    clustering_merge_data = {}
    clustering_merge_data['scores'] = {}
    clustering_merge_data['cell_types'] = {}
    cluster_dif_list = []
    for cluster_id in adata.obs[cluster_type].unique():
        cluster_score_list = []
        cluster_merge_clusters_scores = {}
        cluster_merge_clusters_scores['markers'] = {}
        rank_df = sc.get.rank_genes_groups_df(adata, group=cluster_id)
        for marker_id in rank_df['names'].unique():
            curr_score = rank_df[rank_df['names'] == marker_id]['scores'].values[0]
            cluster_merge_clusters_scores['markers'][marker_id] = float(curr_score)
            cluster_score_list.append(curr_score)
        cluster_merge_clusters_scores['score_max'] = float(max(cluster_score_list))
        cluster_merge_clusters_scores['score_min'] = float(min(cluster_score_list))
        cluster_merge_clusters_scores['score_dif'] = cluster_merge_clusters_scores['score_max'] - cluster_merge_clusters_scores['score_min']
        cluster_dif_list.append(cluster_merge_clusters_scores['score_dif'])
        cluster_merge_clusters_scores['q75'] = float((cluster_merge_clusters_scores['score_max'] - cluster_merge_clusters_scores['score_min']) * 75 / 100 + cluster_merge_clusters_scores['score_min'])
        cluster_merge_clusters_scores['q25'] = float((cluster_merge_clusters_scores['score_max'] - cluster_merge_clusters_scores['score_min']) * 25 / 100 + cluster_merge_clusters_scores['score_min'])
        clustering_merge_data['scores'][cluster_id] = cluster_merge_clusters_scores
        clustering_merge_data['cell_types'][cluster_id] = []
    clustering_merge_data['dif_median'] = float(statistics.median(cluster_dif_list))

    cell_types = pd.read_csv(data_folder + '/cell_types.csv')
    for index, row in cell_types.iterrows():
        for cluster_id in clustering_merge_data['scores']:
            cell_type_prob = check_cell_type(row, cluster_id, clustering_merge_data)
            if cell_type_prob is not None:
                if clustering_merge_data['scores'][cluster_id]['score_dif'] < clustering_merge_data['dif_median']:
                    cell_type_prob = cell_type_prob * clustering_merge_data['scores'][cluster_id]['score_dif'] / clustering_merge_data['dif_median']
                clustering_merge_data['cell_types'][cluster_id].append({ 'cell_type' : row['cell_group'] + '.' + row['cell_type'] + '.' + row['cell_subtype'], 'prob' : cell_type_prob})
    clustering_merge_data['candidates'] = {}
    adata.obs[cluster_type + "_ref"] = adata.obs[cluster_type].astype(str)
    adata.obs[cluster_type + "_ref_p"] = adata.obs[cluster_type].astype(str)
    ordered_cluster_keys = list(clustering_merge_data['cell_types'])
    ordered_cluster_keys.sort()
    for cluster_id in ordered_cluster_keys:
        best_candidate = None
        for curr_cell_type in clustering_merge_data['cell_types'][cluster_id]:
            curr_cell_type['prob'] = curr_cell_type['prob'] / len(clustering_merge_data['cell_types'][cluster_id])
            if best_candidate is None or best_candidate['prob'] < curr_cell_type['prob']:
                best_candidate = { 'cell_type': curr_cell_type['cell_type'], 'prob' : curr_cell_type['prob'], 'real_confidence' : '{:.1%}'.format((curr_cell_type['prob'] * len(clustering_merge_data['cell_types'][cluster_id])) / 100.0) }
        clustering_merge_data['candidates'][cluster_id] = best_candidate
        adata.obs.loc[adata.obs[cluster_type + "_ref"] == cluster_id, cluster_type + "_ref"] = best_candidate['cell_type']
        adata.obs.loc[adata.obs[cluster_type + "_ref_p"] == cluster_id, cluster_type + "_ref_p"] = best_candidate['real_confidence'][:-1]

    with open(data_folder + '/analysis/downstream/cell_types_result_' + cluster_type + '.json', 'w') as outfile:
        json.dump(clustering_merge_data, outfile, indent = 4)


#Function to perform different cluster methods
def clustering(df_norm, markers):
    #We create an anndata object with the data so we can proceed with the clustering analysis using scanpy and squidpy libraries
    adata = sc.AnnData(df_norm.loc[:, markers])
    adata.obs_names = 'cell_id_' + df_norm['cell_id'].astype(str)
    adata.var_names = markers

    adata.obs["id"] = np.array(df_norm.loc[:, 'cell_id'])
    adata.obs["size"] = np.array(df_norm.loc[:, 'size'])
    adata.obs["x"] = np.array(df_norm.loc[:, 'x'])
    adata.obs["y"] = np.array(df_norm.loc[:, 'y'])

    adata.obsm["spatial"] = np.array(df_norm.loc[:, ['x', 'y']])
    print(">>> Anndata object created =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #We take the chance to show spatially the intensities of every marker
    if batch_corr == '' and len(df_norm.index) <= max_samples:
        for marker in markers:
            plt.figure()
            sc.pl.spatial(adata, color=marker, cmap='viridis', spot_size=20, show=False, save='_spatial_' + marker)
            plt.clf()
            plt.close()
    else:
        print(">>> Dataset too big to create spatial plots per marker =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    #We calculate PCA, neighbors and UMAP for the anndata
    sc.pp.pca(adata)
    adata.obsm['X_pca'] = np.nan_to_num(adata.obsm['X_pca'], copy=False)
    print(">>> PCA calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    num_neighbors = int(max(5, 15 * min(1, max_samples / len(df_norm.index))))
    print(">>> n_neighbors set to",num_neighbors,"=", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    sc.pp.neighbors(adata, n_neighbors=num_neighbors)
    print(">>> Neighbors graph calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    sc.tl.umap(adata)
    print(">>> UMAP calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    plt.figure()
    sc.pl.umap(adata, show=False, save='_base')
    plt.clf()
    plt.close()

    #standard_embedding = umap.UMAP(random_state=0,low_memory=True).fit_transform(adata.obsm['X_pca'])
    #plt.scatter(standard_embedding[:, 0], standard_embedding[:, 1], s=0.1, cmap='Spectral');
    #plt.savefig('umaptest.jpg')
    #sys.exit(0)

    if (leiden == 'yes'):
        #We calculate leiden cluster
        sc.tl.leiden(adata)
        print(">>> Leiden cluster calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

        #We print the complete leiden cluster and all related information
        calculate_cluster_info(adata, "leiden")

        if refine_clusters == "yes":
            try:
                refine_clustering(adata, 'leiden')
                calculate_cluster_info(adata, "leiden_ref")
            except Exception as e:
                print(">>> Failed at refining leiden cluster =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    if (kmeans == 'yes'):
        if (elbow == 'yes'):
            if (len(df_norm.index) <= max_samples):
                #We calculate all kmeans clusters with k 1 to 20 so we can show the elbow method plots
                print(">>> Performing kmeans elbow method =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
                distortions = []
                inertias = []
                mapping1 = {}
                mapping2 = {}
                K = range(1, 20)

                for k in K:
                    kmeanModel = KMeans(n_clusters=k).fit(adata.obsm['X_pca'])
                    kmeanModel.fit(adata.obsm['X_pca'])

                    distortions.append(sum(np.min(cdist(adata.obsm['X_pca'], kmeanModel.cluster_centers_,
                                        'euclidean'), axis=1)) / adata.obsm['X_pca'].shape[0])
                    inertias.append(kmeanModel.inertia_)

                    mapping1[k] = sum(np.min(cdist(adata.obsm['X_pca'], kmeanModel.cluster_centers_,
                                   'euclidean'), axis=1)) / adata.obsm['X_pca'].shape[0]
                    mapping2[k] = kmeanModel.inertia_
                    print(">>> Kmeans cluster calculated with k",k," =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

                plt.figure()
                plt.plot(K, distortions, 'bx-')
                plt.xlabel('Values of K')
                plt.ylabel('Distortion')
                plt.title('The Elbow Method using Distortion')
                plt.savefig(data_folder + '/analysis/downstream/elbow_distortion.jpg')
                plt.clf()
                plt.close()

                plt.figure()
                plt.plot(K, inertias, 'bx-')
                plt.xlabel('Values of K')
                plt.ylabel('Inertia')
                plt.title('The Elbow Method using Inertia')
                plt.savefig(data_folder + '/analysis/downstream/elbow_inertia.jpg')
                plt.clf()
                plt.close()
                print(">>> Elbow method calculated =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
            else:
                print(">>> Dataset too big to perform elbow method =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

        #We print the complete spatial kmeans cluster and all related information. Please note that the k used is the one passad as parameter (or 10 by default)
        kmeans_cluster = KMeans(n_clusters=k_clusters, random_state=0).fit(adata.obsm['X_pca'])
        adata.obs['kmeans'] = kmeans_cluster.labels_.astype(str)

        #We print the complete kmeans cluster and all related information
        calculate_cluster_info(adata, "kmeans")

        if refine_clusters == "yes":
            try:
                refine_clustering(adata, 'kmeans')
                calculate_cluster_info(adata, "kmeans_ref")
            except Exception as e:
                print(">>> Failed at refining kmeans cluster =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    if (leiden == 'yes' or kmeans == 'yes'):
        df = pd.read_csv(data_folder + '/analysis/cell_data.csv')
        if (leiden == 'yes'):
            df['leiden'] = df['cell_id'].map(adata.obs.set_index('id')['leiden']).astype(str)
            df['leiden'] = df['leiden'].fillna('')
            if refine_clusters == "yes":
                df['leiden_ref'] = df['cell_id'].map(adata.obs.set_index('id')['leiden_ref']).astype(str)
                df['leiden_ref'] = df['leiden_ref'].fillna('')
                df['leiden_ref_p'] = df['cell_id'].map(adata.obs.set_index('id')['leiden_ref_p']).astype(str)
                df['leiden_ref_p'] = df['leiden_ref_p'].fillna('')
            leiden_color_list = adata.uns['leiden_colors']
            df['leiden_color'] = df.apply(lambda row: generate_leiden_color(row['leiden'], leiden_color_list), axis=1)
            df_norm['leiden'] = df_norm['cell_id'].map(df.set_index('cell_id')['leiden']).astype(str)
            df_norm['leiden_color'] = df_norm['cell_id'].map(df.set_index('cell_id')['leiden_color']).astype(str)
        if (kmeans == 'yes'):
            #We add to the original cell segmentation csv file the calculated kmeans group for each cell
            df['kmeans'] = df['cell_id'].map(adata.obs.set_index('id')['kmeans']).astype(str)
            df['kmeans'] = df['kmeans'].fillna('')
            if refine_clusters == "yes":
                df['kmeans_ref'] = df['cell_id'].map(adata.obs.set_index('id')['kmeans_ref']).astype(str)
                df['kmeans_ref'] = df['kmeans_ref'].fillna('')
                df['kmeans_ref_p'] = df['cell_id'].map(adata.obs.set_index('id')['kmeans_ref_p']).astype(str)
                df['kmeans_ref_p'] = df['kmeans_ref_p'].fillna('')
            kmeans_color_list = adata.uns['kmeans_colors']
            df['kmeans_color'] = df.apply(lambda row: generate_leiden_color(row['kmeans'], kmeans_color_list), axis=1)
            df_norm['kmeans'] = df_norm['cell_id'].map(df.set_index('cell_id')['kmeans']).astype(str)
            df_norm['kmeans_color'] = df_norm['cell_id'].map(df.set_index('cell_id')['kmeans_color']).astype(str)

            df_corr = pd.concat([df_norm[markers], pd.get_dummies(df_norm['kmeans'], prefix='clusterK')], axis=1).corr()

            plt.figure()
            fig1, ax1 = plt.subplots(figsize=(7,5))
            sns_heatmap = sns.heatmap(df_corr, annot=True, annot_kws={"fontsize":5}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = False, linewidths=.1, cbar=False, ax=ax1)
            plt.savefig(data_folder + '/analysis/downstream/kmeans_clusters_correlation_heatmap.jpg')
            plt.clf()
            plt.close()

        df.to_csv(data_folder + '/analysis/cell_data.csv', index=False)
        df_norm.to_csv(data_folder + '/analysis/downstream/cell_data_norm.csv', index=False)
        adata.write(data_folder + '/analysis/downstream/anndata.h5ad')

#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
       print('analysis.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-image_size=<optional, one-side approximate resolution> : example -> -image_size=1000\n\t-analysis_markers=<optional, list of present specific markers to analyze> : example -> -analysis_markers=AMY2A,SST,GORASP2\n\t-cellsize_max=<optional, percentage of biggest cells to remove> : example -> -cellsize_max=5\n\t-cellsize_min=<optional, percentage of smallest cells to remove> : example -> -cellsize_min=5\n\t-custom_filter=<yes or no to apply custom Cell Profiling lab\'s biomarkers filtering> : example -> -custom_filter=yes\n\t-log_norm=<yes or no to apply log n + 1 normalization> : example -> -log_norm=yes\n\t-std_norm=<yes or no to apply 0 to 1 re-scale normalization> : example -> -std_norm=yes\n\t-quantile_norm=<yes or no to apply quantile normalization> : example -> -quantile_norm=yes\n\t-batch_corr=<optional, name of the column in cell_data.csv to perform batch correction by> : example -> batch_id\n\t-use_bin=<optional, yes or no to use binarized columns for clustering> : example -> -use_bin=no\n\t-leiden=<optional, yes or no to perform leiden clustering> : example -> -leiden=yes\n\t-kmeans=<optional, yes or no to perform kmeans clustering> : example -> -kmeans=yes\n\t-elbow=<optional, yes or no to show elbow analysis for kmeans> : example -> -elbow=yes\n\t-k_clusters=<optional, force k number of cluster in kmeans> : example -> -k_clusters=10\n\t-refine_clusters=<optional, yes or no to refine cluster results> : example -> -refine_clusters=yes\n\t-neigh_cluster_id=<optional, cluster column name to base the neighborhood analysis on> : example -> -neigh_cluster_id=kmeans_ref', flush=True)
       sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print('analysis.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-image_size=<optional, one-side approximate resolution> : example -> -image_size=1000\n\t-analysis_markers=<optional, list of present specific markers to analyze> : example -> -analysis_markers=AMY2A,SST,GORASP2\n\t-cellsize_max=<optional, percentage of biggest cells to remove> : example -> -cellsize_max=5\n\t-cellsize_min=<optional, percentage of smallest cells to remove> : example -> -cellsize_min=5\n\t-custom_filter=<yes or no to apply custom Cell Profiling lab\'s biomarkers filtering> : example -> -custom_filter=yes\n\t-log_norm=<yes or no to apply log n + 1 normalization> : example -> -log_norm=yes\n\t-std_norm=<yes or no to apply 0 to 1 re-scale normalization> : example -> -std_norm=yes\n\t-quantile_norm=<yes or no to apply quantile normalization> : example -> -quantile_norm=yes\n\t-batch_corr=<optional, name of the column in cell_data.csv to perform batch correction by> : example -> batch_id\n\t-use_bin=<optional, yes or no to use binarized columns for clustering> : example -> -use_bin=no\n\t-leiden=<optional, yes or no to perform leiden clustering> : example -> -leiden=yes\n\t-kmeans=<optional, yes or no to perform kmeans clustering> : example -> -kmeans=yes\n\t-elbow=<optional, yes or no to show elbow analysis for kmeans> : example -> -elbow=yes\n\t-k_clusters=<optional, force k number of cluster in kmeans> : example -> -k_clusters=10\n\t-refine_clusters=<optional, yes or no to refine cluster results> : example -> -refine_clusters=yes\n\t-neigh_cluster_id=<optional, cluster column name to base the neighborhood analysis on> : example -> -neigh_cluster_id=kmeans_ref', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-image_size='):
                global image_size
                image_size = int(arg[12:])
            elif arg.startswith('-analysis_markers='):
                global analysis_markers
                analysis_markers = arg[18:].split(",")
            elif arg.startswith('-cellsize_max='):
                global cellsize_max
                cellsize_max = int(arg[14:]) / 100.0
            elif arg.startswith('-cellsize_min='):
                global cellsize_min
                cellsize_min = int(arg[14:]) / 100.0
            elif arg.startswith('-custom_filter='):
                global custom_filter
                custom_filter = arg[15:]
            elif arg.startswith('-log_norm='):
                global log_norm
                log_norm = arg[10:]
            elif arg.startswith('-std_norm='):
                global std_norm
                std_norm = arg[10:]
            elif arg.startswith('-quantile_norm='):
                global quantile_norm
                quantile_norm = arg[15:]
            elif arg.startswith('-batch_corr='):
                global batch_corr
                batch_corr = arg[12:]
            elif arg.startswith('-use_bin='):
                global use_bin
                use_bin = arg[9:]
            elif arg.startswith('-leiden='):
                global leiden
                leiden = arg[8:]
            elif arg.startswith('-kmeans='):
                global kmeans
                kmeans = arg[8:]
            elif arg.startswith('-elbow='):
                global elbow
                elbow = arg[7:]
            elif arg.startswith('-k_clusters='):
                global k_clusters
                k_clusters = int(arg[12:])
            elif arg.startswith('-refine_clusters='):
                global refine_clusters
                refine_clusters = arg[17:]
            elif arg.startswith('-neigh_cluster_id='):
                global neigh_cluster_id
                neigh_cluster_id = arg[18:]


if __name__ =='__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()

    print(">>> Start time analysis =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    try:
        os.mkdir(data_folder + '/analysis/downstream')
    except OSError as error:
        print('>>> analysis/downstream folder already exists, overwriting results', flush=True)

    #Saving general settings for libraries
    sc.settings.figdir= data_folder + '/analysis/downstream'
    sc.settings.set_figure_params(format='jpg',figsize=(image_size / 100, image_size / 100))
    plt.rcParams['figure.dpi'] = 200
    sns.set(font_scale=0.6)

    df_norm, markers = data_calculations()

    if use_bin == 'yes':
        markers = [marker + '_bin' for marker in markers]

    clustering(df_norm, markers)

    print(">>> End time analysis =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
