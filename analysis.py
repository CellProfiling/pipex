import sys
import os
import argparse
import math
from collections import OrderedDict

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import datetime
from pipex_utils import log, sanitize_marker_list, validate_marker_columns
import numpy as np
import PIL
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
import qnorm
from combat.pycombat import pycombat
import json
import random


data_folder = os.environ.get('PIPEX_DATA')
image_size = 1000
analysis_markers = []
use_bin = []
cellsize_max = 0
cellsize_min = 0
minmax_norm = "no"
z_norm = "no"
log_norm = "no"
quantile_norm = "no"
batch_corr = ""
leiden = "no"
leiden_res = 0.5
kmeans = "no"
k_estimation = "no"
k_clusters = 10
refine_clusters = "no"
neigh_cluster_id = ""
neigh_k_values = [1, 5, 10]
neigh_density_threshold = 0.05

max_samples = 200000


#Function to perform all data filtering, normalization and derived calculations
def data_calculations():
    #Reading the cell segmentation csv file
    df_norm = pd.read_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'))

    markers = analysis_markers

    if len(use_bin) > 0:
        markers = [m + suffix for m in markers for suffix in use_bin]

    validate_marker_columns(df_norm, markers)

    log(f"List of markers to analyze {markers}")

    for marker in markers:
        df_norm[marker] = pd.to_numeric(df_norm[marker]).fillna(0)

    #We filter biggest and smallest cells
    filter_set = set()
    cellsize_df = df_norm['size']

    if cellsize_min > 0 or cellsize_max > 0:
        size_quantiles = cellsize_df.quantile([0.0000000001 + cellsize_min, 0.9999999999 - cellsize_max])

        filter_set.update(df_norm[df_norm['size'] > size_quantiles[0.9999999999 - cellsize_max]].index.values.tolist())
        filter_set.update(df_norm[df_norm['size'] < size_quantiles[0.0000000001 + cellsize_min]].index.values.tolist())

    df_norm.drop(index=filter_set, axis=0, inplace=True)

    if cellsize_min > 0 or cellsize_max > 0:
        log("Cells filtered")

    #We normalize all the markers through min-max
    for marker in markers:
        if log_norm == 'yes':
            df_norm[marker] = np.log1p(df_norm[marker])
        if z_norm == 'yes':
            std = df_norm[marker].std()
            if std > 0:
                df_norm[marker] = (df_norm[marker] - df_norm[marker].mean()) / std
            else:
                df_norm[marker] = 0.0
        if minmax_norm == 'yes':
            marker_min = df_norm[marker].min()
            marker_max = df_norm[marker].max()
            if marker_max > marker_min:
                df_norm[marker] = df_norm[marker].apply(lambda x: (x - marker_min) / (marker_max - marker_min))
            else:
                df_norm[marker] = 0.0

    if log_norm == 'yes' or z_norm == 'yes' or minmax_norm == 'yes':
        log("Markers normalized")

    #ComBat batch correction
    if batch_corr != '':
        batch = df_norm[batch_corr].tolist()
        marker_floors = df_norm[markers].min()
        corrected = pycombat(df_norm[markers].transpose(), batch).transpose()
        df_norm[markers] = corrected.clip(lower=marker_floors, axis=1)

        log("ComBat batch correction performed")

    #Quantile normalization
    if quantile_norm == 'yes':
        df_norm[markers] = qnorm.quantile_normalize(df_norm[markers].transpose()).transpose()

        log("Quantile normalization performed")

    df_norm.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'), index=False)

    #We calculate and plot the correlations between all markers
    heatmap_fontsize = max(3, int(80 / len(markers)))
    df_corr = df_norm[markers].corr()

    fig1, ax1 = plt.subplots(figsize=(image_size / 100,image_size / 140))
    sns.heatmap(df_corr, annot=True, annot_kws={"fontsize":heatmap_fontsize}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = True, linewidths=.1, cbar=True, ax=ax1)
    fig1.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'correlation_heatmap_pearson.jpg'))
    plt.close(fig1)

    for linkage in ('average', 'complete'):
        sns_clustermap = sns.clustermap(df_corr, figsize=(image_size / 100,image_size / 140), method=linkage)
        sns_clustermap.figure.savefig(os.path.join(data_folder, 'analysis', 'downstream', f'correlation_dendrogram_pearson_{linkage}.jpg'))
        plt.close(sns_clustermap.figure)

    #We calculate and plot the spearman correlations between all markers
    df_corr_spearman = df_norm[markers].corr(method='spearman')

    fig1s, ax1s = plt.subplots(figsize=(image_size / 100,image_size / 140))
    sns.heatmap(df_corr_spearman, annot=True, annot_kws={"fontsize":heatmap_fontsize}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = True, linewidths=.1, cbar=True, ax=ax1s)
    fig1s.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'correlation_heatmap_spearman.jpg'))
    plt.close(fig1s)

    for linkage in ('average', 'complete'):
        sns_clustermap_spearman = sns.clustermap(df_corr_spearman, figsize=(image_size / 100,image_size / 140), method=linkage)
        sns_clustermap_spearman.figure.savefig(os.path.join(data_folder, 'analysis', 'downstream', f'correlation_dendrogram_spearman_{linkage}.jpg'))
        plt.close(sns_clustermap_spearman.figure)

    columns = ['marker', 'num_cells', 'percent_cells_above_mean', 'mean', 'median', 'std', 'q10', 'q25', 'q50', 'q75', 'q90'] + markers

    #We create an extra pandas dataframe with global information about each marker. This will be a new output in PIPEX analysis step as a csv file called 'cell_data_markers.csv'
    df_ext = pd.DataFrame(columns=columns)

    x_tiles = math.ceil(df_norm['x'].max() / 1844)
    y_tiles = math.ceil(df_norm['y'].max() / 1844)
    x_bins = pd.cut(df_norm['x'], bins=x_tiles, labels=False)
    y_bins = pd.cut(df_norm['y'], bins=y_tiles, labels=False)

    for marker in markers:
        marker_df = df_norm[marker]
        quantiles = marker_df.quantile([.10, .25, .50, .75, .90])

        tile_sums = df_norm.groupby([y_bins, x_bins], observed=False)[marker].sum() / (1844 ** 2)
        qif = np.zeros((y_tiles, x_tiles))
        for (i, j), val in tile_sums.items():
            qif[int(i)][int(j)] = val

        row = {'marker': marker,
               'num_cells': len(marker_df),
               'percent_cells_above_mean': (marker_df > marker_df.mean()).sum() / len(marker_df) * 100,
               'mean': marker_df.mean(),
               'median': marker_df.median(),
               'std': marker_df.std(),
               'q10': quantiles[.10],
               'q25': quantiles[.25],
               'q50': quantiles[.50],
               'q75': quantiles[.75],
               'q90': quantiles[.90]}

        for i in range(y_tiles):
            for j in range(x_tiles):
                row['QIF_' + str(i) + '_' + str(j)] = qif[i][j]

        df_ext = pd.concat([df_ext, pd.DataFrame([row])], ignore_index=True)

    for marker in markers:
        df_ext[marker] = df_corr[marker].values

    #We generate a boxplot with each marker showing the full intensity distribution across all cells
    fig2, ax2 = plt.subplots(figsize=(image_size / 100,image_size / 140))
    sns_boxplot = sns.boxplot(x = "variable", y = "value", data = pd.melt(df_norm[markers].reset_index(drop=True)), showfliers = False, color="skyblue", ax=ax2)
    sns_boxplot.set_xticklabels(sns_boxplot.get_xticklabels(), rotation = 45)
    sns_boxplot.set(xlabel=None)
    sns_boxplot.set(ylabel=None)
    fig2.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'markers_boxplot.jpg'))
    plt.close(fig2)

    df_ext.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_markers.csv'), index=False)
    log("Markers information calculated")

    del df_corr
    del df_ext

    fill_surface_html_template(markers, df_norm)
    log("Markers surface plot generated")

    return df_norm, markers


#Function pass data to the surface html template
def fill_surface_html_template(markers, df_norm):
    with open(os.path.join(os.path.dirname(__file__), "markers_surface.html"), "r") as html_template:
        html_content = html_template.read()
    markers_formatted = "["
    for marker in markers:
        markers_formatted = markers_formatted + '\"' + marker + '\",'
    markers_formatted = markers_formatted[:-1] + "]"
    html_content = html_content.replace("$$$MARKERS$$$", markers_formatted)
    min_dim = min(min(df_norm['x'].values), min(df_norm['y'].values))
    max_dim = max(max(df_norm['x'].values), max(df_norm['y'].values))
    factor = math.ceil(max_dim / 100)
    ticktext = [*range(min_dim, max_dim, factor * 20)]
    ticktext_formatted = "["
    for tick in ticktext:
        ticktext_formatted = ticktext_formatted + str(tick) + ','
    ticktext_formatted = ticktext_formatted[:-1] + "]"
    html_content = html_content.replace("$$$TICKS$$$", ticktext_formatted)
    df_surface = df_norm
    if minmax_norm != 'yes':
        df_surface = df_norm.copy()
        for marker in markers:
            marker_min = df_surface[marker].min()
            marker_max = df_surface[marker].max()
            df_surface[marker] = df_surface[marker].apply(lambda x: (x - marker_min) / (marker_max - marker_min))
    z = []
    for marker in markers:
        z.append([])
    for x in range(min_dim, max_dim, factor):
        for m in range(len(markers)):
            z_row = []
            for y in range(min_dim, max_dim, factor):
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
    f = open(os.path.join(data_folder, 'analysis', 'downstream', 'markers_surface.html'), 'w')
    f.write(html_content)
    f.close()


#Function to generate a random RGB packed integer
def random_rgb_color(seed):
    random.seed(seed)
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)

    return (r << 16) + (g << 8) + b


#Function to convert regular RGB to packed integer
def generate_cluster_color(cluster_id, cluster_color_list):
    try:
        rgb = PIL.ImageColor.getcolor(cluster_color_list[int(cluster_id)], "RGB")
        return (rgb[0] << 16) + (rgb[1] << 8) + rgb[2]
    except Exception as e:
        return ''


def check_cell_type_threshold(curr_rule, curr_percentile):
    if curr_rule == 'high':
        score = curr_percentile
        level = "high" if curr_percentile >= 75 else ("medium" if curr_percentile >= 25 else "low")
    elif curr_rule == 'low':
        score = 100 - curr_percentile
        level = "low" if curr_percentile <= 25 else ("medium" if curr_percentile <= 75 else "high")
    else:
        score = 100 - 2 * abs(curr_percentile - 50)
        level = "medium" if 25 <= curr_percentile <= 75 else ("high" if curr_percentile > 75 else "low")
    return score, level


def check_cell_type(row, cluster_id, clustering_merge_data, rank_filter):
    all_scores = list(clustering_merge_data['scores'][cluster_id]['markers'].values())
    score_list = [s for s in all_scores if s >= 0] if rank_filter == "positive_only" else all_scores
    accumulated_score = 0
    rule_match = {}
    num_rules = 1
    while ('marker' + str(num_rules)) in row and not pd.isnull(row['marker' + str(num_rules)]):
        curr_marker = row['marker' + str(num_rules)]
        curr_rule = row['rule' + str(num_rules)]
        num_rules = num_rules + 1
        if curr_marker in clustering_merge_data['scores'][cluster_id]['markers']:
            curr_score = clustering_merge_data['scores'][cluster_id]['markers'][curr_marker]
            if not score_list or (rank_filter == "positive_only" and curr_score < 0):
                curr_percentile = 0.0
            else:
                curr_percentile = float(np.mean(np.array(score_list) <= curr_score) * 100)
            rule_score, marker_level = check_cell_type_threshold(curr_rule, curr_percentile)
            accumulated_score += rule_score
            rule_match[curr_marker] = marker_level
        elif rank_filter != "none":
            return None, None
    if accumulated_score > 0:
        return accumulated_score / (num_rules - 1), rule_match

    return None, None


def calculate_cluster_info(adata, cluster_type, markers):
    sc.settings.set_figure_params(format='jpg', figsize=(image_size / 100, image_size / 100))
    sc.pl.umap(adata, color=[cluster_type], show=False, save='_' + cluster_type)
    sc.pl.spatial(adata, color=cluster_type, spot_size=20, show=False, save='_spatial_' + cluster_type)

    try:
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_type)
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_type, method="single", show=False,
                               save='nhood_enrichment_' + cluster_type + '.jpg')
    except Exception as e:
        log('Neighborhood calculations failed for cluster ' + cluster_type + ': ' + str(e))

    try:
        sq.gr.interaction_matrix(adata, cluster_key=cluster_type)
        sq.pl.interaction_matrix(adata, cluster_key=cluster_type, show=False, save='interaction_matrix_' + cluster_type + '.jpg')
    except Exception as e:
        log('Interaction matrix analysis failed for cluster ' + cluster_type + ': ' + str(e))

    try:
        sc.tl.rank_genes_groups(adata, cluster_type, method='t-test')
        sc.pl.rank_genes_groups(adata, n_genes=len(markers), sharey=False, show=False, save='')
    except Exception as e:
        log('Rank genes groups analysis failed for cluster ' + cluster_type)

    sc.pl.heatmap(adata, markers, groupby=cluster_type, swap_axes=True, cmap='viridis', dendrogram=False, show=False,
                  save='_' + cluster_type)


def _sort_json_keys(obj):
    if isinstance(obj, dict):
        def _key(k):
            try:
                return (0, int(k), '')
            except (ValueError, TypeError):
                return (1, 0, str(k))
        return {k: _sort_json_keys(v) for k, v in sorted(obj.items(), key=lambda x: _key(x[0]))}
    if isinstance(obj, list):
        return [_sort_json_keys(i) for i in obj]
    return obj


def refine_clustering(adata, cluster_type, curr_ref_id, cell_types_ref):
    clustering_merge_data = {}
    clustering_merge_data['scores'] = {}
    clustering_merge_data['cell_types'] = {}
    if adata.uns.get('rank_genes_groups', {}).get('params', {}).get('groupby') != cluster_type:
        sc.tl.rank_genes_groups(adata, cluster_type, method='t-test')
    for cluster_id in adata.obs[cluster_type].unique():
        cluster_score_list = []
        cluster_merge_clusters_scores = {}
        cluster_merge_clusters_scores['markers'] = {}
        rank_df = sc.get.rank_genes_groups_df(adata, group=cluster_id)
        marker_score_map = rank_df.set_index('names')['scores'].to_dict()
        for marker_id, curr_score in marker_score_map.items():
            cluster_merge_clusters_scores['markers'][marker_id] = float(curr_score)
            cluster_score_list.append(curr_score)

        cluster_merge_clusters_scores['rank_filter'] = {}

        cluster_merge_clusters_scores['rank_filter']['all'] = {}
        cluster_merge_clusters_scores['rank_filter']['all']['score_max'] = float(max(cluster_score_list))
        cluster_merge_clusters_scores['rank_filter']['all']['score_min'] = float(min(cluster_score_list))
        cluster_merge_clusters_scores['rank_filter']['all']['score_dif'] = cluster_merge_clusters_scores['rank_filter']['all']['score_max'] - cluster_merge_clusters_scores['rank_filter']['all']['score_min']
        cluster_merge_clusters_scores['rank_filter']['all']['q75'] = float(np.percentile(cluster_score_list, 75))
        cluster_merge_clusters_scores['rank_filter']['all']['q25'] = float(np.percentile(cluster_score_list, 25))

        cluster_score_list_positive = [x for x in cluster_score_list if x >= 0]
        if len(cluster_score_list_positive) > 0:
            cluster_merge_clusters_scores['rank_filter']['positive_only'] = {}
            cluster_merge_clusters_scores['rank_filter']['positive_only']['score_max'] = float(max(cluster_score_list_positive))
            cluster_merge_clusters_scores['rank_filter']['positive_only']['score_min'] = float(min(cluster_score_list_positive))
            cluster_merge_clusters_scores['rank_filter']['positive_only']['score_dif'] = cluster_merge_clusters_scores['rank_filter']['positive_only']['score_max'] - cluster_merge_clusters_scores['rank_filter']['positive_only']['score_min']
            cluster_merge_clusters_scores['rank_filter']['positive_only']['q75'] = float(np.percentile(cluster_score_list_positive, 75))
            cluster_merge_clusters_scores['rank_filter']['positive_only']['q25'] = float(np.percentile(cluster_score_list_positive, 25))

        clustering_merge_data['scores'][cluster_id] = cluster_merge_clusters_scores
        clustering_merge_data['cell_types'][cluster_id] = []

    for index, row in cell_types_ref.iterrows():
        for cluster_id in clustering_merge_data['scores']:
            if row['rank_filter'] != "positive_only" or "positive_only" in clustering_merge_data['scores'][cluster_id]['rank_filter']:
                cell_type_prob, marker_ranks = check_cell_type(row, cluster_id, clustering_merge_data, row['rank_filter'])
                if cell_type_prob is not None:
                    cell_type_name = '.'.join(part for part in [row['cell_group'], row['cell_type'], row['cell_subtype']] if str(part).strip() not in ('', 'nan'))
                    curr_final_merging_data = {'cell_type': cell_type_name, 'prob': cell_type_prob, 'rank_filter': row['rank_filter'], 'confidence_threshold': row['min_confidence']}
                    curr_final_merging_data['marker_ranks'] = ','.join(k + ':' + v for k, v in marker_ranks.items())
                    clustering_merge_data['cell_types'][cluster_id].append(curr_final_merging_data)

    clustering_merge_data['candidates'] = {}
    adata.obs[cluster_type + "_ref" + curr_ref_id] = adata.obs[cluster_type].astype(str)
    adata.obs[cluster_type + "_ref" + curr_ref_id + "_p"] = adata.obs[cluster_type].astype(str)
    ordered_cluster_keys = list(clustering_merge_data['cell_types'])
    ordered_cluster_keys.sort()
    for cluster_id in ordered_cluster_keys:
        best_candidate = None
        best_real_confidence = 0
        for curr_cell_type in clustering_merge_data['cell_types'][cluster_id]:
            if curr_cell_type['prob'] >= int(curr_cell_type['confidence_threshold']):
                best_candidate = { 'cell_type': curr_cell_type['cell_type'], 'prob' : curr_cell_type['prob']  / 100 }
                best_real_confidence = curr_cell_type['prob']
                break

        if best_real_confidence > 0:
            clustering_merge_data['candidates'][cluster_id] = best_candidate
            adata.obs.loc[adata.obs[cluster_type + "_ref" + curr_ref_id] == cluster_id, cluster_type + "_ref" + curr_ref_id] = best_candidate['cell_type']
            adata.obs.loc[adata.obs[cluster_type + "_ref" + curr_ref_id + "_p"] == cluster_id, cluster_type + "_ref" + curr_ref_id + "_p"] = '{:.1%}'.format(best_candidate['prob']) #best_candidate['real_confidence'][:-1]

    with open(os.path.join(data_folder, 'analysis', 'downstream', 'cell_types_result_' + cluster_type + curr_ref_id + '.json'), 'w') as outfile:
        json.dump(_sort_json_keys(clustering_merge_data), outfile, indent = 4)


#Function to perform different cluster methods
def clustering(df_norm, markers):
    adata = sc.AnnData(df_norm[markers])
    adata.obs_names = ('cell_id_' + df_norm['cell_id'].astype(str)).values

    adata.obs["id"] = df_norm['cell_id'].values
    adata.obs["size"] = df_norm['size'].values

    adata.obsm["spatial"] = df_norm[['x', 'y']].values
    log("Anndata object created")

    #We take the chance to show spatially the intensities of every marker
    if batch_corr == '' and len(df_norm.index) <= max_samples:
        for marker in markers:
            fig, ax = plt.subplots(figsize=(image_size / 100, image_size / 100))
            scatter = ax.scatter(df_norm['x'], df_norm['y'], c=df_norm[marker], cmap='viridis', s=1, rasterized=True)
            plt.colorbar(scatter, ax=ax)
            ax.set_title(marker)
            ax.invert_yaxis()
            ax.set_aspect('equal')
            fig.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'spatial_' + marker + '.jpg'))
            plt.close(fig)
        log("Spatial plots per marker created")
    else:
        log("Dataset too big to create spatial plots per marker")

    #We calculate PCA, neighbors and UMAP for the anndata
    sc.pp.pca(adata, n_comps=min(len(markers), 50, adata.n_obs - 1, adata.n_vars - 1))

    pca_loadings = adata.varm['PCs']
    loadings_df = pd.DataFrame(pca_loadings, index=adata.var_names, columns=[f'PC{i + 1}' for i in range(pca_loadings.shape[1])])
    loadings_df.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'PCA_loadings.csv'))
    variance_explained = adata.uns['pca']['variance']
    variance_ratio = adata.uns['pca']['variance_ratio']
    variance_df = pd.DataFrame({
        'Variance': variance_explained,
        'Variance Ratio': variance_ratio,
        'Cumulative Variance Ratio': np.cumsum(variance_ratio)
    }, index=[f'PC{i + 1}' for i in range(len(variance_explained))])
    variance_df.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'PCA_variance.csv'))

    nan_count = np.isnan(adata.obsm['X_pca']).sum()
    if nan_count > 0:
        log(f"WARNING: {nan_count} NaN values found in PCA output — likely caused by zero-variance markers. Replacing with 0.")
    adata.obsm['X_pca'] = np.nan_to_num(adata.obsm['X_pca'], copy=False)
    log("PCA calculated")

    num_neighbors = int(max(5, 15 * min(1, max_samples / len(df_norm.index))))
    log(f"n_neighbors set to {num_neighbors}")
    sc.pp.neighbors(adata, n_neighbors=num_neighbors)
    log("Neighbors graph calculated")

    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=num_neighbors)
    log("Spatial neighbors graph calculated")

    sc.tl.umap(adata)
    log("UMAP calculated")
    sc.pl.umap(adata, show=False, save='_base')

    cell_types = None
    if refine_clusters == "yes":
        try:
            _ct = pd.read_csv(os.path.join(data_folder, 'cell_types.csv'))
            _ct['ref_id'] = _ct['ref_id'].astype(str)
            if set(['ref_id', 'cell_group', 'cell_type', 'cell_subtype', 'rank_filter', 'min_confidence']).issubset(set(_ct.columns.tolist())):
                cell_types = _ct
            else:
                log("cell_types.csv is malformed")
        except Exception as e:
            print(e)
            log("Failed to read cell_types.csv")

    if leiden == 'yes':
        #We calculate leiden cluster
        sc.tl.leiden(adata, resolution=leiden_res, random_state=0)
        log("Leiden cluster calculated")

        #We print the complete leiden cluster and all related information
        calculate_cluster_info(adata, "leiden", markers)

        if refine_clusters == "yes" and cell_types is not None:
            try:
                for curr_ref_id in cell_types["ref_id"].unique():
                    refine_clustering(adata, 'leiden', curr_ref_id, cell_types[cell_types["ref_id"] == curr_ref_id])
                    calculate_cluster_info(adata, "leiden_ref" + curr_ref_id, markers)
            except Exception as e:
                print(e)
                log("Failed at refining leiden cluster")

    if kmeans == 'yes':
        if k_estimation == 'yes':
            if len(df_norm.index) <= max_samples:
                log("Performing kmeans k estimation")
                distortions = []
                inertias = []
                silhouette_scores = []
                K = range(1, 20)
                pca_data = adata.obsm['X_pca']
                silhouette_sample_size = min(10000, len(df_norm.index))

                for k in K:
                    kmeanModel = KMeans(n_clusters=k, random_state=0).fit(pca_data)
                    distortions.append(sum(np.min(cdist(pca_data, kmeanModel.cluster_centers_,
                                        'euclidean'), axis=1)) / pca_data.shape[0])
                    inertias.append(kmeanModel.inertia_)
                    if k >= 2:
                        silhouette_scores.append(silhouette_score(pca_data, kmeanModel.labels_, sample_size=silhouette_sample_size, random_state=0))
                    log(f"Kmeans k estimation calculated for k {k}")

                plt.figure()
                plt.plot(K, distortions, 'bx-')
                plt.xlabel('Values of K')
                plt.ylabel('Distortion')
                plt.title('The Elbow Method using Distortion')
                plt.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'elbow_distortion.jpg'))
                plt.clf()
                plt.close()

                plt.figure()
                plt.plot(K, inertias, 'bx-')
                plt.xlabel('Values of K')
                plt.ylabel('Inertia')
                plt.title('The Elbow Method using Inertia')
                plt.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'elbow_inertia.jpg'))
                plt.clf()
                plt.close()

                plt.figure()
                plt.plot(range(2, 20), silhouette_scores, 'gx-')
                plt.xlabel('Values of K')
                plt.ylabel('Silhouette Score')
                plt.title('Silhouette Score per K (higher is better)')
                plt.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'elbow_silhouette.jpg'))
                plt.clf()
                plt.close()
                log("Kmeans k estimation calculated")
            else:
                log("Dataset too big to perform kmeans k estimation")

        #We print the complete spatial kmeans cluster and all related information. Please note that the k used is the one passad as parameter (or 10 by default)
        kmeans_cluster = KMeans(n_clusters=k_clusters, random_state=0).fit(adata.obsm['X_pca'])
        adata.obs['kmeans'] = kmeans_cluster.labels_.astype(str)

        log("Kmeans cluster calculated")

        #We print the complete kmeans cluster and all related information
        calculate_cluster_info(adata, "kmeans", markers)

        if refine_clusters == "yes" and cell_types is not None:
            try:
                for curr_ref_id in cell_types["ref_id"].unique():
                    refine_clustering(adata, 'kmeans', curr_ref_id, cell_types[cell_types["ref_id"] == curr_ref_id])
                    calculate_cluster_info(adata, "kmeans_ref" + curr_ref_id, markers)
            except Exception as e:
                print(e)
                log("Failed at refining kmeans cluster")

    if neigh_cluster_id != "":
        if len(df_norm.index) > max_samples:
            log("Dataset too large for neighborhood analysis, skipping")
        else:
            if neigh_cluster_id not in adata.obs:
                adata.obs[neigh_cluster_id] = df_norm[neigh_cluster_id].astype('category')
            try:
                sq.gr.centrality_scores(adata, neigh_cluster_id)
                sq.pl.centrality_scores(adata, neigh_cluster_id, save=(neigh_cluster_id + "_centrality_scores.jpg"))
                log("Neighborhood centrality scores calculated")
            except Exception as e:
                log("Neighborhood analysis failed: " + str(e))
            neighborhood_cell_type_analysis(adata, neigh_cluster_id, neigh_k_values, neigh_density_threshold, data_folder, image_size)

    if leiden == 'yes' or kmeans == 'yes':
        df = pd.read_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'))
        if leiden == 'yes':
            obs_by_id = adata.obs.set_index('id')
            df['leiden'] = df['cell_id'].map(obs_by_id['leiden']).astype(str)
            df['leiden'] = df['leiden'].fillna('')
            leiden_color_list = adata.uns['leiden_colors']
            df['leiden_color'] = df.apply(lambda row: generate_cluster_color(row['leiden'], leiden_color_list), axis=1)
            if refine_clusters == "yes" and cell_types is not None:
                list_cell_types = cell_types["ref_id"].unique()
                if len(list_cell_types) > 1:
                    df['leiden_ref_merged'] = ''
                    df['leiden_ref_merged_color'] = 0
                for curr_ref_id in list_cell_types:
                    df['leiden_ref' + curr_ref_id] = df['leiden']
                    df['leiden_ref' + curr_ref_id] = df['cell_id'].map(obs_by_id['leiden_ref' + curr_ref_id]).astype(str)
                    df['leiden_ref' + curr_ref_id + '_p'] = 0
                    df['leiden_ref' + curr_ref_id + '_p'] = df['cell_id'].map(obs_by_id['leiden_ref' + curr_ref_id + '_p']).astype(str)
                    df['leiden_ref' + curr_ref_id + '_color'] = df['leiden_color']
                    for leiden_ref_id in adata.obs['leiden_ref' + curr_ref_id].unique():
                        df.loc[df['leiden_ref' + curr_ref_id] == leiden_ref_id, "leiden_ref" + curr_ref_id + "_color"] = df.loc[df['leiden_ref' + curr_ref_id] == leiden_ref_id, "leiden_ref" + curr_ref_id + "_color"].values[0]
                    if len(list_cell_types) > 1:
                        df['leiden_ref_merged'] = df['leiden_ref_merged'] + '-' + df['leiden_ref' + curr_ref_id].astype(str)
                if len(list_cell_types) > 1:
                    list_cluster_ids = list(df['leiden_ref_merged'].unique())
                    df['leiden_ref_merged_color'] = df.apply(lambda row: random_rgb_color(list_cluster_ids.index(row['leiden_ref_merged'])), axis=1)
            df_norm['leiden'] = df_norm['cell_id'].map(df.set_index('cell_id')['leiden']).astype(str)
            df_norm['leiden_color'] = df_norm['cell_id'].map(df.set_index('cell_id')['leiden_color']).astype(str)

            df_corr = pd.concat([df_norm[markers], pd.get_dummies(df_norm['leiden'], prefix='clusterL')], axis=1).corr()

            fig1, ax1 = plt.subplots(figsize=(image_size / 100,image_size / 140))
            leiden_heatmap_fontsize = max(3, int(80 / (len(markers) + len(adata.obs['leiden'].unique()))))
            sns.heatmap(df_corr, annot=True, annot_kws={"fontsize":leiden_heatmap_fontsize}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = True, linewidths=.1, cbar=True, ax=ax1)
            fig1.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'leiden_clusters_correlation_heatmap.jpg'))
            plt.close(fig1)

        if kmeans == 'yes':
            #We add to the original cell segmentation csv file the calculated kmeans group for each cell
            obs_by_id = adata.obs.set_index('id')
            df['kmeans'] = df['cell_id'].map(obs_by_id['kmeans']).astype(str)
            df['kmeans'] = df['kmeans'].fillna('')
            kmeans_color_list = adata.uns['kmeans_colors']
            df['kmeans_color'] = df.apply(lambda row: generate_cluster_color(row['kmeans'], kmeans_color_list), axis=1)
            if refine_clusters == "yes" and cell_types is not None:
                list_cell_types = cell_types["ref_id"].unique()
                if len(list_cell_types) > 1:
                    df['kmeans_ref_merged'] = ''
                    df['kmeans_ref_merged_color'] = 0
                for curr_ref_id in list_cell_types:
                    df['kmeans_ref' + curr_ref_id] = df['kmeans']
                    df['kmeans_ref' + curr_ref_id] = df['cell_id'].map(obs_by_id['kmeans_ref' + curr_ref_id]).astype(str)
                    df['kmeans_ref' + curr_ref_id + '_p'] = 0
                    df['kmeans_ref' + curr_ref_id + '_p'] = df['cell_id'].map(obs_by_id['kmeans_ref' + curr_ref_id + '_p']).astype(str)
                    df['kmeans_ref' + curr_ref_id + '_color'] = df['kmeans_color']
                    for kmeans_ref_id in adata.obs['kmeans_ref' + curr_ref_id].unique():
                        df.loc[df['kmeans_ref' + curr_ref_id] == kmeans_ref_id, "kmeans_ref" + curr_ref_id + "_color"] = df.loc[df['kmeans_ref' + curr_ref_id] == kmeans_ref_id, "kmeans_ref" + curr_ref_id + "_color"].values[0]
                    if len(list_cell_types) > 1:
                        df['kmeans_ref_merged'] = df['kmeans_ref_merged'] + '-' + df['kmeans_ref' + curr_ref_id].astype(str)
                if len(list_cell_types) > 1:
                    list_cluster_ids = list(df['kmeans_ref_merged'].unique())
                    df['kmeans_ref_merged_color'] = df.apply(lambda row: random_rgb_color(list_cluster_ids.index(row['kmeans_ref_merged'])), axis=1)
            df_norm['kmeans'] = df_norm['cell_id'].map(df.set_index('cell_id')['kmeans']).astype(str)
            df_norm['kmeans_color'] = df_norm['cell_id'].map(df.set_index('cell_id')['kmeans_color']).astype(str)

            df_corr = pd.concat([df_norm[markers], pd.get_dummies(df_norm['kmeans'], prefix='clusterK')], axis=1).corr()

            fig1, ax1 = plt.subplots(figsize=(image_size / 100,image_size / 140))
            kmeans_heatmap_fontsize = max(3, int(80 / (len(markers) + k_clusters)))
            sns.heatmap(df_corr, annot=True, annot_kws={"fontsize":kmeans_heatmap_fontsize}, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1, center = 0, square = True, linewidths=.1, cbar=True, ax=ax1)
            fig1.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'kmeans_clusters_correlation_heatmap.jpg'))
            plt.close(fig1)

        df.to_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'), index=False)
        df_norm.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'), index=False)
        adata.write(os.path.join(data_folder, 'analysis', 'downstream', 'anndata.h5ad'))

def neighborhood_cell_type_analysis(adata, neigh_cluster_id, k_values, density_threshold, data_folder, image_size):
    k_values = sorted(set(k_values))[:3]
    cell_types = adata.obs[neigh_cluster_id].astype(str).values
    unique_types = sorted(set(cell_types), key=lambda x: int(x) if x.lstrip('-').isdigit() else x)
    n_types = len(unique_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    cell_type_idx = np.array([type_to_idx[t] for t in cell_types])
    coords = adata.obsm["spatial"]

    tree = KDTree(coords)
    k_max = min(max(k_values), len(coords) - 1)
    k_values = [k for k in k_values if k <= k_max]
    _, all_idxs = tree.query(coords, k=k_max + 1)

    # --- Neighbor composition matrices ---
    compositions = {}
    for k in k_values:
        neighbor_type_idx = cell_type_idx[all_idxs[:, 1:k + 1]]
        comp = np.zeros((n_types, n_types))
        for ki in range(k):
            np.add.at(comp, (cell_type_idx, neighbor_type_idx[:, ki]), 1)
        row_sums = comp.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        compositions[k] = comp / row_sums * 100

    fontsize = max(4, int(60 / n_types))
    n_plots = len(k_values)
    base = image_size / 100

    # Heatmap
    fig, axes = plt.subplots(1, n_plots, figsize=(base * n_plots, base))
    if n_plots == 1:
        axes = [axes]
    for ax, k in zip(axes, k_values):
        df_comp = pd.DataFrame(compositions[k], index=unique_types, columns=unique_types)
        sns.heatmap(df_comp, annot=True, fmt='.1f', cmap='YlOrRd', vmin=0, vmax=100, ax=ax,
                    annot_kws={"fontsize": fontsize}, linewidths=0.3)
        ax.set_title(f'k={k} neighbors (%)', fontsize=fontsize + 2)
        ax.set_xlabel('Neighbor cell type', fontsize=fontsize)
        ax.set_ylabel('Cell type', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
    fig.tight_layout()
    fig.savefig(os.path.join(data_folder, 'analysis', 'downstream',
                             neigh_cluster_id + '_neighborhood_analysis_cell_types_heatmap.jpg'), dpi=150)
    plt.close(fig)
    log("Neighborhood composition heatmap saved")

    # Stacked bar
    colors = [plt.cm.tab20(i / 20) for i in range(min(n_types, 20))]
    fig, axes = plt.subplots(1, n_plots, figsize=(base * n_plots, base))
    if n_plots == 1:
        axes = [axes]
    for ax, k in zip(axes, k_values):
        df_comp = pd.DataFrame(compositions[k], index=unique_types, columns=unique_types)
        bottom = np.zeros(n_types)
        for j, col_type in enumerate(unique_types):
            ax.bar(unique_types, df_comp[col_type].values, bottom=bottom,
                   label=col_type, color=colors[j % len(colors)])
            bottom += df_comp[col_type].values
        ax.set_title(f'k={k} neighbors (%)', fontsize=fontsize + 2)
        ax.set_xlabel('Cell type', fontsize=fontsize)
        ax.set_ylabel('% of neighbors', fontsize=fontsize)
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', rotation=45, labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower right', fontsize=fontsize, title='Neighbor type',
               title_fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(os.path.join(data_folder, 'analysis', 'downstream',
                             neigh_cluster_id + '_neighborhood_analysis_cell_types_stackbar.jpg'), dpi=150)
    plt.close(fig)
    log("Neighborhood composition stacked bar saved")

    # --- Spatial distribution classification ---
    global_median_nn = np.median(tree.query(coords, k=2)[0][:, 1])
    dbscan_eps = global_median_nn * 3

    results = []
    for cell_type in unique_types:
        mask = cell_types == cell_type
        type_coords = coords[mask]
        n_cells = len(type_coords)
        if n_cells < 10:
            results.append({'cell_type': cell_type, 'n_cells': n_cells,
                             'distribution': 'too few cells', 'n_clusters': 0,
                             'avg_cells_per_cluster': 0})
            continue
        labels = DBSCAN(eps=dbscan_eps, min_samples=5).fit(type_coords).labels_
        noise_fraction = (labels == -1).sum() / n_cells
        clustered = labels[labels != -1]
        n_clusters = len(set(clustered)) if len(clustered) > 0 else 0
        avg_cluster_size = len(clustered) / n_clusters if n_clusters > 0 else 0
        if noise_fraction > 0.5:
            distribution = 'scattered sparsely'
        elif avg_cluster_size >= density_threshold * n_cells:
            distribution = 'clustered densely'
        else:
            distribution = 'clustered sparsely'
        results.append({'cell_type': cell_type, 'n_cells': n_cells, 'distribution': distribution,
                         'n_clusters': n_clusters, 'avg_cells_per_cluster': round(avg_cluster_size, 1)})
        log(f"Cell type '{cell_type}': {distribution} ({n_clusters} clusters, avg {round(avg_cluster_size, 1)} cells/cluster)")

    pd.DataFrame(results).to_csv(
        os.path.join(data_folder, 'analysis', 'downstream', neigh_cluster_id + '_spatial_distribution.csv'),
        index=False)
    log("Spatial distribution classification saved")


#Function to handle the command line parameters passed
def options(argv):
    converted = ['--' + a[1:] if a.startswith('-') and not a.startswith('--') else a for a in argv]
    parser = argparse.ArgumentParser(prog='analysis.py')
    parser.add_argument('--data', default=os.environ.get('PIPEX_DATA'),
        help='path to images folder : example -> -data=/lab/projectX/images')
    parser.add_argument('--image_size', type=int, default=1000,
        help='one-side approximate resolution : example -> -image_size=1000')
    parser.add_argument('--analysis_markers', type=lambda s: [x.strip() for x in s.split(',')], default=[],
        help='list of markers to analyze : example -> -analysis_markers=AMY2A,SST,GORASP2')
    parser.add_argument('--cellsize_max', type=lambda s: float(s) / 100.0, default=0,
        help='percentage of biggest cells to remove : example -> -cellsize_max=5')
    parser.add_argument('--cellsize_min', type=lambda s: float(s) / 100.0, default=0,
        help='percentage of smallest cells to remove : example -> -cellsize_min=5')
    parser.add_argument('--log_norm', choices=['yes', 'no'], default='no',
        help='apply log n+1 normalization : example -> -log_norm=yes')
    parser.add_argument('--z_norm', choices=['yes', 'no'], default='no',
        help='apply z normalization : example -> -z_norm=yes')
    parser.add_argument('--minmax_norm', choices=['yes', 'no'], default='no',
        help='apply 0 to 1 re-scale normalization : example -> -minmax_norm=yes')
    parser.add_argument('--quantile_norm', choices=['yes', 'no'], default='no',
        help='apply quantile normalization : example -> -quantile_norm=yes')
    parser.add_argument('--batch_corr', default='',
        help='column in cell_data.csv for batch correction : example -> -batch_corr=batch_id')
    parser.add_argument('--use_bin', type=lambda s: [x.strip() for x in s.split(',') if x.strip()], default=[],
        help='comma-separated suffixes for marker input columns : example -> -use_bin=_local_90')
    parser.add_argument('--leiden', choices=['yes', 'no'], default='no',
        help='perform leiden clustering : example -> -leiden=yes')
    parser.add_argument('--leiden_res', type=float, default=0.5,
        help='leiden resolution, higher values produce more clusters : example -> -leiden_res=0.5')
    parser.add_argument('--kmeans', choices=['yes', 'no'], default='no',
        help='perform kmeans clustering : example -> -kmeans=yes')
    parser.add_argument('--k_estimation', choices=['yes', 'no'], default='no',
        help='show k estimation analysis for kmeans : example -> -k_estimation=yes')
    parser.add_argument('--k_clusters', type=int, default=10,
        help='force k number of clusters in kmeans : example -> -k_clusters=10')
    parser.add_argument('--refine_clusters', choices=['yes', 'no'], default='no',
        help='refine cluster results : example -> -refine_clusters=yes')
    parser.add_argument('--neigh_cluster_id', default='',
        help='cluster column for neighborhood analysis : example -> -neigh_cluster_id=kmeans')
    parser.add_argument('--neigh_k_values', default=[1, 5, 10],
        type=lambda s: sorted(set(int(x.strip()) for x in s.split(',')))[:3],
        help='up to 3 k values for neighbor composition (max 3) : example -> -neigh_k_values=1,5,10')
    parser.add_argument('--neigh_density_threshold', type=float, default=0.05,
        help='fraction of cell type size to distinguish sparse vs dense clusters : example -> -neigh_density_threshold=0.05')
    if not argv:
        parser.print_help()
        sys.exit()
    return parser.parse_args(converted)


if __name__ =='__main__':
    args = options(sys.argv[1:])
    data_folder = args.data
    image_size = args.image_size
    analysis_markers = sanitize_marker_list(args.analysis_markers)
    cellsize_max = args.cellsize_max
    cellsize_min = args.cellsize_min
    log_norm = args.log_norm
    z_norm = args.z_norm
    minmax_norm = args.minmax_norm
    quantile_norm = args.quantile_norm
    batch_corr = args.batch_corr
    use_bin = args.use_bin
    leiden = args.leiden
    leiden_res = args.leiden_res
    kmeans = args.kmeans
    k_estimation = args.k_estimation
    k_clusters = args.k_clusters
    refine_clusters = args.refine_clusters
    neigh_cluster_id = args.neigh_cluster_id
    neigh_k_values = args.neigh_k_values
    neigh_density_threshold = args.neigh_density_threshold

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()
    with open(os.path.join(data_folder, 'log_settings_analysis.txt'), 'w+', encoding='utf-8') as f:
        f.write(">>> Start time analysis = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))
        f.close()

    log("Start time analysis")

    os.makedirs(os.path.join(data_folder, 'analysis', 'downstream'), exist_ok=True)

    #Saving general settings for libraries
    sc.settings.figdir= os.path.join(data_folder, 'analysis', 'downstream')
    sc.settings.set_figure_params(format='jpg',figsize=(image_size / 100, image_size / 100))
    plt.rcParams['figure.dpi'] = 200
    sns.set(font_scale=0.6)

    df_norm, markers = data_calculations()

    clustering(df_norm, markers)

    log("End time analysis")
