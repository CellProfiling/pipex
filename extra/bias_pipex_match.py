import pandas as pd

# Reading the cell segmentation csv file
df_pipex = pd.read_csv('/Users/R_S4_pip.csv', sep=";")
df_bias = pd.read_csv('/Users/R_S4_fm.csv', sep=";")

has_cluster_id = False
if 'Cluster' in df_pipex:
    has_cluster_id = True

pipex_ref_size = 0
pipex_ref_x = 0
pipex_ref_y = 0
bias_ref_size = 0
bias_ref_x = 0
bias_ref_y = 0

value_counts = df_pipex['size'].value_counts().to_frame()
for index, row in df_pipex.iterrows():
    if round(row['size']) in value_counts.index and value_counts.loc[round(row['size'])]['size'] == 1:
        pipex_ref_size = row['size']
        pipex_ref_x = row['x']
        pipex_ref_y = row['y']
        break

row_bias = df_bias[df_bias['CELL AREA'] == pipex_ref_size]
bias_ref_size = row_bias.iloc[0]['CELL AREA']
bias_ref_x = row_bias.iloc[0]['CELL METCENTER-X']
bias_ref_y = row_bias.iloc[0]['CELL METCENTER-Y']

diff_x = bias_ref_x - pipex_ref_x
diff_y = bias_ref_y - pipex_ref_y

num_correct = 0
num_error = 0
for index, row in df_pipex.iterrows():
    curr_size = round(row['size'])
    curr_x = row['x']
    curr_y = row['y']
    index_bias = df_bias.index[
        (round(df_bias['CELL AREA']) == curr_size) & (abs(df_bias['CELL METCENTER-X'] - (curr_x + diff_x)) < 10) & (
                    abs(df_bias['CELL METCENTER-Y'] - (curr_y + diff_y)) < 10)].tolist()
    if (len(index_bias) == 1):
        df_bias.at[index_bias[0], 'PIPEX_CELL_ID'] = row['cell_id']
        if has_cluster_id:
            df_bias.at[index_bias[0], 'PIPEX_CLUSTER_ID'] = row['Cluster']

        num_correct = num_correct + 1
    else:
        num_error = num_error + 1

print("correct:", num_correct, ", error:", num_error)
df_bias.to_csv('./allfeaturesexp_link.csv', index=False, sep=';', line_terminator=';\n\n')

