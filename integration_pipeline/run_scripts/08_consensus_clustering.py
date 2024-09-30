
print(" Import")
# Import relevant modules
import pandas as pd
import numpy as np
import scanpy as sc
#import louvain
import igraph
#import re
import os
import sys
import json
import gc

import consensus_clustering

print(" Read parameters")

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# test:
json_file= open("/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params/leiden_params_52f4889cb1f95ef0d2b97d97560d2316.json")
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
batch=parameter_dict["batch_var"]
data_filepath_full = parameter_dict["harmonization_folder_path"]+parameter_dict["new_name_suffix"]+".h5ad"
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
target_clusterN = parameter_dict["target_clusterN"]
job_id = parameter_dict["job_id"]
start_res = parameter_dict["start_res"]
end_res = parameter_dict["end_res"]
step_size = parameter_dict["step_size"]
include_low_res = parameter_dict["include_low_res"]
n_neighbors = parameter_dict["k_param"]
snn_name = "SNN_"+parameter_dict["integration_name"]
min_cells_valid = parameter_dict["min_cells_valid"]
additional_clustering_suffix = parameter_dict["additional_clustering_suffix"]
key_prefix = str(parameter_dict["clustering_key_name"]) # "leiden_clusters"

#define resolution range, hardcoded atm
#resolutions = [round(x*step_size,3) for x in range(int(1/step_size)*start_res,int(1/step_size)*end_res+1)]
#if(include_low_res):
  #low_res_list = [0.005,0.05,0.5]
  #low_res_list = [0.001,0.005,0.01,0.05,0.1,0.175,0.25,0.5,0.75]
#  low_res_list = [0.001,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.175,0.25,0.5,0.75]
  #low_res_list = [0.001,0.005,0.0075,0.01,0.05,0.1,0.175,0.25,0.5,0.75,3,5,6,7,15,16,28,37,38,45,50]
#  resolutions = low_res_list + resolutions

# https://lhqing.github.io/ALLCools/cell_level/step_by_step/5kb/03-Clustering.html
# set params for ConsensusClustering
n_neighbors = 25
metric = 'euclidean'
min_cluster_size = 20
consensus_rate = 0.5
leiden_repeats = 100
leiden_resolution = 5
random_state = 0
n_jobs = 40
train_frac = 0.5
train_max_n = 1000
max_iter = 50
  
print(" Read anndata from "+data_filepath_full)
# read adata
adata = sc.read_h5ad(data_filepath_full)

#ensure there are no bytestrings 
str_df = adata.obs
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
str_df = str_df.set_index('Cell_ID',drop=False)
adata.obs = str_df
# for features:
str_df = adata.var
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
if 'features' in str_df.columns:
  str_df = str_df.set_index('features',drop=False)
adata.var = str_df

print(adata.obsm.keys())
embedding = "X_"+parameter_dict["integration_name"] # X_scvi

# initialize ConsensusClustering
cc = consensus_clustering.ConsensusClustering(model=None,
                         n_neighbors=n_neighbors,
                         metric=metric,
                         min_cluster_size=min_cluster_size,
                         leiden_repeats=leiden_repeats,
                         leiden_resolution=leiden_resolution,
                         consensus_rate=consensus_rate,
                         random_state=random_state,
                         train_frac=train_frac,
                         train_max_n=train_max_n,
                         max_iter=max_iter,
                         n_jobs=n_jobs)
# run predict ConsensusClustering with multiple leiden clusterings
cc.fit_predict(adata.obsm[embedding])

adata.obs[clustering_name] = cc.label

cc.save("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_2/consensus_clustering_test/nonNeurons_cc_50it_test.lib")

from joblib import load
cc = load("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_2/consensus_clustering_test/nonNeurons_cc_test.lib")

from sklearn.cluster import DBSCAN

leiden_matrix = cc.leiden_result_df

db = DBSCAN(eps=0.3, min_samples=min_cluster_size, metric = "hamming").fit(leiden_matrix)
np.unique(db.labels_,return_counts=True)


from sklearn.metrics import adjusted_rand_score
import itertools
# Get the number of clusterings
num_clusterings = leiden_matrix.shape[1]
# Initialize an empty dictionary to store ARI scores for each pair of clusterings
ari_scores_dict = {}
# Iterate through each pair of clusterings
for idx1, idx2 in itertools.combinations(range(num_clusterings), 2):
    # Print statement when idx1 changes
    if idx1 != previous_idx1:
        print(f"Current idx1: {idx1}")
        previous_idx1 = idx1
    # Get the predicted labels for clustering 1 and clustering 2
    predicted_labels_1 = leiden_matrix.iloc[:, idx1]
    predicted_labels_2 = leiden_matrix.iloc[:, idx2]
    # Calculate ARI between clustering 1 and clustering 2
    ari = adjusted_rand_score(predicted_labels_1, predicted_labels_2)
    # Add the ARI score to the dictionary
    key = (idx1, idx2)
    ari_scores_dict[key] = ari

# Convert the dictionary to a list of tuples containing index columns and ARI scores
ari_tuples = [(pair[0], pair[1], ari) for pair, ari in ari_scores_dict.items()]
# Create a pandas DataFrame from the list of tuples
ari_df = pd.DataFrame(ari_tuples, columns=['clustering1_idx', 'clustering2_idx', 'ARI'])
