
print(" Import")
# Import relevant modules
import warnings
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import igraph

import joblib
import leidenalg

from natsort import natsorted
from scanpy.neighbors import Neighbors
#import re
import os
import sys
import json
import gc

print(" Read parameters")

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# test:
#json_file= open("/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params/leiden_params_52f4889cb1f95ef0d2b97d97560d2316.json")
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to read the code)
data_filepath_full = parameter_dict["harmonization_folder_path"]+parameter_dict["new_name_suffix"]+".h5ad"
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
job_id = parameter_dict["job_id"]
n_neighbors = parameter_dict["k_param"]
#snn_name = "SNN_"+parameter_dict["integration_name"]
min_cells_valid = parameter_dict["min_cells_valid"]
leiden_repeats = parameter_dict["leiden_repeats"]
leiden_resolution = parameter_dict["leiden_resolution"]
n_jobs = parameter_dict["n_cores"]
additional_clustering_suffix = parameter_dict["additional_clustering_suffix"]
key_prefix = str(parameter_dict["clustering_key_name"]) # "leiden_clusters"

min_cluster_size = min_cells_valid
# leiden_repeats = 100
# leiden_resolution = 5
# n_jobs = 40

# load anndata
# TODO
print(" Read igraph from "+data_filepath_full)
#https://igraph.org/python/doc/api/igraph.Graph.html#Read_Adjacency
knn_filename= results_path+new_name_suffix+"_knn_"+str(n_neighbors)+".graphml"
c_igraph = igraph.read(filename = knn_filename,format="graphml")
#c_igraph = neighbors_object.to_igraph()

# run multi leiden clustering:

# define leiden runner helper
def _leiden_runner(g, random_states, partition_type, **partition_kwargs):
    """
    Run leiden on a graph and return the partition.

    The leiden clustering repeated len(random_states) times with different random states,
    return all clusters as a pd.DataFrame.
    """
    results = []
    for seed in random_states:
        part = leidenalg.find_partition(g, partition_type, seed=seed, **partition_kwargs)
        groups = np.array(part.membership)
        groups = pd.Categorical(
            values=groups.astype("U"),
            categories=natsorted(np.unique(groups).astype("U")),
        )
        results.append(groups)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        result_df = pd.DataFrame(results, columns=random_states)
    return result_df

# generate n different seeds for each single leiden partition
np.random.seed(global_seed)
random_states = np.random.choice(range(99999), size=leiden_repeats, replace=False)
step = max(int(leiden_repeats / n_jobs), 10)
random_state_chunks = [random_states[i : min(i + step, leiden_repeats)] for i in range(0, leiden_repeats, step)]

# kwargs parameters
partition_type=None
partition_kwargs=None
use_weights=True
n_iterations=-1
results = []

# repeat over n runs
results = []
print("Running resolution"+str(leiden_resolution))
print(f"Repeating leiden clustering {leiden_repeats} times")
with ProcessPoolExecutor(max_workers=n_jobs) as executor:
    future_dict = {}
    for random_state_chunk in random_state_chunks:
        # flip to the default partition type if not over writen by the user
        if partition_type is None:
            partition_type = leidenalg.RBConfigurationVertexPartition
        # prepare find_partition arguments as a dictionary, appending to whatever the user provided
        # it needs to be this way as this allows for the accounting of a None resolution
        # (in the case of a partition variant that doesn't take it on input)
        partition_kwargs = {}
        if use_weights:
            partition_kwargs["weights"] = np.array(c_igraph.es["weight"]).astype(np.float64)
        partition_kwargs["n_iterations"] = n_iterations
        partition_kwargs["resolution_parameter"] = leiden_resolution
        # clustering proper
        future = executor.submit(
            _leiden_runner,
            g=c_igraph,
            random_states=random_state_chunk,
            partition_type=partition_type,
            **partition_kwargs,
        )
        future_dict[future] = random_state_chunks
    for future in as_completed(future_dict):
        _ = future_dict[future]
        try:
            data = future.result()
            results.append(data)
        except Exception as exc:
            print(f"_leiden_runner generated an exception: {exc}")
            raise exc
total_result = pd.concat(results, axis=1, sort=True)
cluster_count = total_result.apply(lambda i: i.unique().size)

print(
    f"Found {cluster_count.min()} - {cluster_count.max()} clusters, "
    f"mean {cluster_count.mean():.1f}, std {cluster_count.std():.2f}"
)


# Completed leiden clustering
print(" Leiden clustering complete. Saving results")
# save result as normal tsv
#total_result.set_index(adata.obs['Cell_ID'], inplace=True)
# TODO: Update cfile name properly !!!!
filename_leiden = results_path+new_name_suffix+"_leiden_"+additional_clustering_suffix+".tsv"
total_result.to_csv(filename_leiden,sep="\t",index=False,header=True)


# unfortunately this package throws an error - so for now I use an R implementation in a separate step
#print(" Running consensus clustering")
# array with each row being one clustering !
# labels =  total_result
# labels = np.transpose(labels)
# label_ce = ensembleclustering.cluster_ensembles(labels,solver = 'cspa',verbose=True)
# print(" Saving consensus results")

print(" Complete")
