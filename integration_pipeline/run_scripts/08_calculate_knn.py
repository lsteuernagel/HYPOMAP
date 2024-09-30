
print(" Import")
# Import relevant modules
import pandas as pd
import numpy as np
import scanpy as sc
import igraph
#import re
import os
import sys
import json
import gc

print(" Read parameters")

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
data_filepath_full = parameter_dict["harmonization_folder_path"]+parameter_dict["new_name_suffix"]+".h5ad"
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
job_id = parameter_dict["job_id"]
n_neighbors = parameter_dict["k_param"]
metric = parameter_dict["dist_type"] # metric = 'cosine'
#snn_name = "SNN_"+parameter_dict["integration_name"] # not used currently
min_cells_valid = parameter_dict["min_cells_valid"]
#additional_clustering_suffix = parameter_dict["additional_clustering_suffix"]
key_prefix = str(parameter_dict["clustering_key_name"]) # "leiden_clusters"

# load anndata
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
#embedding = "X_scvi"

# run neighbors
print(" Run neighbors via sc.pp.neighbors")
print(adata.obsm.keys())

npcs = adata.obsm[embedding].shape[1]
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=npcs,use_rep=embedding,random_state=global_seed)
neighbors_object = sc.Neighbors(adata)

# run multi leiden clustering:
# convert neighbors to igraph
c_igraph = neighbors_object.to_igraph()

# save
# https://igraph.org/python/doc/api/igraph.Graph.html#write_adjacency
knn_filename= results_path+new_name_suffix+"_knn_"+str(n_neighbors)+".graphml"
print("Saving igraph to "+knn_filename)
c_igraph.write(f = knn_filename, format="graphml")

print(" Neighbor graph construction complete")
