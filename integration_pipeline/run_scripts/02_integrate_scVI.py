# Run python scVI
# runs different parameters on one assay+featureset

# Import relevant modules
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import igraph
import os
import sys
import json
import scvi
import gc
from os import listdir
from os.path import isfile, join

print(scvi.__version__)

print("Read parameters")
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
# general params from dict (could also save these lines, but I think this way it is easier to digest)
job_id=parameter_dict["job_id"]
batch_var=parameter_dict["batch_var"]
data_filepath_full = parameter_dict["harmonization_folder_path"]+parameter_dict["new_name_suffix"]+".h5ad"
feature_set_file = parameter_dict["feature_set_file"]
results_path = parameter_dict["harmonization_folder_path"]
new_name_suffix = parameter_dict["new_name_suffix"]
global_seed = parameter_dict["global_seed"]
#hvgs_set_name = parameter_dict["hvgs_set_name"]
use_cuda = parameter_dict["use_cuda"]
categorical_covariates = parameter_dict["categorical_covariates"]
continuous_covariates = parameter_dict["continuous_covariates"]

scvi.settings.seed = global_seed

# read RNA assay (RNA assay because scVI needs raw data)
print("Read anndata")
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
# additionally requires the 

# read json into dictionary
json_file = open(feature_set_file)
hvgs = json.load(json_file)
#hvgs = hvg_dict[hvgs_set_name] # add in hvgs variable

# scVI needs raw data and so we subset X to .raw and then relevant features
print("Copy raw into .X")
adata.X = adata.raw.X.copy()
adata = adata[:, hvgs]

# clean up
gc.collect()

## Run scVI
# https://www.scvi-tools.org/en/stable/user_guide/notebooks/harmonization.html
print("Preparing scVI")
# if only 1 variable --> make list

#length_cov = len(categorical_covariates)+len(continuous_covariates)
if len(categorical_covariates) == 0:
    categorical_covariates = None
else:
  if not isinstance(categorical_covariates, list):
    categorical_covariates_save=categorical_covariates
    categorical_covariates = list()
    categorical_covariates.insert(0,categorical_covariates_save)
  
if len(continuous_covariates) == 0:
    continuous_covariates = None
else:
  if not isinstance(continuous_covariates, list):
    continuous_covariates_save=continuous_covariates
    continuous_covariates = list()
    continuous_covariates.insert(0,continuous_covariates_save)
  
print("categorical_covariate_keys: "+str(categorical_covariates))
print("continuous_covariate_keys: "+str(continuous_covariates))
print("batch_key: "+str(batch_var))
adata_scvi = adata.copy()
# setup for scvi
scvi.model.SCVI.setup_anndata(adata_scvi, batch_key = batch_var)#,categorical_covariate_keys=categorical_covariates,continuous_covariate_keys=continuous_covariates)

# for all parameter combinations
print("Run scVI")
  # set up VAE
vae = scvi.model.SCVI(adata_scvi,
    n_layers=int(parameter_dict['n_layers']),
    n_latent=int(parameter_dict['n_latent']),
    n_hidden=int(parameter_dict['n_hidden']),
    dropout_rate=float(parameter_dict['dropout_rate']),
    dispersion = str(parameter_dict['dispersion']),
    gene_likelihood = str(parameter_dict['gene_likelihood'])
    #use_cuda=use_cuda
)
# train
vae.train(max_epochs=int(parameter_dict['max_epochs']), early_stopping = bool(parameter_dict['early_stopping']), use_gpu = use_cuda)
# get result
adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()
output = pd.DataFrame(adata_scvi.obsm["X_scVI"])
output = output.set_index(adata_scvi.obs_names)
output2 = output.set_axis(["scVI_" + str(s) for s in output.axes[1].to_list()], axis=1, inplace=False)
# save
output2.to_csv(results_path+new_name_suffix+"_scVI_reduction.txt", sep='\t',index=True)
               
print("Store corrected counts")
# get corrected counts result
# https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.get_normalized_expression.html#scvi.model.SCVI.get_normalized_expression
norm_counts = vae.get_normalized_expression(library_size=10e4)
norm_counts.to_csv(results_path+new_name_suffix+"_scVI_corrected.txt", sep='\t',index=True)
    
print("Store trained model")
# https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.save.html#scvi.model.SCVI.save
vae.save(results_path+new_name_suffix+"_scVI_model", overwrite=True, save_anndata=False)

# clean up
gc.collect()
                  
## End of for
print("Finalized scVI runs")















