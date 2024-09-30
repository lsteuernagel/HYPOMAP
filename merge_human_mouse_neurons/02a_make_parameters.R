## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/hypothalamus_neurons_cross_species.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_downsampled_example.rds"#
param_list$data_filepath_full = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/hypothalamus_neurons_cross_species.h5ad"
param_list$new_name_suffix = "cross_species_neurons_scvi_1"#

# signature for evaluation
param_list$genes_to_exclude_file = "data/human_features_exclude_list_2.json"

param_list$job_id = "job"

# general
param_list$n_cores = 56
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$sample_column = "Sample_ID"
param_list$batch_var = "Sample_ID"
param_list$feature_set_size = 2000
param_list$feature_set_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/feature_set_cross_species.json"
param_list$assay_name = "RNA"
param_list$integration_name = "scvi"

# scvi integration:
param_list$categorical_covariates = c("species","Dataset") # param_list$batch_var#c("Dataset",param_list$batch_var)
param_list$continuous_covariates =character(0)
param_list$n_layers = 2
param_list$n_latent = 80
param_list$n_hidden = 256
param_list$dropout_rate = 0.1
param_list$max_epochs = 500
param_list$early_stopping = FALSE
param_list$dispersion = "gene"
param_list$gene_likelihood = "zinb"
param_list$use_cuda =FALSE

## general harmonization
param_list$k_param = 30
param_list$dist_type="cosine"

# save
scUtils::writeList_to_JSON(param_list,filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/parameters_cross_species_neurons_v2.json")


