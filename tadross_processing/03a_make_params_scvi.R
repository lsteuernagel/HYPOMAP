## This script creates the json with general parameters --> make otehr jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()
# opts
opts <- workflow_options(project = "HuHy", out_path = "test/")
# basic_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/"
basic_path = opts$data_path


# define these files:
param_list$input_seurat = paste0(basic_path,"human_hypo_raw_filtered.rds")
param_list$merged_file_name = paste0(basic_path,"human_nucseq_processed")# "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/human_nucseq_processed"
param_list$project_name = "HumanNucSeq"#

# files that should not change:
#param_list$qc_path = paste0(basic_path,"quality_control/")
param_list$feature_set_file = paste0(basic_path,"feature_set.json")
param_list$genes_to_exclude_file = "features_exclude_list2.json"
param_list$scvi_script_file = "03_run_scvi.py"
param_list$project_path = basic_path
#param_list$utility_functions_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/sc_seq_recipe/process_scSeq_scvi.py"

# processing
param_list$n_cores = 56
param_list$feature_set_size = opts$nfeatures
param_list$global_seed = opts$seed
param_list$min_cells_sample = 100
param_list$sample_column = "sample"
param_list$k_param = 30
param_list$dist_type="cosine"

# general scvi
param_list$batch_var = "sample" # Provide NULL
param_list$feature_set_size = 3000
param_list$assay_name = "RNA"
param_list$integration_name = "scvi"

# scvi integration:
param_list$categorical_covariates =character(0) # param_list$batch_var#c("Dataset",param_list$batch_var)
param_list$continuous_covariates =character(0)
param_list$n_layers = 2
param_list$n_latent = 80
param_list$n_hidden = 256
param_list$dropout_rate = 0.1
param_list$max_epochs = 250 # set to 200-300
param_list$early_stopping = FALSE
param_list$dispersion = "gene"
param_list$gene_likelihood = "zinb"
param_list$use_cuda =FALSE

# save
scUtils::writeList_to_JSON(param_list,filename = paste0("parameters_scvi.json"))
message("Saving to: ", paste0("parameters_scvi.json"))

