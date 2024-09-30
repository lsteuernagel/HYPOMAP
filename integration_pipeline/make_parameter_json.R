## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$harmonization_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_harmonization/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_merge/human_hypothalamus_merged.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_downsampled_example.rds"#
param_list$new_name_suffix = "human_hypo"#
param_list$additional_clustering_suffix = "_basic"

# signature for evaluation
param_list$genes_to_exclude_file = "data/human_features_exclude_list_2.json"

# general
param_list$n_cores = 56
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$sample_column = "Sample_ID"
param_list$batch_var = "Sample_ID"
param_list$feature_set_size = 2500
param_list$feature_set_file = paste0(param_list$harmonization_folder_path,"feature_set.json")
param_list$assay_name = "RNA"
param_list$integration_name = "scvi"

# scvi integration:
param_list$categorical_covariates =character(0) # param_list$batch_var#c("Dataset",param_list$batch_var)
param_list$continuous_covariates =character(0)
param_list$n_layers = 2
param_list$n_latent = 80
param_list$n_hidden = 256
param_list$dropout_rate = 0.1
param_list$max_epochs = 100
param_list$early_stopping = FALSE
param_list$dispersion = "gene"
param_list$gene_likelihood = "zinb"
param_list$use_cuda =FALSE

## general harmonization
param_list$k_param = 30
param_list$dist_type="cosine"

## initial clustering
param_list$target_clusterN_initial = 500
param_list$start_res_initial = 12
param_list$end_res_initial = 15
param_list$step_size_initial = 1
param_list$include_low_res_initial = TRUE
param_list$clustering_key_name_initial = "leiden_clusters"


#### Partly changed from here on for v2 params
## full clustering
param_list$leiden_repeats = 100
param_list$target_clusterN = 800
param_list$resolutions_basic = c(c(0.001,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.175,0.25,0.5,0.75),seq(1,3,0.5))
param_list$resolutions_neurons = c(c(0.001,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.175,0.25,0.5,0.75),seq(1,3,0.5),seq(4,50,2))
param_list$min_cells_valid = 10

# mrtree
param_list$clusters_for_mrtree_file = "mrtree_input_labels.txt"
param_list$use_recon_labelmat = TRUE # avoids skipping the lowest level of cluster when building the matrix !
param_list$specificity_base = 0.001
param_list$n_cores_markers = 4
# 
# # pruning:
param_list$min_cells = 20 # if below --> merge with neighbor
param_list$min_specificity = 1.5 # min specificity for a sibling marker to count
param_list$max_pvalue_prune = 0.001 # max pvalue for a sibling marker to count
param_list$min_sibling_markers = 15 # how many sibling markers are required to not merge
param_list$min_prune_level = 4 # highest level is 2 (1 does not exist because level is based on destination node)
param_list$start_nodes_pruning_markers = c("all") # use this when there are multiple marker tables after splitting the marker detection
param_list$old_prefix = "K"
param_list$new_prefix = "C"

# basic marker detection
param_list$basic_marker_filename = "_initial_markers"
param_list$assay_markers ="RNA"
param_list$assay_slot = "data"
param_list$test.use = "wilcox-stratified" # either "wilcox-stratified" or a basic Seurat FindMarkers test
param_list$logfc.threshold = 0.3
param_list$min.pct = 0.1
param_list$min.diff.pct = 0.05
param_list$max.cells.per.ident = 20000
param_list$min.cells.feature = 30 # compared to full data ! sure ?
param_list$min.cells.group =  10
param_list$base = 2
param_list$only.pos = TRUE
param_list$batch_var = param_list$batch_var
param_list$downsample_cells = 50000
param_list$add_batch_as_latent = FALSE

# # annotation
# param_list$start_nodes_annotation_markers = c("C2-1","C2-2") # can also just be = start_nodes_pruning_markers  # use this when there are multiple marker tables after splitting the marker detection

# save
scUtils::writeList_to_JSON(param_list,filename = "integration_pipeline/parameters/parameters_human_hypo_v2.json")
