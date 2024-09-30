##########
### Load parameters and packages
##########

message(Sys.time(),": Starting fast marker detection .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)
library(parallel)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/stratified_wilcoxon_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to excludes
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat -- after annotation o have the properly cleaned clusters!
if(is.null(parameter_list$seurat_object_markers )){
  message("Cannot find parameter_list$seurat_object_markers. Defaulting to object from: ",paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated",".rds"))
  harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated",".rds"))
}else{
  message("Reading seurat from: ",parameter_list$seurat_object_markers)
  harmonized_seurat_object = readRDS(parameter_list$seurat_object_markers)
}

# load clusters
#hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$additional_clustering_suffix,"_leiden_clustering.txt"),data.table = F)
# addd
# temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,hypoMap_test_initial_leiden_clustering[,],by="Cell_ID")
# rownames(temp_meta) = temp_meta$Cell_ID
# harmonized_seurat_object@meta.data = temp_meta

if(is.null(parameter_list$cluster_column_markers )){
  message("Cannot find parameter_list$cluster_column_markers. Defaulting to 'leiden_clusters_12'")
  cluster_column = "leiden_clusters_12"
  # use column with highest resolution --> hard coded at the moment !
  # cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[ncol(hypoMap_test_initial_leiden_clustering)]
  # cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[4]
}else{
  message("Using cluster column: ",parameter_list$cluster_column_markers)
  cluster_column = parameter_list$cluster_column_markers
}


# define output file name
output_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$basic_marker_filename,".txt")

# markers
assay_markers = parameter_list$assay_markers
assay_slot =parameter_list$assay_slot
test.use = parameter_list$test.use
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
use_stratified=parameter_list$use_stratified
batch_var = parameter_list$batch_var
downsample_cells = parameter_list$downsample_cells
if(is.null(base)){base = 2}
if(is.null(downsample_cells)){downsample_cells = 50000}

##########
### downsample for marker detection
##########


if(downsample_cells < nrow(harmonized_seurat_object@meta.data)){
  
  message("-----",Sys.time(),": Downsample to ",downsample_cells," cells seurat for cluster detection")
  
  downsampled_meta = scUtils::downsample_balanced_iterative(metadata = harmonized_seurat_object@meta.data,target_sample_size = downsample_cells,predictor_var = cluster_column,stepsize = 1000,global_seed = parameter_list$global_seed )
  min.cells.feature_use = max(3,floor(min.cells.feature * (downsample_cells / nrow(harmonized_seurat_object@meta.data))))
  harmonized_seurat_object_small = subset(harmonized_seurat_object,cells=downsampled_meta$Cell_ID)
  
}else{
  
  message("-----",Sys.time(),": Using full data for cluster detection")
  
  harmonized_seurat_object_small = harmonized_seurat_object
  min.cells.feature_use = min.cells.feature
}

message("min.cells.feature_use: ",min.cells.feature_use)

##########
### Calculate markers
##########

message(Sys.time(),": Start marker detection" )
Idents(harmonized_seurat_object) = cluster_column # set ident!!

# genes to test
if(ncol(harmonized_seurat_object_small@assays[['RNA']]@data) > 30000){
  set.seed(parameter_list$global_seed)
  subset = sample(colnames(harmonized_seurat_object_small@assays[['RNA']]@data),size = 30000)
  gene_expr_dataset = harmonized_seurat_object_small@assays[['RNA']]@data[,subset]
}else{
  gene_expr_dataset = harmonized_seurat_object_small@assays[['RNA']]@data
}
gene_expr_dataset[gene_expr_dataset != 0] <- 1 # set to 1 for occ
min.cells.feature_local = max(3,floor(min.cells.feature * (30000 / nrow(harmonized_seurat_object@meta.data))))
gene_sums = Matrix::rowSums(gene_expr_dataset)
rm(gene_expr_dataset)
genes_to_include = names(gene_sums)[gene_sums > min.cells.feature_local ]
genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
message("Testing ",length(genes_to_include)," genes as cluster markers")
counter=1

# simple for loop
all_clusters = unique(harmonized_seurat_object_small@meta.data[,cluster_column])
message("Using van elteren Test for batch-aware marker detection for",length(all_clusters)," clusters.")
all_markers_list = list()
counter = 1
Idents(harmonized_seurat_object_small) = cluster_column

library(foreach)
library(doParallel)
doParallel::registerDoParallel( min(6,parameter_list$n_cores)) # depends on he available RAM

# run in parallel: I pass errors to be able to set the names of the result list afterwards
all_markers_list <- foreach::foreach(current_cluster = all_clusters,.errorhandling = "remove",.verbose=FALSE) %dopar% {
  # FindAll stratified and use batch_var
  # see stratified_wilcoxon_functions.R
  temp_res =FindMarkers2.Seurat(object = harmonized_seurat_object_small,
                                ident.1 = current_cluster,
                                assay = "RNA",
                                logfc.threshold = logfc.threshold,
                                slot = "data",
                                test.use = "VE",
                                min.pct = min.pct,
                                min.diff.pct = min.diff.pct,
                                max.cells.per.ident=max.cells.per.ident,
                                min.cells.feature = min.cells.feature_use,
                                min.cells.group = min.cells.group,
                                return.thresh = 1,
                                base = 2,
                                only.pos = TRUE,
                                latent.vars = batch_var,
                                genre = "locally-best")
  temp_res$cluster = current_cluster
  temp_res$gene = rownames(temp_res)
  colnames(temp_res)[colnames(temp_res)=="avg_logFC"] = "avg_log2FC"
  gc()
  temp_res
}

gc()

all_markers = as.data.frame(do.call(rbind, all_markers_list))
all_markers$specificity = ((all_markers$pct.1 + 0.001)/ (all_markers$pct.2+ 0.001)) * all_markers$avg_log2FC
all_markers = all_markers[all_markers$p_val < 0.05,]

#all_markers$cluster = stringr::str_remove(all_markers$cluster,pattern = "c_")

all_markers$cluster = stringr::str_remove(all_markers$cluster,pattern = "c_")

##########
### Save result
##########

message("-----",Sys.time(),": Saving marker results")

# save objects as txt
data.table::fwrite(all_markers,file=output_file_name,sep="\t")

message(Sys.time(),": Completed initial marker detection." )
