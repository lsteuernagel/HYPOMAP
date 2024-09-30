##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree construction .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/hbgf_function.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
#param_file="/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/consensus_params_1_73964fb8b28d48e277467d8ffb9eceab.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### # load clusters & Consensus clustering
##########

leiden_resolutions_to_run = parameter_list$leiden_resolutions_to_run

for(leiden_res in leiden_resolutions_to_run){
  
  message(Sys.time(),": Consensus clustering with HBGF on resolution ",leiden_res)
  # define right file name
  additional_clustering_suffix = paste0("res_",leiden_res,"_repeat_",parameter_list$leiden_repeat)
  clustering_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_",additional_clustering_suffix,".tsv")
  # load clusters
  
  # ensure that cells are in right order!
  cluster_matrix_for_consensus = data.table::fread(clustering_file_name,data.table = FALSE,header = TRUE)
  
  # run consensus clustering
  
  if(sum(apply(cluster_matrix_for_consensus,2,function(x){length(unique(x))})) / ncol(cluster_matrix_for_consensus) < 1.5){ # if a majority are the same cluster just return that
    consensus_clusters = data.frame(consensus = rep(1, nrow(cluster_matrix_for_consensus)))
    colnames(consensus_clusters) = additional_clustering_suffix
  }else{
    hbgf_result = HBGF(ClusterEnsembles = cluster_matrix_for_consensus,nstart=50,optimalk = round(median(apply(cluster_matrix_for_consensus,2,max)+1)),useMiniBatch = TRUE) # +1 because python starts with cluster 0
    consensus_clusters = data.frame(consensus = hbgf_result$Clust$Clusters)
    colnames(consensus_clusters) = additional_clustering_suffix
  }
  
  # save result
  dir.create(paste0(parameter_list$harmonization_folder_path,"consensus/"),showWarnings = FALSE)
  consensus_file_name =  paste0(parameter_list$harmonization_folder_path,"consensus/",parameter_list$new_name_suffix,"_consensus_",additional_clustering_suffix,".tsv")
  data.table::fwrite(consensus_clusters,file = consensus_file_name,col.names = TRUE)
  message(Sys.time(),": Saved results")
}

message(Sys.time(),": Complete")
# 
# human_hypo_AstroEpendymal = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/AstroEpendymal/human_hypo_AstroEpendymal.rds")
# human_hypo_AstroEpendymal@meta.data$consensus_clusters = consensus_clusters
# human_hypo_AstroEpendymal@meta.data = cbind(human_hypo_AstroEpendymal@meta.data,cluster_matrix_for_consensus)
# DimPlot(human_hypo_AstroEpendymal,group.by = "46149",reduction = "umap_scvi_AstroEpendymal")

