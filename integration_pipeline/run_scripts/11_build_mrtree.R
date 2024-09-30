##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree construction .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1] # Neuron: "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/mrtree_construction_params_e5ee4c97254bd2cc12632d19f93d3a9f.json" ##ASTRO: param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/mrtree_construction_params_dfb490dc05d9d11822cdb93e6a0fa876.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load clusters

##########
### load clusters
##########

# either provide a file with all clusterings:
# ensure that cells are in right order!
if(parameter_list$load_from_file){ # load_from_file = FALSE
  message("Reading clusters from file: ",parameter_list$clusters_for_mrtree_file)
  cluster_matrix_for_mrtree = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$clusters_for_mrtree_file),data.table = F)
}else{
  # or merge from multiple files (consensus clusterings from multiple jobs)
  message("Reading clusters from folder and merge: ",parameter_list$clustering_folder)
  all_cluster_files = list.files(parameter_list$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(parameter_list$harmonization_folder_path,"consensus/")
  all_clusterings = lapply(paste0(parameter_list$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
  cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)
}
message("parameter_list$use_recon_labelmat: ",parameter_list$use_recon_labelmat)

##########
### Run mrtree
##########

message(Sys.time(),": Prepare clustering" )

# ensure that we have good numerics
cluster_matrix_for_mrtree = as.matrix(cluster_matrix_for_mrtree)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.character)
cluster_matrix_for_mrtree <- apply(cluster_matrix_for_mrtree, 2, as.numeric)
#print(head(cluster_matrix_for_mrtree))

# add to seurat --> allows to also specify e.g. the anno column for 
harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data,cluster_matrix_for_mrtree)

# input clusters
if(is.null(parameter_list$clusterings_to_use)){
  message(" Using all clusterings")
  input_clusters = cluster_matrix_for_mrtree
}else{
  if(any(! parameter_list$clusterings_to_use %in% colnames(harmonized_seurat_object@meta.data))){
    stop("Cannot find all clusterings_to_use in colnames(harmonized_seurat_object@meta.data). Stopping.")
  }else{
    message(" Using clusterings: ",paste0(parameter_list$clusterings_to_use,collapse = " | "))
    input_clusters = harmonized_seurat_object@meta.data[,parameter_list$clusterings_to_use]
  }
}

prefix = "c_"
colnames(input_clusters) = paste0(prefix,colnames(input_clusters))

non_num_indices = which(!sapply(input_clusters,is.numeric))
for(non_num_idx in non_num_indices){
  input_clusters[,non_num_idx] = as.numeric(as.factor(input_clusters[,non_num_idx]))
}

message(Sys.time(),": Build mrtree using matrix with: ",dim(input_clusters)[1]," cells and ",dim(input_clusters)[2]," cluster levels." )
# feed into mrtree
# function from 'mrtree_functions.R' which is copied from original repo
mrtree_res <- mrtree(input_clusters,
                     prefix = prefix,
                     suffix = NULL,
                     max.k = Inf,
                     consensus = FALSE,
                     sample.weighted = FALSE, # parameter_list$mr_tree_weighted,
                     augment.path = FALSE,
                     verbose = TRUE,
                     n.cores = parameter_list$n_cores)

saveRDS(mrtree_res,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_mrtree_result_raw",".rds"))

message(Sys.time(),": Saved raw result of length ",length(mrtree_res)," to: ",paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_mrtree_result_raw",".rds"))

message(Sys.time(),": Create mrtree ouput list" )

# make labelmat with colnames included
# using the recon mat keeps the lowest level ! - but labelmat.mrtree avoids duplicated levels if clusters do not change!
# I am instead now emulating the old code version by binding the last clustering level to labelmat.mrtree 
if(parameter_list$use_recon_labelmat){
  # labelmat=mrtree_res$labelmat.recon
  # Ks.recon = apply(labelmat, 2, function(y) length(unique(y)))
  # unique.idx = 1:length(Ks.recon) # don#t remove last col
  # colnames(labelmat) = paste0("K", Ks.recon[unique.idx])
  # new version:
  labelmat=mrtree_res$labelmat.mrtree
  if(length(unique(labelmat[,ncol(labelmat)])) != length(unique(mrtree_res$labelmat.recon[,ncol(mrtree_res$labelmat.recon)]))){
    labelmat=cbind(labelmat,mrtree_res$labelmat.recon[,ncol(mrtree_res$labelmat.recon)])
    colnames(labelmat)[ncol(labelmat)] = paste0("K", length(unique(mrtree_res$labelmat.recon[,ncol(mrtree_res$labelmat.recon)])))
  }
}else{
  labelmat=mrtree_res$labelmat.mrtree
}

# build dataframe with labels
n=nrow(labelmat)
backup = colnames(labelmat)
labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
colnames(labelmat)=backup
df = as.data.frame(unique(labelmat), stringsAsFactors = F)

# save in data.tree format
require(data.tree)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

# export edgelist  and nodelist from data.tree
edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
#edges$from = stringr::str_split_i(edges$from,pattern = "/",i = -1)
#edges$to = stringr::str_split_i(edges$to,pattern = "/",i = -1)
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c("all",FALSE,1,5,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# make cluster object (list)
cluster_object = list(labelmat = labelmat,
                      edgelist = edges ,
                      nodelist = nodes,
                      data_tree = tree.datatree,
                      mrtree_output = mrtree_res)
#seurat_object_harmonized@misc$mrtree_clustering_results = cluster_object

##########
###Save
##########

# save object as rds
saveRDS(cluster_object,paste0(parameter_list$harmonization_folder_path,parameter_list$clustering_key_name,"_mrtree_results",".rds"))

# bla
message(Sys.time(),": Complete." )