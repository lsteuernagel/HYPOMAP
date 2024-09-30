##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree merging .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1] # param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### Load subset data
##########

subset_names = parameter_list$subset_names #c("HypoNeuron","NonNeuron","Astrocytes","Oligodendrocytes")

additional_offset = 0 # set to 0 for C0 to C6 and to 1 for C1 to C7

#subset_names = c("Astrocytes","NonNeuron","Oligodendrocytes")
cluster_number_offset = rep(0,10)
seurat_object_list = list()
comparisons_all_list = list()
comparisons_siblings_list = list()
edglist_list = list()

counter = 0

for( subname in subset_names){
  counter=counter+1
  message(">>>>>>> ",subname)
  
  param_set = parameter_list
  # load subset specifc
  param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
  param_set$new_name_suffix = paste0(param_set$new_name_suffix,"_",subname)
  param_set$clustering_key_name = paste0(subname,"_clusters")
  
  # set up
  param_set$n_cores_markers = 4
  param_set$marker_suffix = paste0("pruned_",param_set$prune_round)
  param_set$start_node = "all" # "all" for everything
  param_set$mrtree_result_file = paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_pruned_",param_set$prune_round,"_mrtree_clustering_results",".rds")
  #  paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_pruned_",pr,"_mrtree_clustering_results",".rds")
  # paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_pruned_mrtree_clustering_results",".rds")
  # load mrtree clustering
  message(param_set$mrtree_result_file)
  #if(is.null(param_set$mrtree_result_file)){
  mrtree_result = readRDS(param_set$mrtree_result_file)
  message("read from list: ",param_set$mrtree_result_file)
  # }else{
  #   stop("Cannot find param_set$mrtree_result_file")
  #   #mrtree_result = readRDS(param_set$mrtree_result_file)
  #   #message("read from file: ",param_set$mrtree_result_file)
  # }
  
  # load object
  message("Load seurat")
  seurat_object = readRDS(paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,".rds"))
  
  ## get labelmat
  labelmat_neurons = mrtree_result$labelmat
  print(colnames(labelmat_neurons))
  
  # rename levels in neurons
  message("Rename clustering levels")
  labelmat_neurons_new = labelmat_neurons
  updated_names_list = list()
  print(cluster_number_offset)
  for(i in 1:ncol(labelmat_neurons_new)){
    old_col = labelmat_neurons_new[,i]
    all_numbers = stringr::str_extract(string = old_col,pattern = "-[0-9]+") %>% stringr::str_remove(pattern = "-") %>% as.numeric()
    all_numbers = all_numbers + cluster_number_offset[i] 
    new_col = paste0("C",i+additional_offset,"-",all_numbers)
    #new_col = stringr::str_replace(string = old_col,pattern = "C[0-9]+-",paste0("C",i,"-"))
    labelmat_neurons_new[,i] = new_col
    updated_names_list[[i]] = data.frame(old_name = old_col, new_name = new_col)
    cluster_number_offset[i] = max(all_numbers)
  }
  updated_neuron_cluster_names = as.data.frame(do.call(rbind,updated_names_list)) %>% dplyr::distinct(old_name,new_name)
  colnames(labelmat_neurons_new) = paste0("C",1:ncol(labelmat_neurons_new)+additional_offset)
  
  ### clean up metadata and add
  seurat_object@meta.data = seurat_object@meta.data[,c(1:16,which(colnames(seurat_object@meta.data) %in% c("celltype_annotation","celltype_status","leiden_clusters_12_simplified")))]
  # add CO column: 
  seurat_object@meta.data$C0 =paste0("C0-",counter)
  # add other metadata
  seurat_object@meta.data = cbind(seurat_object@meta.data,labelmat_neurons_new)
  
  seurat_object_list[[subname]] = seurat_object
  
  ## load neuron marker genes
  additional_suffix = ""
  # comparisons_all_neurons = data.table::fread(paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_",param_set$start_node,"_markers_all_",param_set$marker_suffix,"_",additional_suffix,".tsv"),data.table = F)
  # comparisons_siblings_neurons = data.table::fread(paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_",param_set$start_node,"_markers_siblings_",param_set$marker_suffix,"_",additional_suffix,".tsv"),data.table = F)
  
  comparisons_all_neurons = data.table::fread(file = paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_",param_set$start_node,"_markers_all_",param_set$marker_suffix,"_",param_set$additional_suffix,".tsv"),data.table = F)
  comparisons_siblings_neurons  =  data.table::fread(file = paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_",param_set$start_node,"_markers_siblings_",param_set$marker_suffix,"_",param_set$additional_suffix,".tsv"),data.table = F)
  
  ## add new labels to markers
  comparisons_all_neurons_new = comparisons_all_neurons %>% 
    dplyr::left_join(updated_neuron_cluster_names,by=c("parent"="old_name")) %>%
    dplyr::rename(parent_new = new_name) %>%
    dplyr::left_join(updated_neuron_cluster_names,by=c("cluster_id"="old_name")) %>%
    dplyr::select(cluster = new_name,gene ,comparison,parent = parent_new,specificity,avg_log2FC,pct.1,pct.2,effect_size,p_val_adj,p_val)
  comparisons_all_neurons_new$parent[is.na(comparisons_all_neurons_new$parent)] = paste0("C0-",counter)
  
  comparisons_siblings_neurons_new = comparisons_siblings_neurons %>% 
    dplyr::left_join(updated_neuron_cluster_names,by=c("parent"="old_name")) %>%
    dplyr::rename(parent_new = new_name) %>%
    dplyr::left_join(updated_neuron_cluster_names,by=c("cluster_id"="old_name")) %>%
    dplyr::select(cluster = new_name,gene ,comparison,parent = parent_new,specificity,avg_log2FC,pct.1,pct.2,effect_size,p_val_adj,p_val)
  comparisons_siblings_neurons_new$parent[is.na(comparisons_siblings_neurons_new$parent)] = paste0("C0-",counter)
  # add to list
  comparisons_all_list[[subname]] = comparisons_all_neurons_new
  comparisons_siblings_list[[subname]] = comparisons_siblings_neurons_new
  
  ### make a new edgelist for the tree
  neuron_edgelist_old = mrtree_result$edgelist
  neuron_edgelist_new = neuron_edgelist_old %>% dplyr::left_join(updated_neuron_cluster_names,by=c("from"="old_name")) %>% 
    dplyr::rename(new_from = new_name) %>%
    dplyr::left_join(updated_neuron_cluster_names,by=c("to"="old_name")) %>%
    dplyr::rename(new_to = new_name)
  neuron_edgelist_new$new_from[is.na(neuron_edgelist_new$new_from)] = paste0("C0-",counter)
  neuron_edgelist_new =  neuron_edgelist_new %>% dplyr::select(from = new_from, to = new_to,isLeaf,level,count,totalCount,height)
  neuron_edgelist_new$level = neuron_edgelist_new$level+1
  neuron_edgelist_new = rbind(data.frame(from = "all", to = paste0("C",additional_offset,"-",counter),isLeaf=FALSE,level=2,count=nrow(neuron_edgelist_old[neuron_edgelist_old$from == "all",]),totalCount=nrow(neuron_edgelist_old),height=neuron_edgelist_old$height[1]+1),neuron_edgelist_new)
  
  edglist_list[[subname]] = neuron_edgelist_new 
  message("Completed this subset")
}

##########
### Merge data
##########

message(Sys.time(),": Merge objects " )

remerged_seurat_object = merge(seurat_object_list[[1]],seurat_object_list[2:length(seurat_object_list)],merge.dr = c("scvi","umap_scvi"))

comparisons_all = as.data.frame(do.call(rbind,comparisons_all_list))
comparisons_siblings = as.data.frame(do.call(rbind,comparisons_siblings_list))
updated_edgelist = as.data.frame(do.call(rbind,edglist_list)) %>% dplyr::arrange((from)) %>% dplyr::arrange(level) %>% dplyr::select(-height) 

### update the non neuron level 4(5)-6 with placeholders
temp_labelmat=remerged_seurat_object@meta.data[,grepl("C[0-9]+",colnames(remerged_seurat_object@meta.data))]
for(j in 2:ncol(temp_labelmat)){
  temp_labelmat[is.na(temp_labelmat[,j]),j] = temp_labelmat[is.na(temp_labelmat[,j]),j-1]
}
remerged_seurat_object@meta.data[,grepl("C[0-9]+",colnames(remerged_seurat_object@meta.data))] = temp_labelmat

##########
### Calculate a new shared UMAP
##########

# run umap and save model
message(Sys.time(),": Build UMAP on merged object with ",param_set$k_param," n.neighbors ..." )
remerged_seurat_object = RunUMAP(remerged_seurat_object,
                                 reduction = param_set$integration_name,
                                 seed.use= param_set$global_seed,
                                 dims=1:ncol(remerged_seurat_object@reductions[[param_set$integration_name]]@cell.embeddings),
                                 reduction.name=paste0("umap_",param_set$integration_name,"_hypo"),
                                 reduction.key = paste0("umap_",param_set$integration_name,"_hypo"),
                                 verbose=F,
                                 n.neighbors = param_set$k_param,
                                 return.model = TRUE)


##########
### Export objects
##########

message(Sys.time(),": Export merged objects " )
# read list again
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

parameter_list$new_name_suffix = "human_hypo_combined"
parameter_list$harmonization_folder_path = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix ,"/")
dir.create(parameter_list$harmonization_folder_path,showWarnings = F)

dummy=matrix(data = as.numeric())
remerged_seurat_object@assays[["RNA"]]@var.features = character()
remerged_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

## add to object
remerged_seurat_object@misc$edgelist = updated_edgelist

# save file
data.table::fwrite(updated_edgelist,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_edgelist_mrtree.txt"),sep="\t")

# save parameter_list$new_name_suffix = "human_hypo_combined"
data.table::fwrite(comparisons_all,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_comparisons_all_updated.txt"),sep="\t")
data.table::fwrite(comparisons_siblings,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_comparisons_siblings_updated.txt"),sep="\t")

# make file name
subset_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix)

# rds
saveRDS(remerged_seurat_object,paste0(subset_file_name,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = remerged_seurat_object,filename = paste0(subset_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(subset_file_name,".h5seurat"), dest =  paste0(subset_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)

message(Sys.time(),": Complete." )