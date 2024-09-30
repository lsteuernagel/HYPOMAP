##########
### Load parameters and packages
##########

message(Sys.time(),": Starting fast marker detection .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")
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

# load mrtree clustering
message(parameter_list$mrtree_result_file)
if(is.null(parameter_list$mrtree_result_file)){
  mrtree_result = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$clustering_key_name,"_mrtree_results",".rds"))
  message("read from list: ",paste0(parameter_list$harmonization_folder_path,parameter_list$clustering_key_name,"_mrtree_results",".rds"))
}else{
  mrtree_result = readRDS(parameter_list$mrtree_result_file)
  message("read from file: ",parameter_list$mrtree_result_file)
}

# load object
neuron_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# important: start node
start_node = "all"
# start_node = parameter_list$start_node
# if(is.null(start_node)){
#   start_node = "all"
# }

# markers
n_cores_markers = parameter_list$n_cores_markers
assay_markers = parameter_list$assay_markers
assay_slot =parameter_list$assay_slot
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
batch_var = parameter_list$batch_var
test.use = parameter_list$test.use
specificity_base = parameter_list$specificity_base
message("test.use ",test.use," ... ",parameter_list$test.use)
add_batch_as_latent = parameter_list$add_batch_as_latent
additional_suffix = parameter_list$additional_suffix

if(is.null(additional_suffix)){additional_suffix = ""}
if(is.null(base)){base = 2}
if(is.null(add_batch_as_latent)){add_batch_as_latent = FALSE}
if(! test.use %in% c('negbinom', 'poisson', 'MAST', "LR") ){
  add_batch_as_latent = FALSE
}

message("additional_suffix: ",additional_suffix)
message("add_batch_as_latent: ",add_batch_as_latent)

##########
### Prepare edgelist and labelmat
##########

# get key objects
edgelist = mrtree_result$edgelist[,1:2]
message("loaded edglist with ",nrow(edgelist)," rows.")
labelmat = mrtree_result$labelmat
message("loaded labelmat with ",nrow(labelmat)," rows.")

message("Object cells: ",ncol(neuron_seurat_object),"  Labelmat cells: ",nrow(labelmat))

# if start node if specified:
# find subtree and subset edgelist:

# walk through tree:
all_children = scUtils::find_children(nodes = start_node, edges = edgelist[,1:2])
if(length(all_children) < 3){
  message("Cannot find enough children for starting node ",start_node," . Using all nodes in tree instead")
  all_children =  unique(edges[,2])
}

# subset
edgelist = edgelist[edgelist$to %in% all_children,]

# whichc ells to include when downsampleing for feature abundance check
if(start_node != "all"){
  cells_to_check = rownames(neuron_seurat_object@meta.data)[labelmat[,which(apply(labelmat,2,function(x,target){target %in% x},target=start_node))] == start_node]
}else{
  cells_to_check = rownames(neuron_seurat_object@meta.data)
}

##########
### Calculate markers between leaf-siblings ?  and merge ?
##########

## subset features to be tested;
if(ncol(neuron_seurat_object@assays[['RNA']]@data) > 30000){
  set.seed(parameter_list$global_seed)
  subset = sample(colnames(neuron_seurat_object@assays[['RNA']]@data[,cells_to_check]),size = 30000)
  gene_expr_dataset = neuron_seurat_object@assays[['RNA']]@data[,subset]
}else{
  gene_expr_dataset = neuron_seurat_object@assays[['RNA']]@data
}
gene_expr_dataset[gene_expr_dataset != 0] <- 1 # set to 1 for occ
gene_sums = Matrix::rowSums(gene_expr_dataset)
rm(gene_expr_dataset)
genes_to_include = names(gene_sums)[gene_sums > min.cells.feature ]
genes_to_include = genes_to_include[! genes_to_include %in% features_exclude_list]
message("Testing ",length(genes_to_include)," genes as cluster markers")

if(add_batch_as_latent){
  latent_vars_seurat = batch_var
  message("Adding batch_var: ",latent_vars_seurat," as latent var to Seurat based test.")
}else{
  latent_vars_seurat = NULL
}

# run marker detection on tree using findMarkers_tree2
message(Sys.time(),": Start marker detection on mrtree" )
all_markers_mrtree_list = findMarkers_tree2(neuron_seurat_object,
                                            edgelist=edgelist[,1:2],
                                            labelmat = labelmat,
                                            n_cores = n_cores_markers,
                                            assay=assay_markers,
                                            slot=assay_slot,
                                            test.use = test.use,
                                            batch_var=batch_var,
                                            latent_vars_seurat = latent_vars_seurat,
                                            genes_to_include = genes_to_include,
                                            comparison_types = c("All","Sibling"),
                                            logfc.threshold = logfc.threshold,
                                            min.pct = min.pct,
                                            min.diff.pct = min.diff.pct,
                                            max.cells.per.ident=max.cells.per.ident,
                                            min.cells.feature = min.cells.feature,
                                            min.cells.group = min.cells.group,
                                            base = base,
                                            only.pos = c("All" = TRUE,"Sibling" =FALSE))

message("names(all_markers_mrtree_list): ",names(all_markers_mrtree_list))

# extract results from list
comparisons_siblings = as.data.frame(all_markers_mrtree_list$comparisons_Sibling)
comparisons_all = as.data.frame(all_markers_mrtree_list$comparisons_All)

# save objects as txt - just to be safe
data.table::fwrite(comparisons_siblings,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,"_",additional_suffix,".tsv"),sep="\t")
data.table::fwrite(comparisons_all,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,"_",additional_suffix,".tsv"),sep="\t")


#calculate specificity --- assumes base 2!
# for siblings
comparisons_siblings = comparisons_siblings %>%
  dplyr::mutate(specificity = case_when(avg_log2FC >= 0 ~ ((pct.1+specificity_base) / (pct.2+specificity_base)) * avg_log2FC,
                                        avg_log2FC < 0 ~ ((pct.2+specificity_base) / (pct.1+specificity_base)) * avg_log2FC,
                                        TRUE ~ NA_real_))
# for all
comparisons_all = comparisons_all %>%
  dplyr::mutate(specificity = case_when(avg_log2FC >= 0 ~ ((pct.1+specificity_base) / (pct.2+specificity_base)) * avg_log2FC,
                                        avg_log2FC < 0 ~ ((pct.2+specificity_base) / (pct.1+specificity_base)) * avg_log2FC,
                                        TRUE ~ NA_real_))

##########
### Save result
##########

# save objects as txt
data.table::fwrite(comparisons_siblings,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,"_",additional_suffix,".tsv"),sep="\t")
data.table::fwrite(comparisons_all,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,"_",additional_suffix,".tsv"),sep="\t")

message(Sys.time(),": Finalized marker detection." )

# end

