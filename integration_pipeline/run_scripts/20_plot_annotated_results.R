##########
### Load parameters and packages
##########

message(Sys.time(),": Starting annotated merged object plotting .." )

message(" Load parameters and packages ")

library(Seurat)
library(tidyverse)
library(ggtree)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

command_args<-commandArgs(TRUE)
param_file = command_args[1]
#param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### Load combined
##########

message(Sys.time(),": Load .." )

plot_path = paste0(parameter_list$harmonization_folder_path,"annotation_plots/")
dir.create(plot_path)

human_hypo_combined = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
# other data
human_hypo_combined_edgelist = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","edgelist_mrtree",".txt"),data.table = F)
#human_hypo_combined_edgelist = rbind(data.frame(from = "all",to=unique(human_hypo_combined_edgelist$from[grepl("C0-",human_hypo_combined_edgelist$from)]),isLeaf=FALSE,level),human_hypo_combined_edgelist)

human_hypo_combined_markers_all = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","comparisons_all_updated",".txt"),data.table = F)
human_hypo_combined_markers_sibling = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","comparisons_siblings_updated",".txt"),data.table = F)

mrtree_result_plotting_labelmat = human_hypo_combined@meta.data[,grepl("C[0-9]+",colnames(human_hypo_combined@meta.data))]

##########
### mrtree
##########

message(Sys.time(),": Plot .." )

# make anno_df
anno_df =human_hypo_combined@misc$edgelist %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})

edgelist = human_hypo_combined@misc$edgelist

##### plot cluster tree:
source("paper_figures/tree_plotting_functions.R")
tree_color = "grey70"

leaf_levels = c(4,5,6)
for(leaf_level in leaf_levels){
  
  
  circular_tree = plot_cluster_tree(edgelist = edgelist,
                                    leaf_level=leaf_level,
                                    anno_df = anno_df ,
                                    metadata=human_hypo_combined@meta.data,
                                    label_size = 2.5, 
                                    label_size_tip = 1.75,
                                    show_genes = TRUE,
                                    vjust_label = -0.25,
                                    annotate_reverse =F,
                                    edge_color = tree_color, 
                                    node_color = tree_color)
  circular_tree = ggtree::rotate_tree(circular_tree, -90)
  circular_tree
  
  #circular_tree
  ggsave(filename = paste0(plot_path,"circular_tree_anno_",leaf_level,"_",parameter_list$new_name_suffix,".pdf"), plot = circular_tree, "pdf",dpi=400,width=300,height = 300,units="mm")
  
  ## add a heatmap
  
  # make data for third heatmap with regions
  edgelist2 = edgelist[edgelist$level <= leaf_level, ]
  leaf_level_column = gsub("-[0-9]+","",edgelist2$to[length(edgelist2$to)])
  heatmap_data3 = human_hypo_combined@meta.data %>% dplyr::select(Cell_ID,celltype_annotation,!!sym(leaf_level_column)) %>%
    dplyr::group_by(!!sym(leaf_level_column),celltype_annotation) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
    dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
    dplyr::distinct(!!sym(leaf_level_column),celltype_annotation,.keep_all=TRUE) %>% as.data.frame()
  heatmap_data3$celltype_annotation[heatmap_data3$celltype_annotation == "P2RX2_OTP"] = "Neurons"
  heatmap_matrix3 = as.matrix(heatmap_data3[,"celltype_annotation",drop=F])
  rownames(heatmap_matrix3) = heatmap_data3[,leaf_level_column]
  #colnames(heatmap_matrix3) = "R"
  
  circular_tree_heat = hypoMapUtils::add_heatmap(circular_tree,heatmap_matrix = heatmap_matrix3,heatmap_colors = getOkabeItoPalette(16),matrix_width=0.05,matrix_offset = 0.12)
  circular_tree_heat = circular_tree_heat + scale_fill_manual(values = getOkabeItoPalette(length(unique(heatmap_data3$celltype_annotation))))+ theme(legend.text=element_text(size=15))
  circular_tree_heat
  
  ggsave(filename = paste0(plot_path,"circular_tree_anno_celltypes_",leaf_level,"_",parameter_list$new_name_suffix,".pdf"), plot = circular_tree_heat, "pdf",dpi=400,width=300,height = 300,units="mm")
  
}

##########
### UMAP plots
##########

# plot per pruned level
cluster_cols = paste0("C",0:4,"_named")
cluster_cols = c(cluster_cols)
umap_name = "umap_scvi_hypo" 
for(col in cluster_cols){
  if(length(unique(human_hypo_combined@meta.data[,col])) < 30){
    label.size = 4.5
  }else if(length(unique(human_hypo_combined@meta.data[,col])) < 150){
    label.size = 2
  }else{
    label.size = 1.5
  }
  clusterplot = DimPlot(human_hypo_combined,group.by = col,raster.dpi = c(2048,2048),pt.size = 1.6,shuffle = TRUE,label = TRUE,raster = TRUE,reduction = umap_name,label.size = label.size)+NoLegend()
  # export plot
  ggsave(filename = paste0(plot_path,parameter_list$new_name_suffix,"_",col,"_umap.pdf"), plot = clusterplot, "pdf",dpi=400,width=200,height = 200,units="mm")
}

message(Sys.time(),": Complete" )
