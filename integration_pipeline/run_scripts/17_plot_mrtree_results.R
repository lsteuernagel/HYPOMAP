
##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree plotting .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

command_args<-commandArgs(TRUE)
param_file = command_args[1]
#param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### For all subsets
##########

subset_names = parameter_list$subset_names#c("HypoNeuron","NonNeuron","Astrocytes","Oligodendrocytes")

for( subname in subset_names){
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
  
  plot_path = paste0(param_set$harmonization_folder_path,"annotation_plots/")
  dir.create(plot_path)
  
  ##########
  ### Load mrtree results
  ##########
  
  # mrtree input
  param_set$clustering_folder = paste0(param_set$harmonization_folder_path,"consensus/") 
  #cluster_matrix_for_mrtree = data.table::fread(paste0(param_set$harmonization_folder_path,param_set$clusters_for_mrtree_file),data.table = F)
  param_set$load_from_file = FALSE
  if(param_set$load_from_file){ # load_from_file = FALSE
    message("Reading clusters from file: ",param_set$clusters_for_mrtree_file)
    cluster_matrix_for_mrtree = data.table::fread(paste0(param_set$harmonization_folder_path,param_set$clusters_for_mrtree_file),data.table = F)
  }else{
    # or merge from multiple files (consensus clusterings from multiple jobs)
    message("Reading clusters from folder and merge: ",param_set$clustering_folder)
    all_cluster_files = list.files(param_set$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(param_set$harmonization_folder_path,"consensus/")
    all_clusterings = lapply(paste0(param_set$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
    cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)
  }
  message("read 1")
  # load mrtree clustering
  mrtree_result = readRDS(paste0(param_set$harmonization_folder_path,param_set$clustering_key_name,"_mrtree_results",".rds"))
  pruned_mrtree_clustering_results = readRDS(param_set$mrtree_result_file)
  #  "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization/human_hypo_neuron_clusters_pruned_mrtree_clustering_results.rds"
  
  mrtree_result_plotting = pruned_mrtree_clustering_results
  mrtree_result_plotting_labelmat = mrtree_result_plotting$labelmat
  
  # check
  table_df = rbind(
    apply(cluster_matrix_for_mrtree,2,function(x){length(unique(x))}),
    apply(mrtree_result$labelmat,2,function(x){length(unique(x))}),
    apply(mrtree_result_plotting_labelmat,2,function(x){length(unique(x))})
  ) %>% as.data.frame()
  colnames(table_df) = paste0("Level_",1:ncol(table_df))
  rownames(table_df) = c("Leiden","Mrtree","PrunedFinal")
  
  write.table(x = table_df,file = paste0(plot_path,"cluster_overview_",subname,".txt"),sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
  
  ##########
  ### plot mrtree results
  ##########
  
  library(hypoMapUtils)
  library(magrittr)
  library(ggplot2)
  library(scales)
  library(treeio)
  library(tidytree)
  library(ggtree)
  
  #########
  edge_color="black"
    
  node_color="black"
    
  mrtree_result_plotting_labelmat = mrtree_result_plotting$labelmat
  anno_df = mrtree_result_plotting_labelmat %>% as.data.frame() %>% tidyr::gather(key="clusterlevel",value="cluster_id")
  
  leaf_level = 6
  edgelist = mrtree_result_plotting$edgelist
  # reduce edgelist to certain level and from and to cols
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level<=as.numeric(leaf_level),1:2]
  edgelist = edgelist[edgelist$to %in% anno_df$cluster_id,] # remove edges/nodes that are not part of anno_df
  
  ## convert to treedata
  # only take
  tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
  tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
  tree_data_tibble$nodesize = 1 # default node size
  tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edgelist) # count children number
  
  # convert back to treedata
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
  #plot circular tree
  label_size = 3
  if(nrow(edgelist[edgelist$level==as.numeric(leaf_level),]) > 50){
    label_size_tip = 2
  }else{
    label_size_tip = label_size 
  }
  vjust_label = -0.5
  circular_tree =ggtree(tree_data,layout = 'circular', branch.length='none',color=edge_color)+
    geom_nodepoint(aes(subset = n_children > 1),color=node_color) + 
    geom_nodelab(ggplot2::aes(x = branch, label = label), size = label_size, vjust = vjust_label, color = "darkred")+ #+
    geom_tiplab(ggplot2::aes(x = branch, label = label), size = label_size_tip, vjust = vjust_label, color = "darkred")
  
  #circular_tree
  ggsave(filename = paste0(plot_path,"circular_tree_",subname,".pdf"), plot = circular_tree, "pdf",dpi=400,width=300,height = 300,units="mm")
  
  ##########
  ### plot mrtree results umap
  ##########
  
  # load object
  seurat_object = readRDS(paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,".rds"))
  
  # add 
  # seurat_object@meta.data = cbind(seurat_object@meta.data, mrtree_result$labelmat)
  seurat_object@meta.data = cbind(seurat_object@meta.data, pruned_mrtree_clustering_results$labelmat)
  #colnames(mrtree_result_plotting_labelmat)
  
  # plot per pruned level
  cluster_cols = colnames(mrtree_result_plotting_labelmat)
  umap_name = paste0("umap_scvi_",subname)
  #umap_name = "umap_scvi_neurons" # comments this out !
  for(col in cluster_cols){
    if(length(unique(seurat_object@meta.data[,col])) < 30){
      label.size = 4
    }else{
      label.size = 2
    }
    clusterplot = DimPlot(seurat_object,group.by = col,raster.dpi = c(2048,2048),pt.size = 1.6,shuffle = TRUE,label = TRUE,raster = TRUE,reduction = umap_name,label.size = label.size)+NoLegend()
    # export plot
    ggsave(filename = paste0(plot_path,subname,"_",col,"_umap.pdf"), plot = clusterplot, "pdf",dpi=400,width=200,height = 200,units="mm")
  }
  
  #pruned_markers = data.table::fread(paste0(param_set$harmonization_folder_path,param_set$new_name_suffix,"_",param_set$clustering_key_name,"_all_markers_all_pruned_.tsv"),data.table = F)

}

##########
### Repeat for combined
##########

param_set = parameter_list
subname = "human_hypo_combined"

# load subset specifc
param_set$harmonization_folder_path = paste0(param_set$harmonization_folder_path,"/",subname,"/")
param_set$new_name_suffix = subname

plot_path = paste0(param_set$harmonization_folder_path,"annotation_plots/")
dir.create(plot_path)

human_hypo_combined = readRDS(paste0(param_set$harmonization_folder_path,subname,".rds"))
# other data
human_hypo_combined_edgelist = data.table::fread(paste0(param_set$harmonization_folder_path,subname,"_","edgelist_mrtree",".txt"),data.table = F)
#human_hypo_combined_edgelist = rbind(data.frame(from = "all",to=unique(human_hypo_combined_edgelist$from[grepl("C0-",human_hypo_combined_edgelist$from)]),isLeaf=FALSE,level),human_hypo_combined_edgelist)

human_hypo_combined_markers_all = data.table::fread(paste0(param_set$harmonization_folder_path,subname,"_","comparisons_all_updated",".txt"),data.table = F)
human_hypo_combined_markers_sibling = data.table::fread(paste0(param_set$harmonization_folder_path,subname,"_","comparisons_siblings_updated",".txt"),data.table = F)

mrtree_result_plotting_labelmat = human_hypo_combined@meta.data[,grepl("C[0-9]+",colnames(human_hypo_combined@meta.data))]

##########
### mrtree
##########

anno_df = mrtree_result_plotting_labelmat %>% as.data.frame() %>% tidyr::gather(key="clusterlevel",value="cluster_id")

leaf_level = 6
edgelist = human_hypo_combined_edgelist
# reduce edgelist to certain level and from and to cols
edgelist$level = as.numeric(edgelist$level)
leaf_level_column = gsub("-[0-9]+","",edgelist$to[length(edgelist$to)])
edgelist = edgelist[edgelist$level<=as.numeric(leaf_level),1:2]
#edgelist = edgelist[edgelist$to %in% anno_df$cluster_id,] # remove edges/nodes that are not part of anno_df

## convert to treedata
# only take
tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
tree_data_tibble$nodesize = 1 # default node size
tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edgelist) # count children number

# convert back to treedata
tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
#plot circular tree
label_size = 3
if(nrow(edgelist[edgelist$level==as.numeric(leaf_level),]) > 50){
  label_size_tip = 2
}else{
  label_size_tip = label_size 
}
vjust_label = -0.5
circular_tree =ggtree(tree_data,layout = 'circular', branch.length='none',color=edge_color)+
  geom_nodepoint(aes(subset = n_children > 1),color=node_color) + 
  geom_nodelab(ggplot2::aes(x = branch, label = label), size = label_size, vjust = vjust_label, color = "darkred")+ #+
  geom_tiplab(ggplot2::aes(x = branch, label = label), size = label_size_tip, vjust = vjust_label, color = "darkred")

circular_tree
#circular_tree
ggsave(filename = paste0(plot_path,"circular_tree_",subname,".pdf"), plot = circular_tree, "pdf",dpi=400,width=300,height = 300,units="mm")

## add a heatmap

# make data for third heatmap with regions

heatmap_data3 = human_hypo_combined@meta.data %>% dplyr::select(Cell_ID,celltype_annotation,!!sym(leaf_level_column)) %>%
  dplyr::group_by(!!sym(leaf_level_column),celltype_annotation) %>% dplyr::count() %>% dplyr::group_by(!!sym(leaf_level_column)) %>%
  dplyr::top_n(n = 1,wt = n) %>% ungroup() %>%
  dplyr::distinct(!!sym(leaf_level_column),celltype_annotation,.keep_all=TRUE) %>% as.data.frame()
heatmap_data3$celltype_annotation[heatmap_data3$celltype_annotation == "P2RX2_OTP"] = "Neurons"
heatmap_matrix3 = as.matrix(heatmap_data3[,"celltype_annotation",drop=F])
rownames(heatmap_matrix3) = heatmap_data3[,leaf_level_column]
#colnames(heatmap_matrix3) = "R"

circular_tree_heat = hypoMapUtils::add_heatmap(circular_tree,heatmap_matrix = heatmap_matrix3,heatmap_colors = getOkabeItoPalette(16))
circular_tree_heat = circular_tree_heat + scale_fill_manual(values = getOkabeItoPalette(length(unique(heatmap_data3$celltype_annotation))))+ theme(legend.text=element_text(size=15))
circular_tree_heat

ggsave(filename = paste0(plot_path,"circular_tree_celltypes_",subname,".pdf"), plot = circular_tree_heat, "pdf",dpi=400,width=300,height = 300,units="mm")

##########
### UMAP plots
##########

# plot per pruned level
cluster_cols = colnames(mrtree_result_plotting_labelmat)
cluster_cols = c(cluster_cols,"celltype_annotation","celltype_status")
umap_name = "umap_scvi_hypo" 
for(col in cluster_cols){
  if(length(unique(human_hypo_combined@meta.data[,col])) < 30){
    label.size = 4
  }else{
    label.size = 2
  }
  clusterplot = DimPlot(human_hypo_combined,group.by = col,raster.dpi = c(2048,2048),pt.size = 1.6,shuffle = TRUE,label = TRUE,raster = TRUE,reduction = umap_name,label.size = label.size)+NoLegend()
  # export plot
  ggsave(filename = paste0(plot_path,subname,"_",col,"_umap.pdf"), plot = clusterplot, "pdf",dpi=400,width=200,height = 200,units="mm")
}

##########
### pre calculate sample contribution
##########

message("cluster qc")

all_mrtree_cols = colnames(mrtree_result_plotting_labelmat)
all_cluster_stats_list = list()
# make overall per cluster stats
for(mrtree_col in all_mrtree_cols){
  all_cluster_stats = human_hypo_combined@meta.data %>%
    dplyr::group_by(!!sym(mrtree_col)) %>% 
    dplyr::mutate(mean_nCount_RNA = mean(nCount_RNA)) %>%
    dplyr::mutate(mean_nFeature_RNA = mean(nFeature_RNA)) %>%
    dplyr::mutate(mean_percent_mt = mean(percent.mt)) %>%
    dplyr::add_count(name="cells_cluster") %>%
    dplyr::group_by(!!sym(mrtree_col),Dataset) %>% dplyr::add_count(name="cells_cluster_dataset") %>%
    dplyr::group_by(!!sym(mrtree_col),Donor_ID) %>% dplyr::add_count(name="cells_cluster_donorID") %>%
    dplyr::group_by(!!sym(mrtree_col),Sample_ID) %>% dplyr::add_count(name="cells_cluster_sampleID") %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!sym(mrtree_col),cells_cluster,Dataset,cells_cluster_dataset,Donor_ID,cells_cluster_donorID,Sample_ID,cells_cluster_sampleID,mean_nCount_RNA,mean_nFeature_RNA,mean_percent_mt) %>%
    dplyr::mutate(pct_cluster_dataset = round(cells_cluster_dataset / cells_cluster * 100,3),
                  pct_cluster_donorID = round(cells_cluster_donorID / cells_cluster * 100,3),
                  pct_cluster_sampleID = round(cells_cluster_sampleID / cells_cluster * 100,3))
  colnames(all_cluster_stats)[1] = "cluster"
  # get the max Donor
  donor_pct_max = all_cluster_stats %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::distinct(cluster,Donor_ID,pct_cluster_donorID) %>% 
    dplyr::slice_max(order_by = pct_cluster_donorID,n = 1,with_ties = FALSE) %>%
    dplyr::select(cluster,pct_donorMax = pct_cluster_donorID,donorMax = Donor_ID)
  # get the max sample
  sample_pct_max = all_cluster_stats %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::distinct(cluster,Sample_ID,pct_cluster_sampleID) %>% 
    dplyr::slice_max(order_by = pct_cluster_sampleID,n = 1,with_ties = FALSE) %>%
    dplyr::select(cluster,pct_sampleMax = pct_cluster_sampleID,sampleMax = Sample_ID)
  # make a continous average df
  avg_stats = all_cluster_stats %>% 
    dplyr::distinct(cluster,mean_nCount_RNA,mean_nFeature_RNA,mean_percent_mt)
  
  # make the final overview DF
  cluster_qc_overview = all_cluster_stats %>% 
    dplyr::distinct(cluster,cells_cluster) %>%
    dplyr::left_join(donor_pct_max,by="cluster") %>%
    dplyr::left_join(sample_pct_max,by="cluster") 
  
  all_cluster_stats_list[[mrtree_col]] = cluster_qc_overview
}
cluster_qc_overview = do.call(rbind,all_cluster_stats_list) %>% as.data.frame()

write.table(x = cluster_qc_overview,file = paste0(plot_path,"cluster_qc_overview",".txt"),sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)


message(Sys.time(),": Complete" )



