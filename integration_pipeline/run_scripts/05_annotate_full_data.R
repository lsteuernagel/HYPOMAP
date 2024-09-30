##########
### Load parameters and packages
##########

message("-----",Sys.time(),": Load packages and object")

require(tidyverse)
require(Seurat)
require(Matrix)
library(rlang)

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

# load seurat
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))

# load clusters
message("-----",Sys.time(),": Load and clear clusters ")
hypoMap_test_initial_leiden_clustering = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,parameter_list$additional_clustering_suffix,"_leiden_clustering.txt"),data.table = F)

# clean clusters from singlets
# neighbor_nn = as.Neighbor(harmonized_seurat_object@graphs$NN_scvi)
# nn_idx = neighbor_nn@nn.idx
# hypoMap_test_initial_leiden_clustering_clear = apply(hypoMap_test_initial_leiden_clustering[,2:ncol(hypoMap_test_initial_leiden_clustering)],2,clear_clustering,min_cells = parameter_list$min_cells_valid,nn_idx = nn_idx) %>% as.data.frame()
# hypoMap_test_initial_leiden_clustering_clear$Cell_ID = hypoMap_test_initial_leiden_clustering$Cell_ID
# data.table::fwrite(hypoMap_test_initial_leiden_clustering_clear,paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering_clear.txt"),sep="\t")
# hypoMap_test_initial_leiden_clustering_clear = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_clustering_clear.txt"),data.table = F)

hypoMap_test_initial_leiden_clustering_clear = hypoMap_test_initial_leiden_clustering

# add to seurat
temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,hypoMap_test_initial_leiden_clustering_clear,by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta

##########
### Define resolution
##########

# use column with highest resolution
#cluster_column = colnames(hypoMap_test_initial_leiden_clustering)[ncol(hypoMap_test_initial_leiden_clustering)]
cluster_column = "leiden_clusters_12" #colnames(hypoMap_test_initial_leiden_clustering)[5]

##########
### signature enrichments
##########

message("-----",Sys.time(),": Annotate with signatures ")

mouseHypoMap_celltype_signatures = jsonlite::read_json("data/mouseHypoMap_celltype_signaturs.json")
mouseHypoMap_celltype_signatures = sapply(mouseHypoMap_celltype_signatures,unlist)

## manually add some other clusters
mouseHypoMap_celltype_signatures$SHOX2_Thalamus = c("SHOX2","CASQ2","PTPN3","GRID2IP","SYN2","SYT13","ELAVL4")
mouseHypoMap_celltype_signatures$LHX6_CHST9_Unknown = c("RBFOX3","SYN2","SYT13","ELAVL4","LHX6","ANO1","MAF","CHST9","IGF1")
mouseHypoMap_celltype_signatures$LHX6_SCARA5_Unknown = c("RBFOX3","SYN2","SYT13","ELAVL4","LHX6","MAF","SCARA5","IGF1")
mouseHypoMap_celltype_signatures$PKP2_IGF1_Unknown = c("RBFOX3","SYN2","SYT13","ELAVL4","ANO1","MAF","PKP2","CXCL14","IGF1","EGFR","PROX1")
mouseHypoMap_celltype_signatures$GALNT5_PAX6_Unknown = c("RBFOX3","SYN2","GALNT5","PAX6","MSR1","CDH15","NEUROD1","CRTAM")
mouseHypoMap_celltype_signatures$SLC17A7_TBR1_Unknown = c("RBFOX3","SYN2","SYT13","ELAVL4","SLC17A7","TBR1","MCOLN3","GLP2R","SATB2","SCHLAP1","LY86-AS1","NRGN")
mouseHypoMap_celltype_signatures$DMBX1_Midbrain = c("RBFOX3","SYN2","SYT13","ELAVL4","OTX2","DMBX1","OTX2-AS1","GATA3","CASR","CPLX3","LINC01210")

# add signature module scores
harmonized_seurat_object = Seurat::AddModuleScore(harmonized_seurat_object,features = mouseHypoMap_celltype_signatures,name="Signature")
colnames(harmonized_seurat_object@meta.data)[grepl("Signature",colnames(harmonized_seurat_object@meta.data))] = names(mouseHypoMap_celltype_signatures)
# make table with all signature scores
signatures_per_cell = harmonized_seurat_object@meta.data %>%
  dplyr::select(c("Cell_ID",cluster_column, names(mouseHypoMap_celltype_signatures)))
# save
data.table::fwrite(signatures_per_cell,paste0(parameter_list$harmonization_folder_path,"human_nucseq_processed_signatureScores.txt"),sep="\t")

# stats per preliminary cluster
signature_per_cluster = harmonized_seurat_object@meta.data %>%
  dplyr::select(c(cluster_column, names(mouseHypoMap_celltype_signatures)))  %>%
  tidyr::gather(key="celltype",value="score", - !!sym(cluster_column)) %>%
  dplyr::group_by(!!sym(cluster_column),celltype) %>%
  dplyr::summarise(median_score = median(score)) %>%
  dplyr::arrange(desc(median_score)) %>%
  dplyr::ungroup()

signature_per_cluster_top = signature_per_cluster %>%
  dplyr::group_by(!!sym(cluster_column)) %>% dplyr::top_n(n = 1,wt = median_score)

# save as well:
data.table::fwrite(signature_per_cluster_top,paste0(parameter_list$harmonization_folder_path,"human_nucseq_processed_signatures_Top.txt"),sep="\t")
#signature_per_cluster_top = data.table::fread(paste0(parameter_list$harmonization_folder_path,"human_nucseq_processed_signatures_Top.txt"),data.table=F)

# signature_per_cluster_top2_ratio = signature_per_cluster %>%
#   dplyr::group_by(!!sym(low_res_clusters)) %>%
#   dplyr::top_n(n = 2,wt = median_score) %>%
#   dplyr::summarise(ratio = median_score[1] / (median_score[2]+0.02))
#
# signature_per_cluster_top3 = signature_per_cluster %>%
#   dplyr::group_by(!!sym(low_res_clusters)) %>%
#   dplyr::top_n(n = 3,wt = median_score)

message("-----",Sys.time(),": Add signatures to object")

## add to object
signature_per_cluster_top_join = signature_per_cluster_top %>% dplyr::select(celltype_annotation = celltype, !!sym(cluster_column) )
temp_meta = dplyr::left_join(harmonized_seurat_object@meta.data,signature_per_cluster_top_join,by=cluster_column)
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta

##########
### save results
##########

message("-----",Sys.time(),": Saving annotated object")

file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated")

# save data to rds
saveRDS(harmonized_seurat_object,paste0(file_name_prefix,".rds"))

##########
### QC metadata plots
##########
# 
# message("read")
# file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated")
# harmonized_seurat_object = readRDS(paste0(file_name_prefix,".rds"))

message("-----",Sys.time(),": Plot preliminary results")

plot_path = paste0(parameter_list$harmonization_folder_path,"annotation_plots/")
system(paste0("mkdir -p ",paste0(plot_path)))

# make a dataset plot
datasetplot = DimPlot(harmonized_seurat_object,group.by = "Dataset",raster.dpi = c(2048,2048),pt.size = 1.4,shuffle = TRUE,label = FALSE,raster = TRUE,reduction = "umap_scvi")+NoLegend()
# export plot
ggsave(filename = paste0(plot_path,"dataset","_umap.pdf"),plot = datasetplot, "pdf",dpi=400,width=230,height = 200,units="mm")

# make plots for all leiden clusters
leiden_colnames = colnames(hypoMap_test_initial_leiden_clustering)[2:ncol(hypoMap_test_initial_leiden_clustering)]
for(col in leiden_colnames){
  clusterplot = DimPlot(harmonized_seurat_object,group.by = col,raster.dpi = c(2048,2048),pt.size = 1.4,shuffle = TRUE,label = TRUE,raster = TRUE,reduction = "umap_scvi",label.size = 2)+NoLegend()
  # export plot
  ggsave(filename = paste0(plot_path,col,"_umap.pdf"), plot = clusterplot, "pdf",dpi=400,width=200,height = 200,units="mm")
}

# plot annotation
annotation_plot = DimPlot(harmonized_seurat_object,group.by = "celltype_annotation",raster.dpi = c(2048,2048),pt.size = 1.4,shuffle = TRUE,label = TRUE,repel = TRUE,raster = TRUE,reduction = "umap_scvi")+NoLegend()
# export plot
ggsave(filename = paste0(plot_path,"annotation_umap.pdf"),plot = annotation_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

##########
### plot key genes
##########

message("-----",Sys.time(),": Plot key genes ")

keyGenes = sapply(jsonlite::read_json("data/human_key_genes.json"),FUN = unlist)

rasterize_px = 1024
seurat_pt_size = 1.1

# key genes non neurons
p <- FeaturePlot(harmonized_seurat_object,features = keyGenes$keyGenes_nonneuron, combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
KeyGenes_nonNeuron_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 4)
#KeyGenes_nonNeuron_feature_plots

# key genes 1
p <- FeaturePlot(harmonized_seurat_object,features = keyGenes$keyGenes_neuron[1:15], combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
KeyGenes_1_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 5)
#KeyGenes_1_feature_plots

# key genes 2
p <- FeaturePlot(harmonized_seurat_object,features = keyGenes$keyGenes_neuron[16:length(keyGenes$keyGenes_neuron)], combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
KeyGenes_2_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 5)
#KeyGenes_2_feature_plots


# save plots
ggsave(filename = paste0(plot_path,"key_non_neuron_genes.pdf"),
       plot = KeyGenes_nonNeuron_feature_plots, "pdf",dpi=400,width=400,height = 300,units="mm")
ggsave(filename = paste0(plot_path,"key_neuron_genes_1.pdf"),
       plot = KeyGenes_1_feature_plots, "pdf",dpi=400,width=500,height = 300,units="mm")
ggsave(filename = paste0(plot_path,"key_neuron_genes_2.pdf"),
       plot = KeyGenes_2_feature_plots, "pdf",dpi=400,width=500,height = 300,units="mm")


message("-----",Sys.time(),": Completed annotation")

