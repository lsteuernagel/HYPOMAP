
##########
### Load parameters and packages
##########
message("-----",Sys.time(),": Load parameters and packages ")

library(magrittr)
library(scUtils)
library(dplyr)
library(Matrix)
library(Seurat)
source("utility_functions.R")

opts <- workflow_options(project = "HuHy")

# load seurat
human_processed = readRDS(paste0(opts$data_path,"HumanNucSeq_processed.rds"))
human_processed@meta.data$RNA_snn_res.4 = as.character(human_processed@meta.data$RNA_snn_res.4)
Idents(human_processed) = "RNA_snn_res.4"

# if we want to add some seurat reduction
#path_to_seurat_reduction_to_add = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/john_processed/06.scrna.all.data.processed_UMAPonly.txt" # won't work
path_to_seurat_reduction_to_add = NULL

# for doublet flagging
min_ratio = 2

# remove_manual_clusters
#### CAREFUL THIS CAN CHANGE WHEN RE RUNNING 
# leave empty in automized pipeline
remove_manual_clusters =("35")

# make folder for reults
plot_path = paste0(opts$out_path,"05_filtering/")
system(paste0("mkdir -p ",plot_path))

##########
### signature enrichments
##########

message("-----",Sys.time(),": Annotate with signatures ")

mouseHypoMap_celltype_signatures = jsonlite::read_json("data/mouseHypoMap_celltype_signaturs.json")
mouseHypoMap_celltype_signatures = sapply(mouseHypoMap_celltype_signatures,unlist)

# add signature module scores
human_processed = Seurat::AddModuleScore(human_processed,features = mouseHypoMap_celltype_signatures,name="Signature")
colnames(human_processed@meta.data)[grepl("Signature",colnames(human_processed@meta.data))] = names(mouseHypoMap_celltype_signatures)
# make table with all signature scores
signatures_per_cell = human_processed@meta.data %>% 
  dplyr::select(c("Cell_ID","SNN_scvi_res.4", names(mouseHypoMap_celltype_signatures)))
# save
data.table::fwrite(signatures_per_cell,paste0(opts$data_path,"human_nucseq_processed_signatureScores.txt"),sep="\t")


# stats per preliminary cluster
signature_per_cluster = human_processed@meta.data %>% 
  dplyr::select(c("SNN_scvi_res.4", names(mouseHypoMap_celltype_signatures)))  %>%
  tidyr::gather(key="celltype",value="score",-SNN_scvi_res.4) %>%
  dplyr::group_by(SNN_scvi_res.4,celltype) %>%
  dplyr::summarise(median_score = median(score)) %>%
  dplyr::arrange(desc(median_score)) %>%
  dplyr::ungroup()

signature_per_cluster_top = signature_per_cluster %>%   
  dplyr::group_by(SNN_scvi_res.4) %>% dplyr::top_n(n = 1,wt = median_score)

signature_per_cluster_top2_ratio = signature_per_cluster %>%   
  dplyr::group_by(SNN_scvi_res.4) %>% 
  dplyr::top_n(n = 2,wt = median_score) %>%
  dplyr::summarise(ratio = median_score[1] / (median_score[2]+0.02))

signature_per_cluster_top3 = signature_per_cluster %>%   
  dplyr::group_by(SNN_scvi_res.4) %>% 
  dplyr::top_n(n = 3,wt = median_score)

# shox2 pct
shox2_pcts = scUtils::gene_pct_cluster(human_processed,genes = "SHOX2",col_name = "SNN_scvi_res.4",return_long = TRUE)

##########
### Find problematic clusters and more Doublet clusters
##########

message("-----",Sys.time(),": Identiy and remove doublets ")

# Select clusters below X ratio or with high Shox
low_qual_clusters = unique(c(as.character(signature_per_cluster_top2_ratio$SNN_scvi_res.4[signature_per_cluster_top2_ratio$ratio < min_ratio]),shox2_pcts$group[shox2_pcts$pct >= 0.3]))

# find updated names for these clusters
result_vec = character()
for(cluster in low_qual_clusters){
  top_names = signature_per_cluster_top3$celltype[signature_per_cluster_top3$SNN_scvi_res.4 == cluster & signature_per_cluster_top3$median_score > 0.05]
  # if all immune: just set to immune
  if(all(grepl("Immune",top_names))){
    result_vec[cluster] = top_names[1]
  }
  # if endo + mural: rename manually
  if("Endothelial" %in% top_names & "Dcn.Fibroblasts" %in% top_names & "Mural" %in% top_names){
    result_vec[cluster] = "Endo+Mural+Fibro"
  }
  # if two oligo: take larger oligo
  if(all(grepl("Oligo|OPC",top_names))){
    result_vec[cluster] = top_names[1]
  }
  # if neuron / oligo / OPC / astro combo
  # mark as doublet
  if("Neurons" %in% top_names & any(grepl("Oligo",top_names))){
    result_vec[cluster] = "Doublet"
  }
  if("Astrocytes" %in% top_names & any(grepl("Oligo|OPC",top_names))){
    result_vec[cluster] = "Doublet"
  }
  if("Astrocytes" %in% top_names & "Neurons" %in% top_names){
    result_vec[cluster] = "Doublet"
  }
  if(any(grepl("Immune",top_names)) & any(grepl("Oligo|OPC",top_names))){
    result_vec[cluster] = "Doublet"
  }
  # if neuron and OPC --> keep as neuorns
  if("Neurons" %in% top_names & any(grepl("OPC",top_names))){
    result_vec[cluster] = "Unknown_Neurons"
  }
  # rename shox2 neurons
  if(shox2_pcts$pct[shox2_pcts$group == cluster] >= 0.1){
    result_vec[cluster] = "Shox2_Neurons"
  }
  if(shox2_pcts$pct[shox2_pcts$group == cluster] >= 0.1){
    result_vec[cluster] = "Shox2_Neurons"
  }
  
  # if there are no top names --> set to unkowmn
  if(length(top_names) < 1){
    result_vec[cluster] = "Unknown"
  }else{
    # tanys are unclear but if there is a cluster
    if(top_names[1] == "Tanycytes"){
      result_vec[cluster] = top_names[1]
    }
  }
  
}

result_vec[is.na(result_vec)] = "Doublet"
result_df = data.frame(cluster = names(result_vec),updated_annotation = result_vec)

# get updated annotations
updated_annotation = dplyr::left_join(signature_per_cluster_top,result_df,by=c("SNN_scvi_res.4"="cluster"))  %>% as.data.frame()
# take top 1 for non problematic ones
updated_annotation$updated_annotation[is.na(updated_annotation$updated_annotation)] = updated_annotation$celltype[is.na(updated_annotation$updated_annotation)]
# add to metadata
updated_annotation = updated_annotation[,c("SNN_scvi_res.4","updated_annotation")]
tmp_meta = dplyr::left_join(human_processed@meta.data[,!colnames(human_processed@meta.data) %in% c(names(mouseHypoMap_celltype_signatures),"updated_annotation")],updated_annotation,by=c("SNN_scvi_res.4"="SNN_scvi_res.4"))
# add to metadata
rownames(tmp_meta) = tmp_meta$Cell_ID
human_processed@meta.data = tmp_meta

##########
### QC plots
##########

# make plots
clusterplot = DimPlot(human_processed,group.by = "SNN_scvi_res.4",raster.dpi = c(2048,2048),pt.size = 2.2,shuffle = TRUE,label = TRUE,reduction = "umap_scvi")+NoLegend()

annotation_plot = DimPlot(human_processed,group.by = "updated_annotation",raster.dpi = c(2048,2048),pt.size = 2.2,shuffle = TRUE,label = TRUE,repel = TRUE,reduction = "umap_scvi")+NoLegend()
# export plots
ggsave(filename = paste0(plot_path,"clusters_preFiltering.pdf"),
       plot = clusterplot, "pdf",dpi=400,width=200,height = 200,units="mm")
# export plots
ggsave(filename = paste0(plot_path,"annotation_preFiltering.pdf"),
       plot = annotation_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

##########
### remove doublet & problematic clusters
##########


remove_cells = human_processed@meta.data$Cell_ID[human_processed@meta.data$updated_annotation == "Doublet"]
remove_cells = c(remove_cells, human_processed@meta.data$Cell_ID[human_processed@meta.data[,"SNN_scvi_res.4"] %in% remove_manual_clusters])
remove_cells =unique(remove_cells)
keep_cells= human_processed@meta.data$Cell_ID[! human_processed@meta.data$Cell_ID %in% remove_cells]
# make a plot
remove_cells_plot = DimPlot(human_processed,group.by = "SNN_scvi_res.4",raster.dpi = c(2048,2048),pt.size = 2.2,shuffle = TRUE,label = TRUE,reduction = "umap_scvi",cells.highlight = remove_cells)+NoLegend()

# subset object
human_processed = subset(human_processed,cells = keep_cells)

# Want to remove shox as well ?

# plot again
clusterplot2 = DimPlot(human_processed,group.by = "SNN_scvi_res.4",raster.dpi = c(2048,2048),pt.size = 2.2,shuffle = TRUE,label = TRUE,reduction = "umap_scvi")+NoLegend()
annotation_plot2 = DimPlot(human_processed,group.by = "updated_annotation",raster.dpi = c(2048,2048),pt.size = 2.2,shuffle = TRUE,label = TRUE,repel = TRUE,reduction = "umap_scvi")+NoLegend()
# export plots
ggsave(filename = paste0(plot_path,"clusters_postFiltering.pdf"),
       plot = clusterplot2, "pdf",dpi=400,width=200,height = 200,units="mm")
# export plots
ggsave(filename = paste0(plot_path,"annotation_postFiltering.pdf"),
       plot = annotation_plot2, "pdf",dpi=400,width=200,height = 200,units="mm")
# export plots
ggsave(filename = paste0(plot_path,"cells_to_remove.pdf"),
       plot = remove_cells_plot, "pdf",dpi=400,width=200,height = 200,units="mm")


##########
### add seurat integration reduction (umap) for plotting
##########
if(!is.null(path_to_seurat_reduction_to_add)){
  message("-----",Sys.time(),": add seurat integration ")
  seurat_reduction = data.table::fread(path_to_seurat_reduction_to_add,data.table = FALSE)
  rownames(seurat_reduction) = seurat_reduction[,1]
  seurat_reduction = seurat_reduction[,2:ncol(seurat_reduction)]
  
  # make dim red
  dimred <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(seurat_reduction),
    stdev = as.numeric(apply(seurat_reduction, 2, stats::sd)),
    assay = "RNA",
    key = "umap_seurat"
  )
  # add  to object
  seurat_object@reductions[["umap_seurat"]] = dimred
  
  ## add plots with seurat
  
}

##########
### Clean seurat
##########


##########
### save data
##########

message("-----",Sys.time(),": Save data ")


## save filtered rds
saveRDS(human_processed,file = paste0(opts$data_path,"HumanNucSeq_processed_filtered.rds"))

## some additional plots ? 

#FeaturePlot(human_processed,features = "AGRP",order = TRUE,reduction = "umap_scvi")

message("-----",Sys.time(),": Complete ")
