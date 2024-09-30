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

## seurat
human_processed = readRDS(paste0(opts$data_path,"HumanNucSeq_processed_filtered.rds"))

# make folder for reults
plot_path = paste0(opts$out_path,"06_keygenes/")
system(paste0("mkdir -p ",plot_path))

##########
### plot key genes
##########

rasterize_px = 1024
seurat_pt_size = 2.2


# key genes non neurons
p <- FeaturePlot(human_processed,features = KeyGenes_nonNeuron, combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
KeyGenes_nonNeuron_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 4)
#KeyGenes_nonNeuron_feature_plots

# key genes 1
p <- FeaturePlot(human_processed,features = KeyGenes[1:15], combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
KeyGenes_1_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 5)
#KeyGenes_1_feature_plots

# key genes 2
p <- FeaturePlot(human_processed,features = KeyGenes[16:length(KeyGenes)], combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
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


##########
### Pnoc plots
##########

pnoc_genes = c("PNOC","SST","NPY","CRABP1","HTR3B","GAL")

p_noc <- FeaturePlot(human_processed,features = pnoc_genes, combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi")
for(i in 1:length(p_noc)) {
  p_noc[[i]] <- p_noc[[i]] + NoLegend() + NoAxes()
}
pnoc_feature_plots = cowplot::plot_grid(plotlist = p_noc,ncol = 3)
pnoc_feature_plots

ggsave(filename = paste0(plot_path,"pnoc_genes.pdf"),
       plot = pnoc_feature_plots, "pdf",dpi=400,width=300,height = 200,units="mm")

##########
### Pnoc cpmbinations
##########

# pnoc cominations
pnoc_combo_list = list(
  c("PNOC","SST"),
  c("PNOC","CRABP1"),
  c("PNOC","CRABP1","SST"),
  c("PNOC","SST","NPY")
)

for(comb in pnoc_combo_list){
  comb_name = paste0(comb,collapse = "_")
  human_processed@meta.data[,comb_name] = as.numeric(scUtils::CalculateMultScore(human_processed,features = comb))
}

comb_names = sapply(pnoc_combo_list,paste0,collapse = "_")

p_noc_comb <- FeaturePlot(human_processed,features = comb_names, combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction="umap_scvi",cols = c("lightgrey","#c96410"))
for(i in 1:length(p_noc_comb)) {
  p_noc_comb[[i]] <- p_noc_comb[[i]] + NoLegend() + NoAxes()
}
pnoc_comb_feature_plots = cowplot::plot_grid(plotlist = p_noc_comb,ncol = 2)
pnoc_comb_feature_plots

ggsave(filename = paste0(plot_path,"pnoc_comb.pdf"),
       plot = pnoc_comb_feature_plots, "pdf",dpi=400,width=200,height = 200,units="mm")

##########
### Pnoc stats
##########

# per cluster stats
per_cluster_pct = scUtils::gene_pct_cluster(seurat_object = human_processed,genes = c(pnoc_genes,comb_names),col_name = "SNN_scvi_res.4")
per_cluster_pct$cluster = rownames(per_cluster_pct)
per_cluster_pct = per_cluster_pct[,c(ncol(per_cluster_pct),1:(ncol(per_cluster_pct)-1))]
per_cluster_pct = per_cluster_pct %>% dplyr::arrange(desc(PNOC_SST))

# total stats
per_annotation_pct = scUtils::gene_pct_cluster(seurat_object = human_processed,genes = c(pnoc_genes,comb_names),col_name = "updated_annotation")
per_annotation_pct = per_annotation_pct["Neurons",]
per_annotation_pct$cluster = "All neurons"
per_annotation_pct = per_annotation_pct[,c(ncol(per_annotation_pct),1:(ncol(per_annotation_pct)-1))]

per_cluster_pct_final = dplyr::bind_rows(per_annotation_pct,per_cluster_pct %>% dplyr::top_n(n = 10,wt = PNOC_SST))
per_cluster_pct_final[,2:ncol(per_cluster_pct_final)] = apply(per_cluster_pct_final[,2:ncol(per_cluster_pct_final)],2,function(x){round(x*100,2)})

data.table::fwrite(per_cluster_pct_final,file = paste0(plot_path,"pnoc_pcts.txt"),sep="\t")


