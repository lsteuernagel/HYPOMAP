
# sbatch run_slurm.sh ~/Documents/r_scvi_v3.simg 04_integrate_seurat.R

##########
### Load libs
##########

message("-----",Sys.time(),"Load libs")

options(future.globals.maxSize= 200000*1024^2) # in MB

library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(scDblFinder)
library(BiocParallel)
library(scater)
library(scran)
library(qs)
library(outliers)
library(scuttle)
library(future)
# scCustomize
# SCpubr
# COSG

source("utility_functions.R")

##########
### Set opts and load seurat
##########

message("-----",Sys.time()," Load data")

opts <- workflow_options(project = "HuHy", out_path = "test/")
opts$cores <- 56
raw_data_path = opts$data_path
set.seed(opts$seed)
#future::plan("multicore", workers = opts$cores)
future::plan("sequential") # currently not using parallel future because it probably leads to crashes on larger datasets!

scrna <- readRDS(paste0(raw_data_path,"human_hypo_raw_filtered.rds"))

print(scrna)

scrna.list <- SplitObject(scrna, split.by = "batch")
# rm(scrna)
# gc()

##########
### QC code ?
##########

# write_dims(scrnalist = scrna.list,
#            out_filename = "04.batchdims.postQC")
# 
# QCcounts.batchDF <- counts_keygenes(seurat_list = scrna.list,
#                                     gene_list = KeyGenes)
# 
# write.csv(QCcounts.batchDF, paste0(opts$out_path,"results/04.QCcounts.batchDF.csv"))

##########
### Processing batches
##########

message("-----",Sys.time()," Processing batches")

# batch check ----

message("-----",Sys.time()," SCTransforming")

scrna.list <- lapply(X = scrna.list, FUN = function(x) {
  print(x[["batch"]][1,])
  x <- Seurat::SCTransform(x, 
                   vst.flavor = 'v1',
                   vars.to.regress = "percent.mt",
                   verbose = 0.5)
})

scrna.list <- lapply(X = scrna.list, FUN = function(x) {
  x <- RunPCA(x, verbose = FALSE, seed.use = opts$seed)
})

message("-----",Sys.time()," PCA done ")

##########
### Not sure what happens here
##########

# scrna.list <- lapply(X = scrna.list, FUN = function(x) {
#   x <- FindNeighbors(x, dims = opts$batch_nPCS)
# })

# scrna.list <- lapply(X = scrna.list, FUN = function(x) {
#   x <- RunUMAP(x, dims = opts$batch_nPCS, 
#                densmap = TRUE, 
#                dens.lambda = 1,
#                verbose = FALSE, 
#                umap.method = "umap-learn") 
# })

# scrna.list <- lapply(X = scrna.list, FUN = function(x) {
#   x <- FindClusters(x, resolution = opts$batch_RES, algorithm = 4, method = "igraph")
# })

#Test this-----
# scrna.list <- lapply(X = scrna.list, FUN = function(x) {
#   x <- reorder_clusters(object = x, prefix ="SCT", nPCS = opts$batch_nPCS)
# })

# scrna.list <- lapply(X = scrna.list, FUN = function(x) {
#   Idents(x) <- x[['qclust_batch']];x
#   x <- reorder_batch_clusters(object = x, 
#                               nPCS = opts$batch_nPCS, 
#                               assay = "RNA", 
#                               cluster_name = "qclust_batch")
#   
# })

# qsave(scrna.list, paste0(opts$out_path,"results/04.scrna.SCT"), 
#       preset = "balanced", 
#       nthreads = opts$cores)
# 
# message(paste0("Saved ", opts$out_path,"results/04.scrna.SCT"))



#scrna.list <- qread(paste0(opts$out_path,"results/04.scrna.SCT"),  nthreads = opts$cores)

##########
### Batch plots
##########
# 
# pdf(paste0(opts$out_path,"plots/04.Batch_plots.pdf"), width = opts$width, height = opts$height)
# 
# for (name in names(scrna.list)) {
#   
#   nsamples <- length(unique(scrna.list[[name]]@meta.data$sample))
#   num_cols <- ceiling(sqrt(nsamples))
#   
#   top40raw <- head(VariableFeatures(scrna.list[[name]]), 40)
#   
#   top40 <- remove_pseudogenes(top40raw)
#   
#   p1 <- LabelPoints(plot = VariableFeaturePlot(scrna.list[[name]], 
#                                                selection.method = "sct"), 
#                     points = top40, 
#                     repel = TRUE, 
#                     max.overlaps = 30) & 
#     NoLegend()
#   
#   p2a <- Stacked_VlnPlot(scrna.list[[name]],
#                          features = c(markers$neurotypes, 
#                                       BroadMarkers[-1]),
#                          assay = "RNA", 
#                          x_lab_rotate = TRUE, 
#                          color_seed = opts$seed) & 
#     theme(axis.text.x = element_blank())
#   
#   p2b <- Stacked_VlnPlot(scrna.list[[name]],
#                          features = c("nCount_RNA", 
#                                       "nFeature_RNA"),
#                          assay = "RNA", 
#                          x_lab_rotate = TRUE,
#                          color_seed = opts$seed) & 
#     scale_y_continuous(trans = "log2", 
#                        breaks = scales::breaks_log(n = 5, base = 2)) & 
#     theme(axis.text.x = element_blank()) &
#     theme(axis.text.y = element_text(size = 6))
#   
#   p2c <- Stacked_VlnPlot(scrna.list[[name]],
#                          features = c("percent.mt"),
#                          assay = "RNA", 
#                          x_lab_rotate = TRUE,
#                          color_seed = opts$seed) &
#     scale_y_continuous(trans = "log1p", breaks = c(0.2, 1:5, 7, 10, 15)) &
#     theme(axis.text.y = element_text(size = 6))
#   
#   #learn to use the area argument for this as it's more elegant and ?programmable
#   design <- "
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   111111111
#   222222222
#   222222222
#   333333333
# "
#   
#   p2 <- wrap_plots(p2a, p2b, p2c, ncol = 1, design = design)
#   rm(p2a, p2b, p2c)
#   
#   p4 <- do_DimPlot(scrna.list[[name]],
#                    split.by = "sample", 
#                    ncol = num_cols) & 
#     NoLegend() & 
#     theme(aspect.ratio=1)
#   
#   p5 <- wrap_plots(p1, p2) + 
#     plot_annotation(title = paste("Batch", name)) 
#   
#   print(p5)
#   
#   print(p4)
#   
#   #TODO try seurat's CombinePlots function
#   
#   print(
#     FeaturePlot(scrna.list[[name]], 
#                 features = BroadMarkers, 
#                 order = TRUE, 
#                 #              pt.size = opts$SIZE, 
#                 raster = TRUE, split.by = "sample", by.col = FALSE) +
#       plot_annotation(title = paste("Batch", name)) & 
#       theme(aspect.ratio = 1) # & NoLegend()
#   ) 
#   
#   for (sample_name in unique(scrna.list[[name]]$sample)) {
#     print(paste0("Doing heatmap for ", sample_name))
#     
#     markers_cg <- cosg(scrna.list[[name]][,scrna.list[[name]]$sample == sample_name], assay = "RNA")
#     write.csv(markers_cg, paste0(opts$out_path, 
#                                  "results/cosg_markers_RNA_batch_", 
#                                  name, 
#                                  "_",
#                                  sample_name, 
#                                  ".csv"))
#     
#     marker_genes_COSG <- unique(unlist(c("SYT1", "PLP1", "AQP4", "CSF1R", markers_cg$names[1:4,])))
#     
#     scrna <- ScaleData(scrna.list[[name]][,scrna.list[[name]]$sample == sample_name],
#                        features = marker_genes_COSG,
#                        assay = "RNA") 
#     
#     p1 <- DoHeatmap(
#       subset(scrna, downsample = 100),
#       features = marker_genes_COSG,
#       assay = "RNA",
#       slot = "scale.data",
#       angle = 90, 
#       # lines.width = 4,
#       # size = 3, 
#       group.bar.height = 0, 
#       hjust = 0.4,
#       raster = FALSE,
#       label = TRUE) + 
#       NoLegend() # +
#     # theme(axis.text.y = element_text(size = 5)) +
#     # theme(axis.text.x.top = element_text(size = 8))
#     
#     print(p1 + plot_annotation(title = paste("Batch", name, 
#                                              "Sample", sample_name)))
#     
#     rm(scrna, p1)
#   }
# }
# 
# dev.off()
# 
# rm(p2, markers_cg, p4, p5); gc()

#Integration----
#scrna.list <- qread(paste0(opts$out_path,"results/04.scrna.SCT"),  nthreads = opts$cores)

##########
### Integration
##########

message("-----",Sys.time()," SelectIntegrationFeatures ")

Ifeatures <- Seurat::SelectIntegrationFeatures(object.list = scrna.list, nfeatures = opts$nfeatures)

# print("XY genes in Ifeatures")
# Ifeatures[Ifeatures %in% c(.Y_genes_raw$hgnc_symbol, "XIST")]
# #consider removing these genes from the integration list?
# 
# print("integration features in KeyGenes")
# KeyGenes[KeyGenes %in% Ifeatures] 

var_genes <- c(purrr::map(scrna.list, Seurat::VariableFeatures))
# print("variable genes in KeyGenes by sample")
# purrr::map(var_genes, intersect, KeyGenes)

scrna.list <- PrepSCTIntegration(
  object.list = scrna.list, 
  anchor.features = Ifeatures)

### integration------
future::plan("sequential") # currently not using parallel future because it probably leads to crashes on larger datasets!
#future::plan("multicore", workers = 4) 

scrna.list <- lapply(X = scrna.list, 
                     FUN = RunPCA, 
                     features = Ifeatures, 
                     seed.use = opts$seed, 
                     verbose = FALSE)

##TODO, switch to something like
#refs <- which(names(scrna.list) == "A3020")

######### THIS IS PROBLEMATIC !!!!!!!
refs <- c(21, 8, 19)

if (length(scrna.list) >= max(refs)) {
  references <- refs
} else {
  references <- NULL
}
message("-----",Sys.time()," FindIntegrationAnchors ")

huhy.anchors <- FindIntegrationAnchors(
  object.list = scrna.list,
  normalization.method = "SCT",
  anchor.features = Ifeatures,
  reduction = "rpca",
  k.anchor = 10, #alter as needed
  reference = references,
  eps = 1)

rm(scrna.list)
gc()

message("-----",Sys.time()," IntegrateData ")

scrna <- IntegrateData(
  anchorset = huhy.anchors,
  dims = opts$integ_nPCS,
  normalization.method = "SCT",
  eps = 1)

##########
### Integration PCA & UMAP
##########

message("-----",Sys.time()," Run PCA+UMAP")

scrna <- Seurat::RunPCA(scrna, verbose = FALSE,npcs = length(opts$integ_nPCS),assay="SCT",features = Ifeatures)
scrna <- Seurat::RunUMAP(scrna, reduction = "pca", dims = opts$integ_nPCS,return.model =TRUE,assay="SCT",reduction = "pca")

message("-----",Sys.time()," Save data")

saveRDS(scrna,paste0(raw_data_path,"human_hypo_seuratIntegration.rds"))

saveRDS(scrna@reductions[["umap"]],paste0(raw_data_path,"human_hypo_seuratIntegration_umap.rds"))

message("-----",Sys.time()," Finished")
