
# sbatch run_slurm.sh ~/Documents/r_scvi_v3.simg 02_filter_data.R

##########
### load
##########
message("-----",Sys.time(),": load libs and data")

library(magrittr)
library(scUtils)
library(dplyr)
library(Matrix)
library(Seurat)
source("utility_functions.R")

opts <- workflow_options(project = "HuHy", out_path = "test/")
raw_data_path = opts$data_path

seurat_object = readRDS(paste0(raw_data_path,"human_hypo_raw.rds"))

doublet_cluster_tresh = 0.75
max_mt = 10
minUMI = 800

##########
### additionally set outlier column
##########

message("-----",Sys.time(),": add doublet and lowq cells")

cluster_column = "seurat_clusters"
doublet_stats_per_cluster = seurat_object@meta.data %>% dplyr::select(cluster = !!rlang::sym(cluster_column),scDblFinder.class) %>%
  dplyr::group_by(cluster) %>% dplyr::add_count(name="cells_per_cluster") %>% dplyr::filter(scDblFinder.class=="doublet") %>% dplyr::add_count(name="doublet_per_cluster") %>%
  dplyr::distinct(cluster,cells_per_cluster,doublet_per_cluster) %>% dplyr::mutate(doublet_pct = doublet_per_cluster / cells_per_cluster)
# also filter out full cluster with certain pct!
cells_in_doublet_clusters = rownames(seurat_object@meta.data)[seurat_object@meta.data[,cluster_column] %in% doublet_stats_per_cluster$cluster[doublet_stats_per_cluster$doublet_pct >= doublet_cluster_tresh]]
seurat_object@meta.data$outlier[rownames(seurat_object@meta.data)  %in% cells_in_doublet_clusters] =TRUE

#DimPlot(seurat_object,cells.highlight = cells_in_doublet_clusters)

## min counts and max mt
other_lowq_cells = rownames(seurat_object@meta.data)[(seurat_object@meta.data$percent.mt >= max_mt | seurat_object@meta.data$nCount_RNA < minUMI) & !seurat_object@meta.data$outlier ]
seurat_object@meta.data$outlier[rownames(seurat_object@meta.data) %in% other_lowq_cells] =TRUE

# DimPlot(seurat_object,group.by = "outlier",raster.dpi = c(2048,2048),pt.size = 1.8,order = F,shuffle = TRUE)
# table(seurat_object@meta.data$outlier)

##########
### subset
##########

# subset
seurat_object.outliers <- subset(seurat_object, outlier == TRUE)
seurat_object.filtered <- subset(seurat_object, outlier == FALSE)

##########
### save data - just a seurat rds
##########

message("-----",Sys.time(),":  save data")

saveRDS(seurat_object.filtered,paste0(raw_data_path,"human_hypo_raw_filtered.rds"))
saveRDS(seurat_object.outliers,paste0(raw_data_path,"human_hypo_raw_outliers.rds"))
