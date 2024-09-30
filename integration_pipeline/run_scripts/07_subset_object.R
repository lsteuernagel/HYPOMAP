##########
### Load parameters and packages
##########

message(Sys.time(),":  Load parameters and package .." )

require(tidyverse)
require(Seurat)
require(Matrix)
library(parallel)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/stratified_wilcoxon_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1] # param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to excludes
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat -- after annotation o have the properly cleaned clusters!
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated",".rds"))

# path for plots
plot_path = paste0(parameter_list$harmonization_folder_path,"annotation_plots/")
dir.create(plot_path,showWarnings = FALSE)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)
getOkabeItoPalette(3)

##########
### Update annotation
##########

message(Sys.time(),": Add subset annotation to object " )

## merge all oligos
temp_meta = harmonized_seurat_object@meta.data
temp_meta$leiden_clusters_12_simplified = temp_meta$leiden_clusters_12
temp_meta$celltype_annotation[temp_meta$leiden_clusters_12_simplified %in% c("261","185","172","321","142")] = "Neurons"
temp_meta$celltype_annotation[temp_meta$leiden_clusters_12_simplified %in% c("108")] = "P2RX2_OTP"
# non neuron
temp_meta$leiden_clusters_12_simplified[temp_meta$celltype_annotation == "Ermn.Oligodendrocytes" ] = "Oligodendrocytes" #& temp_meta$leiden_clusters_12_simplified != "13"

#  which clusters to remove
temp_meta$celltype_status = "NonNeuron"
temp_meta$celltype_status[temp_meta$celltype_annotation %in%  c("Astrocytes","Ependymal")] = "AstroEpendymal"
temp_meta$celltype_status[temp_meta$celltype_annotation %in%  c("Oligodendrocytes","Ermn.Oligodendrocytes","Gpr17.Oligodendrocytes","OPC")] = "Oligodendrocytes"
temp_meta$celltype_status[temp_meta$celltype_annotation %in%  c("Neurons","P2RX2_OTP")] = "HypoNeuron"
temp_meta$celltype_status[temp_meta$celltype_annotation %in% c("PKP2_IGF1_Unknown","SHOX2_Thalamus","SLC17A7_TBR1_Unknown","GALNT5_PAX6_Unknown","LHX6_CHST9_Unknown", "LHX6_SCARA5_Unknown","DMBX1_Midbrain")] = "NonHypoNeuron"

# load additional problematic clusters:
preliminary_cells_to_remove = data.table::fread(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_2/qc_per_cluster/","preliminary_cells_to_remove.txt"),data.table = F)
temp_meta$celltype_status[temp_meta$Cell_ID %in% preliminary_cells_to_remove$Cell_ID] = "AdditionalRemove"

# add back in
rownames(temp_meta) = temp_meta$Cell_ID
harmonized_seurat_object@meta.data = temp_meta

# plot updated leiden 12
leiden_simplified_plot = DimPlot(harmonized_seurat_object,group.by = "leiden_clusters_12_simplified",label = TRUE,label.size = 2,raster.dpi = c(2048,2048),pt.size = 1.4)+NoLegend()
# export plot
ggsave(filename = paste0(plot_path,"leiden_12_simplified_umap.pdf"),plot = leiden_simplified_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

# plot celltype status
status_plot = DimPlot(harmonized_seurat_object,group.by = "celltype_status",label = TRUE,label.size = 4,cols = getOkabeItoPalette(length(unique(temp_meta$celltype_status))),raster.dpi = c(2048,2048),pt.size = 1.4)+NoLegend()
# export plot
ggsave(filename = paste0(plot_path,"celltype_status_umap.pdf"),plot = status_plot, "pdf",dpi=400,width=200,height = 200,units="mm")

data.table::fwrite(temp_meta[,c("Cell_ID","Sample_ID","Donor_ID","Dataset","celltype_annotation","leiden_clusters_12_simplified","celltype_status")],paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotated","status_meta_data.txt"),sep="\t")
#saveRDS(harmonized_seurat_object,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization/human_hypo_annotated_130123.rds")

##########
### Subset and UMAP
##########

subset_names = c("AstroEpendymal","Oligodendrocytes","HypoNeuron","NonNeuron")

message(Sys.time(),": Subset objects " )
subset_list = list()

for( subname in subset_names){ 
  message(subname)
  subset_list[[subname]] = subset(harmonized_seurat_object,subset =  celltype_status == subname) 
  message(Sys.time(),": Build UMAP on subset with ",parameter_list$k_param," n.neighbors ..." )
  subset_list[[subname]] = Seurat::RunUMAP(subset_list[[subname]],
                                           reduction = parameter_list$integration_name,
                                           seed.use= parameter_list$global_seed,
                                           dims=1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                           reduction.name=paste0("umap_",parameter_list$integration_name,"_",subname),
                                           reduction.key = paste0("umap_",parameter_list$integration_name,"_",subname),
                                           verbose=F,
                                           n.neighbors = parameter_list$k_param,
                                           return.model = TRUE)
  
}

##########
### Export subset objects
##########

message(Sys.time(),": Export subset objects " )


for( subname in subset_names){ 
  message(subname)
  
  dummy=matrix(data = as.numeric())
  subset_list[[subname]]@assays[["RNA"]]@var.features = character()
  subset_list[[subname]]@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay
  
  # make file name
  dir.create(paste0(parameter_list$harmonization_folder_path,"/",subname),showWarnings = FALSE)
  subset_file_name = paste0(parameter_list$harmonization_folder_path,subname,"/",parameter_list$new_name_suffix,"_",subname)
  
  # rds
  saveRDS(subset_list[[subname]],paste0(subset_file_name,".rds"))
  
  # save h5seurat
  SeuratDisk::SaveH5Seurat(object = subset_list[[subname]],filename = paste0(subset_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)
  
  # save to anndata
  SeuratDisk::Convert( paste0(subset_file_name,".h5seurat"), dest =  paste0(subset_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
  
}

message(Sys.time(),": Complete." )
