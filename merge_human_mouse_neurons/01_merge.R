
# sbatch run_slurm.sh ~/Documents/r_scvi_v3_42.simg merge_human_mouse_neurons/01_merge.R

##########
### Load parameters and packages
##########

message(Sys.time(),": human and mouse neuron merging .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)
library(parallel)

seed = 12345

merged_species_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/"
system(paste0("mkdir -p ",paste0(merged_species_path)))

features_exclude_list= jsonlite::read_json("data/human_features_exclude_list_2.json")
features_exclude_list = sapply(features_exclude_list,unlist)

##########
### Load human neuron data
##########

message("-----",Sys.time(),": Load human data")

human_neurons_raw = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_2/HypoNeuron/human_hypo_HypoNeuron.rds")
DefaultAssay(object = human_neurons_raw) <- "RNA"
# add clusters
human_neurons_clustering= readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_2/HypoNeuron/human_hypo_HypoNeuron_HypoNeuron_clusters_pruned_mrtree_clustering_results.rds")
human_neurons_clustering_labelmat = human_neurons_clustering$labelmat
human_neurons_raw@meta.data=cbind(human_neurons_raw@meta.data,human_neurons_clustering_labelmat)

##########
### Load mouse neuron data
##########

message("-----",Sys.time(),": Load mouse data")

hypomap_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2c_final/hypoMap_v2_neurons.rds")


##########
### match mouse neuron genes to human / Make List of 1:1 orthologues ?? 
##########

message("-----",Sys.time(),": Convert mouse data")

# could also consider to reolve 1:many
# using Ensembl Archive Release 101 (August 2020) : Update to human GENCODE 35 --° should correspond to siletti et al who used Gencode v35
hypomap_neurons_human = scUtils::convert_species_seurat(seurat_object = hypomap_neurons,host = "https://aug2020.archive.ensembl.org",verbose = TRUE)

##########
### make smaller human object
##########

message("-----",Sys.time(),": Make clean human object")

remove_samples = names(table(human_neurons_raw@meta.data$Sample_ID))[table(human_neurons_raw@meta.data$Sample_ID) < 100]
valid_cells = human_neurons_raw@meta.data$Cell_ID[! human_neurons_raw@meta.data$Sample_ID %in% remove_samples]
message("Valid cells: ",length(valid_cells)," out of ",length(human_neurons_raw@meta.data$Cell_ID)," total cells.")

valid_genes = intersect(rownames(hypomap_neurons_human),rownames(human_neurons_raw))
human_clean_counts = human_neurons_raw@assays$RNA@counts[valid_genes,valid_cells]
message(length(rownames(human_clean_counts)))

human_neurons = CreateSeuratObject(counts = human_clean_counts,meta.data = human_neurons_raw@meta.data[valid_cells,],project = "human_neurons")
DefaultAssay(object = human_neurons) <- "RNA"
human_neurons <- Seurat::NormalizeData(object = human_neurons,  verbose = F, assay = "RNA")

##########
### Only use shared features in HVG selection
##########

# exclusion
features_exclude_list_long = unique(c(features_exclude_list,setdiff(rownames(hypomap_neurons_human),rownames(human_neurons)),setdiff(rownames(human_neurons),rownames(hypomap_neurons_human))))

# exclude all non TFS
# Homo_sapiens_TF = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/uniprot/Homo_sapiens_TF_AnimalTFDB_4.txt",data.table = F)
# Homo_sapiens_TF_genes = unique(Homo_sapiens_TF$Symbol)
# non_tfs = unique(c(setdiff(rownames(hypomap_neurons_human),Homo_sapiens_TF_genes),setdiff(rownames(human_neurons),Homo_sapiens_TF_genes)))
# features_exclude_list_long = unique(c(features_exclude_list,non_tfs))
# message("Excluding ",length(non_tfs)," non_tfs. features")


message("Excluding ",length(features_exclude_list_long)," total features")

##########
### downsample data
##########

# Not doing this right now. Could do in future

##########
### subset data to TFS
##########

# viable_TFs_human = intersect(Homo_sapiens_TF_genes,rownames(human_neurons@assays$RNA@counts))
# seurat_human_TF = CreateSeuratObject(counts = human_neurons@assays$RNA@counts[viable_TFs_human,])
# seurat_human_TF@meta.data = human_neurons@meta.data
# seurat_human_TF = NormalizeData(seurat_human_TF)
# message("Using ",nrow(seurat_human_TF)," nrow seurat_human_TF")
# 
# viable_TFs_mouse = intersect(Homo_sapiens_TF_genes,rownames(hypomap_neurons_human@assays$RNA@counts))
# seurat_mouse_TF = CreateSeuratObject(counts = hypomap_neurons_human@assays$RNA@counts[viable_TFs_mouse,])
# seurat_mouse_TF@meta.data = hypomap_neurons_human@meta.data
# seurat_mouse_TF = NormalizeData(seurat_mouse_TF)
# message("Using ",nrow(seurat_mouse_TF)," nrow seurat_mouse_TF")

##########
### Run split by mouse sample feature detection
##########

## load mouse exclude feature list
#features_exclude_list_mouse = ....

message(Sys.time(),": Add variable features mouse ")

# normalize data
hypomap_neurons_human <- Seurat::NormalizeData(object = hypomap_neurons_human,  verbose = F, assay = "RNA")

# find HVGs
feature_set_mouse = scUtils::identify_variable_features(seurat_mouse_TF, # hypomap_neurons_human,
                                                        n_hvgs_sizes = 600,
                                                        batch_var = "Batch_ID",
                                                        assay_name = "RNA",
                                                        method = "vst",
                                                        ignore_genes_vector = features_exclude_list_long,
                                                        returnSeurat = FALSE,
                                                        seed = seed)
# save:
scUtils::writeList_to_JSON(feature_set_mouse,filename = paste0(merged_species_path,"feature_set_mouse.json"))


##########
### Run split by human sample feature detection
##########

## load mouse exclude feature list
#features_exclude_list_human = ....

message(Sys.time(),": Add variable features human ")

# normalize data
human_neurons <- Seurat::NormalizeData(object = human_neurons,  verbose = F, assay = "RNA")


# check if some sample_ids are included with too few cells!

# find HVGs
feature_set_human = scUtils::identify_variable_features(seurat_human_TF, #human_neurons,
                                                        n_hvgs_sizes = 600,
                                                        batch_var = "Sample_ID",
                                                        assay_name = "RNA",
                                                        method = "vst",
                                                        ignore_genes_vector = features_exclude_list_long,
                                                        returnSeurat = FALSE,
                                                        seed = seed)
# save:
scUtils::writeList_to_JSON(feature_set_human,filename = paste0(merged_species_path,"feature_set_human.json"))


##########
### Make list of robustly shared features
##########

message(Sys.time(),": Take intersection of HVGs ")

shared_hvgs = intersect(feature_set_mouse,feature_set_human)
message("Found ",length(shared_hvgs)," shared hvgs")
# save:
scUtils::writeList_to_JSON(shared_hvgs,filename = paste0(merged_species_path,"feature_set_cross_species.json"))

##########
### merge datasets
##########

message(Sys.time(),": Merging human and mouse datasets ")

# add species column
hypomap_neurons_human@meta.data$species = "mouse"
human_neurons@meta.data$species = "human"

merged_neuron_dataset = merge( x = human_neurons, y = hypomap_neurons_human)

# valid_genes = intersect(rownames(hypomap_neurons_human),rownames(human_neurons))
# merged_neuron_dataset_clean_counts = merged_neuron_dataset@assays$RNA@counts[valid_genes,]
# 
# merged_neuron_dataset = CreateSeuratObject(counts = merged_neuron_dataset_clean_counts,meta.data = merged_neuron_dataset@meta.data,project = "hypothalamus_neurons")
# human_neurons <- Seurat::NormalizeData(object = human_neurons,  verbose = F, assay = "RNA")

##########
### Export objects of merged_neuron_dataset
##########

message(Sys.time(),": Export objects " )

dummy=matrix(data = as.numeric())
merged_neuron_dataset@assays[["RNA"]]@var.features = character()
merged_neuron_dataset@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

merged_neuron_dataset@misc= list()
merged_neuron_dataset@graphs= list()
merged_neuron_dataset@reductions= list()

# make file name
merged_file_name = paste0(merged_species_path,"hypothalamus_neurons_cross_species")

# rds
saveRDS(merged_neuron_dataset,paste0(merged_file_name,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = merged_neuron_dataset,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(merged_file_name,".h5seurat")))

message(Sys.time(),": Completed." )

##########
### Run integration
##########
# 
# old!! :
# sbatch -J crossspecies_scvi -o /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/crossspecies_scvi_slurm-%j.out -e /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_slurmlogs/crossspecies_scvi_slurm-%j.err integration_pipeline/run_scripts/run_Python_slurm.sh ~/Documents/r_scvi_v3_42.simg scripts/02_integrate_scVI_crossspecies.py /beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/parameters_cross_species_neurons_v2.json


##########
### load and test integration
##########
# 
# hypothalamus_neurons_cross_species = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/hypothalamus_neurons_cross_species.rds")
# # 
# current_embedding = read_embedding("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/cross_species_neurons_scVI_reduction.txt",hypothalamus_neurons_cross_species)
# # make dim red
# dimred <- Seurat::CreateDimReducObject(
#   embeddings = as.matrix(current_embedding),
#   stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
#   assay = "RNA",
#   key = "scvi"
# )
# # add  to object
# hypothalamus_neurons_cross_species@reductions[["scvi"]] = dimred

##########
### UMAP
##########
 
# # run umap and save model
# hypothalamus_neurons_cross_species = RunUMAP(hypothalamus_neurons_cross_species,
#                                              reduction ="scvi",
#                                              seed.use= 123456,
#                                              dims=1:ncol(hypothalamus_neurons_cross_species@reductions[["scvi"]]@cell.embeddings),
#                                              reduction.name=paste0("umap_","scvi"),
#                                              reduction.key = paste0("umap_","scvi"),
#                                              verbose=F,
#                                              n.neighbors = 30,
#                                              return.model = TRUE)
# # palette.colors(palette = "Okabe-Ito")
# p1=DimPlot(hypothalamus_neurons_cross_species,group.by = "species",shuffle = TRUE,raster.dpi = c(1024,1024),cols = c("#E69F00","#009E73" ))+ggtitle(NULL) #+NoLegend()
# p2=FeaturePlot(hypothalamus_neurons_cross_species,features = "GHRH",raster.dpi = c(1024,1024),order = TRUE)
# p2
# cowplot::plot_grid(p1,p2)
# 
# FeaturePlot(hypothalamus_neurons_cross_species,features = "HCRT",split.by = "species",order = TRUE,raster.dpi = c(1024,1024))
# 
# hvgs = unlist(jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/feature_set_cross_species.json"))
# 
# FeaturePlot(harmonized_seurat_object,features = "SLC17A7",raster.dpi = c(1024,1024),order = TRUE)+NoAxes()
# 
# p1 = FeaturePlot(harmonized_seurat_object,features = "GHRH",raster.dpi = c(2048,2048),order = TRUE,pt.size = 1.3)+NoAxes()
# ggsave(filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/dump/human_ghrh.pdf"),
#        plot = p1, "pdf",dpi=300,width=320,height = 300,units="mm")
# 
# p2 = FeaturePlot(hypoMap,features = "Ghrh",raster.dpi = c(2048,2048),order = TRUE,pt.size = 1.3)+NoAxes()
# ggsave(filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/dump/mouse_ghrh.pdf"),
#        plot = p2, "pdf",dpi=300,width=320,height = 300,units="mm")
# 
# 
# 
# 
