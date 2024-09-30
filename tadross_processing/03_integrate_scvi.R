##########
### Load object
##########

# make file name

# sbatch run_slurm.sh ~/Documents/r_scvi_v3.simg 03_integrate_scvi.R

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

opts <- workflow_options(project = "HuHy", out_path = "test/")

# get params-filename from commandline
# command_args<-commandArgs(TRUE)
# param_file = command_args[1]
param_file = "parameters_scvi.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})
#parameter_list$parameter_list$merged_file_name
# parameter_list$merged_file_name = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/human_nucseq_processed"

# read features to excludes
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = unlist(lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}}))
features_exclude_list = toupper(features_exclude_list)

print(head(features_exclude_list))

# load seurat
message("Load from : ",paste0(parameter_list$input_seurat))
seurat_object = readRDS(paste0(parameter_list$input_seurat))

print(seurat_object)

##########
### Run feature detection
##########

message("-----",Sys.time(),": Add variable features ")

# normalize data
seurat_object <- Seurat::NormalizeData(object = seurat_object,  verbose = F, assay = "RNA")

# find HVGs
feature_set = scUtils::identify_variable_features(seurat_object,
                                                  n_hvgs_sizes = parameter_list$feature_set_size,
                                                  batch_var = parameter_list$batch_var,
                                                  assay_name = "RNA",
                                                  method = "vst",
                                                  ignore_genes_vector = features_exclude_list,
                                                  returnSeurat = FALSE,
                                                  seed = parameter_list$global_seed)
#feature_set = seurat_object@misc$var_features[[1]]
seurat_object@assays$RNA@var.features = as.character(feature_set)

# save:
scUtils::writeList_to_JSON(feature_set,filename = paste0(parameter_list$feature_set_file))

##########
### clean object
##########

seurat_object@misc= list()
seurat_object@graphs= list()
#seurat_object@reductions= list()
dummy=matrix(data = as.numeric())
seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

##########
### Export to anndata
##########

message("-----",Sys.time(),": Save object to anndata ..." )

seurat_object@meta.data$Cell_ID = rownames(seurat_object@meta.data)

dummy=matrix(data = as.numeric())
#seurat_object@assays[["RNA"]]@var.features = character()

# save h5seurat
SeuratDisk::SaveH5Seurat(object = seurat_object,filename = paste0(parameter_list$merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(parameter_list$merged_file_name,".h5seurat"), dest =  paste0(parameter_list$merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(parameter_list$merged_file_name,".h5seurat")))

message(" Completed export ")

##########
### Run scvi
##########

message("-----",Sys.time(),": Starting scvi execution script in python")

# use this script
# https://github.sf.mpg.de/lsteuernagel/scHarmonization/blob/main/python/integrate_scVI_v015.py

system(paste0("python3 -u ",parameter_list$scvi_script_file," ",param_file))


##########
### Read scvi result
##########

message("-----",Sys.time(),": Adding integrated scvi embedding ")
# add embedding manually
current_embedding = read_embedding(paste0(parameter_list$project_path,parameter_list$project_name,"_scVI_reduction.txt"),seurat_object)

# make dim red
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(current_embedding),
  stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
  assay = "RNA",
  key = parameter_list$integration_name
)
# add  to object
seurat_object@reductions[[parameter_list$integration_name]] = dimred


##########
### UMAP
##########

# run umap and save model
message("-----",Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors ..." )

seurat_object = Seurat::RunUMAP(seurat_object,
                                reduction = parameter_list$integration_name,
                                seed.use= parameter_list$global_seed,
                                dims=1:ncol(seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                reduction.name=paste0("umap_",parameter_list$integration_name),
                                reduction.key = paste0("umap_",parameter_list$integration_name),
                                verbose=F,
                                n.neighbors = parameter_list$k_param,
                                return.model = TRUE)

##########
### basic SNN:
##########

# run seurat SNN annoy
message("-----",Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors ..." )
seurat_object = Seurat::FindNeighbors(seurat_object,
                                      reduction="scvi",
                                      dims = 1:ncol(seurat_object@reductions[["scvi"]]@cell.embeddings),
                                      k.param = parameter_list$k_param,
                                      nn.method="annoy",
                                      annoy.metric=parameter_list$dist_type,
                                      graph.name = paste0("SNN_","scvi"), verbose=FALSE)

##########
### Preliminary clustering
##########

message("-----",Sys.time(),": Running preliminary louvain clustering")

# find best cluster resolution:

seurat_object = Seurat::FindClusters(seurat_object, resolution = 1:4, 
                     graph.name = paste0("SNN_","scvi"), random.seed = opts$seed, 
                     verbose = FALSE, min_cells = 5,return_seurat = TRUE)

##########
### Save results
##########

message("-----",Sys.time(),": Save results" )

# save RDS
saveRDS(seurat_object,file = paste0(parameter_list$project_path,parameter_list$project_name,"_processed.rds"))


message("-----",Sys.time(),": Complete." )



