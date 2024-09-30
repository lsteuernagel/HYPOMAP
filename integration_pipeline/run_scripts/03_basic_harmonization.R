##########
### Load parameters and packages
##########
message(Sys.time(),": Load parameters and packages ")

library(magrittr)
library(scUtils)
library(Seurat)

source("integration_pipeline/harmonization_functions.R")

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
harmonized_seurat_object = readRDS(paste0(parameter_list$merged_file))


##########
### clean object
##########

harmonized_seurat_object@project.name = parameter_list$new_name_suffix
harmonized_seurat_object@misc= list()
harmonized_seurat_object@graphs= list()
harmonized_seurat_object@reductions= list()
dummy=matrix(data = as.numeric())
harmonized_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

##########
### Raw data & new object
##########

# put scvi imputed into another assay
# message(Sys.time(),": Adding counts as new assay from: ",paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_scVI_corrected.txt"))
#
# corrected_counts_scvi = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_scVI_corrected.txt"),data.table = F)
# rownames(corrected_counts_scvi) = corrected_counts_scvi[,1]
# corrected_counts_scvi = t(as.matrix(corrected_counts_scvi[,2:ncol(corrected_counts_scvi)]))
# scvi_assay = CreateAssayObject(counts=corrected_counts_scvi, min.cells = 0, min.features = 0)
#
# harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]] = scvi_assay
# harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]]@key =paste0(parameter_list$integration_name,"_corrected")

# harmonized_seurat_object@assays$scvi_corrected@counts = harmonized_seurat_object@assays$scvi_corrected@data
# harmonized_seurat_object =NormalizeData(harmonized_seurat_object,assay = "scvi_corrected")

# # add hvgs
# message("Adding variable feature names to assays")
# harmonized_seurat_object@misc$var_features = as.character(unlist(parameter_list$var_features))
# # harmonized_seurat_object@assays[['RNA']]@var.features=parameter_list$var_features
# # harmonized_seurat_object@assays[[paste0(parameter_list$integration_name,'_corrected')]]@var.features=parameter_list$var_features


##########
### Reductions
##########

message(Sys.time(),": Adding integrated embedding and calculating default umap")
# add embedding manually
current_embedding = read_embedding(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_scVI_reduction.txt"),harmonized_seurat_object)

# make dim red
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(current_embedding),
  stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
  assay = "RNA",
  key = parameter_list$integration_name
)
# add  to object
harmonized_seurat_object@reductions[[parameter_list$integration_name]] = dimred


##########
### UMAP
##########

# run umap and save model
message(Sys.time(),": Build UMAP with ",parameter_list$k_param," n.neighbors ..." )
harmonized_seurat_object = RunUMAP(harmonized_seurat_object,
                                   reduction = parameter_list$integration_name,
                                   seed.use= parameter_list$global_seed,
                                   dims=1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                   reduction.name=paste0("umap_",parameter_list$integration_name),
                                   reduction.key = paste0("umap_",parameter_list$integration_name),
                                   verbose=F,
                                   n.neighbors = parameter_list$k_param,
                                   return.model = TRUE)

##########
### basic SNN:
##########

# run seurat SNN annoy
message(Sys.time(),": Build SNN with ",parameter_list$k_param," n.neighbors ..." )
harmonized_seurat_object = FindNeighbors(harmonized_seurat_object,
                                         reduction=parameter_list$integration_name,
                                         dims = 1:ncol(harmonized_seurat_object@reductions[[parameter_list$integration_name]]@cell.embeddings),
                                         k.param = parameter_list$k_param,
                                         nn.method="annoy",
                                         annoy.metric=parameter_list$dist_type,
                                         graph.name = c(paste0("NN_",parameter_list$integration_name),paste0("SNN_",parameter_list$integration_name)), verbose=TRUE)
## TODO:
# need to drop NN if SNN is to be exported via SaveH5Seurat & Convert ?  ---> should be ok!
 # If a nearest-neighbor graph is associated with assay, it will be added to uns/neighbors/distances; if more than one graph is present, then only the last graph according to the index will be added. 


##########
### Save object
##########

message(Sys.time(),": Save objects ..." )

file_name_prefix = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix)

# save data to rds
saveRDS(harmonized_seurat_object,paste0(file_name_prefix,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = harmonized_seurat_object,filename = paste0(file_name_prefix,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(file_name_prefix,".h5seurat"), dest =  paste0(file_name_prefix,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(file_name_prefix,".h5seurat")))

message("Finalized basic seurat object harmonization.")
