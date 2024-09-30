##########
### Load parameters and packages
##########
message(Sys.time(),": Load parameters and packages ")

library(magrittr)
library(scUtils)

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
merged_seurat = readRDS(paste0(parameter_list$merged_file))

#get metadata and save
seurat_metadata = merged_seurat@meta.data
data.table::fwrite(x = seurat_metadata,file = paste0(parameter_list$harmonization_folder_path,paste0(parameter_list$new_name_suffix,"_metadata.txt")),sep = "\t")


##########
### Run feature detection
##########

message(Sys.time(),": Add variable features ")

# normalize data
merged_seurat <- Seurat::NormalizeData(object = merged_seurat,  verbose = F, assay = "RNA")

# find HVGs
feature_set = scUtils::identify_variable_features(merged_seurat,
                                                    n_hvgs_sizes = parameter_list$feature_set_size,
                                                    batch_var = parameter_list$batch_var,
                                                    assay_name = parameter_list$assay_name,
                                                    method = "vst",
                                                    ignore_genes_vector = features_exclude_list,
                                                    returnSeurat = FALSE,
                                                    seed = parameter_list$global_seed)
#feature_set = merged_seurat@misc$var_features[[1]]

# save:
scUtils::writeList_to_JSON(feature_set,filename = paste0(parameter_list$feature_set_file))


##########
### Export to anndata
##########

message(Sys.time(),": Save objects ..." )

dummy=matrix(data = as.numeric())
merged_seurat@assays[["RNA"]]@var.features = character()
merged_seurat@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

merged_seurat@misc= list()
merged_seurat@graphs= list()
merged_seurat@reductions= list()

# make file name
merged_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix)

# save h5seurat
SeuratDisk::SaveH5Seurat(object = merged_seurat,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(merged_file_name,".h5seurat")))

message(" Complete ")
