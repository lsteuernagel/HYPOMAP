
#sbatch -J umap_crossspecies_scvi -o /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/other_logs/umap_crossspecies_scvi_slurm-%j.out -e /beegfs/scratch/bruening_scratch/lsteuernagel/slurm/other_logs/umap_crossspecies_scvi_slurm-%j.err integration_pipeline/run_scripts/run_Rscript_slurm.sh ~/Documents/r_scvi_v3_42.simg merge_human_mouse_neurons/03_calculate_umap.R 


#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons/cross_species_neurons_mgvi_1_MultiGroupVI_mouse_reduction.txt"

# emebeddings_to_load = c("cross_species_neurons_scVI_reduction",
#                         "cross_species_neurons_2_scVI_reduction", 
#                         "cross_species_neurons_3_scVI_reduction", 
#                         "cross_species_neurons_4_scVI_reduction", 
#                         "cross_species_neurons_mgvi_1_MultiGroupVI_shared_reduction" )

emebeddings_to_load = c("cross_species_neurons_scvi_1_scVI_reduction")

# load libs
library(Seurat)
source("utility_functions.R")

# load parameters
message("Load  params and data ...")
param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/cross_species_hypothalamus_neurons_2/parameters_cross_species_neurons_v2.json"
param_list = jsonlite::read_json(param_file)

# make path for umaps
umap_pathname = "evaluated_umaps/"
umap_pathname = paste0(param_list$harmonization_folder_path,umap_pathname)
system(paste0("mkdir -p ",paste0(umap_pathname)))

# load seurat object
merged_seurat_object = readRDS(param_list$merged_file)

message("umaps: ")
# iterate over embeddings
for ( embed in emebeddings_to_load){
  message(embed)
  # find file
  embed_file = c()
  embed_file = list.files(paste0(param_list$harmonization_folder_path),pattern = embed,recursive = TRUE,full.names = TRUE)
  if(length(embed_file) != 1){
    message("Found ",length(embed_file)," files for embed. Please provide exactly one valid file." )
    next
  }
  
  # load
  message("Load ",embed," ...")
  embedding = read_embedding(filename_withpath = embed_file,seurat_object = merged_seurat_object)
  
  #add
  dimred <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(embedding),
    stdev = as.numeric(apply(embedding, 2, stats::sd)),
    assay = "RNA",
    key = embed
  )
  # add
  merged_seurat_object@reductions[[embed]] = dimred
  
  # run UMAP
  message("Run UMAP ...")
  new_name = paste0("umap_",embed)#"scvi"
  merged_seurat_object = Seurat::RunUMAP(merged_seurat_object,
                                         reduction = embed,
                                         seed.use=  param_list$global_seed,
                                         dims=1:ncol(merged_seurat_object@reductions[[embed]]@cell.embeddings),
                                         reduction.name= new_name,
                                         reduction.key = new_name,
                                         verbose=F,
                                         n.neighbors = param_list$k_param,
                                         return.model = TRUE)
  
  message("Save ... ")
  saveRDS(merged_seurat_object@reductions[[new_name]],file = paste0(umap_pathname,new_name,".rds"))
  
}
