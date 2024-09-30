##########
### Load parameters and packages
##########

message(Sys.time(),": Starting final object build .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

command_args<-commandArgs(TRUE)
param_file = command_args[1]
#param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### Load combined
##########

message(Sys.time(),": Load data .." )

plot_path = paste0(parameter_list$harmonization_folder_path,"annotation_plots/")
dir.create(plot_path)

# read seurat:
human_hypo_combined = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
# other data
human_hypo_combined_edgelist = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","edgelist_mrtree",".txt"),data.table = F)

# load marker genes
human_hypo_combined_markers_all = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","comparisons_all_updated",".txt"),data.table = F)
human_hypo_combined_markers_sibling = data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","comparisons_siblings_updated",".txt"),data.table = F)
human_hypo_combined_markers_CO =  data.table::fread(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_markers_C0",".txt"),data.table = F)
# "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3//human_hypo_combined/human_hypo_combined_markers_C0.txt"

# load annotation results
annotated_labelmat_final = data.table::fread(file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","cluster_annotation",".txt"),data.table = F)
annotation_results_list = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_","cluster_annotation_allData",".rds"))


#########
### Add C0 markers & filter a bit
##########

human_hypo_combined_markers_all_final = bind_rows(human_hypo_combined_markers_CO  %>%
                                                    dplyr::mutate(comparison = "All",parent="all") %>%
                                                    dplyr::select(cluster = cluster,gene ,comparison,parent,specificity,avg_log2FC,pct.1,pct.2,effect_size,p_val_adj,p_val),
                                                  human_hypo_combined_markers_all) %>%
  dplyr::filter( p_val_adj < 0.01 & specificity >= 0.33333) # some minimal filtering

human_hypo_combined_markers_sibling_final = human_hypo_combined_markers_sibling %>%
  dplyr::filter( p_val_adj < 0.01 & specificity >= 0.33333) # some minimal filtering

# for genetics slightly more strict cell type signatures for c4:
human_hypo_combined_markers_C4_all_strict = human_hypo_combined_markers_all_final %>%
  dplyr::filter( p_val_adj < 0.0001 & specificity >= 1 & grepl("C4",cluster))

##########
### Add annotations to seurat metadata,marker tables and edgelist
##########

# loop over all rows and add to metadata
tempmeta = human_hypo_combined@meta.data[,1:24]
for(col in 0:4){
  cname = paste0("C",col)
  anno_to_join = annotated_labelmat_final[,c(cname,paste0(cname,"_named"))] %>%
    dplyr::distinct()
  tempmeta = dplyr::left_join(tempmeta,anno_to_join)
}
tempmeta = tempmeta %>% dplyr::select(-leiden_clusters_12_simplified,-celltype_status) %>% as.data.frame()
rownames(tempmeta) = tempmeta$Cell_ID
human_hypo_combined@meta.data = tempmeta

# make long version
annotated_labelmat_final_long_l = list()
for(col in 0:4){
  cname = paste0("C",col)
  anno_to_join = annotated_labelmat_final[,c(cname,paste0(cname,"_named"))] %>%
    dplyr::distinct()
  colnames(anno_to_join) =c("id","name")
  annotated_labelmat_final_long_l[[cname]] = anno_to_join
}
annotated_labelmat_final_long = do.call(rbind,annotated_labelmat_final_long_l) %>% as.data.frame() %>% dplyr::distinct()

# add to edgelist
human_hypo_combined_edgelist_anno = human_hypo_combined_edgelist %>%
  dplyr::left_join(annotated_labelmat_final_long,by=c("from"="id")) %>%
  dplyr::rename(from_named = name) %>%
  dplyr::left_join(annotated_labelmat_final_long,by=c("to"="id")) %>%
  dplyr::rename(to_named = name) 

# add to markers
human_hypo_combined_markers_all_final = human_hypo_combined_markers_all_final %>%
  dplyr::left_join(annotated_labelmat_final_long,by=c("cluster"="id"))  %>%
  dplyr::select(cluster = cluster, name = name,gene ,comparison,parent,specificity,avg_log2FC,pct.1,pct.2,effect_size,p_val_adj,p_val)
human_hypo_combined_markers_sibling_final = human_hypo_combined_markers_sibling_final %>%
  dplyr::left_join(annotated_labelmat_final_long,by=c("cluster"="id"))  %>%
  dplyr::select(cluster = cluster, name = name,gene ,comparison,parent,specificity,avg_log2FC,pct.1,pct.2,effect_size,p_val_adj,p_val)

##########
### Prepare tables with markers used for annotation
##########

colnames(annotation_results_list$c3_level_allMarkers)

c3_markers_anno = annotation_results_list$c3_level_allMarkers %>%
  dplyr::mutate(valid = case_when(
    .default = TRUE,
    grepl("LINC|AL[0-9]{3,}|AP[0-9]{3,}|AC[0-9]{3,}|MIR|-AS|-PS|orf[0-9]{2,}|[0-9]{2,}\\.[0-9]",gene) ~ FALSE,
    gene %in% c("ZNF385B","PLCG2","ELAVL2","SLC44A5","GALNTL6, MGAT4C") ~ FALSE,
    pct.2 > 0.3 ~ FALSE,
    abs(score_siblings_children) > 5 ~ FALSE
  )) %>%
  dplyr::filter(combined_score >= 1) %>%
  dplyr::select(cluster,gene,score_all_self,score_siblings_self,score_siblings_children,combined_score,valid) %>%
  dplyr::arrange(cluster,desc(combined_score)) 

c4_markers_anno = annotation_results_list$c4_level_allMarkers %>%
  dplyr::mutate(valid = case_when(
    .default = TRUE,
    grepl("LINC|AL[0-9]{3,}|AP[0-9]{3,}|AC[0-9]{3,}|MIR|-AS|-PS|orf[0-9]{2,}|[0-9]{2,}\\.[0-9]",gene) ~ FALSE,
    gene %in% c("ZNF385B","PLCG2","ELAVL2","SLC44A5","GALNTL6, MGAT4C") ~ FALSE,
    pct.2 > 0.3 ~ FALSE,
    abs(score_siblings_children) > 5 ~ FALSE
  )) %>%
  dplyr::filter(combined_score >= 1) %>%
  dplyr::select(cluster,gene,score_all_self,score_siblings_self,score_siblings_children,combined_score,valid) %>%
  dplyr::arrange(cluster,desc(combined_score)) 

##########
### Add addditional misc data
##########

## add to object
human_hypo_combined@misc$edgelist = human_hypo_combined_edgelist_anno

##########
### Export objects
##########

message(Sys.time(),": Export merged objects " )
# read list again
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

dir.create(parameter_list$harmonization_folder_path,showWarnings = F)

dummy=matrix(data = as.numeric())
human_hypo_combined@assays[["RNA"]]@var.features = character()
human_hypo_combined@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

# save file
data.table::fwrite(human_hypo_combined_edgelist_anno,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_edgelist_mrtree_annotated.txt"),sep="\t")

# save parameter_list$new_name_suffix = "human_hypo_combined"
data.table::fwrite(human_hypo_combined_markers_all_final,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_comparisons_all_annotated.txt"),sep="\t")
data.table::fwrite(human_hypo_combined_markers_sibling_final,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_comparisons_siblings_annotated.txt"),sep="\t")
# for genetics & anno
data.table::fwrite(human_hypo_combined_markers_C4_all_strict,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_comparisons_C4_all_filtered.txt"),sep="\t")
data.table::fwrite(c3_markers_anno,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotation_genes_C3.txt"),sep="\t")
data.table::fwrite(c4_markers_anno,file=paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_annotation_genes_C4.txt"),sep="\t")

# make file name
subset_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix)
message("Saving to seurat: ",subset_file_name)

# rds
saveRDS(human_hypo_combined,paste0(subset_file_name,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = human_hypo_combined,filename = paste0(subset_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(subset_file_name,".h5seurat"), dest =  paste0(subset_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)

message(Sys.time(),": Complete." )


