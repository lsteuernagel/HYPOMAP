
##########
### Load
##########

message("-----",Sys.time(),": Load parameters and packages ")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(MetaNeighbor)
source("utility_functions.R")
source("integration_pipeline/harmonization_functions.R")
source("merge_human_mouse_neurons/cluster_matching_functions.R")
source("paper_figures/tree_plotting_functions.R")

param_file = "merge_human_mouse_neurons/parameters_cross_species_neurons_scvi_v1.json"
param_list = jsonlite::read_json(param_file)

# make path for umaps
umap_pathname = "evaluated_umaps/"
umap_pathname = paste0(param_list$harmonization_folder_path,umap_pathname)

message("-----",Sys.time(),": Load data ")

# load seurat object
human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"
hypothalamus_neurons_cross_species =  readRDS(paste0(human_hypo_path,"hypothalamus_neurons_cross_species.rds")) 

embed = c( "cross_species_neurons_scvi_1_scVI_reduction" )
embed_short = substr(embed, start = 1, stop = 45)
hypothalamus_neurons_cross_species[[paste0("umap_",embed_short)]] = readRDS(paste0(umap_pathname,"umap_",embed,".rds"))

## mouse
hypoMapUtils::load_required_files()

hypoMap_annos = hypoMap@misc$annotation_result
hypoMap_edgelist =  hypoMap@misc$clustering_edgelist %>% 
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,from_named = clean_names_withID),by=c("from"="cluster_id")) %>%
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,to_named = clean_names_withID),by=c("to"="cluster_id"))


##########
### Summarized experiment
##########

library(SummarizedExperiment)
library(Matrix)
library(SingleCellExperiment)
# SummarizedExperiment object containing gene-by-sample expression matrix.
mn_data = SingleCellExperiment::SingleCellExperiment(assays = as.matrix(t(embedding),"sparseMatrix"))
scvi_names= colnames(embedding)

##########
### Add clustering results from human neurons using embedding -- we skip this for now
##########

# hypothalamus_neurons_cross_species@meta.data$combined_clusters = hypothalamus_neurons_cross_species@meta.data$C465_named
# hypothalamus_neurons_cross_species@meta.data$combined_clusters[is.na(hypothalamus_neurons_cross_species@meta.data$combined_clusters)] = hypothalamus_neurons_cross_species@meta.data$C4_named[is.na(hypothalamus_neurons_cross_species@meta.data$combined_clusters)]
# 
# library(MetaNeighbor)
# 
# celltype_NV = MetaNeighborUS(var_genes = scvi_names,
#                              dat = mn_data,
#                              study_id = hypothalamus_neurons_cross_species@meta.data$species,
#                              cell_type = hypothalamus_neurons_cross_species@meta.data$combined_clusters,
#                              fast_version = TRUE)


##########
### Add clustering results from human neurons using genes
##########

mn_data_genes = SingleCellExperiment::SingleCellExperiment(assays = hypothalamus_neurons_cross_species@assays$RNA@data)
feature_set_cross_species = unlist(jsonlite::read_json("merge_human_mouse_neurons/feature_set_cross_species.json"))

hypothalamus_neurons_cross_species@meta.data$combined_clusters = hypothalamus_neurons_cross_species@meta.data$C465_named
hypothalamus_neurons_cross_species@meta.data$combined_clusters[is.na(hypothalamus_neurons_cross_species@meta.data$combined_clusters)] = hypothalamus_neurons_cross_species@meta.data$C4_named[is.na(hypothalamus_neurons_cross_species@meta.data$combined_clusters)]

celltype_NV_genes = MetaNeighborUS(var_genes = feature_set_cross_species,
                                   dat = mn_data_genes,
                                   study_id = hypothalamus_neurons_cross_species@meta.data$species,
                                   cell_type = hypothalamus_neurons_cross_species@meta.data$combined_clusters,
                                   fast_version = TRUE)

##########
### Filter Metaneighbor to final matching clusters
##########

main_filtering_threshold = 0.9

# topHitsByStudy
top_hits_study_genes = topHitsByStudy(
  auroc = celltype_NV_genes,
  threshold = main_filtering_threshold,
  n_digits = 2,
  collapse_duplicates = TRUE
)

# for genes similar format with edgelist and removing prefixes:
top_hits_study_genes_renamed = top_hits_study_genes %>%
  dplyr::mutate(from = stringr::str_remove(`Study_ID|Celltype_1`,"human\\|"),to = stringr::str_remove(`Study_ID|Celltype_2`,"mouse\\|")) %>%
  dplyr::mutate(from_save = from, to_save = to ) %>%
  dplyr::mutate(from = case_when(
    grepl("mouse",from) ~ to_save,
    .default = from_save)) %>%
  dplyr::mutate(to = case_when(
    grepl("human",to) ~ from_save,
    .default = to_save)) %>%
  dplyr::mutate(from = stringr::str_remove(from,"human\\|"),to = stringr::str_remove(to,"mouse\\|")) %>%
  dplyr::select(from, to, similarity = AUROC, Match_type )

top_hits_study_genes_reciprocal = top_hits_study_genes_renamed %>% dplyr::filter(Match_type == "Reciprocal_top_hit")

# prune_similarity_network genes
similarity_edgelist_genes = celltype_NV_genes %>% as.data.frame() %>%
  dplyr::mutate(from = rownames(celltype_NV_genes)) %>%
  tidyr::gather(key="to",value =similarity,-from) %>%
  dplyr::filter(grepl("human",from) & grepl("mouse",to)) %>%
  dplyr::mutate(from = stringr::str_remove(from,"human\\|"),to = stringr::str_remove(to,"mouse\\|")) %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::filter(grepl("C3",from)) %>% dplyr::select(from_parent = from_named,to_named),by=c("from"="to_named")) %>%
  dplyr::left_join(hypoMap_edgelist %>% dplyr::filter(grepl("C286",from_named)) %>% dplyr::select(to_parent = from_named,to_named),by=c("to"="to_named"))

pruned_neighbors_genes= prune_similarity_network(similarity_edgelist_genes, adjust_similarity = FALSE,similarity_threshold = main_filtering_threshold,  min_similarity_siblings = main_filtering_threshold,min_similarity = main_filtering_threshold)
pruned_neighbors_genes_adj= prune_similarity_network(similarity_edgelist_genes, adjust_similarity = TRUE,similarity_threshold = main_filtering_threshold,  min_similarity_siblings = main_filtering_threshold,min_similarity = main_filtering_threshold)

##########
### Correlate clusters between mouse & human C4
##########

colnames(hypothalamus_neurons_cross_species@meta.data)

# get list with all clusters
cluster_col = "C4_named"
cells_cluster_list_human = split(hypothalamus_neurons_cross_species@meta.data$Cell_ID,f = hypothalamus_neurons_cross_species@meta.data[,cluster_col])
cells_cluster_list_human = cells_cluster_list_human[sapply(cells_cluster_list_human,length) >= 30 ]

# calc avg scvi per cluster
mean_scvi_cluster_list_human = lapply(cells_cluster_list_human, function(cells,expr_matrix){
  if(length(cells)>1){
    return(colMeans(expr_matrix[cells,]))
  }else{
    return(expr_matrix[cells,])
  }
},expr_matrix = hypothalamus_neurons_cross_species@reductions[[embed]]@cell.embeddings)

mean_human_cluster = do.call(cbind,mean_scvi_cluster_list_human)

## repeat with hypoMap
# get list with all clusters
cluster_col = "C465_named"
cells_cluster_list_mouse = split(hypothalamus_neurons_cross_species@meta.data$Cell_ID,f = hypothalamus_neurons_cross_species@meta.data[,cluster_col])
cells_cluster_list_mouse = cells_cluster_list_mouse[sapply(cells_cluster_list_mouse,length) >= 30 ]

# calc avg scvi per cluster
mean_scvi_cluster_list_mouse = lapply(cells_cluster_list_mouse, function(cells,expr_matrix){
  if(length(cells)>1){
    return(colMeans(expr_matrix[cells,]))
  }else{
    return(expr_matrix[cells,])
  }
},expr_matrix = hypothalamus_neurons_cross_species@reductions[[embed]]@cell.embeddings)

mean_scvi_cluster_mouse = do.call(cbind,mean_scvi_cluster_list_mouse)

### correlate
cor_clusters = cor(x=mean_human_cluster,y=mean_scvi_cluster_mouse)
# cor_clusters = as.matrix(pdist::pdist(X=t(mean_scvi_cluster),Y=t(mean_scvi_cluster_hypomap)))
# rownames(cor_clusters) = colnames(mean_scvi_cluster)
# colnames(cor_clusters) = colnames(mean_scvi_cluster_hypomap)
cor_clusters_df_c4 = cor_clusters %>% as.data.frame() %>% dplyr::mutate(human_clusters = rownames(cor_clusters))

heatmap_df_long_c4 = cor_clusters_df_c4 %>% tidyr::gather(-human_clusters,key="hypomap_clusters",value="cor")# %>% dplyr::rename(time_step=time_from_intervention_int)
heatmap_df_long_c4_corrected = heatmap_df_long_c4 %>% dplyr::group_by(hypomap_clusters)  %>% dplyr::mutate(corrected_cor = cor - mean(cor)) %>%
  dplyr::group_by(human_clusters)  %>%  dplyr::mutate(human_zscore = (corrected_cor - mean(corrected_cor) ) / sd(corrected_cor))  # %>


##########
### Prune correlation network to final matches
##########

# need: 
# heatmap_df_long_c4_corrected from 05_correlate...

# Input: edgelist with correlatios
cor_edgelist = heatmap_df_long_c4_corrected %>%
  dplyr::select(from = human_clusters, to = hypomap_clusters, similarity = cor)  %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::filter(grepl("C3",from)) %>% dplyr::select(from_parent = from_named,to_named),by=c("from"="to_named")) %>%
  dplyr::left_join(hypoMap_edgelist %>% dplyr::filter(grepl("C286",from_named)) %>% dplyr::select(to_parent = from_named,to_named),by=c("to"="to_named"))

pruned_cor_edgelist= prune_similarity_network(cor_edgelist, adjust_similarity = FALSE,similarity_threshold = 0.7,  min_similarity_siblings = 0.8,min_similarity = 0.7)
pruned_cor_edgelist_adj= prune_similarity_network(cor_edgelist, adjust_similarity = TRUE,similarity_threshold = 0.7,  min_similarity_siblings = 0.8,min_similarity = 0.7)


##########
### Upset Comparisons
##########

## Aim: compare scvi_metaneighbor, scvi_metaneighbor_reciprocal, scvi_cor, genes_metaneighbor , genes_metaneighbor_reciprocal

## make edge id and show as upset plot
all_edgelists = list(
  "MetaNeighborTopHits" = top_hits_study_genes_renamed,
  "MetaNeighborReciprocal" = top_hits_study_genes_reciprocal,
  #pruned_neighbors_genes = pruned_neighbors_genes,
  "MetaNeighborPruned" = pruned_neighbors_genes_adj,
  #pruned_cor_edgelist = pruned_cor_edgelist,
  "scviCorPruned" =  pruned_cor_edgelist_adj
)
all_edgelists = sapply(all_edgelists,function(df){
  df$edge_id = paste0(df$from,"_",df$to)
  df
})

all_edgelists_ids = sapply(all_edgelists,function(df){df$edge_id},simplify = F)

## plot upset !!
upset_comparison = UpSetR::upset(data = UpSetR::fromList(all_edgelists_ids),nsets = 10,nintersects = 50,order.by="freq",text.scale = 3)


# some more intersection comparisons
intersect_from_list = function(input){
  # this code is adapted from : UpSetR::fromList
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  rownames(data) <- elements # added this to include rownames!
  return(data)
}

all_intersects = intersect_from_list(all_edgelists_ids)

nrow(all_intersects[all_intersects$MetaNeighborReciprocal == 1 & all_intersects$scviCorPruned == 1,]) / nrow(all_intersects[all_intersects$MetaNeighborReciprocal == 1 , ])
nrow(all_intersects[all_intersects$MetaNeighborReciprocal == 1 , ])
nrow(all_intersects[all_intersects$scviCorPruned == 1 , ])
nrow(all_intersects[all_intersects$scviCorPruned == 1 , ])-nrow(all_intersects[all_intersects$MetaNeighborReciprocal == 1 & all_intersects$scviCorPruned == 1,]) 

#pruned_cor_edgelist_adj_SPECIFIC = pruned_cor_edgelist_adj[pruned_cor_edgelist_adj$edgename %in% rownames(all_intersects[rowSums(all_intersects) == 1 & all_intersects[,"pruned_cor_edgelist_adj"] == 1,]),]

##########
### Visual example comparisons
##########
# ,"PVH-CRH/TRH"="C3-113"
validation_scenarios = c("POMC" = "C3-129","AGRP/SST"="C3-123","PMCH"="C2-43","PVH-AVP"="C3-106","PPP1R17"="C3-126","FOXD2"="C3-61",
                         "SCN-VIP"="C3-42","HDC"="C2-45","LHX9"="C3-117","GHRH/GAL"="C2-22","Post-LMX1A"="C3-84")

get_validation_mapping = function(validation_scenarios,matching_df,full_edgelist){
  all_scenarios = list()
  for(i in 1:length(validation_scenarios)){
    scenario_name = names(validation_scenarios)[i]
    scenario_node = validation_scenarios[i]
    all_scenario_nodes = scUtils::find_children(scenario_node,edges = full_edgelist)
    all_scenario_nodes_named = full_edgelist$to_named[full_edgelist$to %in% all_scenario_nodes]
    mapped_to = matching_df[matching_df$from %in% all_scenario_nodes_named,] 
    additional_mapped_from = matching_df[!matching_df$from %in% all_scenario_nodes_named & matching_df$to %in% mapped_to$to,] 
    scenario = dplyr::bind_rows(mapped_to,additional_mapped_from) %>%
      dplyr::mutate(source = case_when(from %in% all_scenario_nodes_named ~ "scenario",
                                       .default = "additional")) %>%
      as.data.frame()
    if(nrow(scenario) > 0){
      scenario$scenario_name = scenario_name
    }
    all_scenarios[[scenario_name]] = scenario
  }
  all_scenarios 
}

names(all_edgelists)

versions_to_compare = names(all_edgelists)#c("top_hits_study_genes_renamed","pruned_neighbors_genes_adj","pruned_neighbors_scvi_adj","pruned_cor_edgelist_adj")
all_validation_mappings  = list()
for(v in versions_to_compare){
  all_validation_mappings[[v]] = get_validation_mapping(validation_scenarios,matching_df =  all_edgelists[[v]],full_edgelist=combined_edgelist_mrtree)
}

all_scenario_plots = list()
for(i in 1:length(validation_scenarios)){
  scenario_name = names(validation_scenarios)[i]
  message(">>>>",scenario_name)
  scenario_node = validation_scenarios[i]
  scenario_plots = list()
  for(j in 1:length(all_validation_mappings)){
    message("   >>",names(all_validation_mappings)[j])
    matched_clusters = all_validation_mappings[[j]][[scenario_name]]
    if(nrow(matched_clusters)>0){
      scenario_plots[[j]] = plot_tree_comparison(matched_clusters,edgelist_human = combined_edgelist_mrtree, edgelist_mouse = hypoMap_edgelist ,general_label_size = 5,min_similarity=0.7,color_order_regex = "C4",start_level = "C2")+ggtitle(names(all_validation_mappings)[j])
    }else{
      scenario_plots[[j]] = NULL 
    }
  }
  all_scenario_plots[[scenario_name]] = cowplot::plot_grid(plotlist = scenario_plots,nrow = 2)
}


# plot results
names(validation_scenarios)
all_scenario_plots$POMC
all_scenario_plots$`AGRP/SST`
all_scenario_plots$PMCH
all_scenario_plots$PPP1R17
all_scenario_plots$FOXD2
all_scenario_plots$HDC
all_scenario_plots$`SCN-VIP`

##########
### Export 
##########

output_dir = "merge_human_mouse_neurons/evaluation/"
dir.create(output_dir)

# export the updset plot
pdf(file =  paste0(output_dir,"upset_comparison.pdf"), width = 20,height = 12) # The height of the plot in inches
upset_comparison
dev.off()

# export some select cell type plots
for(i in 1:length(all_scenario_plots)){
  sc_name = names(all_scenario_plots)[i]
  sc_name = gsub("\\/|-","",sc_name)
  message(sc_name)
  ggsave(filename = paste0(output_dir,"celltype_comp_plot_",sc_name,".pdf"),
         plot = all_scenario_plots[[i]], "pdf",dpi=400,width=500,height = 350,units="mm")
  
}
# export matching results for all comapred options
for(i in 1:length(all_edgelists)){
  match_name = names(all_edgelists)[i]
  data.table::fwrite(all_edgelists[[i]],paste0(output_dir,"matched_clusters_",match_name,".tsv"),sep="\t")
}

