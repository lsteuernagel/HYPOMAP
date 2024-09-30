##########
### Load data
##########

figure_dir = "paper_figures/revision_figures/"
dir.create(figure_dir,showWarnings = FALSE)

# source
source("paper_figures/tree_plotting_functions.R")
source("paper_figures/celltype_heatmaps_functions.R")

# load files
human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"
library(Seurat)
library(tidyverse)
library(ggtree)

human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP.rds"))

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

# load cell type matching
# load matching result
matched_clusters_final = data.table::fread("merge_human_mouse_neurons/matched_clusters_scviCorPruned.tsv",data.table = F)

## mouse
hypoMap = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap.rds") # requires mouse HypoMap from: https://doi.org/10.17863/CAM.87955

##########
### Campbell comparison
##########

# I do only C286

# get C465 to Campbell & subset to cell types matched to Campbell
cluster_matching = hypoMap@meta.data %>%
  dplyr::filter(Dataset == "CampbellDropseq") %>%
  dplyr::group_by(Author_CellType,C286_named) %>% #C465_named
  dplyr::count() %>%
  dplyr::filter(n >= 10)

# map to human clusters
cluster_matching = dplyr::left_join(cluster_matching,matched_clusters_final[,c("from","to_parent","similarity")],by=c("C286_named"="to_parent")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(C286_named) %>%
  dplyr::slice_max(similarity,n = 1,with_ties = F) %>%
  dplyr::rename(human_cluster = from,mouse_cluster = C286_named, Campbell_cluster = Author_CellType) %>%
  dplyr::filter(grepl("GABA|GLU",mouse_cluster))

# add human cell numbers
human_cellcounts = human_hypo_combined@meta.data %>% dplyr::group_by(C4_named) %>% dplyr::count(name="nhuman")
cluster_matching = dplyr::left_join(cluster_matching,human_cellcounts,by=c("human_cluster"="C4_named"))

cluster_matching$Campbell_cluster[is.na(cluster_matching$Campbell_cluster)] = "NA"
cluster_matching$human_cluster[is.na(cluster_matching$human_cluster)] = "NA"

# maybe reduced to C286

# make 3 level sankey



# compare with ARC clustering



# Maybe pull out ARH clusters and re-run ?

##########
### Campbell comparison
##########

library(networkD3)

## using fun code
nlink=40
sankey_edges = cluster_matching %>% dplyr::filter(n >= nlink)
sankey_edges$Campbell_cluster[sankey_edges$Campbell_cluster=="NA"] = paste0("NA_campbell")
sankey_edges$mouse_cluster[sankey_edges$mouse_cluster=="NA"] = paste0("NA_mm")
sankey_edges$human_cluster[sankey_edges$human_cluster=="NA"] = paste0("NA_hs")

sankey_nodes1 = data.frame(name = sankey_edges$Campbell_cluster, group = "Campbell_cluster")
sankey_nodes2 = data.frame(name = sankey_edges$mouse_cluster, group = "mouse_cluster")
sankey_nodes3 = data.frame(name = sankey_edges$human_cluster, group = "human_cluster")

sankey_nodes <- rbind(rbind(sankey_nodes3,sankey_nodes2),sankey_nodes1) %>% 
  dplyr::distinct(name,group) %>% as.data.frame()

sankey_edges1 = sankey_edges %>% dplyr::select(NAMEsource = human_cluster, NAMEtarget = mouse_cluster, value = nhuman)
sankey_edges1$value[is.na(sankey_edges1$value)]=1#1000
sankey_edges1$value[sankey_edges1$value > 1000]=1000
sankey_edges1$value = sankey_edges1$value / sum(sankey_edges1$value)
sankey_edges2 = sankey_edges %>% dplyr::select(NAMEsource = mouse_cluster, NAMEtarget = Campbell_cluster, value = n)
sankey_edges2$value = sankey_edges2$value / sum(sankey_edges2$value)
sankey_edges_b = rbind(sankey_edges1,sankey_edges2)

sankey_edges_b$IDsource <- match(sankey_edges_b$NAMEsource, sankey_nodes$name) - 1
sankey_edges_b$IDtarget <- match(sankey_edges_b$NAMEtarget, sankey_nodes$name) - 1
sankey_edges_b = as.data.frame(sankey_edges_b)

# save files
data.table::fwrite(sankey_edges_b,paste0("data/campbell_comparison_sankey_links_",nlink,".tsv"),sep="\t")
data.table::fwrite(sankey_nodes,paste0("data/campbell_comparison_sankey_nodes",nlink,".tsv"),sep="\t")

#plot
text_size = 13
p <- networkD3::sankeyNetwork(Links = sankey_edges_b, Nodes = sankey_nodes, 
                              Source = "IDsource", Target = "IDtarget", Value = "value", 
                              NodeID = "name", NodeGroup = "group" , sinksRight = FALSE, fontSize = text_size)
p

#, colourScale = networkD3::JS(paste0("d3.scaleOrdinal() .domain([\"clustering_1\",\"clustering_2\",\"clustering_1_added\",\"clustering_2_added\"]) .range([\"",  col1, "\",\"", col2, "\",\"", color_vector_lighter[1],  "\",\"", color_vector_lighter[2], "\"]);")), 

# # save sankeys
library(webshot)
# https://stackoverflow.com/questions/65158327/how-can-i-save-a-networkd3sankeynetwork-into-a-static-image-automatically-vi
save_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/dump/"
networkD3::saveNetwork(p, paste0(save_path,"sankey_clusters_campbell_human_",nlink,".html"))
# convert it
# need: webshot::install_phantomjs()
webshot::webshot(paste0(save_path,"sankey_clusters_campbell_human_",nlink,".html"),file=paste0(save_path,"sankey_clusters_campbell_human_",nlink,".png"), vwidth = 1000, vheight = 900)
webshot::webshot(paste0(save_path,"sankey_clusters_campbell_human_",nlink,".html"),file=paste0(save_path,"sankey_clusters_campbell_human_",nlink,".pdf"), vwidth = 1000, vheight = 900)

# I used a local script to export with webshot !!
