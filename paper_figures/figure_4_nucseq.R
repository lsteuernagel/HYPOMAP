
##########
###  load data
##########

library(tidyverse)
library(ggplot2)
library(ggh4x)
library(ggtree)
source("utility_functions.R")
source("integration_pipeline/harmonization_functions.R")
source("merge_human_mouse_neurons/cluster_matching_functions.R")
source("paper_figures/tree_plotting_functions.R")

human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"

# load matching result
matched_clusters_final = data.table::fread("merge_human_mouse_neurons/matched_clusters_scviCorPruned.tsv",data.table = F)

# save
figure_dir = "paper_figures/revision_figures/"
dir.create(figure_dir,showWarnings = F)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

##########
###  load scseq objects
##########

# load seurat object
if(!exists("hypothalamus_neurons_cross_species", envir = .GlobalEnv )){
  hypothalamus_neurons_cross_species = readRDS(paste0(human_hypo_path,"hypothalamus_neurons_cross_species.rds")) 
}

## mouse
hypoMap = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap.rds") # requires mouse HypoMap from: https://doi.org/10.17863/CAM.87955

hypoMap_annos = hypoMap@misc$annotation_result
hypoMap_edgelist =  hypoMap@misc$clustering_edgelist %>% 
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,from_named = clean_names_withID),by=c("from"="cluster_id")) %>%
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,to_named = clean_names_withID),by=c("to"="cluster_id"))

### Add clustering results from human neurons using the final object
if(!exists("human_hypo_combined", envir = .GlobalEnv )){
  human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP.rds"))
}
## load markers
human_marker_genes = data.table::fread("paper_figures/revision_figures/suppl_tables/marker_genes_all_top100.txt",data.table = F)

# edgelist
combined_edgelist_mrtree = human_hypo_combined@misc$edgelist # data.table::fread(paste0(human_hypo_path,"human_hypo_combined_edgelist_mrtree_annotated.txt"),data.table = F)

##########
### Prepare for Circular Tree ----- Figure 2
##########

# Run Figure1_nucseq to load all relevant data and update meta.data

# need:
# TODO NEEED COLOR ORDER
# color_order from figure 1 script

edgelist = human_hypo_combined@misc$edgelist

# anno_df
# make anno_df
anno_df =human_hypo_combined@misc$edgelist %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})
# manually shorten CO level
anno_df$first_cluster_name[anno_df$first_cluster_name == "AstroEpendymal"] = "AstroEpen"
anno_df$first_cluster_name[anno_df$first_cluster_name == "Oligodendrocytes"] = "Oligo"
anno_df$first_cluster_name[anno_df$first_cluster_name == "Immune/Vascular"] = "NN"
anno_df$first_cluster_name = stringr::str_remove(anno_df$first_cluster_name,pattern = "Oligo-")

# reorder first level but keep C0 as is
tmp = anno_df[anno_df$cluster_id %in% color_order$cluster_id,]
tmp_c0  = anno_df[grepl("C0",anno_df$cluster_id ),]
anno_df[anno_df$cluster_id %in% color_order$cluster_id,]=tmp[match(color_order$cluster_id,tmp$cluster_id),]
anno_df[grepl("C0",anno_df$cluster_id ),]=tmp_c0[match(c("C0-3","C0-1","C0-2","C0-4"),tmp_c0$cluster_id),]
anno_df$order_plot = 1:nrow(anno_df)
#anno_df

##########
### Regions
##########

color_value_df =data.table::fread("data/region_anno_revision/spatial_clusters_colours_OkabeItoPalette.txt",data.table = F)
color_value_vector = color_value_df$col
names(color_value_vector) = color_value_df$V1

region_mapping = data.table::fread("data/region_anno_revision/nucseq_clusters_region_assignment.txt",data.table = F)
region_mapping = dplyr::left_join(region_mapping,color_value_df,by=c("region_threshold_manual"="V1"))

region_mapping_forOverview = region_mapping %>% 
  dplyr::left_join(anno_df,by=c("C4"="cluster_id")) %>%
  dplyr::select(cluster_name,region = region_threshold_manual, region_color = col,C4,C3)
region_mapping_forOverview$region[region_mapping_forOverview$region == "Ventricular"]="Vent"
region_mapping_forOverview$region[region_mapping_forOverview$region == "Periventricular"]="Perivent"

names(color_value_vector)[names(color_value_vector)=="Ventricular"]="Vent"
names(color_value_vector)[names(color_value_vector)=="Periventricular"]="Perivent"


##########
### Overview UMAPs ----- Figure 1
##########

names(hypothalamus_neurons_cross_species@reductions)
colnames(hypothalamus_neurons_cross_species@meta.data)

# by dataset
species_umap = DimPlot(hypothalamus_neurons_cross_species,group.by = "species",reduction = "umap_cross_species_neurons_scvi_1_scVI_reduction",
                       pt.size = 1.5,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(4),label=F,label.size = 2.5,shuffle = TRUE)+
  NoAxes()+ggtitle(NULL)+NoAxes()#+NoLegend()
species_umap
## save
ggsave(filename = paste0(figure_dir,"EXT_06_b_umap_species.pdf"),plot = species_umap, "pdf",dpi=400,width=320,height = 300,units="mm")

## human clusters
human_umap = DimPlot(hypothalamus_neurons_cross_species,group.by = "C4",reduction = "umap_cross_species_neurons_scvi_1_scVI_reduction",
                     pt.size = 2,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(length(unique(hypothalamus_neurons_cross_species@meta.data$C4))),na.value = "grey80",order=TRUE,label=F,label.size = 2.5,shuffle = F)+
  NoAxes()+ggtitle(NULL)+NoLegend()
human_umap

ggsave(filename = paste0(figure_dir,"EXT_06_c_umap_species_humanclusters.pdf"),
       plot = human_umap, "pdf",dpi=400,width=300,height = 300,units="mm")

## mouse clusters
mouse_umap = DimPlot(hypothalamus_neurons_cross_species,group.by = "C465_named",reduction = "umap_cross_species_neurons_scvi_1_scVI_reduction",
                     pt.size = 2,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(length(unique(hypothalamus_neurons_cross_species@meta.data$C465_named))),na.value = "grey80",order=TRUE,label=F,label.size = 2.5,shuffle = F)+
  NoAxes()+ggtitle(NULL)+NoLegend()
mouse_umap

ggsave(filename = paste0(figure_dir,"EXT_06_d_umap_species_mouseclusters.pdf"),
       plot = mouse_umap, "pdf",dpi=400,width=300,height = 300,units="mm")


##########
### General mapping statistics
##########

# no, 1:1, 1:many

relationship_data2 = matched_clusters_final %>% 
  as.data.frame() %>%
  dplyr::mutate(edge_type = case_when(
    from_occ > 1 & to_occ == 1 ~ "1:n",
    from_occ > 1 & to_occ > 1 ~ "PROBLEM", # special case
    from_occ == 1 & to_occ > 1 ~ "m:1", # means multiple human match to one mouse
    from_occ == 1 & to_occ == 1 ~ "1:1",
    TRUE ~ "other"
  ))

# numbers
#length(unique())

##########
### General mapping statistics 
##########

# make type stat for each human (from) node
relationship_data_stat_from = relationship_data2 %>%
  dplyr::group_by(from,edge_type) %>%
  dplyr::count() %>%
  dplyr::distinct(from,edge_type,n)

## make also a version with missing:
relationship_data_stat_from2 = data.frame(cluster_name = unique(human_hypo_combined@meta.data$C4_named)) %>%
  dplyr::left_join(relationship_data_stat_from,by=c("cluster_name"="from")) %>%
  dplyr::filter(grepl("GABA|GLU",cluster_name))
length(unique(relationship_data_stat_from2$cluster_name))

# %>%
# dplyr::left_join(anno_df,by=c("cluster_id"="cluster_id"))
relationship_data_stat_from2$n[is.na(relationship_data_stat_from2$n)] = 0
relationship_data_stat_from2$edge_type[relationship_data_stat_from2$n == 0] = "human only"

## look into Mid-2 -- Mediobasal neurons
mapping_stats_2=relationship_data_stat_from2[grepl("Mid-2",relationship_data_stat_from2$cluster_name),]
table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)])
round(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)]) / sum(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)])) * 100,3)
# DMH PPP1R17 neurons
# HDC and related WNT7B neurons 

# mid 1
mapping_stats_1=relationship_data_stat_from2[grepl("Mid-1",relationship_data_stat_from2$cluster_name),]
table(relationship_data_stat_from2$edge_type[grepl("Mid-1",relationship_data_stat_from2$cluster_name)])
round(table(relationship_data_stat_from2$edge_type[grepl("Mid-1",relationship_data_stat_from2$cluster_name)]) / sum(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)])) * 100,3)
# some DMH (Hmx2/hmx3 neurons)
# Lhx1/Lhx6 neurons
# TRPC4 NTS is the only cell type with some SHOX2 expression remaining ...
# VIP neurons are found in both species, here clster matching likely mismatched clusters across species, partly owing to different patterns of tf expression such as ZIC1,2,4,5

round(table(relationship_data_stat_from2$edge_type[grepl("Mid-3",relationship_data_stat_from2$cluster_name)]) / sum(table(relationship_data_stat_from2$edge_type[grepl("Mid-3",relationship_data_stat_from2$cluster_name)])) * 100,3)
mapping_stats_3=relationship_data_stat_from2[grepl("Mid-3",relationship_data_stat_from2$cluster_name),]

table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)])
round(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)]) / sum(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)])) * 100,3)
## all
table(relationship_data_stat_from2$edge_type)
round(table(relationship_data_stat_from2$edge_type) / sum(table(relationship_data_stat_from2$edge_type)) * 100,3)
length(unique(relationship_data_stat_from2$cluster_name))
#length(unique(relationship_data_stat_from2$))

## other norm ????
round(table(relationship_data_stat_from2$edge_type[grepl("Mid-2",relationship_data_stat_from2$cluster_name)]) / length(unique(human_hypo_combined@meta.data$C4_named)[grepl("Mid-2",unique(human_hypo_combined@meta.data$C4_named))]) * 100,3)
round(table(relationship_data_stat_from2$edge_type) / length(unique(human_hypo_combined@meta.data$C4_named)[grepl("GABA|GLU",unique(human_hypo_combined@meta.data$C4_named))]) * 100,3)

# make overview for barplot
relationship_data_all = relationship_data2  %>%
  dplyr::distinct(from,to,edge_type,.keep_all = TRUE) %>%
  dplyr::group_by(edge_type)%>%
  dplyr::count()  

# count missing ones
toadd = data.frame(edge_type=c("only mouse","only human"),n = c(length(setdiff(unique(hypoMap@meta.data$C465_named[grepl("GABA|GLU",hypoMap@meta.data$C465_named)]),matched_clusters_final$to)),length(setdiff(unique(human_hypo_combined@meta.data$C4_named)[grepl("GABA|GLU",unique(human_hypo_combined@meta.data$C4_named))],matched_clusters_final$from))))
relationship_data_all = bind_rows(relationship_data_all,toadd)
relationship_data_all
# numers
# length(unique(cor_edgelist3$to))
# length(unique(cor_edgelist3$from))
# length(setdiff(unique(cor_edgelist3$to),matched_clusters_final$to)) / length(unique(cor_edgelist3$to))*100
# length(setdiff(unique(cor_edgelist3$from),matched_clusters_final$from)) / length(unique(cor_edgelist3$from))*100

length(unique(human_hypo_combined@meta.data$C4_named[grepl("GABA|GLU",human_hypo_combined@meta.data$C4_named)]))

setdiff(unique(hypoMap@meta.data$C465_named[grepl("GABA|GLU",hypoMap@meta.data$C465_named)]),matched_clusters_final$to)

length(find_children(nodes = "C0-3",edges = combined_edgelist_mrtree)[grepl("C4",find_children(nodes = "C0-3",edges = combined_edgelist_mrtree))])
length(combined_edgelist_mrtree$to[combined_edgelist_mrtree$level==6]) # 413 in combined neurons !

relationship_data_all$edge_type = factor(relationship_data_all$edge_type,levels = c("only human","m:1", "1:1","1:n","only mouse"))
# plot as abrplot or similar
text_size = 25
barplot = ggplot(relationship_data_all,aes(x=edge_type,y=n,fill = edge_type))+geom_col()+
  scale_fill_manual(values=c("#E69F00","#e0c68d","#0072B2","#659487","#009E73"))+ # ,"#9c9c9c"
  geom_text(aes(label = n),nudge_y = 5,size=10)+
  theme_light()+xlab(NULL)+ylab("Type of match")+NoLegend()+
  theme(text = element_text(size=text_size),panel.border = element_blank(),axis.line=element_line(color="grey80",linewidth = 0.5))# ,panel.background = element_rect(fill = "white"),axis.line=element_line(color="grey40",size = 1),legend.key=element_rect(fill="white"))+

barplot

sum(relationship_data_all$n)
#184+102+106

ggsave(filename = paste0(figure_dir,"EXT_06_x_barplot_species_matchtype.pdf"),
       plot = barplot, "pdf",dpi=200,width=280,height = 200,units="mm")

##########
### Plot on UMAP
##########

# make type stat for each mouse (to) node
relationship_data_stat_to = relationship_data2 %>%
  dplyr::group_by(to,edge_type) %>%
  dplyr::count()  %>%
  dplyr::distinct(to,edge_type,n)

# make overview for barplot
relationship_data_all = relationship_data2  %>%
  dplyr::distinct(from,to,edge_type,.keep_all = TRUE) %>%
  dplyr::group_by(edge_type)%>%
  dplyr::count()  

# plot on UMAP -- use per cluster value
anno_crossspecies = hypothalamus_neurons_cross_species@meta.data %>% dplyr::select(Cell_ID,C465_named,C4_named)
anno_crossspecies$C4_named[is.na(anno_crossspecies$C4_named)] = ""
anno_crossspecies$C465_named[is.na(anno_crossspecies$C465_named)] = ""
anno_crossspecies$anno_both = anno_crossspecies$C465_named
anno_crossspecies$anno_both[anno_crossspecies$anno_both == ""] = anno_crossspecies$C4_named[anno_crossspecies$anno_both == ""]
#anno_crossspecies$anno_both[is.na(anno_crossspecies$anno_both)] = anno_crossspecies$C4[is.na(anno_crossspecies$anno_both)]

#length(unique(cor_edgelist3$to))+length(unique(cor_edgelist3$from))
relationship_data_stat_both = bind_rows(relationship_data_stat_from %>% dplyr::rename(anno_both = from ) , relationship_data_stat_to %>% dplyr::rename(anno_both = to )) %>% as.data.frame()
anno_crossspecies = dplyr::left_join(anno_crossspecies, relationship_data_stat_both,by=c("anno_both"="anno_both"))
anno_crossspecies$edge_type[is.na(anno_crossspecies$edge_type) & anno_crossspecies$C4_named == ""] = "only mouse"
anno_crossspecies$edge_type[is.na(anno_crossspecies$edge_type) & anno_crossspecies$C465_named == ""] = "only human"

# add via assignment 
hypothalamus_neurons_cross_species@meta.data$edge_type = anno_crossspecies$edge_type
hypothalamus_neurons_cross_species@meta.data$edge_type = factor(hypothalamus_neurons_cross_species@meta.data$edge_type,levels = c("only human","m:1", "1:1","1:n","only mouse"))
## human clusters
edgetype_umap = DimPlot(hypothalamus_neurons_cross_species,group.by = "edge_type",reduction = "umap_cross_species_neurons_scvi_1_scVI_reduction",
                        pt.size = 1.5,raster.dpi = c(2048,2048),cols = c("#E69F00","#e0c68d","#0072B2","#659487","#009E73"),order=F,label=F,label.size = 2.5,shuffle = TRUE)+
  NoAxes()+ggtitle(NULL)#+NoLegend()

edgetype_umap

ggsave(filename = paste0(figure_dir,"EXT_06_x_umap_species_edgetype.pdf"),
       plot = edgetype_umap, "pdf",dpi=400,width=330,height = 300,units="mm")


##########
### Tree plot
##########

# neuron edgelist
# all_neurons = scUtils::find_children(nodes = "C0-1",edges = combined_edgelist_mrtree)
# combined_edgelist_mrtree_neurons = combined_edgelist_mrtree[combined_edgelist_mrtree$from %in% all_neurons | combined_edgelist_mrtree$to %in% all_neurons,]

##### plot cluster tree:
tree_color = "grey70"

##### plot cluster tree:
tree_color = "grey70"
label_color = "black"
leaf_level = 6
reduction_constant_last_level = 0.3

circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=human_hypo_combined@meta.data,
                                  #      branch_length_base = 0.1,
                                  label_size = 4, 
                                  label_size_tip = 3,
                                  show_genes = FALSE,
                                  show_genes_tip = FALSE,
                                  vjust_label = -0.25,
                                  nudge_x_tip = -0.15,
                                  annotate_reverse =F,
                                  edge_color = tree_color, 
                                  node_color = tree_color,
                                  label_color = label_color,
                                  label_lowest = FALSE,
                                  order_column = "order_plot"
)
circular_tree = ggtree::rotate_tree(circular_tree, -90)

circular_tree

##########
### Add bg colors
##########

tree_data =circular_tree$data

### cluster coloring 
data_bg_fill <- color_order
data_bg_fill$node =tree_data$node[match(data_bg_fill$cluster_id,tree_data$label)]
data_bg_fill$fill = data_bg_fill$color
circular_tree2 <- circular_tree +
  ggtree::geom_hilight(
    data = data_bg_fill,
    mapping = aes(
      node = node,
      fill = I(fill)
    ),
    align = "none",
    alpha = 0.4,
    extendto = 4
  )

# move recatngles as first layer so that tree is plotted on top
circular_tree2$layers = c(circular_tree2$layers[[length(circular_tree2$layers)]],circular_tree2$layers[1:(length(circular_tree2$layers)-1)])

# adjust xmin and xmax so that it looks better -- somehow need to invert here !
# circular_tree2$layers[[1]]$data$xmin = 3.6#circular_tree2$layers[[1]]$data$xmin -1
# circular_tree2$layers[[1]]$data$xmax = 0.5
circular_tree2$layers[[1]]$data$xmin = 0.5
save = circular_tree2$layers[[1]]$data$xmax
circular_tree2$layers[[1]]$data$xmax = circular_tree2$layers[[1]]$data$xmin
circular_tree2$layers[[1]]$data$xmin = save

circular_tree2

##########
### Add heatmap with type
##########

leaf_level_column = "C4"

heatmap_data1 = anno_crossspecies %>%
  dplyr::filter(C4_named != "") %>%
  dplyr::distinct(C4_named, edge_type) %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::select(to_named,from = to),by=c("C4_named"="to_named")) %>%
  dplyr::select(from, type = edge_type)

heatmap_matrix1 = as.matrix(heatmap_data1[,2,drop=FALSE])
rownames(heatmap_matrix1) = heatmap_data1[,"from"]
#heatmap_matrix2[is.na(heatmap_matrix2)] = 0


# plot tree with heatmap 1
bg_col = "grey90"

circular_tree_heat = hypoMapUtils::add_heatmap(circular_tree=circular_tree2,
                                               heatmap_matrix = heatmap_matrix1,
                                               heatmap_colors=c(bg_col,"darkred"),
                                               scale_limits = c(0,1),
                                               heatmap_colnames =FALSE, 
                                               legend_title = "Type",
                                               matrix_offset = 0.1,
                                               matrix_width =0.3,
                                               colnames_angle=0,
                                               legend_text_size = 5,
                                               hjust_colnames=1,
                                               na_color = "white",
                                               heatmap_text_size=2)+
  scale_fill_manual(values = c("only human" = "#E69F00","m:1" = "#e0c68d", "1:1" = "#0072B2","1:n" = "#659487","only mouse" = "#009E73" ),na.value = "white")
#scale_fill_gradient2(low = "orange",mid = "yellow",high = "darkgreen",na.value = "white",midpoint = 0.8,limits=c(0.6,1),name = "Species Cor")

circular_tree_heat$labels$fill <- "Type"
#circular_tree_heat
#circular_tree_heat

##########
### Add heatmap
##########

leaf_level_column = "C4"

#  USE matched_clusters_final !!!!!

#all_edges_filt2_nfix = all_edges_filt2[all_edges_filt2$from %in% c("C6-122","C6-123","C6-124","C6-125"),]

heatmap_data2 = matched_clusters_final %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::select(to_named,from_id = to),by=c("from"="to_named")) %>% 
  dplyr::group_by(from) %>% 
  dplyr::arrange(desc(similarity)) %>%
  dplyr::mutate(within_id = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::select(from = from_id, within_id, similarity) %>%
  dplyr::filter(within_id <= 5) %>%
  dplyr::mutate(similarity = round(similarity,3)) %>%
  tidyr::spread(key=within_id,value = similarity) %>%
  as.data.frame()

heatmap_matrix2 = as.matrix(heatmap_data2[,2:6,drop=FALSE])
#colnames(heatmap_matrix2) = "%"
rownames(heatmap_matrix2) = heatmap_data2[,"from"]
#heatmap_matrix2[is.na(heatmap_matrix2)] = 0


# plot tree with heatmap 1
bg_col = "grey90"

circular_tree_heat2 = hypoMapUtils::add_heatmap(circular_tree=circular_tree_heat,
                                                heatmap_matrix = heatmap_matrix2,
                                                heatmap_colors=c(bg_col,"darkred"),
                                                scale_limits = c(0,1),
                                                heatmap_colnames =FALSE, 
                                                legend_title = "Pct Dataset",
                                                matrix_offset = 2,
                                                matrix_width = 1.5,
                                                colnames_angle=0,
                                                legend_text_size = 5,
                                                hjust_colnames=1,
                                                na_color = "white",
                                                heatmap_text_size=2)+
  scale_fill_gradient2(low = "orange",mid = "yellow",high = "darkgreen",na.value = "white",midpoint = 0.85,limits=c(0.7,1),name = "Species Cor")

circular_tree_heat2 = circular_tree_heat2+theme(legend.text = element_text(size=12))

circular_tree_heat2


##########
### Save tree
##########

# shorter was with matrix_width=1 for second heatmap !
## save
ggsave(filename = paste0(figure_dir,"FIG_03_A_circular_tree_mouseCorrelations_short.pdf"),  plot = circular_tree_heat2, "pdf",dpi=600,width=400,height = 400,units="mm")


##########
### POMC plot
##########

matched_clusters_pomc = matched_clusters_final[grepl("POMC",matched_clusters_final$from),]
pomc_comparison_plot = plot_tree_comparison(matched_clusters_pomc,edgelist_human = combined_edgelist_mrtree, edgelist_mouse = hypoMap_edgelist ,include_extra_level_human =TRUE,general_label_size = 5,min_similarity=0.7,color_order_regex = "C4",start_level = "C2")+
  ggtitle(NULL)

pomc_comparison_plot

ggsave(filename = paste0(figure_dir,"FIG_03_d_pomc_phylotree_crossspecies.pdf"),
       plot = pomc_comparison_plot, "pdf",dpi=300,width=300,height = 150,units="mm")

##########
### GPCR plots
##########

## see the separate script !



