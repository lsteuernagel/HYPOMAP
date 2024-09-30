##########
### Load data
##########

figure_dir = "paper_figures/revision_figures/"
dir.create(figure_dir,showWarnings = FALSE)

# source
source("paper_figures/tree_plotting_functions.R")
source("paper_figures/celltype_heatmaps_functions.R")

# load files
human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/" # change to local path
library(Seurat)
library(tidyverse)
library(ggtree)

human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP.rds"))
combined_edgelist_mrtree = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_edgelist_mrtree_annotated.txt"),data.table = F) # available from paper supplementary tables
comparisons_all_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_all_annotated.txt"),data.table = F)
comparisons_siblings_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_siblings_annotated.txt"),data.table = F)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

##########
### UMAP with annotation
##########

colnames(human_hypo_combined@meta.data)

# human_hypo_combined@meta.data$celltype_annotation[ human_hypo_combined@meta.data$celltype_annotation == "P2RX2_OTP"] = "Neurons"
# human_hypo_combined@meta.data$celltype_annotation[ human_hypo_combined@meta.data$celltype_annotation == "Tmem119.Immune"] = "Microglia"
length(unique(human_hypo_combined@meta.data$C1_named))
DimPlot(human_hypo_combined,group.by = "C1_named",reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(13),label = TRUE,label.size = 7,repel = TRUE)+NoAxes()+NoLegend()

DimPlot(human_hypo_combined,group.by = "C2_named",reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(52),label = TRUE,label.size = 7,repel = TRUE)+NoAxes()+NoLegend()

# FeaturePlot(human_hypo_combined,"SIM1",reduction = "umap",order = TRUE,raster.dpi = c(1536,1536),pt.size = 1.5)+NoAxes()

##########
### Prepare final annotations
##########

## make annotation data frame
color_order = combined_edgelist_mrtree %>% dplyr::filter(grepl("C1",to)) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
#anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
color_order$first_cluster_name = sapply(color_order$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})

## prepare colors and order
color_order = color_order[,c("cluster_id","cluster_name","first_cluster_name")]

# automatically assign colors:
color_order$color = getOkabeItoPalette(nrow(color_order))

# manually assign colors:
color_order$color = NA
color_order$color[color_order$first_cluster_name=="Astrocytes"] = "#D06B53"
color_order$color[color_order$first_cluster_name=="Ependymal"] = "#D55E00"
color_order$color[color_order$first_cluster_name=="Oligo-Precursor"] = "#9EA974"
color_order$color[color_order$first_cluster_name=="Oligo-Mature"] = "#6A6858"
color_order$color[color_order$first_cluster_name=="Vascular"] = "#F0E442"
color_order$color[color_order$first_cluster_name=="Immune"] = "#CC79A7"

color_order$color[color_order$first_cluster_name=="Pre-1"] = "#E69F00"
color_order$color[color_order$first_cluster_name=="Pre-2"] ="#0072B2"
color_order$color[color_order$first_cluster_name=="Mid-1"] = "#56B4E9" #
color_order$color[color_order$first_cluster_name=="Mid-2"] = "#009E73"
color_order$color[color_order$first_cluster_name=="Mid-3"] = "#2AA9AD" # 
color_order$color[color_order$first_cluster_name=="Post-1"] = "#77C15A" # 
color_order$color[color_order$first_cluster_name=="Post-2"] = "#78AB79"

# setdiff(unique(color_order$color),getOkabeItoPalette(nrow(color_order)))
# setdiff(getOkabeItoPalette(nrow(color_order)),unique(color_order$color))
# make factor

color_order$final_name = factor(color_order$first_cluster_name,levels = c("Vascular","Immune","Oligo-Precursor","Oligo-Mature","Astrocytes","Ependymal","Pre-1", "Pre-2","Mid-1", "Mid-2", "Mid-3","Post-1", "Post-2"))
color_order = color_order %>% dplyr::arrange((final_name)) %>% as.data.frame()
color_order_vec = unlist(split(color_order$color,color_order$cluster_name))
#human_hypo_combined@meta.data$final_name = factor(human_hypo_combined@meta.data$final_name,levels = color_order$final_name)

DimPlot(human_hypo_combined,group.by = "C1_named",reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(13),label = TRUE,label.size = 7,repel = TRUE)+
  scale_color_manual(values = color_order_vec)+
  NoAxes()+NoLegend()


##########
### Base Tree plot
##########

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
anno_df

edgelist = human_hypo_combined@misc$edgelist

##### plot cluster tree:
tree_color = "grey70"
label_color = "black"
leaf_level = 5
reduction_constant_last_level = 0.3

circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=leaf_level,
                                  anno_df = anno_df ,
                                  metadata=human_hypo_combined@meta.data,
                                  #      branch_length_base = 0.1,
                                  label_size = 4, 
                                  label_size_tip = 3,
                                  show_genes = TRUE,
                                  show_genes_tip = TRUE,
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

# circular_tree$data$branch = circular_tree$data$branch*0.2
# circular_tree$data$x = circular_tree$data$x*0.2

# make branches shorter
circular_tree$data$clusterlevel[is.na(circular_tree$data$clusterlevel)] = "all"

# clevel="C1"
# circular_tree$data$branch[circular_tree$data$clusterlevel==clevel] = circular_tree$data$branch[circular_tree$data$clusterlevel==clevel] - reduction_constant_last_level
# circular_tree$data$x[circular_tree$data$clusterlevel==clevel] = circular_tree$data$x[circular_tree$data$clusterlevel==clevel]- reduction_constant_last_level
# clevel="C2"
# circular_tree$data$branch[circular_tree$data$clusterlevel==clevel] = circular_tree$data$branch[circular_tree$data$clusterlevel==clevel]- reduction_constant_last_level
# circular_tree$data$x[circular_tree$data$clusterlevel==clevel] = circular_tree$data$x[circular_tree$data$clusterlevel==clevel]- reduction_constant_last_level

# also tip
circular_tree$data$branch[circular_tree$data$isTip] = circular_tree$data$branch[circular_tree$data$isTip]-reduction_constant_last_level
circular_tree$data$x[circular_tree$data$isTip] = circular_tree$data$x[circular_tree$data$isTip]-reduction_constant_last_level


# plot
circular_tree


##########
### Add backgound fill 
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
### Add heatmaps
##########

leaf_level_column = "C3"
# make data for second heatmap with n cells
# heatmap_data2 = human_hypo_combined@meta.data %>% dplyr::select(Cell_ID,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column)) %>% #dplyr::filter(predicted_Campbell!="NA") 
#   dplyr::count(name = "ncells")  %>% dplyr::ungroup()  %>% dplyr::mutate(ncells_pct = ncells / sum(ncells)*100)  %>% as.data.frame()
heatmap_data2 = human_hypo_combined@meta.data %>% dplyr::select(Cell_ID,!!sym(leaf_level_column),Dataset) %>% 
  dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::add_count(name = "cluster_ncells")  %>% 
  dplyr::group_by(!!sym(leaf_level_column),Dataset) %>% dplyr::add_count(name = "ncells")  %>% 
  dplyr::ungroup() %>% dplyr::distinct(!!sym(leaf_level_column),Dataset,.keep_all = TRUE) %>% 
  dplyr::mutate(ncells_pct = round(ncells / cluster_ncells *100,2))  %>% 
  dplyr::select(!!sym(leaf_level_column),Dataset,ncells_pct) %>% tidyr::spread(key = Dataset,value = ncells_pct) %>% as.data.frame() 

heatmap_matrix2 = as.matrix(heatmap_data2[,c("Tadross","Siletti"),drop=FALSE])
#colnames(heatmap_matrix2) = "%"
rownames(heatmap_matrix2) = heatmap_data2[,leaf_level_column]
heatmap_matrix2[is.na(heatmap_matrix2)] = 0


# plot tree with heatmap 1
bg_col = "grey90"

circular_tree_heat = hypoMapUtils::add_heatmap(circular_tree=circular_tree2,
                                               heatmap_matrix = heatmap_matrix2,
                                               heatmap_colors=c(bg_col,"darkred"),
                                               scale_limits = c(0,100),
                                               heatmap_colnames =TRUE, 
                                               legend_title = "Pct Dataset",
                                               matrix_offset = 0.05,
                                               matrix_width =0.05,
                                               colnames_angle=0,
                                               legend_text_size = 5,
                                               hjust_colnames=1,
                                               na_color = "white",
                                               heatmap_text_size=2)
circular_tree_heat

#circular_tree_heat + geom_hilight(node=358, fill='purple', color='white', alpha=0.5,linewidth=0)
# need correct node:
node_to_highlight = tree_data$node[tree_data$label == "C2-36"]
circular_tree_heat_highlight = circular_tree_heat +  geom_cladelab(node=node_to_highlight, label="", align=TRUE, offset = 0.01, textcolor='blue', barcolor='blue',barsize = 1.3)

circular_tree_heat_highlight$labels$fill <- "Pct"
circular_tree_heat_highlight

## save
ggsave(filename = paste0(figure_dir,"FIG_01_b_circular_tree_heatmap.pdf"),
       plot = circular_tree_heat_highlight, "pdf",dpi=600,width=400,height = 400,units="mm")

#### source data tree
# figure2_tree_heatmap_sourcedata = heatmap_data %>% bind_cols(heatmap_data2 %>% dplyr::select(-C185))  %>% bind_cols(heatmap_data3 %>% dplyr::select(-C185,-n))
# data.table::fwrite(figure2_tree_heatmap_sourcedata,paste0(results_path_figure2,"source_figure2_a_tree_heatmap.txt"),sep="\t")
# data.table::fwrite(hypoMap_v2_seurat@misc$clustering_edgelist,paste0(results_path_figure2,"source_figure2_a_tree_edgelist.txt"),sep="\t")
# data.table::fwrite(anno_df,paste0(results_path_figure2,"source_figure2_a_tree_annotations.txt"),sep="\t")

##########
### Corresponding UMAP plot
##########

#plot
#human_hypo_combined@meta.data$final_name = factor(human_hypo_combined@meta.data$final_name,levels = color_order$final_name)
# cols = color_order$color,
final_anno_umap = DimPlot(human_hypo_combined,group.by = "C1_named",label = TRUE,reduction = "umap",raster.dpi = c(2048,2048),pt.size = 1.7,shuffle = TRUE,repel = TRUE,label.size = 5)+
  scale_color_manual(values = color_order_vec)+NoAxes()+NoLegend()+ggtitle(NULL)

final_anno_umap

## save
ggsave(filename = paste0(figure_dir,"FIG_01_a_umap_annotated_celltypes.pdf"),
       plot = final_anno_umap, "pdf",dpi=400,width=300,height = 300,units="mm")


## plot with no prefix
target_level = "C1"
target_short = paste0(target_level,"_noPrefix")
color_order_vec_short = color_order_vec
names(color_order_vec_short) =  stringr::str_remove(names(color_order_vec_short),pattern = paste0(target_level,"-[0-9]+"))
human_hypo_combined@meta.data[,target_short] = stringr::str_remove(human_hypo_combined@meta.data[,paste0(target_level,"_named")],pattern = paste0(target_level,"-[0-9]+"))
final_anno_umap = DimPlot(human_hypo_combined,group.by = target_short,label = TRUE,cols = color_order$color,reduction = "umap",raster.dpi = c(2048,2048),pt.size = 1.7,shuffle = TRUE,repel = TRUE,label.size = 5)+
  scale_color_manual(values = color_order_vec_short)+NoAxes()+NoLegend()+ggtitle(NULL)+NoAxes()+NoLegend()+ggtitle(NULL)

final_anno_umap

## save
ggsave(filename = paste0(figure_dir,"FIG_01_a_umap_annotated_celltypes_shortLabels.pdf"),
       plot = final_anno_umap, "pdf",dpi=400,width=300,height = 300,units="mm")


#plot without labels
#human_hypo_combined@meta.data$final_name = factor(human_hypo_combined@meta.data$final_name,levels = color_order$final_name)
final_anno_umap = DimPlot(human_hypo_combined,group.by = "C1_named",label = F,cols = color_order$color,reduction = "umap",raster.dpi = c(2048,2048),pt.size = 1.7,shuffle = TRUE,repel = TRUE,label.size = 5)+
  scale_color_manual(values = color_order_vec)+NoAxes()+NoLegend()+ggtitle(NULL)+NoAxes()+NoLegend()+ggtitle(NULL)

final_anno_umap

## save
ggsave(filename = paste0(figure_dir,"FIG_01_a_umap_annotated_celltypes_noLabels.pdf"),
       plot = final_anno_umap, "pdf",dpi=400,width=300,height = 300,units="mm")


##########
### PVH Corresponding UMAP plot
##########
# source
source("paper_figures/tree_plotting_functions.R")
source("paper_figures/celltype_heatmaps_functions.R")
# 
# figure_dir = "paper_figures/final_figures/pvh_plots/"
# dir.create(figure_dir)

ngenes_factor= 17

## PVH cluster
# C2-36 would be old ones
node_to_use = "C2-36"

remove_name_part =  paste0("C[0-4]-[0-9]+ ",gsub("C[0-4]-[0-9]+ ","",anno_df$cluster_name[anno_df$cluster_id == node_to_use])," ") #paste0(anno_df$cluster_name[anno_df$cluster_id == node_to_use]," ")
additional_genes = c()#c("SIM1")#c("SIM1")
#additional_genes = c("SST","OTP")
additional_nodes = character(0)

############## avg.exp capped

## make avg.exp plots
current_plots_expr  = make_cluster_minitree_helper(node_to_use=node_to_use,keep_quantile = 0.975,hjust_colnames=0.3,heatmap_text_size = 4,gene_col = "black",
                                                   additional_nodes = additional_nodes,remove_name_part=remove_name_part,additional_genes=additional_genes,vjust_label = -0.5,hjust_label = 0.35,
                                                   label_size=5,label_size_tip=5,tree_nodesize = 5, human_hypo_combined = human_hypo_combined,combined_edgelist_mrtree = combined_edgelist_mrtree,
                                                   expr_type = "avg.exp",avg_exp_cap = 2,use_custom_scale = F,hilight_col="blue") # pct.exp or avg.exp or avg.exp.scaled

# show
current_plots_expr$tree_heatmap#+theme(plot.margin = margin(2,2,10,2, "cm"))
current_plots_expr$cluster_umap
current_plots_expr$overview_umap

# save plots
filename = paste0("FIG_01_c_",gsub("-","_",node_to_use),"_","avgExpCapped","_")
ggsave(filename = paste0(figure_dir,filename,"tree_heatmap",".pdf"),plot = current_plots_expr$tree_heatmap, "pdf",dpi=300,width=200+ngenes_factor*nrow(current_plots_expr$heatmap_data),height = 250,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"cluster_umap",".pdf"),plot = current_plots_expr$cluster_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"overview_umap",".pdf"),plot = current_plots_expr$overview_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
data.table::fwrite(current_plots_expr$heatmap_data, paste0(figure_dir,filename,"heatmap_data",".txt"),sep = "\t")


