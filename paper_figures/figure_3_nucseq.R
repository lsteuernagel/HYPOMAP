##########
### Load data
##########

figure_dir = "paper_figures/revision_figures/"
dir.create(figure_dir)

# load files
human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"
library(Seurat)
library(tidyverse)
library(ggtree)

human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP.rds"))
combined_edgelist_mrtree = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_edgelist_mrtree_annotated.txt"),data.table = F)
comparisons_all_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_all_annotated.txt"),data.table = F)
comparisons_siblings_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_siblings_annotated.txt"),data.table = F)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

# source
source("paper_figures/tree_plotting_functions.R")
source("paper_figures/celltype_heatmaps_functions.R")

##########
### VMH small tree plot
##########

### VMH neurons
node_to_use = "C2-40"
remove_name_part =  paste0("C[0-4]-[0-9]+ ",gsub("C[0-4]-[0-9]+ ","",anno_df$cluster_name[anno_df$cluster_id == node_to_use])," ") #paste0(anno_df$cluster_name[anno_df$cluster_id == node_to_use]," ")
additional_genes = c("NR5A1","FEZF1")
additional_nodes =c()
ngenes_factor= 14
############## avg.exp capped

## make avg.exp plots
current_plots_expr  = make_cluster_minitree_helper(node_to_use=node_to_use,keep_quantile = 0.975,hjust_colnames=0.3,heatmap_text_size = 4,
                                                   additional_nodes = additional_nodes,remove_name_part=remove_name_part,additional_genes=additional_genes,vjust_label = -0.5,hjust_label = 0.35,
                                                   label_size=5,label_size_tip=5,tree_nodesize = 5, human_hypo_combined = human_hypo_combined,combined_edgelist_mrtree = combined_edgelist_mrtree,
                                                   expr_type = "avg.exp",avg_exp_cap = 2) # pct.exp or avg.exp or avg.exp.scaled

# show
current_plots_expr$tree_heatmap#+theme(plot.margin = margin(2,2,10,2, "cm"))
current_plots_expr$cluster_umap
current_plots_expr$overview_umap

# save plots
filename = paste0("FIG_02_d_","VMH","_","avgExpCapped","_")
ggsave(filename = paste0(figure_dir,filename,"tree_heatmap",".pdf"),plot = current_plots_expr$tree_heatmap, "pdf",dpi=300,width=200+ngenes_factor*nrow(current_plots_expr$heatmap_data),height = 250,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"cluster_umap",".pdf"),plot = current_plots_expr$cluster_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"overview_umap",".pdf"),plot = current_plots_expr$overview_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
data.table::fwrite(current_plots_expr$heatmap_data, paste0(figure_dir,filename,"heatmap_data",".txt"),sep = "\t")


##########
### VMH small tree plot JUST C3
##########

### VMH neurons
node_to_use = "C2-40"
remove_name_part = ""  #paste0("C[0-4]-[0-9]+ ",gsub("C[0-4]-[0-9]+ ","",anno_df$cluster_name[anno_df$cluster_id == node_to_use])," ") #paste0(anno_df$cluster_name[anno_df$cluster_id == node_to_use]," ")
additional_genes = c("FEZF1","NR5A1","BDNF","ESR1","NOS1","ADCYAP1")
additional_nodes=c()

############## avg.exp capped

## make avg.exp plots
vmh_c3_plots  = make_cluster_minitree_helper(node_to_use=node_to_use,keep_quantile = 0.975,hjust_colnames=0.3,heatmap_text_size = 6,
                                             additional_nodes = c(),remove_name_part=remove_name_part,
                                             additional_genes=additional_genes,vjust_label = -0.5,hjust_label = 0.35,
                                             label_size=6,label_size_tip=6,tree_nodesize = 6, 
                                             human_hypo_combined = human_hypo_combined,combined_edgelist_mrtree = combined_edgelist_mrtree,
                                             expr_type = "avg.exp",avg_exp_cap = 2,hilight_col = "blue",
                                             heatmap_offset_base = 3, matrix_width_per_gene = 0.3,
                                             leaf_level=5,final_name_column="C3_named",seurat_ident_column="C3") # pct.exp or avg.exp or avg.exp.scaled

# show
vmh_c3_plots$tree_heatmap#+theme(plot.margin = margin(2,2,10,2, "cm"))
vmh_c3_plots$cluster_umap
vmh_c3_plots$overview_umap

# save plots
ngenes_factor= 90
filename = paste0("FIG_02_d_","VMH_C3","_","avgExpCapped","_")
ggsave(filename = paste0(figure_dir,filename,"tree_heatmap",".pdf"),plot = vmh_c3_plots$tree_heatmap, "pdf",dpi=300,width=200+ngenes_factor*nrow(vmh_c3_plots$heatmap_data),height = 250,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"cluster_umap",".pdf"),plot = vmh_c3_plots$cluster_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
ggsave(filename = paste0(figure_dir,filename,"overview_umap",".pdf"),plot = vmh_c3_plots$overview_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
data.table::fwrite(vmh_c3_plots$heatmap_data, paste0(figure_dir,filename,"heatmap_data",".txt"),sep = "\t")


