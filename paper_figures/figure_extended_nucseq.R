
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

# load matching result
matched_clusters_final = data.table::fread("merge_human_mouse_neurons/matched_clusters_scviCorPruned.tsv",data.table = F)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

# make anno_df
anno_df =human_hypo_combined@misc$edgelist %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})

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
### nucseq QC and scvi Figure
##########

## see plotx exports from pipeline and scIntegration for this !

complete_evaluation_results_human = data.table::fread("data/complete_evaluation_results.txt",data.table = F)

##########
### Plot ASW by hyperparameters
##########

## prepare df for boxplots
plot_df = complete_evaluation_results_human[complete_evaluation_results_human$method=="scVI",]
plot_df$ndim= factor(plot_df$ndim,levels = c("50","80"))
plot_df$max_epochs= factor(plot_df$max_epochs,levels = as.character(sort(as.numeric(unique(plot_df$max_epochs)))))
plot_df$dropout_rate= factor(plot_df$dropout_rate,levels = as.character(sort(as.numeric(unique(plot_df$dropout_rate)))))
plot_df$features_ngenes= factor(plot_df$features_ngenes,levels = as.character(sort(as.numeric(unique(plot_df$features_ngenes)))))

okabe_ito_colors = 4 # 3 would work but 4 gives nicer ones
dodge_width = 0.75
text_size = 20

## plot ndim
p_asw_ndim = ggplot(plot_df, aes(x=ndim, y=asw)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers)) +
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

## plot epochs
p_asw_epochs = ggplot(plot_df, aes(x=max_epochs, y=asw)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers)) +
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

## plot features_ngenes
p_asw_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=asw)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers)) +
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

asw_combined_plot = cowplot::plot_grid(p_asw_ndim,p_asw_epochs,p_asw_nfeatures,ncol = 3)
asw_combined_plot
# save plot
# ggsave(filename = paste0(figure_dir,"evaluation_asw_by_hyperparam.pdf"),
#        plot = asw_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

##########
### Plot mean_knn_purity by hyperparameters
##########

## plot ndim
p_mean_knn_purity_ndim = ggplot(plot_df, aes(x=ndim, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

## plot epochs
p_mean_knn_purity_epochs = ggplot(plot_df, aes(x=max_epochs, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

## plot features_ngenes
p_mean_knn_purity_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

mean_knn_purity_combined_plot = cowplot::plot_grid(p_mean_knn_purity_ndim,p_mean_knn_purity_epochs,p_mean_knn_purity_nfeatures,ncol = 3)
mean_knn_purity_combined_plot
# save plot
# ggsave(filename = paste0(figure_dir,"evaluation_purity_by_hyperparam.pdf"),
#        plot = mean_knn_purity_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

##########
### Plot mixingknn by hyperparameters
##########

## plot ndim
p_mixingknn_ndim = ggplot(plot_df, aes(x=ndim, y=mixingknn)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

## plot epochs
p_mixingknn_epochs = ggplot(plot_df, aes(x=max_epochs, y=mixingknn)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

# ## plot dropout_rate
# p_mixingknn_dropout = ggplot(plot_df, aes(x=dropout_rate, y=mixingknn)) +
#   geom_boxplot(aes(fill=as.factor(nlayers))) +
#   geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))

## plot features_ngenes
p_mixingknn_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mixingknn)) +
  geom_boxplot(aes(fill=as.factor(nlayers))) +
  geom_point(position=position_dodge(width=dodge_width),aes(group=nlayers))+
  scale_fill_manual(values=getOkabeItoPalette(okabe_ito_colors))+
  theme_minimal()+theme(
    text = element_text(size=text_size),
    panel.background = element_rect(color=NA),
    panel.border = element_blank(),axis.line=element_line(color="grey40",linewidth = 0.5))

mixingknn_combined_plot = cowplot::plot_grid(p_mixingknn_ndim,p_mixingknn_epochs,p_mixingknn_nfeatures,ncol = 3)
mixingknn_combined_plot
# save plot
#ggsave(filename = paste0(evaluation_plot_path,"evaluation_mixing_by_hyperparam.pdf"),
#       plot = mixingknn_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

##########
### Re arrange plots by hyperparameter
##########

# features
features_combined_plot = cowplot::plot_grid(p_mixingknn_nfeatures,p_mean_knn_purity_nfeatures,p_asw_nfeatures,ncol = 3)
features_combined_plot
# save plot
ggsave(filename = paste0(figure_dir,"EXT_01_evaluation_metrics_by_features.pdf"),
       plot = features_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

# epochs
epochs_combined_plot = cowplot::plot_grid(p_mixingknn_epochs,p_mean_knn_purity_epochs,p_asw_epochs,ncol = 3)
epochs_combined_plot
# save plot
ggsave(filename = paste0(figure_dir,"EXT_01_evaluation_metrics_by_epochs.pdf"),
       plot = epochs_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

# features
ndims_combined_plot = cowplot::plot_grid(p_mixingknn_ndim,p_mean_knn_purity_ndim,p_asw_ndim,ncol = 3)
ndims_combined_plot
# save plot
ggsave(filename = paste0(figure_dir,"EXT_01_evaluation_metrics_by_ndims.pdf"),
       plot = ndims_combined_plot, "pdf",dpi=400,width=650,height = 200,units="mm")

##########
### Plot by dataset
##########

dataset_umap = DimPlot(human_hypo_combined,group.by = "Dataset",reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols =c("#D55E00","#56B4E9"),label=F,label.size = 2.5,shuffle = TRUE)+
  NoAxes()+ggtitle(NULL)#+NoLegend()
dataset_umap

## save
ggsave(filename = paste0(figure_dir,"EXT_02_umap_datasets.pdf"),
       plot = dataset_umap, "pdf",dpi=400,width=320,height = 300,units="mm")

## split with custom cols
col_split = "Dataset"
# cols
col_values =c("Siletti" = "#D55E00","Tadross" = "#56B4E9")
plist = list()
for(col_level in unique(human_hypo_combined@meta.data[,col_split])){
  human_hypo_combined@meta.data[,"plotcol"] = human_hypo_combined@meta.data[,col_split]
  human_hypo_combined@meta.data[human_hypo_combined@meta.data[,"plotcol"] != col_level,"plotcol"] = NA
  level_umap = DimPlot(human_hypo_combined,group.by = "plotcol",reduction = "umap",cols = col_values[col_level],pt.size = 1.5,raster.dpi = c(2048,2048),na.value = "grey80",label=F,label.size = 2.5,shuffle = F,order = TRUE)+
    NoAxes()+ggtitle(col_level)+NoLegend()
  temp_data = level_umap$data
  temp_data[,"plotcol"] = factor(temp_data[,"plotcol"],levels = c(col_level,NA))
  temp_data = temp_data %>% dplyr::arrange(!is.na(plotcol), plotcol)
  rownames(temp_data) = temp_data$Cell_ID
  level_umap$data = temp_data
  plist[[col_level]] = level_umap
}
split_umap = cowplot::plot_grid(plotlist = plist,nrow=1)
split_umap
## save
ggsave(filename = paste0(figure_dir,"EXT_02_umap_datasets_split.pdf"),
       plot = split_umap, "pdf",dpi=400,width=620,height = 300,units="mm")

##########
### Plot by Sex
##########

sex_umap = DimPlot(human_hypo_combined,group.by = "sex",reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols =c("darkred","darkblue"),label=F,label.size = 2.5,shuffle = TRUE)+
  NoAxes()+ggtitle(NULL)#+NoLegend()
sex_umap

## save
ggsave(filename = paste0(figure_dir,"EXT_02_umap_sex.pdf"),
       plot = sex_umap, "pdf",dpi=400,width=320,height = 300,units="mm")

## split with custom cols
col_split = "sex"
# cols
col_values =c("Male" = "darkblue","Female" = "darkred")
plist = list()
for(col_level in unique(human_hypo_combined@meta.data[,col_split])){
  human_hypo_combined@meta.data[,"plotcol"] = human_hypo_combined@meta.data[,col_split]
  human_hypo_combined@meta.data[human_hypo_combined@meta.data[,"plotcol"] != col_level,"plotcol"] = NA
  level_umap = DimPlot(human_hypo_combined,group.by = "plotcol",reduction = "umap",cols = col_values[col_level],pt.size = 1.5,raster.dpi = c(2048,2048),na.value = "grey80",label=F,label.size = 2.5,shuffle = F,order = TRUE)+
    NoAxes()+ggtitle(col_level)+NoLegend()
  temp_data = level_umap$data
  temp_data[,"plotcol"] = factor(temp_data[,"plotcol"],levels = c(col_level,NA))
  temp_data = temp_data %>% dplyr::arrange(!is.na(plotcol), plotcol)
  rownames(temp_data) = temp_data$Cell_ID
  level_umap$data = temp_data
  plist[[col_level]] = level_umap
}
split_umap = cowplot::plot_grid(plotlist = plist,nrow=1)
split_umap
## save
ggsave(filename = paste0(figure_dir,"EXT_02_umap_sex_split.pdf"),
       plot = split_umap, "pdf",dpi=400,width=620,height = 300,units="mm")


##########
### plot C2 overview and heatmap
##########

###
target_level="C2"
target_named="C2_named"
remove_prefix = FALSE


# make version without prefix 
if(remove_prefix){
  target_short = paste0(target_level,"_noPrefix")
  human_hypo_combined@meta.data[,target_short] = stringr::str_remove(human_hypo_combined@meta.data[,target_named],pattern = paste0(target_level,"-[0-9]+"))
  target_named = target_short
}

### umap plot

c2_umap = DimPlot(human_hypo_combined,group.by = target_named,reduction = "umap",pt.size = 1.5,raster.dpi = c(2048,2048),cols = getOkabeItoPalette(53),label=TRUE,label.size = 4,repel = TRUE)+
  NoAxes()+ggtitle(NULL)+NoLegend()

c2_umap

# select gene 

## make annotation data frame
anno_df_dotplot = anno_df


cluster_to_include = anno_df_dotplot$cluster_id[anno_df_dotplot$clusterlevel==target_level]
#neuron_node_ids = find_children(nodes = "C0-3",edges = combined_edgelist_mrtree)
#cluster_to_include = cluster_to_include[cluster_to_include %in% neuron_node_ids]

anno_df_dotplot = anno_df[anno_df$cluster_id %in% cluster_to_include,]
anno_df_dotplot$gene_dotplot = anno_df_dotplot$first_cluster_name
# use names:
anno_df_dotplot$gene_dotplot[! anno_df_dotplot$gene_dotplot %in% rownames(human_hypo_combined@assays$RNA@counts)] = NA

#use just the top feature ?
# features_dotplot = comparisons_all_updated %>% dplyr::filter(specificity > 2 & cluster %in% cluster_to_include) %>%
#   dplyr::filter(!grepl("LINC|AC|AL|\\-AS",gene)) %>%
#   dplyr::arrange(desc(specificity)) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 1,wt = specificity)
features_dotplot = comparisons_all_updated %>% dplyr::filter(specificity > 1 & cluster %in% cluster_to_include & pct.2 < 0.3 & pct.1 > 0.3) %>%
  dplyr::filter(!grepl("LINC|AC|AL|\\-AS",gene)) %>%
  dplyr::arrange(desc(specificity)) %>% 
  dplyr::group_by(cluster) %>% dplyr::top_n(n = 10,wt = specificity) %>% # take top N
  dplyr::mutate(rank_within = row_number()) %>% dplyr::arrange(rank_within) %>%
  dplyr::ungroup() %>% dplyr::distinct(gene,.keep_all = TRUE) %>% # then only use unique genes
  dplyr::group_by(cluster) %>% dplyr::top_n(n = 1,wt = specificity)%>%
  dplyr::arrange(desc(specificity))  # then reduce to top 1 (thereby taking 2nd or 3rd place gene if first is already included via another cluster

anno_df_dotplot = dplyr::left_join(anno_df_dotplot,features_dotplot %>% dplyr::select(gene,cluster),by=c("cluster_id"="cluster"))
anno_df_dotplot$gene_dotplot[is.na(anno_df_dotplot$gene_dotplot)] = anno_df_dotplot$gene[is.na(anno_df_dotplot$gene_dotplot)]

anno_df_dotplot = anno_df_dotplot %>% dplyr::distinct(cluster_name,.keep_all = TRUE)

# reorder anno_df_dotplot
unordered = unique(human_hypo_combined@meta.data[,target_level])
unordered = as.numeric(stringr::str_remove(stringr::str_extract(unordered,"-[0-9]+"),"-")) # order by number --> then it also is ordered by tree
names(unordered) = unique(human_hypo_combined@meta.data[,target_level])
rownames(anno_df_dotplot) = anno_df_dotplot$cluster_id
anno_df_dotplot = anno_df_dotplot[match(stringr::str_extract(names(sort(unordered)),"C2-[0-9]+"),rownames(anno_df_dotplot)),] 
anno_df_dotplot = anno_df_dotplot[!is.na(anno_df_dotplot$cluster_id),]
if(grepl("noPrefix",target_named)){
  anno_df_dotplot$factor_levels = stringr::str_remove( anno_df_dotplot$cluster_name,pattern = paste0(target_level,"-[0-9]+"))
}else{
  anno_df_dotplot$factor_levels = anno_df_dotplot$cluster_name
}

# need to reorder factor level in seurat
target_col_named =target_named
human_hypo_combined@meta.data[,target_col_named] = factor(as.character(human_hypo_combined@meta.data[,target_col_named]),levels = anno_df_dotplot$factor_levels)

# plot
cols_for_feature_plot=c("lightgrey", "blue")
Idents(human_hypo_combined) = target_named
dotplot_celltypes = Seurat::DotPlot(object = human_hypo_combined,
                                    features = unique(anno_df_dotplot$gene_dotplot[!is.na(anno_df_dotplot$gene_dotplot)]),
                                    #idents= as.character(cluster_to_include),
                                    scale = FALSE,
                                    cluster.idents = F)
# prettify
dotplot_celltypes2 = dotplot_celltypes + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  #  scale_size_continuous(limits = c(0,100)) +
  scale_color_gradient(low = cols_for_feature_plot[1],high = cols_for_feature_plot[2],limits = c(0,3), oob = scales::squish)
dotplot_celltypes2

## save
ggsave(filename = paste0(figure_dir,"EXT_02_umap_",target_named,"_celltypes.pdf"),
       plot = c2_umap, "pdf",dpi=400,width=300,height = 300,units="mm")

ggsave(filename = paste0(figure_dir,"EXT_02_dotplot_",target_named,"_celltypes.pdf"),
       plot = dotplot_celltypes2, "pdf",dpi=400,width=500,height = 400,units="mm")


##########
### Plot top genes as umaps
##########

## make feature plot with multiple genes:
rasterize_px = 1024
seurat_pt_size = 1
gene_set_to_plot = c("MEIS2","LHX6","FEZF1","TBX3","SIX3","OTP","SIM1","FOXB1")
p <- FeaturePlot(human_hypo_combined,features = gene_set_to_plot,keep.scale = "all", combine = FALSE,raster = TRUE,order=TRUE,raster.dpi = c(rasterize_px,rasterize_px),pt.size = seurat_pt_size,reduction = "umap")
for(i in 1:length(p)) {
  legend_plot = p[[i]]
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
  #p[[i]] <- scUtils::rasterize_ggplot(p[[i]],pixel_raster = 2048,pointsize = 1.8)
}
combined_feature_plots = cowplot::plot_grid(plotlist = p,ncol = 2)
combined_feature_plots

ggsave(filename = paste0(figure_dir,"EXT_03_umaps_tf_features.pdf"),
       plot = combined_feature_plots, "pdf",dpi=400,width=200,height = 400,units="mm")

# TODO: need to fix namespace problem and export the legend again
# Convert to a ggplot and print
leg <- ggpubr::get_legend(legend_plot)
leg <- ggpubr::as_ggplot(leg)
leg
ggsave(filename = paste0(figure_dir,"EXT_03_umaps_tf_features_legend.pdf"),
       plot = leg, "pdf",dpi=100,width=20,height = 40,units="mm")




##########
### Figure Extended GLP1R
##########

source("paper_figures/celltype_gene_overview_function.R")

other_genes_glp1r = c("POMC","SST","CALCR","TBX19","LEPR","TRH","AVP","SIM1","GIPR","SLC32A1")

incretin_overview_data = celltype_overview_data(
  main_gene = c("GLP1R"), # ,"GIPR"
  other_genes =other_genes_glp1r,
  human_hypo_combined = human_hypo_combined,
  hypoMap = hypoMap,region_mapping = region_mapping_forOverview,
  human_marker_genes = human_marker_genes,
  q_prob = 0.95 ,
  ignore_mouse_expression=TRUE,
  min_expr_value=0.05
)

#cutoff = 0.1

# incretin_overview_data_filtered = incretin_overview_data$gene_pct_anno[incretin_overview_data$gene_pct_anno$GLP1R > cutoff , ]
# 
# incretin_overview_plot = celltype_overview_plot(
#   gene_pct_anno_filter = incretin_overview_data_filtered,#incretin_overview_data$gene_pct_anno_filter,
#   combined_edgelist_mrtree= combined_edgelist_mrtree,
#   main_gene = c("GLP1R"),#,"GIPR"),
#   other_genes = c("POMC","SST","CALCR","TBX19","LEPR","TRH","AVP","SIM1","LMX1A","GIPR","SLC32A1"),
#   anno_df= anno_df,
#   highlight_leaves = c()#c("C6-278","C6-209","C6-261","C6-198","C6-413", "C6-414","C6-396"),
# )
# 
# incretin_overview_plot

incretin_overview_plot = celltype_overview_plot(
  gene_pct_anno_filter = incretin_overview_data$gene_pct_anno_filter,
  combined_edgelist_mrtree= combined_edgelist_mrtree,
  main_gene = c("GLP1R"),#,"GIPR"),
  other_genes = other_genes_glp1r,
  anno_df= anno_df,
  highlight_leaves = c()#c("C6-278","C6-209","C6-261","C6-198","C6-413", "C6-414","C6-396"),
)

incretin_overview_plot

# save plot
width=480
height = 300
ggsave(filename = paste0(figure_dir,"FIG_05_c_incretin_overview_plot.pdf"),plot = incretin_overview_plot, "pdf",dpi=400,width=width,height = height,units="mm")


## save table
data.table::fwrite(incretin_overview_data$gene_pct_anno,paste0(figure_dir,"incretin_overview_table_full.txt"),sep="\t")
data.table::fwrite(incretin_overview_data$gene_pct_anno_filter,paste0(figure_dir,"incretin_overview_table_filtered.txt"),sep="\t")

other_genes_both = c("POMC","SST","CALCR","TBX19","LEPR","TRH","AVP","SIM1","LMX1A","SLC32A1")

## version with both
incretin_overview_data2 = celltype_overview_data(
  main_gene = c("GLP1R","GIPR"), # ,
  other_genes = other_genes_both,
  human_hypo_combined = human_hypo_combined,
  hypoMap = hypoMap,
  region_mapping = region_mapping_forOverview,
  human_marker_genes = human_marker_genes,
  q_prob = 0.95 ,
  ignore_mouse_expression=TRUE,
  min_expr_value=0.05
)

## save table
data.table::fwrite(incretin_overview_data$gene_pct_anno,paste0(figure_dir,"incretin_overview_table_both_full.txt"),sep="\t")
data.table::fwrite(incretin_overview_data2$gene_pct_anno_filter,paste0(figure_dir,"incretin_overview_table_both_filtered.txt"),sep="\t")

##########
### Figure EXT Gipr only
##########

source("merge_human_mouse_neurons/celltype_gene_overview_function.R")

other_genes_gipr = c("PITX2","CCK","VIP","LEPR","CALCR","AVP","SIM2","SIM1","LMX1A","SLC32A1")

gipr_overview_data = celltype_overview_data(
  main_gene = c("GIPR"), # ,"GIPR"
  other_genes =other_genes_gipr,
  human_hypo_combined = human_hypo_combined,
  hypoMap = hypoMap,region_mapping = region_mapping_forOverview,
  human_marker_genes = human_marker_genes,
  q_prob = 0.95 ,
  ignore_mouse_expression=TRUE,
  min_expr_value=0.05
)

gipr_overview_plot = celltype_overview_plot(
  gene_pct_anno_filter = gipr_overview_data$gene_pct_anno_filter,
  combined_edgelist_mrtree= combined_edgelist_mrtree,
  main_gene = c("GIPR"),#,"GIPR"),
  other_genes = other_genes_gipr,
  anno_df= anno_df,
  highlight_leaves = c()#c("C6-278","C6-209","C6-261","C6-198","C6-413", "C6-414","C6-396"),
)

gipr_overview_plot

# save plot
width=500
height = 300
ggsave(filename = paste0(figure_dir,"EXT_08_a_incretin_overview_plot.pdf"),plot = gipr_overview_plot, "pdf",dpi=400,width=width,height = height,units="mm")


##########
### Figure Extended MC4R/MC4R
##########


mcr_overview_data = celltype_overview_data(
  main_gene = c("MC4R","MC3R"),
  other_genes = c("CHAT","NTS","GAL","KISS1","TAC2","TAC3","LEPR","CRABP1","SIM1","TRH","SLC32A1"),
  human_hypo_combined = human_hypo_combined,
  hypoMap = hypoMap,
  region_mapping = region_mapping_forOverview,
  human_marker_genes = human_marker_genes,
  q_prob = 0.95 ,
  ignore_mouse_expression=TRUE,
  min_expr_value=0.051
)

# mcr_overview_data_filtered = mcr_overview_data$gene_pct_anno[mcr_overview_data$gene_pct_anno$MC4R > cutoff | 
#                                                                mcr_overview_data$gene_pct_anno$MC3R > cutoff, ]
# 
# mcr_overview_plot = celltype_overview_plot(
#   gene_pct_anno_filter = mcr_overview_data_filtered,#incretin_overview_data$gene_pct_anno_filter,
#   combined_edgelist_mrtree= combined_edgelist_mrtree,
#   main_gene = c("MC4R","MC3R"),
#   other_genes = c("CHAT","NTS","GAL","KISS1","LEPR","TAC2","TAC3","CRABP1","SIM1","TRH","SLC32A1"),
#   anno_df= anno_df,
#   highlight_leaves = c()#c("C6-278","C6-209","C6-261","C6-198","C6-413", "C6-414","C6-396"),
# )

mcr_overview_plot = celltype_overview_plot(
  gene_pct_anno_filter = mcr_overview_data$gene_pct_anno_filter,
  combined_edgelist_mrtree= combined_edgelist_mrtree,
  main_gene = c("MC4R","MC3R"),
  other_genes = c("CHAT","NTS","GAL","KISS1","TAC2","TAC3","LEPR","CRABP1","SIM1","TRH","SLC32A1"),
  anno_df= anno_df,
  highlight_leaves = c()#c("C6-278","C6-209","C6-261","C6-198","C6-413", "C6-414","C6-396"),
)

mcr_overview_plot

# save plot
width=480
height = 300
ggsave(filename = paste0(figure_dir,"FIG_04_c_melanocortin_overview_plot.pdf"),plot = mcr_overview_plot, "pdf",dpi=400,width=width,height = height,units="mm")

## save table
data.table::fwrite(mcr_overview_data$gene_pct_anno,paste0(figure_dir,"melanocortin_overview_table_full.txt"),sep="\t")
data.table::fwrite(mcr_overview_data$gene_pct_anno_filter,paste0(figure_dir,"melanocortin_overview_table_filtered.txt"),sep="\t")

##########
### Density plot
##########

threshold = 0.1

gene = "GLP1R"
glp1r_density = ggplot(incretin_overview_data$gene_pct_anno ,aes_string(x=gene))+
  geom_density(size=1, color = "lightblue", fill="lightblue")+
  theme_light()+
  theme(text=element_text(size=25))+
  geom_vline(xintercept = 0.1,color="grey40", linetype="dashed", size=0.8)+
  xlab(paste0(gene," Average Expression"))
gene = "GIPR"
gipr_density = ggplot(incretin_overview_data$gene_pct_anno ,aes_string(x=gene))+
  geom_density(size=1, color = "lightblue", fill="lightblue")+
  theme_light()+
  theme(text=element_text(size=25))+
  geom_vline(xintercept = 0.1,color="grey40", linetype="dashed", size=0.8)+
  xlab(paste0(gene," Average Expression"))
gene = "MC4R"
mc4r_density = ggplot(mcr_overview_data$gene_pct_anno ,aes_string(x=gene))+
  geom_density(size=1, color = "lightblue", fill="lightblue")+
  theme_light()+
  theme(text=element_text(size=25))+
  geom_vline(xintercept = 0.1,color="grey40", linetype="dashed", size=0.8)+
  xlab(paste0(gene," Average Expression"))
gene = "MC3R"
mc3r_density = ggplot(mcr_overview_data$gene_pct_anno ,aes_string(x=gene))+
  geom_density(size=1, color = "lightblue", fill="lightblue")+
  theme_light()+
  theme(text=element_text(size=25))+
  geom_vline(xintercept = 0.1,color="grey40", linetype="dashed", size=0.8)+
  xlab(paste0(gene," Average Expression"))

density_plots = cowplot::plot_grid(glp1r_density,gipr_density,mc4r_density,mc3r_density,nrow = 1)
density_plots

ggsave(filename = paste0(figure_dir,"EXT_07_d_pomc_phylotree_crossspecies.pdf"), plot = density_plots, "pdf",dpi=300,width=600,height = 150,units="mm")





