##########
### Load spatial data
##########

# load the spatial cluster to grouped simplified anno mapping
avg_abundance_C3_ungrouped = data.table::fread("data/region_anno_revision/average_abundance_per_ungrouped_cluster/all_clusters_avg_abundance_C3_regions_NN12_res_0.5_named_median_240516.txt",data.table = F)
avg_abundance_C3_grouped = data.table::fread("data/region_anno_revision/average_abundance_per_grouped_cluster/all_clusters_avg_abundance_C3_regions_NN12_res_0.5_grouped_median_240516.txt",data.table = F)

# make a grouped to ungrouped df
region_grouping = data.frame(region_name_ungrouped = avg_abundance_C3_ungrouped[,1])
region_grouping$region_name_grouped = stringr::str_remove(region_grouping$region_name_ungrouped,pattern = " [0-9]+")

# load the color scheme for the grouped simplified anno
clusters_grouped_colors = data.table::fread("data/region_anno_revision/spatial_clusters_colours.txt",data.table = F)
region_color_mapping = as.character(clusters_grouped_colors$col)
names(region_color_mapping) = as.character(clusters_grouped_colors$V1)

# add to above also
region_grouping = dplyr::left_join(region_grouping,clusters_grouped_colors,by=c("region_name_grouped"="V1"))
region_color_mapping_ungrp = as.character(region_grouping$col)
names(region_color_mapping_ungrp) = as.character(region_grouping$region_name_ungrouped)


# load the avg c4 abundances per spatial cluster
avg_abundances_c4 = data.table::fread("data/region_anno_revision/average_abundance_per_ungrouped_cluster/all_clusters_avg_abundance_C4_regions_NN12_res_0.5_named_mean_240516.txt",data.table = F)
avg_abundances_c3 = data.table::fread("data/region_anno_revision/average_abundance_per_ungrouped_cluster/all_clusters_avg_abundance_C3_regions_NN12_res_0.5_named_mean_240516.txt",data.table = F)
avg_abundances_c2 = data.table::fread("data/region_anno_revision/average_abundance_per_ungrouped_cluster/all_clusters_avg_abundance_C2_regions_NN12_res_0.5_named_mean_240516.txt",data.table = F)
# make long format
avg_abundances_c4_long = avg_abundances_c4 %>% tidyr::gather(key = "cluster",value="abundance",-V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))
avg_abundances_c3_long = avg_abundances_c3 %>% tidyr::gather(key = "cluster",value="abundance",-V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))
avg_abundances_c2_long = avg_abundances_c2 %>% tidyr::gather(key = "cluster",value="abundance",-V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))

all_abundances = list(
  C2 = avg_abundances_c2_long,
  C3 = avg_abundances_c3_long,
  C4 = avg_abundances_c4_long
)

human_hypo_combined@meta.data = human_hypo_combined@meta.data[,1:27] 

##########
### Check assigned region cluster by taking max
##########


all_abundances_max_stats =list()

for(i in 1:length(all_abundances)){
  target_level = names(all_abundances)[i]
  
  avg_abundances_long_stats = all_abundances[[i]] %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(median_ab = median(abundance),
                  max_ab = max(abundance), 
                  max_median_ratio = max_ab / median_ab, 
                  current_max_ratio = abundance / max_ab) %>%
    dplyr::slice_max(order_by = abundance,n = 1)
  
  all_abundances_max_stats[[target_level]] = avg_abundances_long_stats
  
  add_regions = avg_abundances_long_stats %>% 
    dplyr::left_join(region_grouping[,1:2],by=c("region"="region_name_ungrouped")) %>%
    dplyr::select(cluster,region_name_grouped) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(cluster = stringr::str_replace(cluster,"\\.","-"))
  colnames(add_regions) = c(target_level,paste0(target_level,"_region_max"))
  temp_meta = human_hypo_combined@meta.data %>% dplyr::left_join(add_regions) %>% as.data.frame()
  rownames(temp_meta) = temp_meta$Cell_ID
  human_hypo_combined@meta.data = temp_meta
}

# and make plots
colnames( human_hypo_combined@meta.data)
p2= DimPlot(human_hypo_combined,group.by = c("C2_region_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p3 = DimPlot(human_hypo_combined,group.by = c("C3_region_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p4 = DimPlot(human_hypo_combined,group.by = c("C4_region_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
cowplot::plot_grid(plotlist = list(p2,p3,p4),nrow = 1)

#
barplot_list = list()
for(i in 1:length(all_abundances_max_stats)){
  barplot_list[[ names(all_abundances_max_stats)[i]]] = ggplot(all_abundances_max_stats[[i]],aes(region,fill=region))+geom_bar()+
    ggtitle(names(all_abundances_max_stats)[i])+theme_light()+theme(text=element_text(size=15),axis.text.x = element_text(size=10,angle = 90))+NoLegend()+ 
    scale_fill_manual(values = region_color_mapping_ungrp)
}
cowplot::plot_grid(plotlist = barplot_list,nrow = 1)

##########
### Make ridgeplot
##########

ridgeplot_list = list()
for(i in 1:length(all_abundances)){
  ridgeplot_list[[ names(all_abundances)[i]]] = ggplot2::ggplot(data = all_abundances[[i]],mapping = aes(x=log(abundance),y=region,fill=region)) + 
    ggridges::geom_density_ridges(alpha=0.3)+
    ggtitle(names(all_abundances)[i])+theme_light()+theme(text=element_text(size=15))+NoLegend()+ 
    scale_fill_manual(values = region_color_mapping_ungrp)
}
cowplot::plot_grid(plotlist = ridgeplot_list,nrow = 1)

#+NoLegend()



##########
### Check assigned region cluster by taking max after adjusting
##########


all_abundances_adj_max_stats =list()

for(i in 1:length(all_abundances)){
  target_level = names(all_abundances)[i]
  
  # with adjust
  avg_abundances_long_stats = all_abundances[[i]] %>%
    dplyr::group_by(region) %>%
    dplyr::mutate(abundance_adj = abundance - median(abundance)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(median_ab = median(abundance_adj),max_ab = max(abundance_adj), max_median_ratio = max_ab / median_ab, current_max_ratio = abundance_adj / max_ab) %>%
    dplyr::slice_max(order_by = abundance_adj,n = 1)
  
  all_abundances_adj_max_stats[[target_level]] = avg_abundances_long_stats
  
  add_regions = avg_abundances_long_stats %>% 
    dplyr::left_join(region_grouping[,1:2],by=c("region"="region_name_ungrouped")) %>%
    dplyr::select(cluster,region_name_grouped) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(cluster = stringr::str_replace(cluster,"\\.","-"))
  colnames(add_regions) = c(target_level,paste0(target_level,"_region_adj_max"))
  temp_meta = human_hypo_combined@meta.data %>% dplyr::left_join(add_regions) %>% as.data.frame()
  rownames(temp_meta) = temp_meta$Cell_ID
  human_hypo_combined@meta.data = temp_meta
}

# and make plots
colnames( human_hypo_combined@meta.data)
p2= DimPlot(human_hypo_combined,group.by = c("C2_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p3 = DimPlot(human_hypo_combined,group.by = c("C3_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p4 = DimPlot(human_hypo_combined,group.by = c("C4_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
cowplot::plot_grid(plotlist = list(p2,p3,p4),nrow = 1)

#
barplot_list = list()
for(i in 1:length(all_abundances_adj_max_stats)){
  barplot_list[[ names(all_abundances_adj_max_stats)[i]]] = ggplot(all_abundances_adj_max_stats[[i]],aes(region,fill=region))+geom_bar()+
    ggtitle(names(all_abundances_adj_max_stats)[i])+theme_light()+theme(text=element_text(size=15),axis.text.x = element_text(size=10,angle = 90))+NoLegend()+ 
    scale_fill_manual(values = region_color_mapping_ungrp)
}
cowplot::plot_grid(plotlist = barplot_list,nrow = 1)


##########
### Check assigned region cluster by taking max after adjusting + filter for MAD
##########


all_abundances_adj_max_mad_stats =list()
human_hypo_combined@meta.data= human_hypo_combined@meta.data[,1:27]
for(i in 1:length(all_abundances)){
  target_level = names(all_abundances)[i]
  print(target_level)
  # with adjust
  avg_abundances_long_stats = all_abundances[[i]] %>%
    dplyr::group_by(region) %>%
    dplyr::mutate(abundance_adj = abundance - median(abundance)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(median_ab = median(abundance_adj),max_ab = max(abundance_adj), max_median_ratio = max_ab / median_ab, current_max_ratio = abundance_adj / max_ab, abundance_mad = mad(abundance_adj), mad_x = abundance_adj / abundance_mad) %>%
    dplyr::slice_max(order_by = abundance_adj,n = 1)# %>%
   # dplyr::filter( mad_x >= 5)
  
  all_abundances_adj_max_mad_stats[[target_level]] = avg_abundances_long_stats
  
  add_regions = avg_abundances_long_stats %>% 
    dplyr::left_join(region_grouping[,1:2],by=c("region"="region_name_ungrouped")) %>%
    dplyr::select(cluster,region_name_grouped) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(cluster = stringr::str_replace(cluster,"\\.","-"))
  colnames(add_regions) = c(target_level,paste0(target_level,"_region_adj_max"))
  temp_meta = human_hypo_combined@meta.data %>% dplyr::left_join(add_regions) %>% as.data.frame()
  rownames(temp_meta) = temp_meta$Cell_ID
  human_hypo_combined@meta.data = temp_meta
}

# and make plots
colnames( human_hypo_combined@meta.data)
sapply(all_abundances_adj_max_mad_stats,nrow)
p2= DimPlot(human_hypo_combined,group.by = c("C2_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p3 = DimPlot(human_hypo_combined,group.by = c("C3_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
p4 = DimPlot(human_hypo_combined,group.by = c("C4_region_adj_max"),reduction = "umap_scvi_hypo",label = F,shuffle=TRUE,raster.dpi = c(1536,1536),pt.size = 1.3) + scale_color_manual(values = region_color_mapping)#+NoLegend()
cowplot::plot_grid(plotlist = list(p2,p3,p4),nrow = 1)

#
barplot_list = list()
for(i in 1:length(all_abundances_adj_max_stats)){
  barplot_list[[ names(all_abundances_adj_max_stats)[i]]] = ggplot(all_abundances_adj_max_stats[[i]],aes(region,fill=region))+geom_bar()+
    ggtitle(names(all_abundances_adj_max_stats)[i])+theme_light()+theme(text=element_text(size=15),axis.text.x = element_text(size=10,angle = 90))+NoLegend()+ 
    scale_fill_manual(values = region_color_mapping_ungrp)
}
cowplot::plot_grid(plotlist = barplot_list,nrow = 1)

##########
### make table for print
##########

all_t= data.frame()
for(i in 1:length(all_abundances_max_stats)){
  t=all_abundances_max_stats[[i]]
  t2 = t[,c("cluster","region")]
  t2$level = names(all_abundances_max_stats)[i]
  t2$type="max"
  all_t = bind_rows(all_t,t2)
}
for(i in 1:length(all_abundances_adj_max_stats)){
  t=all_abundances_adj_max_stats[[i]]
  t2 = t[,c("cluster","region")]
  t2$level = names(all_abundances_adj_max_stats)[i]
  t2$type="max_adj"
  all_t = bind_rows(all_t,t2)
}
all_t = all_t %>% tidyr::spread(key = "type",value = "region")


for(i in 1:length(all_abundances_adj_max_mad_stats)){
  t=all_abundances_adj_max_mad_stats[[i]]
  t2 = t[,c("cluster","region")]
  t2$level = names(all_abundances_adj_max_mad_stats)[i]
  t2$type="max_adj"
  all_t = bind_rows(all_t,t2)
}
all_t = all_t %>% tidyr::spread(key = "type",value = "region")



data.table::fwrite(all_t,"data/region_anno_revision/preliminary_assignments_v1.tsv",sep="\t")

##########
### Georgie
##########

regions_zscore_based = readxl::read_xlsx("data/region_anno_revision/all_clusters_avg_abundance_C4_regions_NN12_res_0.5_named_median_240516_zscoretest.xlsx")
colnames(regions_zscore_based)[1] ="region"
regions_zscore_based_long = regions_zscore_based %>% tidyr::gather(key = "cluster",value = "value"-region)


##########
### Georgie
##########

all_abundances_adj_max_mad_stats$C3
all_abundances_adj_max_mad_stats$C3$cluster_id = gsub("\\.","-",all_abundances_adj_max_mad_stats$C3$cluster)
threshold=5
all_c3_clusters = data.frame(cluster_id = combined_edgelist_mrtree$to[grepl("C3-",combined_edgelist_mrtree$to)]) %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::select(to,cluster_name = to_named),by=c("cluster_id"="to")) %>%
  dplyr::left_join(all_abundances_adj_max_mad_stats$C3[,c("cluster_id","region","mad_x","abundance")],by="cluster_id") %>%
  dplyr::left_join(region_grouping %>% dplyr::rename(region_final = region_name_grouped),by=c("region"="region_name_ungrouped")) %>%
    dplyr::mutate(region_final = case_when(
      mad_x < threshold ~ NA,
      .default = region_final
    )) %>%
  plyr::mutate(col = case_when(
    mad_x < threshold ~ NA,
    .default = col
  ))
data.table::fwrite(all_c3_clusters,"data/region_anno_revision/mad_assignments_C3.tsv",sep="\t")

all_c4_clusters = data.frame(cluster_id = combined_edgelist_mrtree$to[grepl("C4",combined_edgelist_mrtree$to)])  %>%
  dplyr::left_join(combined_edgelist_mrtree %>% dplyr::select(to,cluster_name = to_named, parent_id = from),by=c("cluster_id"="to")) %>%
  dplyr::left_join(all_c3_clusters %>% dplyr::select(-cluster_name),by=c("parent_id"="cluster_id"))


data.table::fwrite(all_c4_clusters,"data/region_anno_revision/mad_assignments_C4.tsv",sep="\t")


