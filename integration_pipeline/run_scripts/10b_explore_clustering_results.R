param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/mrtree_construction_params_310a3d4d3b60c8b2f05e125f65365740.json"
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### HypoNeuron
##########

subname = "HypoNeuron"
parameter_list$harmonization_folder_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/")
parameter_list$clustering_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/consensus/")
parameter_list$new_name_suffix = paste0("human_hypo_",subname)

### load seurat
if(!exists("neuron_seurat_object")){
  neuron_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
}
### load clusters
#subname = gsub("_clusters","",parameter_list$clustering_key_name)

#  merge from multiple files (consensus clusterings from multiple jobs)
message("Reading clusters from folder and merge: ",parameter_list$clustering_folder)
all_cluster_files = list.files(parameter_list$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(parameter_list$harmonization_folder_path,"consensus/")
all_clusterings = lapply(paste0(parameter_list$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)

### Stats
factor = 2
clusters_stats = as.data.frame(cluster_matrix_for_mrtree) %>% 
  dplyr::mutate(Cell_ID = rownames(cluster_matrix_for_mrtree)) %>%
  tidyr::gather(key="cluster_level",value="cluster_id",-Cell_ID) %>%
  dplyr::group_by(cluster_level,cluster_id) %>%
  dplyr::add_count(name="n_cells") %>%
  dplyr::distinct(cluster_level,cluster_id,n_cells) %>%
  dplyr::group_by(cluster_level) %>%
  dplyr::add_count(name="n_clusters") %>%
  dplyr::distinct(cluster_level,.keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(double_n_clusters = dplyr::lag(n_clusters,n=1)*factor) %>%
  dplyr::mutate(diff_to_lag = abs(n_clusters-double_n_clusters))
clusters_stats$res =as.numeric(stringr::str_remove_all(stringr::str_extract(clusters_stats$cluster_level,"\\_[0-9\\.]+\\_"),"_"))
clusters_stats = clusters_stats %>% dplyr::arrange(res)
clusters_stats$cluster_level = factor(clusters_stats$cluster_level,levels = clusters_stats$cluster_level)
clusters_stats$res = factor(clusters_stats$res,levels = clusters_stats$res)

### Plot

neuron_seurat_object@meta.data = neuron_seurat_object@meta.data[,1:56]
neuron_seurat_object@meta.data = cbind(neuron_seurat_object@meta.data,cluster_matrix_for_mrtree)

cols_to_plot = c("res_0.05_repeat_100","res_0.5_repeat_100","res_4_repeat_100","res_46_repeat_100")
plist=list()
for(col in cols_to_plot){
  plist[[col]] = DimPlot(neuron_seurat_object,group.by = col,reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()
}
cowplot::plot_grid(plotlist = plist,ncol=3)


##########
###  ARI heatmap for neurons
##########

target_res = "res_46_repeat_100"
clustering_file_name = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_leiden_",target_res,".tsv")

# needs above section to get cluster_matrix_for_mrtree
consensus_clusters = cluster_matrix_for_mrtree[,target_res]
# loadd base clusterings
cluster_matrix_for_consensus = data.table::fread(clustering_file_name,data.table = FALSE,header = TRUE)

ari_input = cbind(consensus_orginal = consensus_clusters, cluster_matrix_for_consensus)
ari_mat = matrix(data = 0,nrow = ncol(ari_input),ncol = ncol(ari_input))
for(i in 1:ncol(ari_input)){
  for(j in 1:ncol(ari_input)){
    ari_mat[i,j] = aricode::ARI(c1 = ari_input[,i],c2 = ari_input[,j])
  }
  message("Average ARI of ",colnames(ari_input)[i]," is ", mean(ari_mat[i,]))
}
colnames(ari_mat) = colnames(ari_input)
rownames(ari_mat) = colnames(ari_input)

gplots::heatmap.2(ari_mat[1:ncol(ari_mat),1:ncol(ari_mat)],key=F,scale="none",trace="none")


#Perform hierarchical clustering
colnames(ari_mat)[colnames(ari_mat)=="consensus_orginal"] = "consensus"
rownames(ari_mat)[rownames(ari_mat)=="consensus_orginal"] = "consensus"
dist_mat <- dist(ari_mat)
hc <- hclust(dist_mat)

# Generate dendrograms for rows and columns
row_dendrogram <- as.dendrogram(hc)
# row_dendrogram[[2]] <- rev(row_dendrogram[[2]])
# colden_reordered <- reorder(row_dendrogram, c(10, 1, 1, 100, 300, 200))

leaves <- labels(row_dendrogram)
leave_weights = rep(-100, length(leaves))
names(leave_weights)= leaves
leave_weights[names(leave_weights)=="consensus"]= -10000
row_dendrogram = reorder(row_dendrogram,wts=rev(leave_weights), agglo.FUN = max)

# Plot heatmap with dendrograms
png("results/consensus_cluster_heatmap_res_4_repeat_100.png", width = 1400, height = 1400)
gplots::heatmap.2(
  ari_mat[1:ncol(ari_mat), 1:ncol(ari_mat)],
  key = FALSE,
  scale = "none",
  trace = "none",
  Rowv = row_dendrogram,
  Colv = row_dendrogram
)
dev.off()

##########
### Oligodendrocytes
##########

subname = "Oligodendrocytes"
parameter_list$new_name_suffix = paste0("human_hypo_",subname)
parameter_list$harmonization_folder_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/")
parameter_list$clustering_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/consensus/")

### load seurat
if(!exists("oligo_seurat_object")){
  oligo_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
}

### load clusters
#subname = gsub("_clusters","",parameter_list$clustering_key_name)

#  merge from multiple files (consensus clusterings from multiple jobs)
message("Reading clusters from folder and merge: ",parameter_list$clustering_folder)
all_cluster_files = list.files(parameter_list$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(parameter_list$harmonization_folder_path,"consensus/")
all_clusterings = lapply(paste0(parameter_list$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)

### Stats

clusters_stats = as.data.frame(cluster_matrix_for_mrtree) %>% 
  dplyr::mutate(Cell_ID = rownames(cluster_matrix_for_mrtree)) %>%
  tidyr::gather(key="cluster_level",value="cluster_id",-Cell_ID) %>%
  dplyr::group_by(cluster_level,cluster_id) %>%
  dplyr::add_count(name="n_cells") %>%
  dplyr::distinct(cluster_level,cluster_id,n_cells) %>%
  dplyr::group_by(cluster_level) %>%
  dplyr::add_count(name="n_clusters") %>%
  dplyr::distinct(cluster_level,.keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(double_n_clusters = dplyr::lag(n_clusters,n=1)*factor) %>%
  dplyr::mutate(diff_to_lag = abs(n_clusters-double_n_clusters))
clusters_stats$res =as.numeric(stringr::str_remove_all(stringr::str_extract(clusters_stats$cluster_level,"\\_[0-9\\.]+\\_"),"_"))
clusters_stats = clusters_stats %>% dplyr::arrange(res)
clusters_stats$cluster_level = factor(clusters_stats$cluster_level,levels = clusters_stats$cluster_level)
clusters_stats$res = factor(clusters_stats$res,levels = clusters_stats$res)

### Plot

oligo_seurat_object@meta.data = oligo_seurat_object@meta.data[,1:56]
oligo_seurat_object@meta.data = cbind(oligo_seurat_object@meta.data,cluster_matrix_for_mrtree)

cols_to_plot = c("res_0.25_repeat_100","res_0.5_repeat_100","res_0.75_repeat_100","res_1_repeat_100","res_1.5_repeat_100","res_2.5_repeat_100") # Oligodendrocytes
cols_to_plot = c("res_0.25_repeat_100","res_0.75_repeat_100","res_1.5_repeat_100") # Oligodendrocytes
plist=list()
for(col in cols_to_plot){
  plist[[col]] = DimPlot(oligo_seurat_object,group.by = col,reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()
}
cowplot::plot_grid(plotlist = plist,ncol=3)

##########
### AstroEpendymal
##########

subname = "AstroEpendymal"
parameter_list$new_name_suffix = paste0("human_hypo_",subname)
parameter_list$harmonization_folder_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/")
parameter_list$clustering_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/consensus/")
parameter_list$new_name_suffix = paste0("human_hypo_",subname)

### load seurat
if(!exists("astro_seurat_object")){
  astro_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
}

### load clusters
#subname = gsub("_clusters","",parameter_list$clustering_key_name)

#  merge from multiple files (consensus clusterings from multiple jobs)
message("Reading clusters from folder and merge: ",parameter_list$clustering_folder)
all_cluster_files = list.files(parameter_list$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(parameter_list$harmonization_folder_path,"consensus/")
all_clusterings = lapply(paste0(parameter_list$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)

### Stats

clusters_stats = as.data.frame(cluster_matrix_for_mrtree) %>% 
  dplyr::mutate(Cell_ID = rownames(cluster_matrix_for_mrtree)) %>%
  tidyr::gather(key="cluster_level",value="cluster_id",-Cell_ID) %>%
  dplyr::group_by(cluster_level,cluster_id) %>%
  dplyr::add_count(name="n_cells") %>%
  dplyr::distinct(cluster_level,cluster_id,n_cells) %>%
  dplyr::group_by(cluster_level) %>%
  dplyr::add_count(name="n_clusters") %>%
  dplyr::distinct(cluster_level,.keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(double_n_clusters = dplyr::lag(n_clusters,n=1)*factor) %>%
  dplyr::mutate(diff_to_lag = abs(n_clusters-double_n_clusters))
clusters_stats$res =as.numeric(stringr::str_remove_all(stringr::str_extract(clusters_stats$cluster_level,"\\_[0-9\\.]+\\_"),"_"))
clusters_stats = clusters_stats %>% dplyr::arrange(res)
clusters_stats$cluster_level = factor(clusters_stats$cluster_level,levels = clusters_stats$cluster_level)
clusters_stats$res = factor(clusters_stats$res,levels = clusters_stats$res)

### Plot

astro_seurat_object@meta.data = astro_seurat_object@meta.data[,1:56]
astro_seurat_object@meta.data = cbind(astro_seurat_object@meta.data,cluster_matrix_for_mrtree)

cols_to_plot = c("celltype_annotation","res_0.25_repeat_100","res_3_repeat_100") # AstroEpendymal
plist=list()
for(col in cols_to_plot){
  plist[[col]] = DimPlot(astro_seurat_object,group.by = col,reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()
}
cowplot::plot_grid(plotlist = plist,ncol=3)

#FeaturePlot(astro_seurat_object,reduction = paste0("umap_scvi_",subname),features = "GJB6",pt.size = 1.5,order = TRUE,raster = TRUE)+NoAxes()
#DimPlot(astro_seurat_object,group.by = "celltype_annotation",reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()

##########
### NonNeuron
##########

subname = "NonNeuron"
parameter_list$new_name_suffix = paste0("human_hypo_",subname)
parameter_list$harmonization_folder_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/")
parameter_list$clustering_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/",subname,"/consensus/")
parameter_list$new_name_suffix = paste0("human_hypo_",subname)

### load seurat
if(!exists("nonneuron_seurat_object")){
  nonneuron_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
}
#  merge from multiple files (consensus clusterings from multiple jobs)
message("Reading clusters from folder and merge: ",parameter_list$clustering_folder)
all_cluster_files = list.files(parameter_list$clustering_folder,pattern = ".tsv|.txt|.csv") # clustering_folder = paste0(parameter_list$harmonization_folder_path,"consensus/")
all_clusterings = lapply(paste0(parameter_list$clustering_folder,all_cluster_files),data.table::fread,data.table=FALSE,header=TRUE)
cluster_matrix_for_mrtree = do.call(cbind,all_clusterings)

### Stats

clusters_stats = as.data.frame(cluster_matrix_for_mrtree) %>% 
  dplyr::mutate(Cell_ID = rownames(cluster_matrix_for_mrtree)) %>%
  tidyr::gather(key="cluster_level",value="cluster_id",-Cell_ID) %>%
  dplyr::group_by(cluster_level,cluster_id) %>%
  dplyr::add_count(name="n_cells") %>%
  dplyr::distinct(cluster_level,cluster_id,n_cells) %>%
  dplyr::group_by(cluster_level) %>%
  dplyr::add_count(name="n_clusters") %>%
  dplyr::distinct(cluster_level,.keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(double_n_clusters = dplyr::lag(n_clusters,n=1)*factor) %>%
  dplyr::mutate(diff_to_lag = abs(n_clusters-double_n_clusters))
clusters_stats$res =as.numeric(stringr::str_remove_all(stringr::str_extract(clusters_stats$cluster_level,"\\_[0-9\\.]+\\_"),"_"))
clusters_stats = clusters_stats %>% dplyr::arrange(res)
clusters_stats$cluster_level = factor(clusters_stats$cluster_level,levels = clusters_stats$cluster_level)
clusters_stats$res = factor(clusters_stats$res,levels = clusters_stats$res)

### Plot

nonneuron_seurat_object@meta.data = nonneuron_seurat_object@meta.data[,1:56]
nonneuron_seurat_object@meta.data = cbind(nonneuron_seurat_object@meta.data,cluster_matrix_for_mrtree)

cols_to_plot = c("res_0.1_repeat_100","res_0.25_repeat_100","res_0.5_repeat_100","res_1_repeat_100","res_2_repeat_100","res_2.5_repeat_100") # NonNeuron
cols_to_plot = c("res_0.01_repeat_100","res_0.1_repeat_100","res_2_repeat_100")
plist=list()
for(col in cols_to_plot){
  plist[[col]] = DimPlot(nonneuron_seurat_object,group.by = col,reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()
}
cowplot::plot_grid(plotlist = plist,ncol=3)

#DimPlot(nonneuron_seurat_object,group.by = "Donor_ID",reduction = paste0("umap_scvi_",subname),label = TRUE,label.size = 1.5,raster = TRUE,pt.size = 1.5)+NoLegend()+NoAxes()


