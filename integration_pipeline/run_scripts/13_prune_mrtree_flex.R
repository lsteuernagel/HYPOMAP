##########
### Load parameters and packages
##########

message(Sys.time(),": Starting mrtree pruning .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1] # param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/mrtree_construction_params_310a3d4d3b60c8b2f05e125f65365740.json" ,,"/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/AstroEpendymalmrtree_pruning_3_params_896fa900f2de2faec4ee2cb6ea813b64.json"
param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/AstroEpendymalmrtree_pruning_5_params_d86674071e0576831712925a0a9fa4c3.json" #. astro
param_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/human_hypo_params_2/HypoNeuronmrtree_pruning_5_params_cf9b44c323b09505e298930d80d9ade8.json" # neuron

# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

#test:
# parameter_list$marker_suffix = "raw"
# parameter_list$start_node = "all"
# parameter_list$merge_sample_based=TRUE
# parameter_list$pct_threshold = 80
# parameter_list$min_cells = 50 # if below --> merge with neighbor
# parameter_list$min_specificity = 1.5 # min specificity for a sibling marker to count
# parameter_list$max_pvalue_prune = 0.001 # max pvalue for a sibling marker to count
# parameter_list$min_sibling_markers = 15 # how many sibling markers are required to not merge

additional_suffix = parameter_list$additional_suffix
if(is.null(additional_suffix)){additional_suffix = ""}

# load object
message(Sys.time(),": Load data .." )
harmonized_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,".rds"))
original_cols = ncol(harmonized_seurat_object@meta.data)

# load mrtree clustering
mrtree_result = readRDS(parameter_list$mrtree_result_file) #readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$clustering_key_name,"_mrtree_results",".rds"))

if(parameter_list$merge_marker_based){
  # load markers 
  markers_comparisons_all = data.table::fread(file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_all_",parameter_list$marker_suffix,"_",parameter_list$additional_suffix,".tsv"),data.table = F)
  markers_comparisons_siblings  =  data.table::fread(file = paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_",parameter_list$clustering_key_name,"_",parameter_list$start_node,"_markers_siblings_",parameter_list$marker_suffix,"_",parameter_list$additional_suffix,".tsv"),data.table = F)
}

##########
### prepare raw edgelists
##########

# get key objects
edgelist = mrtree_result$edgelist
# edgelist$from = stringr::str_split_i(edgelist$from,pattern = "/",i = -1)
# edgelist$to = stringr::str_split_i(edgelist$to,pattern = "/",i = -1)
labelmat = mrtree_result$labelmat

# add to seurat:
harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data[,1:original_cols],labelmat)

##########
### pre calculate sample contribution
##########

all_mrtree_cols = colnames(labelmat)
all_cluster_stats_list = list()
# make overall per cluster stats
for(mrtree_col in all_mrtree_cols){
  all_cluster_stats = harmonized_seurat_object@meta.data %>%
    dplyr::group_by(!!sym(mrtree_col)) %>% 
    dplyr::mutate(mean_nCount_RNA = mean(nCount_RNA)) %>%
    dplyr::mutate(mean_nFeature_RNA = mean(nFeature_RNA)) %>%
    dplyr::mutate(mean_percent_mt = mean(percent.mt)) %>%
    dplyr::add_count(name="cells_cluster") %>%
    dplyr::group_by(!!sym(mrtree_col),Dataset) %>% dplyr::add_count(name="cells_cluster_dataset") %>%
    dplyr::group_by(!!sym(mrtree_col),Donor_ID) %>% dplyr::add_count(name="cells_cluster_donorID") %>%
    dplyr::group_by(!!sym(mrtree_col),Sample_ID) %>% dplyr::add_count(name="cells_cluster_sampleID") %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!sym(mrtree_col),cells_cluster,Dataset,cells_cluster_dataset,Donor_ID,cells_cluster_donorID,Sample_ID,cells_cluster_sampleID,mean_nCount_RNA,mean_nFeature_RNA,mean_percent_mt) %>%
    dplyr::mutate(pct_cluster_dataset = round(cells_cluster_dataset / cells_cluster * 100,3),
                  pct_cluster_donorID = round(cells_cluster_donorID / cells_cluster * 100,3),
                  pct_cluster_sampleID = round(cells_cluster_sampleID / cells_cluster * 100,3))
  colnames(all_cluster_stats)[1] = "cluster"
  # get the max Donor
  donor_pct_max = all_cluster_stats %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::distinct(cluster,Donor_ID,pct_cluster_donorID) %>% 
    dplyr::slice_max(order_by = pct_cluster_donorID,n = 1,with_ties = FALSE) %>%
    dplyr::select(cluster,pct_donorMax = pct_cluster_donorID,donorMax = Donor_ID)
  # get the max sample
  sample_pct_max = all_cluster_stats %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::distinct(cluster,Sample_ID,pct_cluster_sampleID) %>% 
    dplyr::slice_max(order_by = pct_cluster_sampleID,n = 1,with_ties = FALSE) %>%
    dplyr::select(cluster,pct_sampleMax = pct_cluster_sampleID,sampleMax = Sample_ID)
  # make a continous average df
  avg_stats = all_cluster_stats %>% 
    dplyr::distinct(cluster,mean_nCount_RNA,mean_nFeature_RNA,mean_percent_mt)
  
  # make the final overview DF
  cluster_qc_overview = all_cluster_stats %>% 
    dplyr::distinct(cluster,cells_cluster) %>%
    dplyr::left_join(donor_pct_max,by="cluster") %>%
    dplyr::left_join(sample_pct_max,by="cluster") 
  
  all_cluster_stats_list[[mrtree_col]] = cluster_qc_overview
}
cluster_qc_overview = do.call(rbind,all_cluster_stats_list) %>% as.data.frame()

# hist(cluster_qc_overview$pct_donorMax[grepl("K616",cluster_qc_overview$cluster)],breaks = 20)
# 
# hist((cluster_qc_overview$cells_cluster[grepl("K616",cluster_qc_overview$cluster)]),breaks = 40)
##########
### helper function
##########

closest_cluster = function(object, target_cluster, other_clusters,cluster_col = "leiden_clusters_12_simplified", reduction="scvi", metric = "cosine"){
  if(length(other_clusters)>1){
    cells_cluster_list = split(rownames(object@meta.data),f = object@meta.data[,cluster_col])
    cells_cluster_list = cells_cluster_list[names(cells_cluster_list) %in% c(target_cluster, other_clusters)]
    embedding = as.matrix(object@reductions[[reduction]]@cell.embeddings)
    # calc avg scvi per cluster
    mean_embedding_cluster_list = lapply(cells_cluster_list, function(cells,expr_matrix){
      if(length(cells)>1){
        return(colMeans(expr_matrix[cells,]))
      }else{
        return(expr_matrix[cells,])
      }
    },expr_matrix = embedding)
    mean_embedding_cluster = do.call(cbind,mean_embedding_cluster_list)
    distances = as.matrix(dist(t(mean_embedding_cluster)))
    #print(distances)
    distances_to_other = distances[rownames(distances) == target_cluster, colnames(distances) %in% other_clusters]
    closest_cluster = names(which(distances_to_other == min(distances_to_other)))
    if(is.numeric(target_cluster)){closest_cluster = as.numeric(closest_cluster)}
  }else{
    closest_cluster= other_clusters 
  }
  return(closest_cluster)
}

# 
# siblings_to_merge = function(sibling_clusters, markers_comparisons_siblings,max_pvalue_prune,min_specificity,min_sibling_markers){
#   closest_clusters = c()
#   for(c in sibling_clusters){
#     c_markers = markers_comparisons_siblings %>% dplyr::filter(cluster_id == c) %>% 
#       dplyr::arrange(desc(specificity)) %>%
#       dplyr::filter( p_val_adj < max_pvalue_prune & specificity > min_specificity)
#     message(c,": ",nrow(c_markers))
#   }
#   if(nrow(c_markers) < min_sibling_markers){closest_clusters = c(closest_clusters,c)}
#   return(closest_clusters)
# }


##########
### Run mrtree pruning based on markers
##########

# for iterative run test
# edgelist = edgelist_updated
# edgelist = dplyr::left_join(edgelist,mrtree_result$edgelist[,c("to","level")],by=c("to"="to"))
# labelmat = labelmat_updated
# harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data[,1:original_cols],labelmat)

# get key objects
message(Sys.time(),": Prune clusters .." )

# init edgelist and all_nodes
edgelist = edgelist[,c("from","to","level")]
cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
  dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
all_nodes = unique(edgelist$to)

## iterate over all nodes
merge_ncells = vector()
merge_siblingmarkers = vector()
merge_donorPct = vector()
merge_list = list()
for(n in 1:length(all_nodes)){
  
  # get information
  current_node = all_nodes[n]
  message(current_node)
  #  message("At: ",current_node," with ",length(labelmat[labelmat[,current_level] == current_node,current_level])," cells")
  parent_node = edgelist$from[edgelist$to==current_node]
  sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
  current_level = edgelist$clusterlevel[edgelist$to==current_node]
  if(edgelist$level[edgelist$to==current_node] < parameter_list$min_prune_level){next}
  
  ## Check if the parent node was already merged: --> use merge_list
  other_parent_nodes = as.character(unlist(sapply(merge_list,function(x,p){if(p %in% x) return(x[x!=p]) else return(NULL) },p=parent_node)))
  # if yes:
  if(length(other_parent_nodes) > 0){
    # update sibling nodes with children of that node
    sibling_nodes = c(sibling_nodes,edgelist$to[edgelist$from %in% other_parent_nodes & edgelist$to != current_node])
    message("For node ",current_node," an additional parent_node is included: ",other_parent_nodes," as this one will be merged with the orginal parent_node: ",parent_node," .")
    message("Sibling_nodes included for merging: ",sibling_nodes)
    # update parent markers with union of both parents
    if(parameter_list$merge_marker_based){
      parent_global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id %in% c(parent_node,other_parent_nodes)) %>% 
        dplyr::arrange(desc(specificity)) %>%
        dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity) %>%
        dplyr::distinct(gene,.keep_all = TRUE)
    }
  }
  
  # filter siblings by ncells
  # sibling_nodes = sibling_nodes[sibling_nodes %in% edgelist$to[edgelist$ncells >= parameter_list$min_cells]]
  
  # if no sibling nodes: go to next step
  if(length(sibling_nodes) == 0){next}
  
  #### section for ncells of current
  ncells_current_node=edgelist$ncells[edgelist$to==current_node]
  
  #### section for marker checks
  
  if(parameter_list$merge_marker_based){
    # get number of sibling markers:
    sibling_markers = markers_comparisons_siblings %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
      dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)
    
    # get global markers:
    global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == current_node) %>% dplyr::arrange(desc(specificity)) %>%
      dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)
    
    # parent global markers:
    parent_global_markers = markers_comparisons_all %>% dplyr::filter(cluster_id == parent_node) %>% dplyr::arrange(desc(specificity)) %>%
      dplyr::filter( p_val_adj < parameter_list$max_pvalue_prune & specificity > parameter_list$min_specificity)
  }
  
  #### section for sample contribution
  if(parameter_list$merge_sample_based){
    pct_donorMax = cluster_qc_overview$pct_donorMax[cluster_qc_overview$cluster == current_node]
  }
  ## check all potential flags and set merge_node to TRUE
  merge_node = FALSE
  if(ncells_current_node < parameter_list$min_cells){ merge_node = TRUE;merge_ncells = c(merge_ncells,current_node) }
  if(parameter_list$merge_marker_based){ if(nrow(sibling_markers) < parameter_list$min_sibling_markers){merge_node = TRUE;merge_siblingmarkers = c(merge_siblingmarkers,current_node)} }
  if(parameter_list$merge_sample_based){ if(pct_donorMax > parameter_list$pct_threshold ){merge_node = TRUE;merge_donorPct = c(merge_donorPct,current_node)} }
  
  # if any n_sibling_markers are < min_sibling_markers or ncells < min_cells
  if(merge_node){
    
    # if(nrow(parent_global_markers) > 3){
    #   # get the expression of the parent markers in current cluster:
    #   current_gene_expression = Seurat::FetchData(harmonized_seurat_object,vars = unique(parent_global_markers$gene),cells = harmonized_seurat_object@meta.data$Cell_ID[harmonized_seurat_object@meta.data[,current_level] == current_node])
    #   current_gene_expression_mean = colMeans(current_gene_expression)
    #   
    #   # iterate over all sibling markers:
    #   intersection_lengths = c()
    #   coexpr = c()
    #   for(sib in sibling_nodes){
    #     # calc co-expression using marker genes
    #     sibling_gene_expression = Seurat::FetchData(harmonized_seurat_object,vars = unique(parent_global_markers$gene),cells = harmonized_seurat_object@meta.data$Cell_ID[harmonized_seurat_object@meta.data[,current_level] == sib])
    #     sibling_gene_expression_mean = colMeans(sibling_gene_expression)
    #     coexpr[sib] = cor(current_gene_expression_mean,sibling_gene_expression_mean)
    #   }
    #   
    #   # add node with highest number of shared global markers
    #   merge_nodes = names(coexpr)[coexpr == max(coexpr)]
    # }else{
    #   merge_nodes = sibling_nodes[1] # fall back
    # }
    # merge_nodes = siblings_to_merge(sibling_clusters = sibling_nodes, 
    #                                markers_comparisons_siblings = markers_comparisons_siblings,
    #                                max_pvalue_prune = parameter_list$max_pvalue_prune,
    #                                min_specificity = parameter_list$min_specificity ,
    #                                min_sibling_markers = parameter_list$min_sibling_markers)
    
    # use euclidean distance in latent space to define nearest sibling
    merge_nodes = closest_cluster(object = harmonized_seurat_object,
                                  target_cluster = current_node,
                                  other_clusters = sibling_nodes,
                                  cluster_col = current_level,
                                  reduction="scvi")
    

    if(length(merge_nodes) < 1){stop("fqu")}
    message("Merging node ",current_node," (",n,"/",length(all_nodes),") into node(s) ",paste0(merge_nodes,sep=", "))
    all_nodes_to_merge = c(merge_nodes,current_node)
    merge_list[[paste0(all_nodes_to_merge,collapse = "_")]] = all_nodes_to_merge
  }
}


#UpSetR::upset(UpSetR::fromList(list(merge_siblingmarkers=merge_siblingmarkers,merge_donorPct=merge_donorPct,merge_ncells=merge_ncells)),nsets = 5,nintersects = 20,text.scale = 3)

##########
### Create pruned edgelist
##########

message(Sys.time(),": Update labels and save .." )

edgelist_updated = edgelist
labelmat_updated = labelmat
merge_list2 = merge_list

for(i in 1:length(merge_list2)){
  nodes_to_merge = merge_list2[[i]]
  print(nodes_to_merge)
  # only update edgelist and labelmat if there are truly 2 unique labels remaining in the current entry
  if(length(unique(nodes_to_merge)) > 1){
    # use the first node to label
    merge_node = as.character(sort((nodes_to_merge))[1])
    print(paste0(" >> merge to ",merge_node))
    # update edgelist
    edgelist_updated$from[edgelist_updated$from %in% nodes_to_merge] = merge_node
    edgelist_updated$to[edgelist_updated$to %in% nodes_to_merge] = merge_node
    # remove repeated entries
    edgelist_updated = edgelist_updated %>% distinct(from,to,clusterlevel) # merge duplicate rows to one entry and updating the cell total
    # update labelmat
    current_level = edgelist$clusterlevel[edgelist$to==merge_node]
    labelmat_updated[labelmat_updated[,current_level] %in% nodes_to_merge,current_level] = merge_node
    # change all occurrences in merge_list
    merge_list2 = lapply(merge_list2,function(x,remove_nodes,new_node){x[x %in% remove_nodes] = new_node;return(x)},remove_nodes=nodes_to_merge,new_node=merge_node)
    
  }
}

##########
### Update labels in edgelist and labelmat
##########

old_prefix = parameter_list$old_prefix # "K"
new_prefix = parameter_list$new_prefix # "C"
all_cluster_levels_updated = edgelist_updated %>% dplyr::group_by(clusterlevel) %>% dplyr::count()
all_cluster_levels_updated$clusterlevel_new = paste0(all_cluster_levels_updated$clusterlevel %>% stringr::str_extract(old_prefix)  %>% stringr::str_replace(old_prefix, new_prefix),all_cluster_levels_updated$n)
all_cluster_levels_updated$clusterlevel_new = make.unique(all_cluster_levels_updated$clusterlevel_new) # this is just to avoid the code afterwards from breaking -- having two levels with the same number of clusters does not make much sense

# create new labels using the pruned numbers of labels per level
all_rows = list()
for(c in 1:nrow(all_cluster_levels_updated)){
  edgelist_updated_current_level = edgelist_updated[edgelist_updated$clusterlevel == all_cluster_levels_updated$clusterlevel[c],]
  for(i in 1:nrow(edgelist_updated_current_level)){
    new_node_label = data.frame(old_node = edgelist_updated_current_level$to[i],
                                new_node = paste0(all_cluster_levels_updated$clusterlevel_new[c],"-",i),
                                new_cluster_level = all_cluster_levels_updated$clusterlevel_new[c])
    all_rows[[new_node_label$new_node]] = new_node_label
  }
}
new_labels = as.data.frame(do.call(rbind,all_rows))

# make new edgelist
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated,new_labels[,c("old_node","new_node")],by=c("from"="old_node"))
edgelist_updated_new_labels = dplyr::left_join(edgelist_updated_new_labels,new_labels,by=c("to"="old_node"))
edgelist_updated_new_labels = edgelist_updated_new_labels %>% dplyr::select(from = new_node.x, to = new_node.y,clusterlevel = new_cluster_level)
edgelist_updated_new_labels$from[is.na(edgelist_updated_new_labels$from)] = "all"

# make new labelmat
labelmat_updated_new_labels = apply(labelmat_updated,2,function(x,new_labels){
  new_labels$new_node[match(x,new_labels$old_node)]
},new_labels=new_labels)

colnames(labelmat_updated_new_labels) = stringr::str_extract(labelmat_updated_new_labels[1,],pattern = paste0(new_prefix,"[0-9]+"))


##########
### test explore
##########

# ## add post:
harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data[,1:56],labelmat_updated_new_labels)
colnames(labelmat_updated_new_labels)
DimPlot(harmonized_seurat_object,group.by = "C331",reduction = "umap_scvi_HypoNeuron",label = TRUE,label.size = 2.5,raster.dpi = c(1024,1024),pt.size = 2)+NoLegend()+NoAxes()
# 
# add pre:
harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data[,1:56],labelmat)
colnames(labelmat)
DimPlot(harmonized_seurat_object,group.by = "C4",reduction = "umap_scvi_HypoNeuron",label = TRUE,label.size = 2.5,raster.dpi = c(1024,1024),pt.size = 2)+NoLegend()+NoAxes()

# # stats
# target_cluster = "K126-55"
# all_cluster_stats %>% dplyr::filter(cluster == target_cluster) %>% dplyr::distinct(cluster,Donor_ID,cells_cluster_donorID,pct_cluster_donorID) %>% dplyr::arrange(desc(pct_cluster_donorID))
# markers_comparisons_siblings %>% dplyr::filter(cluster_id == target_cluster & specificity > 1.5 & p_val_adj < 0.001) %>% dplyr::arrange(desc(specificity)) %>% dplyr::select(cluster_id,parent,gene,p_val_adj,specificity)
# target_cluster %in% merge_siblingmarkers
# merge_list[grepl(target_cluster,names(merge_list))]
# 
# subset
cluster1="K39-7"
cellsh = harmonized_seurat_object@meta.data$Cell_ID[ harmonized_seurat_object@meta.data$K39 %in% c(cluster1)]
edgelist[edgelist$from == cluster1,]
edges[edges$from == cluster1,]

DimPlot(harmonized_seurat_object,cells.highlight = cellsh,sizes.highlight = 0.1,reduction = "umap_scvi_HypoNeuron",label = F,raster = F)+NoLegend()

subset_seurat = hypoMapUtils::subset_zoom_seurat(harmonized_seurat_object,cells = cellsh,keep_quantile = 0.99,reduction = "umap_scvi_HypoNeuron")

DimPlot(subset_seurat,group.by = "C455",reduction = "umap_scvi_HypoNeuron",label = TRUE,label.size = 4)+NoLegend()
FeaturePlot(subset_seurat,"LHX1",reduction = "umap_scvi_HypoNeuron",order = TRUE)
DimPlot(subset_seurat,group.by = "Donor_ID",reduction = "umap_scvi_HypoNeuron",label = F,label.size = 2.5)#+NoLegend()

FeaturePlot(harmonized_seurat_object,"AASS",reduction = "umap_scvi_HypoNeuron",order = TRUE)

subset_clusters = unique(subset_seurat@meta.data$K654)
UpSetR::upset(UpSetR::fromList(list(merge_siblingmarkers=merge_siblingmarkers[merge_siblingmarkers %in% subset_clusters],
                                    merge_donorPct=merge_donorPct[merge_donorPct %in% subset_clusters],
                                    merge_ncells=merge_ncells[merge_ncells %in% subset_clusters]
                                    )),nsets = 5,nintersects = 20,text.scale = 3)

##########
### Create pruned result version
##########

# save in data.tree format
require(data.tree)
df = as.data.frame(unique(labelmat_updated_new_labels), stringsAsFactors = F)
df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
tree.datatree = data.tree::as.Node(df)

# export edgelist
edges= data.tree::ToDataFrameNetwork(tree.datatree,"isLeaf","level","count","totalCount","height")
nodes = data.frame(id = c("all",as.character(unique(edges$to))),label=c("all",as.character(unique(edges$to))))
nodes = rbind(c(NA,"all",FALSE,1,5,max(edges$height)+1),edges[,2:ncol(edges)]) %>% dplyr::rename(id = to) %>% dplyr::mutate(label = id)

# make cluster object ? and add to misc
updated_cluster_object = list(labelmat = labelmat_updated_new_labels,
                              edgelist = edges ,
                              nodelist = nodes,
                              data_tree = tree.datatree)

#edges_new2 = edges
##########
### Save
##########

message("Reduced edgelist from: ",nrow(edgelist)," to: ",nrow(updated_cluster_object$edgelist)," . Saving result ...")

saveRDS(updated_cluster_object,parameter_list$mrtree_result_file_output)
#saveRDS(updated_cluster_object,paste0(parameter_list$harmonization_folder_path,parameter_list$clustering_key_name,"_pruned_mrtree_clustering_results",".rds"))

message(Sys.time(),": Complete!" )

##########
### Explore new clusters
##########

# harmonized_seurat_object@meta.data = cbind(harmonized_seurat_object@meta.data,labelmat_updated_new_labels)
# colnames(labelmat_updated_new_labels)
# cells_per_cluster_C376 = harmonized_seurat_object@meta.data %>% dplyr::group_by(C376) %>% dplyr::count(name = "ncells")
# cells_per_cluster_C376 = harmonized_seurat_object@meta.data %>% dplyr::group_by(C376) %>% dplyr::count(name = "ncells")
# 
# cells_per_cluster_C376_K588 = harmonized_seurat_object@meta.data %>% dplyr::group_by(C376,K588) %>% dplyr::count(name = "ncells")
# 
# human_hypo_combined@meta.data$C6
