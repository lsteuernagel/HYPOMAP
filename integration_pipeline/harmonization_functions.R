##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){
  
  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames
  if(is.character(current_embedding[,1])){
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      message("Found ",length(rnames)," rows in embedding and ",length(rownames(metadata))," rows in metadata.")
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
}

##########
### clear_clustering
##########

#' Eliminate small clusters from a vector of labels substituting with the labels of NN
#' @param x vector of labels
#' @param min_cells minimum number of cells to keep cluster
#' @param nn_idx matrix of cells x k NN idx --> e.g. output of annoy or rann
#' @return updated vector of labels

clear_clustering = function(x,min_cells,nn_idx){
  x = as.character(x)
  new_x=x
  # which clusters are too small ?
  small_clusters = names(table(x)[table(x) < min_cells])
  # go through small clusters and move cells to laregr clusters based on neighbors
  if(length(small_clusters)>0){
    #message("Removing ",length(small_clusters)," clusters")
    for(i in 1:length(small_clusters)){
      if(i %% 1000 == 0){message(i, " of ",length(small_clusters)," clusters to remove.")}
      current_cluster = small_clusters[i]
      which_idx = which(x == current_cluster)
      # get idx for k NN
      neighbor_idx = nn_idx[which_idx,]
      # substitute with cluster labels
      if(length(which_idx)>1){
        neighbor_clusters = apply(neighbor_idx,1,function(z,cluster_vector){return(cluster_vector[z])},cluster_vector=x)
        # extract that most common label in neighbors
        clusters_vote = apply(neighbor_clusters,2,function(z,exclude_clusters){
          return(names(sort(table(z),decreasing = TRUE))[! names(sort(table(z),decreasing = TRUE)) %in% exclude_clusters][1])
        },exclude_clusters=small_clusters)
      }else{
        neighbor_clusters = x[neighbor_idx]
        clusters_vote = names(sort(table(neighbor_clusters),decreasing = TRUE))[! names(sort(table(neighbor_clusters),decreasing = TRUE)) %in% small_clusters][1]
      }
      # overwrite cluster label
      new_x[which_idx]=clusters_vote
    }
  }
  
  return(new_x)
}


##########
### annotate_tree
##########


#' Add gene based annotation to tree labels
#' @param edgelist edgelist from mrtree corresponding to labelmat
#' @param labelmat dataframe of clusterlabels per cell with one column per level
#' @param markers_comparisons_all cluster markers for all clusters in labelmat. expects format from run_marker_detection.R
#' @param markers_comparisons_siblings cluster markers vs sibling for all clusters in labelmat. expects format from run_marker_detection.R
#' @param manual_names a named vector with names of clusters in edgelist --> will overwrite with manual names
#' @param overwrite_with_manual whether to overwrite the full tree until this point with the manual name or just append similar to the best gene. defaults to FALSE
#' @param manual_exclude_genes a set of genes that should not be included in anno
#' @param max_pval_adj max pavlue allowed for a marker to be considered
#' @param min_combined_score min_specificity allowed for a marker to be considered
#' @param lower_pct_min when calculating the specificty score add this pseudo-count to pct.2 to avoid favouring of highly specific but low abundant genes
#' @param min_cells min_cells
#' @param reverse_order todo: add explanation
#' @return list with annotated_nodes_gene, all_combined_markers, all_combined_markers_full)

annotate_tree = function(edgelist,labelmat,markers_comparisons_all,markers_comparisons_siblings,manual_names=c(),overwrite_with_manual=FALSE,manual_exclude_genes=character(0),
                         max_pval_adj=0.0001, min_combined_score = 5,lower_pct_min=0.1,min_cells=20,max_score_siblings_children = 20,reverse_order = FALSE){
  # init edgelist and all_nodes
  edgelist = edgelist[,c("from","to","level")]
  cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
    dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
  edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
  all_nodes = unique(edgelist$to)
  
  # annotation
  edgelist = edgelist %>% dplyr::arrange(level) # ensure top down order for tree traversal with for loop!+
  
  # list to store intermediate results
  annotation_list = list()
  annotated_nodes_gene=c()
  all_combined_markers = list()
  all_combined_markers_full = list()
  
  message("Running annotation")
  message("Using ",length(manual_names)," manually provided names to overwrite annotation for specific clusters.")
  
  # for each node in edgelist:
  for(n in 1:length(all_nodes)){
    
    # get information
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    children_nodes = scUtils::find_children(current_node,edgelist)
    direct_children_nodes = edgelist$to[edgelist$from==current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]
    
    message(current_node)
    # message("parent: ",parent_node)
    
    # get parent annotation:
    parent_genes = unlist(annotated_nodes_gene[[parent_node[1]]])
    
    # global markers
    markers_all_self = markers_comparisons_all %>% 
      dplyr::filter(cluster == current_node) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    # global markers siblings
    markers_all_siblings = markers_comparisons_all %>% 
      dplyr::filter(cluster %in% sibling_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    
    # sibling markers
    markers_siblings_self = markers_comparisons_siblings %>% 
      dplyr::filter(cluster == current_node) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    
    # VS sibling markers
    markers_siblings_siblings = markers_comparisons_siblings %>% 
      dplyr::filter(cluster %in% sibling_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    # ...
    # sibling markers of children vs all
    markers_all_children = markers_comparisons_all %>% 
      dplyr::filter(cluster %in% direct_children_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_)) #%>%
    # sibling markers between children
    markers_siblings_children = markers_comparisons_siblings %>% 
      dplyr::filter(cluster %in% children_nodes) %>% # direct_children_nodes or children_nodes
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_)) %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(abs(score) == max(abs(score))) %>%
      dplyr::distinct(gene,.keep_all = TRUE)
    
    #potential_descriptive_markers_of_siblings
    # combined scoring
    combined_markers = dplyr::full_join(markers_all_self[,c("gene","score")],markers_siblings_self[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_all_siblings[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_siblings_siblings[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_all_children[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_siblings_children[,c("gene","score")],by="gene") %>%
      dplyr::filter(! gene %in% parent_genes) %>%
      dplyr::filter(! gene %in% manual_exclude_genes) %>%
      dplyr::filter(! grepl("LINC|AL[0-9]+|AP[0-9]+|AC[0-9]+|MIR|-AS|-PS",gene))
    
    
    
    if(nrow(combined_markers) >= 1){
      colnames(combined_markers)[2:ncol(combined_markers)] =c("score_all_self","score_siblings_self","score_all_siblings","score_siblings_siblings","score_all_children","score_siblings_children")
      combined_markers[is.na(combined_markers)] = 0
      if(length(direct_children_nodes) <= 1){ combined_markers$score_all_children = combined_markers$score_all_self}
      combined_markers$combined_score = ((combined_markers$score_all_self+1) * (combined_markers$score_siblings_self+1)) / (abs(combined_markers$score_siblings_children)+1)
      
      combined_markers_full = combined_markers
      
      # make releavnt marker table
      combined_markers = combined_markers %>% 
        dplyr::ungroup() %>%
        dplyr::filter(combined_score > min_combined_score & abs(score_siblings_children) <= max_score_siblings_children) %>%
        dplyr::group_by(gene)  %>%
        dplyr::slice_min(order_by = combined_score,with_ties = FALSE,n = 1) %>%
        # dplyr::filter(combined_score == min(combined_score)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(combined_score))
      
      # add result
      all_combined_markers[[current_node]] = combined_markers
      all_combined_markers_full[[current_node]] = combined_markers_full 
      if(length(sibling_nodes) > 0){
        if(nrow(combined_markers) > 0){
          annotated_nodes_gene[[current_node]] = c(parent_genes,combined_markers$gene[1])
        }else{
          annotated_nodes_gene[[current_node]] = c(parent_genes,"TBD")
        }
      }else{
        annotated_nodes_gene[[current_node]] = c(parent_genes)
      }
    }else{
      annotated_nodes_gene[[current_node]] = c(parent_genes,"TBD")
    }
    if(current_node %in% names(manual_names)){
      annotated_nodes_gene[[current_node]]  = c(parent_genes,manual_names[current_node])
    }
    
    res_list =list(annotated_nodes_gene = annotated_nodes_gene,all_combined_markers = all_combined_markers,all_combined_markers_full = all_combined_markers_full)
    message(paste0(annotated_nodes_gene[[current_node]],collapse = "_"))
    
  }
  message("Returning results")
  return(res_list)
}



##########
### annotate_level
##########


#' Add gene based annotation to tree labels
#' @param edgelist edgelist from mrtree corresponding to labelmat
#' @param labelmat dataframe of clusterlabels per cell with one column per level
#' @param markers_comparisons_all cluster markers for all clusters in labelmat. expects format from run_marker_detection.R
#' @param markers_comparisons_siblings cluster markers vs sibling for all clusters in labelmat. expects format from run_marker_detection.R
#' @param manual_names a named vector with names of clusters in edgelist --> will overwrite with manual names
#' @param overwrite_with_manual whether to overwrite the full tree until this point with the manual name or just append similar to the best gene. defaults to FALSE
#' @param manual_exclude_genes a set of genes that should not be included in anno
#' @param max_pval_adj max pavlue allowed for a marker to be considered
#' @param min_combined_score min_specificity allowed for a marker to be considered
#' @param lower_pct_min when calculating the specificty score add this pseudo-count to pct.2 to avoid favouring of highly specific but low abundant genes
#' @param pct.2_max max pct.2 to choose gene for anno
#' @param min_cells min_cells
#' @param reverse_order todo: add explanation
#' @return list with annotated_nodes_gene, all_combined_markers, all_combined_markers_full)

annotate_level = function(edgelist,labelmat,level,markers_comparisons_all,markers_comparisons_siblings,manual_names=c(),overwrite_with_manual=FALSE,manual_exclude_genes=character(0),
                          max_pval_adj=0.0001, min_combined_score = 5,lower_pct_min=0.1,pct.2_max=1,min_cells=20,max_score_siblings_children = 20,max_score_self=100,max_score_sibling=100,reverse_order = FALSE,parent_names=NULL){
  # init edgelist and all_nodes
  edgelist = edgelist[,c("from","to","level")]
  cluster_levels = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")  %>% dplyr::group_by(cluster) %>%
    dplyr::add_count(name="ncells") %>% dplyr::distinct(clusterlevel,cluster,ncells)
  edgelist = dplyr::left_join(edgelist,cluster_levels,by=c("to"="cluster")) %>% dplyr::arrange(level)
  all_nodes = unique(edgelist$to[grepl(level,edgelist$to)])
  
  # annotation
  edgelist = edgelist %>% dplyr::arrange(level) # ensure top down order for tree traversal with for loop!+
  
  # list to store intermediate results
  annotation_list = list()
  annotated_nodes_gene=c()
  all_combined_markers = list()
  all_combined_markers_full = list()
  
  message("Running annotation")
  message("Using ",length(manual_names)," manually provided names to overwrite annotation for specific clusters.")
  
  # for each node in edgelist:
  for(n in 1:length(all_nodes)){
    
    # get information
    current_node = all_nodes[n]
    parent_node = edgelist$from[edgelist$to==current_node]
    sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
    children_nodes = scUtils::find_children(current_node,edgelist)
    direct_children_nodes = edgelist$to[edgelist$from==current_node]
    current_level = edgelist$clusterlevel[edgelist$to==current_node]
    
    message(current_node)
    # message("parent: ",parent_node)
    
    # get parent annotation:
    if(!is.null(parent_names)){
      parent_genes = unlist(parent_names[[parent_node[1]]])
    }else{
      parent_genes=c() 
    }
    
    # global markers
    markers_all_self = markers_comparisons_all %>% 
      dplyr::filter(cluster == current_node) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    # global markers siblings
    markers_all_siblings = markers_comparisons_all %>% 
      dplyr::filter(cluster %in% sibling_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    
    # sibling markers
    markers_siblings_self = markers_comparisons_siblings %>% 
      dplyr::filter(cluster == current_node) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    
    # VS sibling markers
    markers_siblings_siblings = markers_comparisons_siblings %>% 
      dplyr::filter(cluster %in% sibling_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_))
    # ...
    # sibling markers of children vs all
    markers_all_children = markers_comparisons_all %>% 
      dplyr::filter(cluster %in% direct_children_nodes) %>% 
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_)) #%>%
    # sibling markers between children
    markers_siblings_children = markers_comparisons_siblings %>% 
      dplyr::filter(cluster %in% children_nodes) %>% # direct_children_nodes or children_nodes
      dplyr::arrange(desc(specificity)) %>%
      dplyr::filter(p_val_adj< max_pval_adj ) %>% 
      dplyr::mutate(pct.2_min = pmax(lower_pct_min,pct.2), pct.1_min = pmax(lower_pct_min,pct.1)) %>%
      dplyr::mutate(score = case_when(avg_log2FC >= 0 ~ (pct.1 / pct.2_min) * avg_log2FC,
                                      avg_log2FC < 0 ~ (pct.2 / pct.1_min) * avg_log2FC,
                                      TRUE ~ NA_real_)) %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(abs(score) == max(abs(score))) %>%
      dplyr::distinct(gene,.keep_all = TRUE)
    
    #potential_descriptive_markers_of_siblings
    # combined scoring
    combined_markers = dplyr::full_join(markers_all_self[,c("gene","score","pct.2")],markers_siblings_self[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_all_siblings[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_siblings_siblings[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_all_children[,c("gene","score")],by="gene") %>%
      dplyr::left_join(markers_siblings_children[,c("gene","score")],by="gene") %>%
      dplyr::filter(! gene %in% parent_genes)# %>%
      
    
    if(nrow(combined_markers) >= 1){
      colnames(combined_markers)[2:ncol(combined_markers)] =c("score_all_self","pct.2","score_siblings_self","score_all_siblings","score_siblings_siblings","score_all_children","score_siblings_children")
      combined_markers[is.na(combined_markers)] = 0
      if(length(direct_children_nodes) <= 1){ combined_markers$score_all_children = combined_markers$score_all_self}
      # if(length(sibling_nodes) ==0){ 
      #   combined_markers$score_all_self = combined_markers$score_all_self / 2
      #   combined_markers$score_siblings_self = combined_markers$score_all_self
      # }
      
      ## making scoring:
      
      # expects both numbers to be positive (so that max agle can be 90 degree (1.570796 rad)
      angular_similarity = function(A){
        x1 = acos(sum(A*c(0,1)) / sqrt(sum(A^2)*sum(c(0,1)^2)))
        x2 = acos(sum(A*c(1,0)) / sqrt(sum(A^2)*sum(c(1,0)^2)))
        angular_sim = 1 - round((max(x1,x2)-0.7853982) / (1.570796-0.7853982),4) # -90 degree - 45 degree 
        if(angular_sim < 0){angular_sim=0}
        return(angular_sim)
      }
      # weigh with angular distance
      # if(length(sibling_nodes) ==0){ 
      #   combined_markers$combined_score = (combined_markers$score_all_self+1) / abs(combined_markers$score_siblings_children+1)
      #   combined_markers$score1 = (combined_markers$score_all_self+1)
      #   combined_markers$score2 =  abs(combined_markers$score_siblings_children)+1
      # }else{
      #   combined_markers$score1 = ((apply(combined_markers[,c("score_all_self","score_siblings_self")],1,angular_similarity)*(combined_markers$score_all_self*combined_markers$score_siblings_self))+1) 
      #   combined_markers$score2 = abs(combined_markers$score_siblings_children)+1
      #   combined_markers$combined_score = ((apply(combined_markers[,c("score_all_self","score_siblings_self")],1,angular_similarity)*(combined_markers$score_all_self*combined_markers$score_siblings_self))+1) / (abs(combined_markers$score_siblings_children)+1)
      #   
      # }
      # 
      # min-max normalise the scores
      if(length(sibling_nodes) ==0){ 
        combined_markers$score1 = (pmin(combined_markers$score_all_self,max_score_self) - min(combined_markers$score_all_self)) / (min(max(combined_markers$score_all_self),max_score_self) - min(combined_markers$score_all_self))
        combined_markers$score2 = 0
        combined_markers$score3 =  abs(combined_markers$score_siblings_children)+1
        combined_markers$combined_score = combined_markers$score1 / combined_markers$score3 * max(combined_markers$score_all_self)
      }else{
        combined_markers$score1 = (pmin(combined_markers$score_all_self,max_score_self) - min(combined_markers$score_all_self)) / (min(max(combined_markers$score_all_self),max_score_self) - min(combined_markers$score_all_self))
        combined_markers$score2 = (pmin(combined_markers$score_siblings_self,max_score_sibling) - 0) / (min(max(combined_markers$score_siblings_self),max_score_sibling) - 0) # changed to 0
        combined_markers$score3 = abs(combined_markers$score_siblings_children)+1
        combined_markers$combined_score = combined_markers$score1 * combined_markers$score2 / combined_markers$score3 * max(combined_markers$score_all_self) * max(combined_markers$score_siblings_self)
      }
      #combined_markers$combined_score = ((combined_markers$score_all_self+1) * (combined_markers$score_siblings_self+1)) / (abs(combined_markers$score_siblings_children)+1)
      combined_markers$cluster = current_node
      combined_markers_full = combined_markers
      
      # make releavnt marker table
      combined_markers = combined_markers %>% 
        dplyr::ungroup() %>%
        dplyr::filter(! gene %in% manual_exclude_genes) %>%
        dplyr::filter(! grepl("LINC|AL[0-9]{3,}|AP[0-9]{3,}|AC[0-9]{3,}|MIR|-AS|-PS|orf[0-9]{2,}|[0-9]{2,}\\.[0-9]",gene)) %>%
        dplyr::filter(combined_score > min_combined_score & abs(score_siblings_children) <= max_score_siblings_children & pct.2 < pct.2_max) %>%
        dplyr::group_by(gene)  %>%
        dplyr::slice_min(order_by = combined_score,with_ties = FALSE,n = 1) %>%
        # dplyr::filter(combined_score == min(combined_score)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(combined_score))
      
      
      # add result
      all_combined_markers[[current_node]] = combined_markers
      all_combined_markers_full[[current_node]] = combined_markers_full  %>%
        dplyr::group_by(gene)  %>%
        dplyr::slice_min(order_by = combined_score,with_ties = FALSE,n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(combined_score))
      # if(length(sibling_nodes) > 0){
      if(nrow(combined_markers) > 0){
        annotated_nodes_gene[[current_node]] = c(combined_markers$gene[1])
      }else{
        annotated_nodes_gene[[current_node]] = c("NA")
      }
      #   }else{
      #     annotated_nodes_gene[[current_node]] = c(parent_genes)
      #   }
      # }else{
      #   annotated_nodes_gene[[current_node]] = c(parent_genes,"TBD")
      # }
      if(current_node %in% names(manual_names)){
        annotated_nodes_gene[[current_node]]  = as.character(manual_names[current_node])
      }
      
      res_list =list(annotated_nodes_gene = annotated_nodes_gene,all_combined_markers = all_combined_markers,all_combined_markers_full = all_combined_markers_full)
      message(paste0(annotated_nodes_gene[[current_node]],collapse = "_"))
      
    }
  }
  message("Returning results")
  return(res_list)
}


##########
### findMarkers_tree2
##########

#' Walk through a mrtree / hierchical clustering tree using an edgelist and calculate Marker genes using Seurat's FindMarkers or a stratified version:
#' https://github.com/KChen-lab/stratified-tests-for-seurat
#' Two results: Vs sibling and Vs all . Runs in parallel to be faster
#' @param seurat_object seurat_object to call FindMarkers
#' @param edgelist minimum number of cells to keep cluster
#' @param labelmat
#' @param n_cores doParallel cores
#' @param use_stratified use stratfied
#' @param batch_var batch_var for stratfied mode
#' @param comparison_types 'All' or/and 'Sibling'
#' @param assay
#' @param slot
#' @param ... passed to detection functions
#' @return updated vector of labels

findMarkers_tree2 = function(seurat_object,edgelist,labelmat,n_cores=1,use_stratified=TRUE,test.use="wilcox-stratified",batch_var="Batch_ID",latent_vars_seurat=NULL,genes_to_include=NULL,comparison_types = c("Sibling","All"),only.pos=TRUE,assay="RNA",slot="data",...){
  
  require(doParallel)
  require(tidyselect)
  require(foreach)
  
  # info:
  message(n_cores)
  
  # set edglist
  edgelist = edgelist[,c("from","to")]
  
  # genes
  if(is.null(genes_to_include)){
    genes_to_include = rownames(seurat_object@assays[[assay]]$counts)
  }
  
  # add labelmat to seurat
  if(any(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
    if(!all(colnames(labelmat) %in% colnames(seurat_object@meta.data))){
      stop("Cannot handle duplicated column names in seurat_object@meta.data and labelmat from mrtree. Either include all labelmat columns beforehand or none at all.")
    }
  }else{
    seurat_object@meta.data = cbind(seurat_object@meta.data,labelmat )
  }
  
  # handle only.pos
  if(length(only.pos) > 1){
    message("only.pos: ",only.pos)
    if(length(only.pos) == length(comparison_types) & all(names(only.pos) %in% c("All","Sibling"))){
      message("Using comparison type specific only.pos")
    }else{
      stop("Either provide one boolean value for only.pos or a vector of booleans with the same length as comparison_types") 
    }
  }
  
  # init
  current_node="all"
  colnames(edgelist) = c("from","to")
  label_mat_long = as.data.frame(labelmat) %>% tidyr::pivot_longer(everything(),names_to = "clusterlevel", values_to = "cluster")
  edgelist = dplyr::left_join(edgelist,label_mat_long,by=c("to"="cluster")) %>% dplyr::distinct(from,to,clusterlevel)
  all_nodes = unique(edgelist[,2])
  #all_nodes = all_nodes[!all_nodes %in% c("all","root")]
  
  comparisons_siblings = NULL
  comparisons_all = NULL
  
  registerDoParallel(cores=n_cores)
  
  return_list=list()
  for(comp in comparison_types){
    if(comp == "Sibling"){
      if(length(only.pos) > 1){
        only.pos_local = only.pos["Sibling"]
      }else{
        only.pos_local = only.pos 
      }
      message("Sibling Comparisons with ",test.use, " and only.pos = ",only.pos_local)
    }
    if(comp == "All"){
      if(length(only.pos) > 1){
        only.pos_local = only.pos["All"]
      }else{
        only.pos_local = only.pos 
      }
      message("All Comparisons with ",test.use, " and only.pos = ",only.pos_local)
    }
    comparisons_main <- foreach(n = 1:length(all_nodes),.errorhandling = 'remove', .combine='rbind') %dopar% {
      current_node = all_nodes[n]
      parent_node = edgelist$from[edgelist$to==current_node]
      sibling_nodes = edgelist$to[edgelist$from==parent_node & edgelist$to != current_node]
      current_level = edgelist$clusterlevel[edgelist$to==current_node]
      # set ident to current level!
      Idents(seurat_object) = current_level
      
      #decide what to call against:
      cluster_1 = current_node
      if(comp == "Sibling"){
        cluster_2 = sibling_nodes
      }
      if(comp == "All"){
        cluster_2 = all_nodes[! all_nodes %in% current_node]
        cluster_2 = cluster_2[cluster_2 %in% edgelist$to[edgelist$clusterlevel == current_level]] # need to ensure that only current level is used!
      }
      #calculate markers
      #  current_markers <- tryCatch({
      if(length(cluster_2)>0){
        if(test.use == "wilcox-stratified"){
          current_markers = FindMarkers2.Seurat(object = seurat_object,
                                                ident.1 = cluster_1,
                                                ident.2 = cluster_2,
                                                assay = assay_markers,
                                                features = genes_to_include,
                                                logfc.threshold = logfc.threshold,
                                                slot = assay_slot,
                                                test.use = "VE",
                                                min.pct = min.pct,
                                                min.diff.pct = min.diff.pct,
                                                max.cells.per.ident=max.cells.per.ident,
                                                min.cells.feature = min.cells.feature,
                                                min.cells.group = min.cells.group,
                                                return.thresh = 1,
                                                base = base,
                                                only.pos = only.pos_local,
                                                latent.vars = batch_var,
                                                genre = "locally-best")
          if(base==2){
            colnames(current_markers)[colnames(current_markers)=="avg_logFC"] = "avg_log2FC"
          }
        }else{
          current_markers=Seurat::FindMarkers(object = seurat_object,
                                              ident.1 = cluster_1,
                                              ident.2 = cluster_2,
                                              assay = assay_markers,
                                              logfc.threshold = logfc.threshold,
                                              features = genes_to_include,
                                              slot = assay_slot,
                                              test.use =test.use,
                                              min.pct = min.pct,
                                              latent.vars = latent_vars_seurat,
                                              min.diff.pct = min.diff.pct,
                                              max.cells.per.ident=max.cells.per.ident,
                                              min.cells.feature = min.cells.feature,
                                              min.cells.group = min.cells.group,
                                              base = base,
                                              only.pos = only.pos_local)
        }
        current_markers$gene = rownames(current_markers)
        current_markers$cluster_id = cluster_1
        current_markers$comparison = comp
        current_markers$parent = parent_node
        current_markers
      }else{
        NULL
      }
      # },error=function(cond) {
      #   message(cond)
      #   # Choose a return value in case of error
      #   return(NULL)
      # })
      
      # return
      #current_markers
    }
    return_list[[paste0("comparisons_",comp)]] = comparisons_main
  }
  # return
  # print(return_list[1:4])
  return(return_list)
}
