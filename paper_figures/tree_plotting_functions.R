
## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

##########
### plot_cluster_tree
##########

#' Plot metadata entries in cluster tree using ggtree
#'
#' Plots HypoMap cluster tree
#'
#' TODO: not tested in package yet!
#'
#' See also addHeatmap
#'
#' https://github.com/YuLab-SMU/ggtree/issues/400#issuecomment-845781670
#'
#' @param edgelist TODO
#' @param leaf_level which level to use as leaves ?
#' @param anno_df provide a dataframe which contains mappings of all cluster ids and names in the cluster tree. Defaults to NULL which means automatic construction.
#' @param metadata metadata with cluster ids and names, if anno_df is NULL this has to be provided to construct anno_df
#' @param level_pattern regex for cluster level, if anno_df is NULL this has to be provided to construct anno_df. defaults to 'K0-9+'
#' @param cluster_id_pattern regex for cluster id column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_pruned'
#' @param cluster_name_pattern regex for cluster name column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_named'
#' @param label_size label size on tree. defaults to 2
#' @param label_size_tip label size of lowest level on tree. label_size
#' @param show_genes TRUE
#' @param pruned_edgelist ...
#' @param vjust_label ...
#' @param annotate_reverse ...
#' @param label_lowest description
#' @param na_color ...
#' @param edge_color ...
#' @param node_color ...
#'
#' @return ggtree object or plot
#'
#' @export
#'
#' @import dplyr ggplot2 scales treeio tidytree ggtree
#' @importFrom igraph graph_from_edgelist
#'

plot_cluster_tree = function(edgelist,leaf_level=NULL,anno_df=NULL,metadata=NULL,pruned_edgelist=NULL,level_pattern = "K[0-9]+",cluster_id_pattern = "_pruned",
                             cluster_name_pattern = "_named",branch_length_base = 1,label_size = 2,label_size_tip = label_size, show_genes = TRUE,show_genes_tip=TRUE,vjust_label = -0.5,nudge_x = 0,nudge_x_tip = 0, nudge_y=0,nudge_y_tip = 0,
                             annotate_reverse =TRUE, label_lowest =TRUE,ladderize_right =FALSE,order_column = NULL,na_color = "grey90",edge_color="black",node_color="black",label_color = "darkred"){
  
  # I am loading packages here instead of using them as part of a pachage and playing around with namespaces because for some reason the fucks up the ggtree package
  # for example geom_nodelab behaves different when called via ggtree::geom_nodelab, indicating taht is uses some other version?
  # also ggtree has some bugs, like not working with dplyr > 1.06
  # library(dplyr)
  # library(igraph)
  # library(ggplot2)
  # library(scales)
  # library(treeio)
  # library(tidytree)
  # library(ggtree)
  
  # --> let's see if this works
  
  # check that req columns exist
  if(length(setdiff(c("from","to","level"),colnames(edgelist)))>0){
    
    warning("Error: Wrong edgelist format. Requires columns: from, to, level")
    return(NULL)
  }
  # check that leaf_level is there
  if(!leaf_level %in% edgelist$level){
    warning("Error: leaf_level '",leaf_level,"' cannot be found in level column of edgelist")
    return(NULL)
  }
  
  # construct a dataframe with the required annotations
  if(is.null(anno_df)){
    if(is.null(metadata)){
      warning("Error: Please provide metadata with corresponding cluster ids and names that match provided edgelist")
      return(NULL)
    }
    if(!any(grepl(cluster_id_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_id_pattern '",cluster_id_pattern,"' in provided metadata")
      return(NULL)
    }
    if(!any(grepl(cluster_name_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_name_pattern '",cluster_name_pattern,"' in provided metadata")
      return(NULL)
    }
    # anno_df: cluster id, clean_names, clusterlevel, ncells, first_cluster_name
    pruned_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_id_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = stringr::str_extract(colname,level_pattern)) %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    named_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_name_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = stringr::str_extract(colname,level_pattern))  %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    both_map = dplyr::left_join(pruned_ids,named_ids,by=c("Cell_ID"="Cell_ID","clusterlevel"="clusterlevel")) %>% dplyr::select(cluster_id = id.x,cluster_name = id.y)
    
    # anno_df = neuron_map_seurat@misc$pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
    anno_df =pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
    if(annotate_reverse){
      anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})
    }else{
      anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][length(strsplit(x,"\\.")[[1]])]})
    }
  }else{
    # check that provided anno_df is valid:
    if(length(setdiff(c("cluster_id","clusterlevel","cluster_name","first_cluster_name"),colnames(anno_df)))>0){
      stop("Wrong anno_df format. Required columns: cluster_id, clusterlevel, cluster_name, first_cluster_name")
    }
  }
  
  # if a heatmap matrix is provided, this function tries to infer the leaflevel based on the matrix
  # if(!is.null(heatmap_matrix) & is.null(leaf_level)){
  #   # TODO
  # }
  #
  # reduce edgelist to certain level and from and to cols
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level<=as.numeric(leaf_level),1:2]
  edgelist = edgelist[edgelist$to %in% anno_df$cluster_id,] # remove edges/nodes that are not part of anno_df
  
  ## convert to treedata
  # only take
  tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
  tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
  
  # add labels from annotation_df
  tree_data_tibble = dplyr::left_join(tree_data_tibble,anno_df,by=c("label"="cluster_id"))
  
  # update additional columns
  tree_data_tibble$first_cluster_name[is.na(tree_data_tibble$first_cluster_name)]=""
  tree_data_tibble$nodesize = 1 # default node size
  tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edgelist) # count children number
  tree_data_tibble$n_siblings = sapply(tree_data_tibble$label,function(x,el){ # count siblings
    parent = el$from[el$to==x]
    return(length(el$to[el$from==parent])-1)
  },el=edgelist)
  tree_data_tibble$tip.label =NA
  tree_data_tibble$tip.label[tree_data_tibble$n_children==0] = tree_data_tibble$node[tree_data_tibble$n_children==0] # add tip labels if leaf
  if(label_lowest){
    tree_data_tibble$first_cluster_name[ tree_data_tibble$n_children<2 & is.na(tree_data_tibble$tip.label)] = "" # if only one child node: changeto ""
  }else{
    tree_data_tibble$first_cluster_name[ tree_data_tibble$n_siblings == 0] = "" # if no siblings --> assume that labeled on level above
  }
  tree_data_tibble$branch_length = branch_length_base
  
  # convert back to treedata
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
  
  #plot circular tree
  circular_tree =ggtree(tree_data,layout = 'circular',branch.length = "branch_length",color=edge_color,right = ladderize_right )+ #  branch.length='none',
    geom_nodepoint(aes(subset = n_children > 1),color=node_color)#+geom_tippoint() + layout_circular()
  if(!is.null(order_column)){
    if(order_column %in% colnames(tree_data_tibble)){
      circular_tree = order_tree(tree_object = circular_tree,order_column = "order_plot",verbose = F)
    }else{
      stop("Cannot find order_columnin tree data. Please specifiy a valid column from anno_df or provide nULL to skip reordering") 
    }
  }
  # add genes on edges if necessary
  if(show_genes){
    circular_tree = circular_tree +
      geom_nodelab(aes(x=branch, label=first_cluster_name), size=label_size,vjust=vjust_label, nudge_x = nudge_x, nudge_y= nudge_y, color= label_color)
    if(show_genes_tip){
      circular_tree = circular_tree +
        geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name), size=label_size_tip,vjust=vjust_label, nudge_x = nudge_x_tip, nudge_y= nudge_y_tip,color= label_color)
      
    }
  }
  
  return(circular_tree)
  
}

##########
### plot_mini_tree
##########

#' Plot a smaller rectangular tree using edgelist and anno_df, similar to hypoMapUtils::plot_cluster_tree for the circular dendrogram plot

plot_mini_tree = function(edgelist, leaf_level = NULL, anno_df = NULL, metadata = NULL, layout = "rectangular",
                          pruned_edgelist = NULL, level_pattern = "K[0-9]+", cluster_id_pattern = "_pruned", 
                          cluster_name_pattern = "_named", label_size = 2, label_size_tip = label_size, 
                          show_genes = TRUE, vjust_label = -0.5,hjust_label=0, ignore_direct_parents = TRUE, 
                          na_color = "grey90", edge_color = "black", node_color = "black",node_size = 1,linesize = 1) {
  if (length(setdiff(c("from", "to", "level"), colnames(edgelist))) > 
      0) {
    warning("Error: Wrong edgelist format. Requires columns: from, to, level")
    return(NULL)
  }
  if (!leaf_level %in% edgelist$level) {
    warning("Error: leaf_level '", leaf_level, "' cannot be found in level column of edgelist")
    return(NULL)
  }
  if (length(setdiff(c("cluster_id", "clusterlevel", "cluster_name", 
                       "first_cluster_name"), colnames(anno_df))) > 0) {
    stop("Wrong anno_df format. Required columns: cluster_id, clusterlevel, cluster_name, first_cluster_name")
  }
  
  library(ggtree)
  library(treeio)
  library(igraph)
  edgelist_save = edgelist
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level <= as.numeric(leaf_level), 
                      1:2]
  edgelist = edgelist[edgelist$to %in% anno_df$cluster_id, 
  ]
  tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
  tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
  tree_data_tibble = dplyr::left_join(tree_data_tibble, anno_df, 
                                      by = c(label = "cluster_id"))
  tree_data_tibble$first_cluster_name[is.na(tree_data_tibble$first_cluster_name)] = ""
  tree_data_tibble$nodesize = 1
  tree_data_tibble$cluster_id = tree_data_tibble$label
  tree_data_tibble = dplyr::left_join(tree_data_tibble,edgelist_save[,c("level","to")],by=c("cluster_id"="to"))
  tree_data_tibble$n_children = sapply(tree_data_tibble$label, 
                                       function(x, el) {
                                         length(el$to[el$from == x])
                                       }, el = edgelist)
  tree_data_tibble$n_siblings = sapply(tree_data_tibble$label, 
                                       function(x, el) {
                                         parent = el$from[el$to == x]
                                         return(length(el$to[el$from == parent]) - 1)
                                       }, el = edgelist)
  tree_data_tibble$tip.label = NA
  tree_data_tibble$tip.label[tree_data_tibble$n_children ==   0] = tree_data_tibble$node[tree_data_tibble$n_children ==  0]
  if(ignore_direct_parents){
    tree_data_tibble$first_cluster_name[tree_data_tibble$n_children < 2 & is.na(tree_data_tibble$tip.label)] = ""
  }else{
    print("ab")
    tree_data_tibble$first_cluster_name[tree_data_tibble$n_children == 1 & tree_data_tibble$n_siblings == 0 ] = "" 
  }
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
  # tree
  circular_tree = ggtree(tree_data, layout = layout, branch.length = "level", color = edge_color,size=linesize) + 
    geom_nodepoint(aes(subset = n_children > 1), color = node_color, size = node_size) + 
    geom_tippoint(color = node_color, size = node_size)
  if (show_genes) {
    circular_tree = circular_tree + 
      geom_nodelab(aes(x = branch, label = first_cluster_name), size = label_size, vjust = vjust_label,hjust = hjust_label, color = "darkred") + 
      geom_tiplab(ggplot2::aes(x = branch,  label = first_cluster_name), align=F, size = label_size_tip,  vjust = vjust_label,hjust = hjust_label, color = "darkred")
  }
  return(circular_tree)
}

##########
### plot_mini_tree_ov
##########

#' Extended function to plot a smaller rectangular tree using edgelist and anno_df, similar to hypoMapUtils::plot_cluster_tree for the circular dendrogram plot

plot_mini_tree_ov = function(edgelist, leaf_level = NULL, anno_df = NULL, metadata = NULL, layout = "rectangular",ladderize =TRUE,
                             pruned_edgelist = NULL, level_pattern = "K[0-9]+", cluster_id_pattern = "_pruned", highlight_leaves = character(0),hilight_col = "red",highlight_bg=TRUE,gene_col = "darkred",
                             cluster_name_pattern = "_named", label_size = 2, label_size_tip = label_size, label_size_tip2 = label_size,
                             show_genes = TRUE, vjust_label = -0.5,hjust_label=0, ignore_direct_parents = TRUE, only_relevant_nodepoints =FALSE,reorder_tree=FALSE,order_column = "node",
                             na_color = "grey90", edge_color = "black", node_color = "black",node_size = 1,linesize = 1) {
  if (length(setdiff(c("from", "to", "level"), colnames(edgelist))) > 
      0) {
    warning("Error: Wrong edgelist format. Requires columns: from, to, level")
    return(NULL)
  }
  if (!leaf_level %in% edgelist$level) {
    warning("Error: leaf_level '", leaf_level, "' cannot be found in level column of edgelist")
    return(NULL)
  }
  if (length(setdiff(c("cluster_id", "clusterlevel", "cluster_name", 
                       "first_cluster_name"), colnames(anno_df))) > 0) {
    stop("Wrong anno_df format. Required columns: cluster_id, clusterlevel, cluster_name, first_cluster_name")
  }
  
  library(ggtree)
  library(treeio)
  library(igraph)
  # prepare edgelist via graph_from_edgelist and as.phylo
  edgelist_save = edgelist
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level <= as.numeric(leaf_level), 1:2]
  edgelist = edgelist[edgelist$to %in% anno_df$cluster_id, ]
  tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
  tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
  
  # join anno df
  tree_data_tibble = dplyr::left_join(tree_data_tibble, anno_df, by = c(label = "cluster_id"))
  
  # update & add various columns
  tree_data_tibble$first_cluster_name[is.na(tree_data_tibble$first_cluster_name)] = ""
  tree_data_tibble$nodesize = 1
  tree_data_tibble$cluster_id = tree_data_tibble$label
  tree_data_tibble = dplyr::left_join(tree_data_tibble,edgelist_save[,c("level","to")],by=c("cluster_id"="to"))
  tree_data_tibble$level[is.na(tree_data_tibble$level)] = min(tree_data_tibble$level,na.rm = TRUE)-1
  tree_data_tibble$cluster_id_show = ""
  tree_data_tibble$cluster_id_show[tree_data_tibble$level >= leaf_level] = tree_data_tibble$cluster_id[tree_data_tibble$level >= leaf_level]
  tree_data_tibble$n_children = sapply(tree_data_tibble$label, function(x, el) {
    length(el$to[el$from == x])}, el = edgelist)
  tree_data_tibble$n_siblings = sapply(tree_data_tibble$label, function(x, el) {
    parent = el$from[el$to == x]
    return(length(el$to[el$from == parent]) - 1)}, el = edgelist)
  tree_data_tibble$tip.label = NA
  tree_data_tibble$tip.label[tree_data_tibble$n_children ==   0] = tree_data_tibble$node[tree_data_tibble$n_children ==  0]
  if(ignore_direct_parents){
    tree_data_tibble$first_cluster_name[tree_data_tibble$n_children < 2 & is.na(tree_data_tibble$tip.label)] = ""
  }else{
    #tree_data_tibble$first_cluster_name[tree_data_tibble$n_children == 1 & tree_data_tibble$n_siblings == 0 ] = "" 
  }
  
  ## handle selected
  # add selection
  tree_data_tibble$selected_edges = "not"
  highlight_leaves = highlight_leaves[highlight_leaves %in% edgelist$to]
  if(length(highlight_leaves) > 0){
    if(any(highlight_leaves %in% tree_data_tibble$label)){
      all_selected_nodes = c(highlight_leaves,find_ancestors(highlight_leaves,edgelist[,c("from","to")]))
      tree_data_tibble$selected_edges[tree_data_tibble$label %in% all_selected_nodes] = "selected"
      tree_data_tibble$selected_edges[tree_data_tibble$level <= 3] = "not"
    }
    
  }
  #print(tree_data_tibble)
  # prepare as.treedata
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
  # make tree
  circular_tree = ggtree(tree_data, layout = layout, branch.length = "level", size=linesize,color=edge_color,ladderize = ladderize)+ # aes(color=selected_edges)
    # scale_color_manual(values = c("selected" = hilight_col,"not" = edge_color))+guides(color="none") + 
    # geom_nodepoint(color = node_color, size = node_size) + #aes(subset = n_children > 1) # aes(color=selected_edges)
    geom_tippoint(color = node_color, size = node_size)
  
  # depending on only_relevant_nodepoints show all node points or only those at brnaching or end points
  if(only_relevant_nodepoints){
    circular_tree = circular_tree +  geom_nodepoint(aes(subset = n_children > 1),color = node_color, size = node_size)
  }else{
    circular_tree = circular_tree + geom_nodepoint(color = node_color, size = node_size)
  }
  if (show_genes) {
    circular_tree = circular_tree + 
      geom_nodelab(aes(x = branch, label = first_cluster_name), size = label_size, vjust = vjust_label,hjust = hjust_label, color = gene_col) + 
      geom_tiplab(ggplot2::aes(x = branch,  label = first_cluster_name), align=F, size = label_size_tip,  vjust = vjust_label,hjust = hjust_label, color =gene_col)+
      geom_tiplab(ggplot2::aes(label = cluster_id_show), align=F, size = label_size_tip2,hjust = -0.1, color = "black")
  }
  if(reorder_tree){
    circular_tree = order_tree(circular_tree,order_column = order_column,verbose = FALSE) 
  }
  if(length(highlight_leaves) > 0){
    node_ids_highlight = tree_data_tibble$node[tree_data_tibble$cluster_id %in% highlight_leaves]
    circular_tree = circular_tree + geom_tiplab(ggplot2::aes(label = cluster_id_show,subset = label %in% highlight_leaves), align=F, size = label_size_tip2,hjust = -0.1, color = hilight_col)
    if(highlight_bg){
      circular_tree = circular_tree + geom_hilight(node=node_ids_highlight, fill= hilight_col, type="rect")
    }
  }
  return(circular_tree)
}


##########
### add_heatmap
##########

#' Plot metadata entries in cluster tree using ggtree as heatmaps.
#'
#' Requires an existing tree plot object created with plot_cluster_tree
#'
#' TODO: not tested in package yet!
#'
#' @param circular_tree a tree object
#' @param heatmap_matrix a matrix with rownames corresponding to tip node names
#' @param heatmap_colnames olors to be passed to scale_fill_gradientn
#' @param legend_title "Legend"
#' @param matrix_offset 0.2
#' @param matrix_width 0.2
#' @param colnames_angle 0
#' @param legend_text_size 4
#' @param hjust_colnames colors to be passed to scale_fill_gradientn
#' @param heatmap_colors ...
#' @param scale_limits ...
#' @param heatmap_text_size ...
#' @param na_color ...
#'
#' @return ggtree object or plot
#'
#' @export
#'
#' @import dplyr ggplot2 scales treeio tidytree ggtree ggnewscale
#'

add_heatmap = function(circular_tree, heatmap_matrix, heatmap_colors = c("white", "darkred"), scale_limits = NULL, heatmap_colnames = TRUE, 
                       legend_title = "Legend1", matrix_offset = 0.2, matrix_width = 0.2, 
                       colnames_angle = 0, legend_text_size = 4, hjust_colnames = 0.5, colnames_offset_y = 0,
                       heatmap_text_size = 4, na_color = "grey90") {
  circular_tree_heat <- circular_tree
  circular_tree_heat <- circular_tree_heat + ggnewscale::new_scale_fill()
  if (is.numeric(heatmap_matrix[, 1])) {
    if (is.null(scale_limits)) {
      scale_limits = c(min(heatmap_matrix), max(heatmap_matrix))
    }
    circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix, 
                                   offset = matrix_offset, width = matrix_width, colnames = heatmap_colnames, 
                                   colnames_angle = colnames_angle, colnames_offset_y = colnames_offset_y, 
                                   font.size = heatmap_text_size, hjust = hjust_colnames) + 
      scale_fill_gradientn(colours = heatmap_colors, limits = scale_limits, 
                           oob = squish, na.value = na_color,name=legend_title) + theme(legend.text = element_text(size = legend_text_size))
  }
  else {
    circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix, 
                                   offset = matrix_offset, width = matrix_width, colnames = heatmap_colnames, 
                                   colnames_angle = colnames_angle, colnames_offset_y = colnames_offset_y, 
                                   font.size = heatmap_text_size, hjust = hjust_colnames) + 
      scale_fill_discrete(na.value = na_color,name=legend_title) + theme(legend.text = element_text(size = legend_text_size))
  }
  return(circular_tree_heat)
}

##########
### order_tree
##########

#' Order Tree Nodes
#'
#' This function orders the nodes of a tree based on their vertical position (y-coordinate).
#'
#' @param tree_object A tree object, typically generated using the ggtree package.
#' @param parent The parent node from which to start ordering. Default is 0 (root node). Don't change this !!!
#' @param order_column which column in tree$data to use for bubblesort. Must be present in data and must be comparable by > . Defaults to the node index "node", but can be another numeric or a alpha-numerically comparable string. (e.g. cluster_ids).
#' @param depth_message > for verbose during recursion. Don't change this !!!
#' @param verbose whether to give updates
#'
#' @return A tree object with ordered nodes.
#'
#' @details
#' This function takes a tree object and recursively orders its nodes based on their
#' vertical position (y-coordinate). The ordering is performed using a modified
#' bubble sort algorithm. The function starts from the specified parent node and
#' recursively orders its children.

order_tree = function(tree_object,parent = 0,order_column="node",depth_message="",verbose=FALSE){
  # get tree data
  tree_data = tree_object$data %>% as.data.frame()
  # set parent of root node to 0
  tree_data$parent[tree_data$parent == tree_data$node] = 0
  # get children
  children = tree_data$node[tree_data$parent==parent]
  #message("found ",paste0(children,collapse = ", "))
  # if there is more than one child
  if(length(children) >= 2){
    # 
    if(verbose) message(depth_message,"reordering children of ",tree_data$cluster_name[tree_data$node==parent],":  ",paste0(children,collapse = ", "))
    ## get the y order of children
    children_order = tree_data$y[tree_data$node %in% children]
    names(children_order) = tree_data[tree_data$node %in% children,order_column]
    # this is then the current order by y
    children_to_sort = names(sort(children_order))# as.numeric(names(sort(children_order)))
    if(verbose) message("children to sort: ",paste0(children_to_sort,collapse = ", "))
    ### order children using bubblesort back into a normal numeric order
    # get n
    n = length(children_to_sort)
    swapped = TRUE
    # repeat
    while(swapped){
      swapped = FALSE
      #for i := 1 to n-1 inclusive do
      for(i in 2:n){
        # if this pair is out of order
        if(children_to_sort[i-1] < children_to_sort[i]){ # I changed the direction here from > to < because then it flips in the top-down direction that I want
          # swap them and remember something changed 
          save = children_to_sort[i-1]
          children_to_sort[i-1] = children_to_sort[i]
          children_to_sort[i] = save
          # flip the nodes in the tree object using the bubblesort swap
          if(verbose) message("Flipping: ",children_to_sort[i-1]," and ",children_to_sort[i])
          node_1 = tree_data$node[tree_data[,order_column] == children_to_sort[i-1]] %>% na.omit()# as.numeric(children_to_sort[i-1])
          node_2 = tree_data$node[tree_data[,order_column] == children_to_sort[i]]  %>% na.omit()# as.numeric(children_to_sort[i])
          if(verbose) message("node_1: ",node_1," ;node_2: ",node_2)
          tree_object= ggtree::flip(tree_view = tree_object,node1 = node_1,node2 = node_2)
          # if(verbose) message("swapped!")
          swapped = TRUE
        }
      }
    }
    # then recursively start for all children
    for(child in children){
      depth_message = paste0(depth_message,">")
      tree_object = order_tree(tree_object,parent = child,depth_message = depth_message,verbose=verbose,order_column=order_column) 
    }
  }else if(length(children) == 1){
    if(verbose) message(depth_message," recursively start with ",children)
    # else recursively start with child as new parent
    depth_message = paste0(depth_message,">")
    tree_object = order_tree(tree_object,parent = children,depth_message = depth_message,verbose=verbose,order_column=order_column) 
  }else{
    # do nothing
  }
  return(tree_object) 
}

##########
### plot_tree_comparison
##########

# plot_tree_comparison

plot_tree_comparison = function(matched_clusters,edgelist_human, edgelist_mouse , leaf_level_tree=6,remove_name_part ="",color_order_regex = "C4",start_level = "C2",include_extra_level_human=F,general_label_size = 7,min_similarity=0.7, tree_color="grey40",gene_col="black",node_to_ignore=c()){
  
  ## prepare anno df
  # make anno_df
  anno_df = edgelist_human %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
  anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
  anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})
  
  
  anno_df_mouse = edgelist_mouse %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
  anno_df_mouse$clusterlevel = stringr::str_extract(anno_df_mouse$cluster_id,pattern = "C[0-9]+")
  anno_df_mouse$first_cluster_name = sapply(gsub("C[0-9]+-[0-9]+: ","",anno_df_mouse$cluster_name),function(x){head(strsplit(x,"\\.")[[1]],n=1)})
  
  ### get nodes
  #edgelist_human_filt = c(unique(edgelist_human$from[grepl(start_level,edgelist_human$from)]),find_children(unique(edgelist_human$from[grepl(start_level,edgelist_human$from)]),edgelist_human))
  #edgelist_human_filt = edgelist_human[edgelist_human$from %in% edgelist_human_filt | edgelist_human$isLeaf,]
  primary_from_nodes = edgelist_human$to[edgelist_human$to_named %in% matched_clusters$from]
  all_from_nodes = c(primary_from_nodes,unique(find_ancestors(primary_from_nodes,edgelist_human)))
  
  # find good root
  keep_levels = c(stringr::str_extract(primary_from_nodes,"C[0-9]+"),names(sort(table(stringr::str_extract(all_from_nodes,"C[0-9]+")),decreasing = TRUE))[sort(table(stringr::str_extract(all_from_nodes,"C[0-9]+")),decreasing = TRUE) > 1])
  if(include_extra_level_human){keep_levels = unique(c(keep_levels,stringr::str_extract(unique(find_ancestors("C4-375",combined_edgelist_mrtree)),"C[0-9]+")[1]))}
  edgelist_human$level_name = stringr::str_extract(edgelist_human$to,"C[0-9]+")
  edgelist_human_filt = edgelist_human[edgelist_human$level_name %in% keep_levels,]
  relevantNodes = c(primary_from_nodes,unique(find_ancestors(nodes = primary_from_nodes,edges = edgelist_human_filt)))
  ignore_nodes = c(node_to_ignore,scUtils::find_children(nodes = node_to_ignore,edges = edgelist_human))
  relevantNodes = relevantNodes[!relevantNodes %in% ignore_nodes]
  ### prepare annotation
  nodes_for_umap_annotation = data.frame(cluster_id = relevantNodes)
  nodes_for_umap_annotation = dplyr::left_join(nodes_for_umap_annotation,anno_df[,c("cluster_id","cluster_name","first_cluster_name"),drop=FALSE])
  
  ## prepare colors and order
  color_order = nodes_for_umap_annotation[,c("cluster_id","first_cluster_name","cluster_name"),drop=FALSE]
  color_order = color_order[grepl(color_order_regex,color_order$cluster_id),]
  color_order$color = getOkabeItoPalette(nrow(color_order))
  color_order$cluster_name = factor(color_order$cluster_name,levels = color_order$cluster_name )# make.unique
  color_order$final_name = gsub(remove_name_part,"",color_order$cluster_name)
  
  ### make tree
  message("tree_plot")
  message("relevantNodes: ",paste0(relevantNodes,collapse = " | "))
  subset_edge_list = combined_edgelist_mrtree[combined_edgelist_mrtree$from %in% c(relevantNodes),]
  mini_tree = plot_mini_tree_ov(edgelist = subset_edge_list,
                                leaf_level=leaf_level_tree,
                                #layout="dendrogram",
                                anno_df = anno_df ,
                                metadata=human_hypo_combined@meta.data,
                                label_size = general_label_size, 
                                label_size_tip = general_label_size,
                                gene_col = gene_col,
                                show_genes = TRUE,
                                vjust_label = -0.5,
                                hjust_label = 0,
                                node_size = 5,
                                linesize = 1.5,
                                edge_color = tree_color, 
                                node_color = tree_color,
                                only_relevant_nodepoints =TRUE)
  ###
  mini_tree
  
  ############### Mouse
  
  node_to_use_mouse = matched_clusters$to
  node_to_use_mouse_id = stringr::str_split_i(node_to_use_mouse,pattern = ": ",i = 1)
  relevantNodes_mouse = unique(find_ancestors(nodes = node_to_use_mouse_id,edges = edgelist_mouse))
  mouse_levels_names = stringr::str_extract(relevantNodes_mouse,"C[0-9]+")
  mouse_levels = as.numeric(gsub("C","",mouse_levels_names))
  names(mouse_levels) = mouse_levels_names
  mouse_levels= sort(mouse_levels)
  # print(table(mouse_levels))
  extra_level_to_keep = paste0("C",names(table(mouse_levels))[which(table(mouse_levels) == 1)[length(which(table(mouse_levels) == 1))]])
  other_levels_to_keep = paste0("C",names(table(mouse_levels))[which(table(mouse_levels) > 1)])
  #print(other_levels_to_keep)
  #print(extra_level_to_keep)
  keep_levels = unique(c(stringr::str_extract(node_to_use_mouse_id,"C[0-9]+"),other_levels_to_keep,extra_level_to_keep))
  keep_levels = keep_levels[keep_levels != "C" ]
  edgelist_mouse$level_name = stringr::str_extract(edgelist_mouse$to,"C[0-9]+")
  #print(keep_levels)
  edgelist_mouse_filt = edgelist_mouse[edgelist_mouse$level_name %in% keep_levels,]
  #print(edgelist_mouse_filt)
  relevantNodes_mouse = unique(find_ancestors(nodes = node_to_use_mouse_id,edges = edgelist_mouse_filt))
  ### make tree mouse
  
  tree_color="grey40"
  
  message("tree_plot")
  message("relevantNodesMouse: ",paste0(relevantNodes_mouse,collapse = " | "))
  subset_edge_list_mouse = hypoMap_edgelist[hypoMap_edgelist$from %in% c(relevantNodes_mouse),,drop=FALSE]
  #print(subset_edge_list_mouse)
  # add root if necessary ?
  #root_df = data.frame(from = "rest",to =unique(subset_edge_list_mouse$from[subset_edge_list_mouse$level == min(subset_edge_list_mouse$level)]))
  #root_df$level = min(subset_edge_list_mouse$level) - 1
  #root_df$height = max(subset_edge_list_mouse$height) + 1
  #subset_edge_list_mouse = bind_rows(root_df,subset_edge_list_mouse) %>% as.data.frame()
  
  mini_tree_mouse = plot_mini_tree_ov(edgelist = subset_edge_list_mouse,
                                      leaf_level=max(subset_edge_list_mouse$level),
                                      #layout="dendrogram",
                                      anno_df = anno_df_mouse ,
                                      #   metadata=hypoMap@meta.data,
                                      label_size = general_label_size, 
                                      label_size_tip = general_label_size,
                                      gene_col=gene_col,
                                      show_genes = TRUE,
                                      vjust_label = -0.25,
                                      hjust_label = -1,
                                      node_size = 5,
                                      linesize = 1.5,
                                      edge_color = tree_color,
                                      node_color = tree_color,
                                      only_relevant_nodepoints =TRUE)
  ###
  mini_tree_mouse
  
  ## make a flip (can only do this manually after seeing final result !):
  #mini_tree = flip(mini_tree, 8, 10)
  #mini_tree = flip(mini_tree, 8, 9)
  #mini_tree_mouse = flip(mini_tree_mouse, 12, 13)
  
  data_human <- mini_tree$data
  data_mouse <- mini_tree_mouse$data
  ## reverse x-axis and 
  ## set offset to make the tree on the right-hand side of the first tree
  x_offset = 10
  data_mouse$x <- max(data_mouse$x) - data_mouse$x + max(data_human$x) + x_offset 
  #data_mouse$branch_rev <- max(data_mouse$branch) - data_mouse$branch + max(data_human$x) + 1
  data_mouse$branch_rev <-  max(data_mouse$x) - data_mouse$branch 
  
  # manipulate y of smaller tree
  #data_human$y = (mini_tree$data$y* max(data_mouse$y) / max(mini_tree$data$y) ) - (max(data_mouse$y) / max(mini_tree$data$y)) + 1
  mini_tree$data$y = (mini_tree$data$y* max(data_mouse$y) / max(mini_tree$data$y) ) - (max(data_mouse$y) / max(mini_tree$data$y)) + 1
  
  ### plot new one
  ## NEED TO ADD VIA geom_tree not second ggtree !!
  linesize = 1.5
  pp <- mini_tree +
    geom_tree(data = data_mouse,layout = "rectangular",size = linesize,color=tree_color)  +
    ggplot2::xlim(0, max(data_mouse$x)+1)  + ggplot2::ylim(0, max(data_mouse$y)+1)
  pp
  
  # params for annos
  label_size = general_label_size
  label_size_tip = general_label_size
  vjust_label = -0.4
  hjust_label = 0
  hjust_tiplabel = 0
  static_tiplabel_nudge = 4
  node_size = 3.5
  
  
  #pp + geom_nodepoint(data = data_mouse,color = tree_color, size = node_size)+
  #  geom_nodelab(data = data_mouse,aes(x = branch_rev, label = first_cluster_name), size = label_size, vjust = vjust_label,hjust = hjust_label, color = "darkred")
  
  ## add annos again
  pp2 = pp + 
    # geom_nodepoint(data = data_mouse,color = tree_color, size = node_size)+
    geom_nodepoint(data = data_mouse,aes(subset = n_children > 1),color = tree_color,  size = node_size)+ # only_relevant_nodepoints
    geom_tippoint(data = data_mouse,color = tree_color, size = node_size)+
    geom_nodelab(data = data_mouse,aes(x = branch_rev, label = first_cluster_name), size = label_size, vjust = vjust_label,hjust = hjust_label, color =  gene_col,) + 
    geom_tiplab(data = data_mouse,ggplot2::aes(x = branch_rev,  label = first_cluster_name), align=F, size = label_size_tip,  vjust = vjust_label,hjust = hjust_label, color = gene_col)+
    geom_tiplab(data = data_mouse,ggplot2::aes(x = x- static_tiplabel_nudge,label = cluster_id_show), align=F, size = label_size_tip,hjust = hjust_tiplabel, color = "black")
  
  ## make new plot
  pp2
  
  ## prepare edges
  if(! "edge_id" %in% colnames(matched_clusters) ){
    matched_clusters$edge_id = paste0(matched_clusters$from,"_",matched_clusters$to)
  }
  ## add edge to data
  data_human_new =  mini_tree$data
  data_human_edges = dplyr::left_join(data_human_new,matched_clusters %>%ungroup() %>% dplyr::select(from = from,edge_id, alt_label = to,similarity),by=c("cluster_name"="from"),multiple = "all")
  data_human_edges = data_human_edges[!is.na(data_human_edges$edge_id),]
  data_human_edges$species = "human"
  
  if("label" %in% colnames(data_mouse)){
    #data_mouse = data_mouse %>% dplyr::select(-label) %>% as.data.frame() # fixing a bug where to lebel columns are in there after join
    if(!is.character(data_mouse$label)){
      data_mouse$label = as.character(data_mouse$label )
    }
  }
  data_mouse_edges = dplyr::left_join(data_mouse,matched_clusters %>%ungroup()  %>% dplyr::select(to,edge_id, alt_label = from,similarity) ,by=c("cluster_name"="to"),multiple = "all")
  data_mouse_edges = data_mouse_edges[!is.na(data_mouse_edges$edge_id),]
  data_mouse_edges$species = "mouse"
  
  # print(data_mouse_edges)
  # print(data_human_edges)
  setdiff_edges = setdiff(data_mouse_edges$edge_id,data_human_edges$edge_id)
  dd <- bind_rows(data_human_edges, data_mouse_edges) %>% 
    filter(!edge_id %in% setdiff_edges)
  #dd$x[dd$species=="human"] = dd$x[dd$species=="human"] - 21
  
  pp3 = pp2 + geom_line(aes(x, y, group=edge_id,color=similarity), data=dd,linewidth=2) +
    scale_color_gradient2(low = "orange",mid = "yellow",high = "darkgreen",na.value = "white",midpoint = 0.8,limits=c(min_similarity,1),name = "Similarity")
  pp3
  
  
  ## fancy bg fill
  alphaval = 0.2
  human_root = data_human_new$node[data_human_new$branch.length==0]
  print(human_root)
  mouse_root = data_mouse$node[data_mouse$branch.length==0]
  print(mouse_root)
  pp4 <- pp3 +
    ggtree::geom_hilight(
      data = data_human_new,#[data_human_new$node == human_root,],
      # mapping = aes(
      #   #subset = node == human_root,
      #   node=node
      # ),
      node = human_root,
      fill = "#E69F00",
      align = "none",
      alpha = alphaval,
      extendto = 2
    )+
    ggtree::geom_hilight(
      data = data_mouse,
      # mapping = aes(
      #   #subset = node == mouse_root,
      #   node=node
      # ),
      node = mouse_root,
      fill = "#009E73",
      align = "none",
      alpha = alphaval,
      extendto = 2
    )
  # move layers to background
  total_layers = length(pp4$layers)
  pp4$layers = c(pp4$layers[(total_layers-2):total_layers],pp4$layers[1:(total_layers-3)])
  
  # pp4 <- pp3
  # # move layers to background
  # total_layers = length(pp4$layers)
  # pp4$layers = c(pp4$layers[(total_layers-0):total_layers],pp4$layers[1:(total_layers-1)])
  # 
  #plot
  return(pp4)
  
}

##########
### find_children
##########

## add scUtils::find_children here to avoid dependency

find_children = function(nodes, edges){
  current_children = edges$to[edges$from %in% nodes]
  if (length(current_children) > 0) {
    all_children = c(current_children, find_children(current_children, 
                                                     edges))
  }
  else {
    all_children = current_children
  }
  return(all_children)
}

# find find_ancestors in tree recusrively based on simple edgelist
find_ancestors = function(nodes,edges){
  current_ancestors = edges$from[edges$to %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_ancestors)>0){
    all_ancestors = c(current_ancestors,find_ancestors(current_ancestors,edges))
  }else{
    all_ancestors = current_ancestors
  }
  return(all_ancestors)
}

