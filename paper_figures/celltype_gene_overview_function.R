


##########
### Function that prepares the cell type overview data
##########

#region_mapping = data.table::fread("data/region_anno_revision/mad_assignments_C4.tsv",data.table = F)

celltype_overview_data = function(
    # input clusters & genes
  main_gene ,
  other_genes,
  expr_type = "avg.exp",
  # input data sets
  edgelist = matched_clusters_final,
  human_hypo_combined,
  hypoMap,
  region_mapping = region_mapping,
  human_marker_genes,
  col_name = "C4_named",
  mouse_col_name = "C465_named",
  q_prob = 0.98, # for filtering 
  filter_alternative = TRUE,
  ignore_mouse_expression =FALSE,
  min_expr_value = 0.1
){
  
  require(Seurat)
  # prepare edgelist to only join top cluster
  edgelist_sort = edgelist %>% dplyr::arrange(desc(similarity)) %>%
    dplyr::group_by(from) %>%
    dplyr::slice_max(order_by = similarity,n = 1,with_ties = F)
  
  
  ##########
  ### Basic stats
  ##########
  
  message("get expression")
  ## human expression
  Idents(human_hypo_combined) = col_name
  dotplot= Seurat::DotPlot(human_hypo_combined,features = unique(c(main_gene,other_genes)),scale = FALSE) # with scale FALS avg.exp.scaled contains the normalized counts !!
  dotplot_data = dotplot$data
  
  ## mouse expression
  Idents(hypoMap) = mouse_col_name
  dotplot_mouse= Seurat::DotPlot(hypoMap,features = unique(stringr::str_to_title(c(main_gene,other_genes))),scale = FALSE) # with scale FALS avg.exp.scaled contains the normalized counts !!
  dotplot_mouse_data = dotplot_mouse$data
  
  ## make pct dfs human
  if(expr_type == "pct"){
    human_gene_pct = dotplot_data %>% dplyr::select(group = id, pct = pct.exp, gene= features.plot )
    #human_gene_pct = dplyr::left_join(human_gene_pct,annotated_cluster_overview,by=c("group"="cluster"))
    human_gene_pct_wide = human_gene_pct %>% tidyr::spread(key = "gene",value="pct")
  }else if(expr_type == "avg.exp"){
    human_gene_pct = dotplot_data %>% dplyr::select(group = id, pct = avg.exp.scaled, gene= features.plot )
    #human_gene_pct = dplyr::left_join(human_gene_pct,annotated_cluster_overview,by=c("group"="cluster"))
    human_gene_pct_wide = human_gene_pct %>% tidyr::spread(key = "gene",value="pct")
  }else{
    stop("This expr_type is not implemented!") 
  }
  
  ## make pct dfs mouse
  if(expr_type == "pct"){
    mouse_gene_pct = dotplot_mouse_data %>% dplyr::select(group = id, pct = pct.exp, gene= features.plot )
    mouse_gene_pct_wide = mouse_gene_pct %>% tidyr::spread(key = "gene",value="pct")
  }else if(expr_type == "avg.exp"){
    mouse_gene_pct = dotplot_mouse_data %>% dplyr::select(group = id, pct = avg.exp.scaled, gene= features.plot )
    #  mouse_gene_pct = dplyr::left_join(mouse_gene_pct,annotated_cluster_overview,by=c("group"="cluster"))
    mouse_gene_pct_wide = mouse_gene_pct %>% tidyr::spread(key = "gene",value="pct")
  }else{
    stop("This expr_type is not implemented!") 
  }
  
  ## IMPORTANT: pct is just the name ---> can also be avg.exp !!!!!
  
  ##########
  ### Get gene markers
  ##########
  
  ##### make markers
  message("get marker genes")
  markers_per_cluster = human_marker_genes %>% 
    dplyr::filter(p_val_adj < 0.0001 &  specificity > 3 & !grepl("AC|AL|MIR|-AS|LINC",gene)) %>% 
    dplyr::select(gene,name,specificity) %>%
    #  dplyr::filter(cluster %in% glp1r_pct_per_cluster_top$group ) %>%
    #  dplyr::filter(grepl(stringr::str_extract(col_name,"C[0-9]+"),cluster)) %>%
    dplyr::group_by(name) %>% dplyr::slice_max(order_by = specificity,n = 5) %>%
    mutate(top_markers = paste0(gene, collapse = "|")) %>%
    dplyr::distinct(name,top_markers)
  
  
  ##########
  ### Make combined dataframe
  ##########
  
  message("combine")
  
  gene_pct_anno = dplyr::left_join(human_gene_pct_wide,markers_per_cluster,by=c("group"="name") ) %>% 
    #### HERE avg_abundances_c6_long_stats_summarized
    dplyr::left_join(region_mapping ,by=c("group"="cluster_name")) %>%
    dplyr::left_join(edgelist_sort %>% 
                       dplyr::select(human_cluster = from,mouse_cluster = to,species_cor = similarity,from_occ) %>%
                       dplyr::mutate(species_cor = round(species_cor,3)),
                     by=c("group"="human_cluster") ) %>%
    dplyr::arrange(desc(!!sym(main_gene[1]))) %>%
    dplyr::left_join(mouse_gene_pct_wide,by=c("mouse_cluster"="group")) #%>%
  
  cols_to_round = c(c(main_gene,other_genes),stringr::str_to_title(c(main_gene,other_genes)))
  cols_to_round = cols_to_round[cols_to_round %in% colnames(gene_pct_anno)]
  if(expr_type == "pct"){
    gene_pct_anno[,cols_to_round] = apply(gene_pct_anno[,cols_to_round],2,function(x){round(x,digits = 2)})
  }else{
    gene_pct_anno[,cols_to_round] = apply(gene_pct_anno[,cols_to_round],2,function(x){round(x,digits = 2)})
  }
  colnames(gene_pct_anno)[1] = c("human_cluster")
  
  
  ##########
  ### FILTER
  ##########
  
  message("filter")
  if(ignore_mouse_expression){
    message("ignore_mouse_expression")
    if(length(main_gene) == 1){
      gene_pct_anno_filter = gene_pct_anno[gene_pct_anno[,main_gene] >= max(min_expr_value,quantile(human_gene_pct_wide[,main_gene],probs = q_prob,na.rm = TRUE)),]
      # a bit more complicated when using multiple main genes: need to iterate to find all relevant to rows to filter to  
    }else if(length(main_gene) > 1){
      message("Filtering for multiple main genes !")
      allidx = c()
      for(i in 1:length(main_gene)){
        current_main_gene_human = main_gene[i]
        idx = which(gene_pct_anno[,current_main_gene_human] >= max(min_expr_value,quantile(human_gene_pct_wide[,current_main_gene_human],probs = q_prob,na.rm = TRUE)))
        allidx = c(allidx,idx)
      }
      gene_pct_anno_filter = gene_pct_anno[unique(allidx),]
    }else{
      # some error catching
      stop("Please provide at least one valid main_gene") 
    }
  }else{
    if(!filter_alternative){
      
      # conversion simple by to title:
      main_gene_mouse = stringr::str_to_title(main_gene)
      # simple when there is one maingene:
      if(length(main_gene) == 1){
        gene_pct_anno_filter = gene_pct_anno[gene_pct_anno[,main_gene] >= max(min_expr_value,quantile(human_gene_pct_wide[,main_gene],probs = q_prob,na.rm = TRUE)) |
                                               gene_pct_anno[,main_gene_mouse] >= max(min_expr_value,quantile(mouse_gene_pct_wide[,main_gene_mouse],probs = q_prob,na.rm = TRUE)) ,]
        # a bit more complicated when using multiple main genes: need to iterate to find all relevant to rows to filter to  
      }else if(length(main_gene) > 1){
        message("Filtering for multiple main genes !")
        allidx = c()
        for(i in 1:length(main_gene)){
          current_main_gene_human = main_gene[i]
          current_main_gene_mouse = main_gene_mouse[i]
          idx = which(gene_pct_anno[,current_main_gene_human] >= max(min_expr_value,quantile(human_gene_pct_wide[,current_main_gene_human],probs = q_prob,na.rm = TRUE)) |
                        gene_pct_anno[,current_main_gene_mouse] >= max(min_expr_value,quantile(mouse_gene_pct_wide[,current_main_gene_mouse],probs = q_prob,na.rm = TRUE)))
          allidx = c(allidx,idx)
        }
        gene_pct_anno_filter = gene_pct_anno[unique(allidx),]
      }else{
        # some error catching
        stop("Please provide at least one valid main_gene") 
      }
      # filter alternative uses dotplot data
    }else if(filter_alternative & expr_type == "avg.exp"){
      # conversion simple by to title:
      main_gene_mouse = stringr::str_to_title(main_gene)
      
      # uses dotplot data of pct to select filter cutoff
      human_cutoffs = dotplot_data %>% dplyr::group_by(features.plot) %>% dplyr::summarise(cutoff = quantile(pct.exp,q_prob)) %>% dplyr::filter(features.plot %in% main_gene)
      
      # get list of clusters matching criteria from human
      dotplot_data2 = dplyr::left_join(dotplot_data,human_cutoffs) %>% dplyr::filter(pct.exp >= cutoff)
      human_clusters = unique(as.character(dotplot_data2$id))
      
      # for mouse:
      mouse_cutoffs = dotplot_mouse_data %>% dplyr::group_by(features.plot) %>% dplyr::summarise(cutoff = quantile(pct.exp,q_prob)) %>% dplyr::filter(features.plot %in% main_gene_mouse)
      
      # get list of clusters matching criteria from mouse
      dotplot_mouse_data2 = dplyr::left_join(dotplot_mouse_data,mouse_cutoffs) %>% dplyr::filter(pct.exp >= cutoff)
      mouse_clusters = unique(as.character(dotplot_mouse_data2$id))
      
      ## filter main df
      gene_pct_anno_filter = gene_pct_anno[gene_pct_anno$human_cluster %in% human_clusters | gene_pct_anno$mouse_cluster %in% mouse_clusters,]
      
    }else{
      stop("Select other filtering combination") 
    }
  }
  return(list(gene_pct_anno_filter = gene_pct_anno_filter, gene_pct_anno = gene_pct_anno))
}


##########
### Function that makes the cell typer overview tree plot
##########


celltype_overview_plot = function(
    # target
  # input clusters & genes
  main_gene ,
  other_genes,
  highlight_leaves=c(),
  # input
  gene_pct_anno_filter,
  combined_edgelist_mrtree,#
  anno_df,
  # main plotting params
  col_name = "C4_named",
  show_intermediate = FALSE,
  addsuffix="_resub",
  root_node_use = "C0-3",
  second_scale = TRUE, # receptos with their own scale
  # other plotting params
  one_gene_width = 0.06, # per gene width for heatmaps
  offset_elements = 1, # additional offset to space elemnts of plot
  heatmap_text_size = 4.5, # gene names text size
  legend_text_size = 6, # ???
  remove_irrelevevant_c3 = TRUE, ## whether to remove any non relevant cluster on level C4/C3
  tree_label_size = 3, # size of gene labels on tree
  tree_tip_label_size = 3.5, # size of cluster id label 
  expr_type = "avg.exp", # or pct
  min_expr_value = 0.05, # 0r 2.5% --> set a minimum of what is considered expressed for color cutoff to full grey and for filtering cell types !!
  tree_color="grey40",
  text_color_edges = "black",
  nbreaks_scale = 5,
  leaf_level = 6,
  label_textanno_size = 4,
  mouse_label_offset = 28
){
  
  require(scales)
  # join cluster_id 
  gene_pct_anno_filter = gene_pct_anno_filter  %>% dplyr::left_join(anno_df[,c("cluster_id","cluster_name")],by=c("human_cluster"="cluster_name"))
  
  ##########
  ### Make a mini tree
  ##########
  
  use_edgelist_mrtree = combined_edgelist_mrtree #%>%
  #dplyr::rename(from_id = from,to_id = to) %>%
  #dplyr::rename(from=from_named,to=to_named)
  
  all_root_children = find_children(nodes = root_node_use,edges = use_edgelist_mrtree)
  use_edgelist_mrtree = use_edgelist_mrtree[use_edgelist_mrtree$to %in% all_root_children,]
  
  # define which nodes to use -- only down to C4 at the moment
  relevantNodes = c(root_node_use,find_children(nodes = root_node_use,edges = use_edgelist_mrtree))
  all_c2_nodes =  relevantNodes[grepl("C1",relevantNodes)]
  relevantNodes = relevantNodes[grepl("C0|C1|C2",relevantNodes)] # "C0|C1|C2|C3|C4"
  
  # get top clusters and their ancestors and additionally include
  top_gene_clusters = anno_df$cluster_id[anno_df$cluster_name %in% gene_pct_anno_filter$human_cluster] #[gene_pct_anno[,main_gene] > 30]
  top_gene_clusters = top_gene_clusters[top_gene_clusters %in% find_children(nodes = root_node_use,edges = use_edgelist_mrtree)] # only those that are in neurons
  top_gene_clusters_ancestors = find_ancestors(nodes = top_gene_clusters,edges = use_edgelist_mrtree)
  relevantNodes = unique(c(relevantNodes,top_gene_clusters_ancestors,top_gene_clusters))
  
  ## delele C4 nodes of C2 nodes that are not in ancestors # changed to C3 of C1
  if(remove_irrelevevant_c3){
    not_used_C2_nodes = all_c2_nodes[!all_c2_nodes %in% top_gene_clusters_ancestors[grepl("C1",top_gene_clusters_ancestors)]]
    not_used_C2_nodes_C4_children = find_children(nodes = not_used_C2_nodes,edges = use_edgelist_mrtree)
    not_used_C2_nodes_C4_children = not_used_C2_nodes_C4_children[grepl("C3",not_used_C2_nodes_C4_children)]
    relevantNodes = relevantNodes[! relevantNodes %in% not_used_C2_nodes_C4_children]
  }
  
  ### make tree - rectangle
  #subset_edge_list = use_edgelist_mrtree[use_edgelist_mrtree$from %in% c(relevantNodes),]
  subset_edge_list = use_edgelist_mrtree[use_edgelist_mrtree$from %in% c(relevantNodes) & use_edgelist_mrtree$to %in% c(relevantNodes),]
  anno_df$first_cluster_name[anno_df$first_cluster_name == "TBD"] = ""
  
  # change order in list: --> only works when node is used
  subset_edge_list_save = subset_edge_list
  subset_edge_list = bind_rows(subset_edge_list_save[1:2,],subset_edge_list_save[5,]) %>% 
    bind_rows(subset_edge_list_save[3:4,]) %>%
    bind_rows(subset_edge_list_save[6:nrow(subset_edge_list_save),]) 
  
  ## make modified anno_df to enforce order via cluster_id
  anno_df_mod = anno_df
  #anno_df_mod$cluster_order = seq_along(anno_df_mod$cluster_id)
  anno_df_mod$cluster_order = 1:nrow(anno_df_mod)
  # anno_df_mod$cluster_order[anno_df_mod$cluster_id=="C1-5"] = 7
  # anno_df_mod$cluster_order[anno_df_mod$cluster_id=="C1-3"] = 8
  # anno_df_mod$cluster_order[anno_df_mod$cluster_id=="C1-4"] = 9
  # anno_df_mod$cluster_order[anno_df_mod$cluster_id=="C2-8"] = 31
  # anno_df_mod$cluster_order[anno_df_mod$cluster_id=="C2-7"] = 32
  
  
  
  set.seed(42)
  celltype_tree_plot = plot_mini_tree_ov(edgelist = subset_edge_list,#subset_edge_list[1:10,],
                                         leaf_level=leaf_level,#length(unique(stringr::str_extract(relevantNodes,"C[0-9]+")))+1+as.numeric(stringr::str_remove(stringr::str_extract(node_to_use,"C[0-9]+"),"C")),
                                         #layout="dendrogram",
                                         anno_df = anno_df_mod ,
                                         metadata=human_hypo_combined@meta.data,
                                         label_size = tree_label_size, 
                                         label_size_tip = tree_label_size,
                                         label_size_tip2 = tree_tip_label_size,
                                         gene_col = text_color_edges,
                                         show_genes = TRUE,
                                         vjust_label = -0.1,
                                         hjust_label = 0.5,
                                         node_size = 1,
                                         linesize = 0.5,
                                         ignore_direct_parents = FALSE,
                                         edge_color = tree_color, 
                                         node_color = tree_color,
                                         ladderize = TRUE,
                                         reorder_tree = TRUE,
                                         order_column = "cluster_order",# TRUE,
                                         highlight_leaves = highlight_leaves,
                                         highlight_bg =FALSE)
  ## order treeso that it always is the same ---> also done within above function if  reorder_tree = TRUE
  
  if(show_intermediate){celltype_tree_plot}
  
  current_offset = 5
  
  ##########
  ### prepare all heatmap data
  ##########
  
  min_pct_to_show = 1
  
  ## add main_gene
  heatmap_matrix1 = gene_pct_anno_filter[,main_gene,drop=F]
  rownames(heatmap_matrix1) = gene_pct_anno_filter$cluster_id
  #heatmap_matrix1[heatmap_matrix1 < min_pct_to_show] = NA
  
  #### second set
  heatmap_matrix2 = gene_pct_anno_filter[,other_genes[other_genes %in% colnames(gene_pct_anno_filter)],drop=F]
  rownames(heatmap_matrix2) = gene_pct_anno_filter$cluster_id
  #heatmap_matrix2[heatmap_matrix2 < min_pct_to_show] = NA
  
  #### main gene mouse
  heatmap_matrix3 = gene_pct_anno_filter[,stringr::str_to_title(main_gene),drop=F]
  rownames(heatmap_matrix3) = gene_pct_anno_filter$cluster_id
  #heatmap_matrix3[heatmap_matrix3 < min_pct_to_show] = NA
  
  #### other genes mouse
  heatmap_matrix4 = gene_pct_anno_filter[,stringr::str_to_title(other_genes)[stringr::str_to_title(other_genes) %in% colnames(gene_pct_anno_filter)],drop=F]
  rownames(heatmap_matrix4) = gene_pct_anno_filter$cluster_id
  #heatmap_matrix4[heatmap_matrix4 < min_pct_to_show] = NA
  
  ## SET SCALING AND LEGEND NAMES
  if(expr_type == "pct"){
    legend_name = "Pct"
    scale_limits = limits=c(0,50)
    scale_midpoint = min_expr_value # 2.5%
  }else{
    legend_name = "Expr"
    scale_limits = limits=c(0,1)
    scale_midpoint = min_expr_value # 0.05 or so
  }
  
  ##########
  ### Plot potential region
  ##########
  
  # get data from existing plot --> need for a and y !
  celltype_tree_plot_dat = celltype_tree_plot$data
  
  # dataframe to use as data in additional plot
  region_anno = gene_pct_anno_filter[,c("human_cluster","cluster_id","region")]# %>% dplyr::left_join(anno_df[,c("cluster_id","cluster_name")],by=c("human_cluster"="cluster_name"))
  region_anno$label = region_anno$cluster_id
  region_anno = dplyr::left_join(region_anno,celltype_tree_plot_dat[,c("label","x","y")],by="label")
  
  # make shorter
  region_anno$region_short = stringr::str_sub(region_anno$region, start = 1L, end = 12)
  #names(color_value_vector) = stringr::str_sub(names(color_value_vector), start = 1L, end = 20)
  
  region_anno$x = region_anno$x + current_offset
  region_anno$regions_label = region_anno$region_short
  region_anno$regions_label[is.na(region_anno$regions_label)] = "          "
  
  # make mapping vec
  tmp1 = gene_pct_anno_filter%>% dplyr::distinct(region,region_color)
  region_color_mapping = tmp1$region_color
  names(region_color_mapping) = tmp1$region
  region_color_mapping = region_color_mapping[region_color_mapping!=""]
  # add geom_label_to plot
  # celltype_tree_plot = celltype_tree_plot +
  #   ggnewscale::new_scale_fill() +
  #   geom_label(mapping = aes(x=x,y=y,group=label,label =regions_label,  fill=region_short),
  #              data = region_anno,
  #              color="white",
  #              size = label_textanno_size,label.padding= unit(0.05, "lines"),
  #              hjust = 0,
  #              inherit.aes=F)+
  #   ggplot2::scale_fill_manual(values = region_color_mapping,na.value = "white",guide= "none")
  # ggplot2::scale_fill_manual(values = getOkabeItoPalette(length(unique(region_anno$region_short))),na.value = "white",guide= "none")
  
  ### comment out:
  celltype_tree_plot = celltype_tree_plot +
    ggnewscale::new_scale_fill() +
    geom_text(mapping = aes(x=x,y=y,group=label,label =regions_label,  color=region_short),
              data = region_anno,
              size = label_textanno_size,label.padding= unit(0.05, "lines"),
              hjust = 0,
              fontface="bold",
              inherit.aes=F)+
    ggplot2::scale_color_manual(values = region_color_mapping,na.value = "white",guide= "none")
  # 
  if(show_intermediate){celltype_tree_plot}
  # 
  current_offset = current_offset + 5
  
  ##########
  ### Plot heatmap 2
  ##########
  
  celltype_tree_plot_1 = add_heatmap(circular_tree=celltype_tree_plot,
                                     heatmap_matrix = heatmap_matrix2,
                                     heatmap_colors=c("grey90","#E69F00"),
                                     scale_limits = scale_limits,
                                     heatmap_colnames =TRUE, 
                                     legend_title = paste0(legend_name," Human"),
                                     matrix_offset = current_offset,
                                     matrix_width =one_gene_width*ncol(heatmap_matrix2),
                                     colnames_angle=90,
                                     legend_text_size = legend_text_size,
                                     hjust_colnames=0,
                                     colnames_offset_y = 0.5,
                                     na_color = "white",
                                     heatmap_text_size=heatmap_text_size)+
    scale_fill_gradient2(low = "black",mid = "grey95",high = "#E69F00",na.value = "white",midpoint = scale_midpoint,limits=scale_limits,name = paste0(legend_name," Human"),n.breaks=nbreaks_scale, oob=squish)
  if(show_intermediate){celltype_tree_plot_1}
  
  ### this seems to be the one:
  current_offset = (max(celltype_tree_plot_1$data$x)*ncol(heatmap_matrix2)*one_gene_width)+current_offset+offset_elements
  
  ### old:
  #current_offset = current_offset + 2.5 * ncol(heatmap_matrix2)
  
  ##########
  ### Plot heatmap 1
  ##########
  
  # optionally use a second scale
  if(second_scale){
    scale_limits_2nd = c(0,max(heatmap_matrix1,na.rm = TRUE))
    scale_midpoint_2nd = (scale_midpoint / scale_limits[2])*scale_limits_2nd[2]
    legend_name_2nd = paste0(legend_name," Rec Human")
  }else{
    scale_limits_2nd = scale_limits  
    scale_midpoint_2nd = scale_midpoint
    legend_name_2nd = paste0(legend_name," Human")
  }
  
  library(scales)
  celltype_tree_plot_2 = add_heatmap(circular_tree=celltype_tree_plot_1,
                                     heatmap_matrix = heatmap_matrix1,
                                     heatmap_colors=c("grey90","#E69F00"), # "darkred"
                                     scale_limits = scale_limits,
                                     heatmap_colnames =TRUE, 
                                     legend_title = paste0(legend_name," Human"),
                                     matrix_offset = current_offset,
                                     matrix_width =one_gene_width* length(main_gene),
                                     colnames_angle=90,
                                     legend_text_size = legend_text_size,
                                     hjust_colnames=0,
                                     colnames_offset_y = 0.5,
                                     na_color = "white",
                                     heatmap_text_size=heatmap_text_size)+
    scale_fill_gradient2(low = "black",mid = "grey95",high = "#E69F00",na.value = "white",midpoint = scale_midpoint_2nd,limits=scale_limits_2nd,name = legend_name_2nd,n.breaks=nbreaks_scale, oob=squish)
  if(show_intermediate){celltype_tree_plot_2}
  
  current_offset = (max(celltype_tree_plot_1$data$x)*ncol(heatmap_matrix1)*one_gene_width)+current_offset+offset_elements
  
  ##########
  ### Plot heatmap 3 (mouse)
  ##########
  
  # optionally use a second scale
  if(second_scale){
    scale_limits_2nd = c(0,max(heatmap_matrix3,na.rm = TRUE))
    scale_midpoint_2nd = (scale_midpoint / scale_limits[2])*scale_limits_2nd[2]
    legend_name_2nd = paste0(legend_name," Rec Mouse")
  }else{
    scale_limits_2nd = scale_limits  
    scale_midpoint_2nd = scale_midpoint
    legend_name_2nd = paste0(legend_name," Mouse")
  }
  
  celltype_tree_plot_3 = add_heatmap(circular_tree=celltype_tree_plot_2,
                                     heatmap_matrix = heatmap_matrix3,
                                     heatmap_colors=c("grey90","#009E73"), # "darkblue"
                                     scale_limits = scale_limits,
                                     heatmap_colnames =TRUE, 
                                     legend_title = paste0(legend_name," Mouse"),
                                     matrix_offset = current_offset,
                                     matrix_width =one_gene_width* length(main_gene),
                                     colnames_angle=90,
                                     legend_text_size = legend_text_size,
                                     hjust_colnames=0,
                                     colnames_offset_y = 0.5,
                                     na_color = "white",
                                     heatmap_text_size=heatmap_text_size)+
    scale_fill_gradient2(low = "black",mid = "grey95",high = "#009E73",na.value = "white",midpoint = scale_midpoint_2nd,limits = scale_limits_2nd,name = legend_name_2nd,n.breaks=nbreaks_scale ,oob=squish) # ,breaks = seq(scale_limits_2nd[1], scale_limits_2nd[2], length.out = nbreaks_scale + 1)
  if(show_intermediate){celltype_tree_plot_3}
  
  current_offset = (max(celltype_tree_plot_1$data$x)*ncol(heatmap_matrix3)*one_gene_width)+current_offset+offset_elements
  
  ##########
  ### Plot heatmap 4 (mouse)
  ##########
  
  celltype_tree_plot_4 = add_heatmap(circular_tree=celltype_tree_plot_3,
                                     heatmap_matrix = heatmap_matrix4,
                                     heatmap_colors=c("grey90","#009E73"),
                                     scale_limits = scale_limits,
                                     heatmap_colnames =TRUE, 
                                     legend_title = paste0(legend_name," Mouse"),
                                     matrix_offset = current_offset,
                                     matrix_width =one_gene_width*ncol(heatmap_matrix4),
                                     colnames_angle=90,
                                     legend_text_size = legend_text_size,
                                     hjust_colnames=0,
                                     colnames_offset_y = 0.5,
                                     na_color = "white",
                                     heatmap_text_size=heatmap_text_size)+
    scale_fill_gradient2(low = "black",mid = "grey95",high = "#009E73",na.value = "white",midpoint = scale_midpoint,limits = scale_limits,name = paste0(legend_name," Mouse"),n.breaks=nbreaks_scale, oob=squish)
  if(show_intermediate){celltype_tree_plot_4}
  
  current_offset = (max(celltype_tree_plot_1$data$x)*ncol(heatmap_matrix4)*one_gene_width)+current_offset+offset_elements+1.5
  
  # increase legend text
  celltype_tree_plot_4 = celltype_tree_plot_4+theme(legend.text = element_text(size=10))
  
  ##########
  ### Plot mouse clusters dots
  ##########
  
  # get data from existing plot --> need for a and y !
  celltype_tree_plot_dat = celltype_tree_plot$data
  
  # dataframe to use as data in additional plot
  mouse_mapping_occ = gene_pct_anno_filter[,c("human_cluster","mouse_cluster","from_occ","species_cor")]
  mouse_mapping_occ$label = stringr::str_extract(mouse_mapping_occ$human_cluster,"C[0-9]-[0-9]+")
 # print(head(celltype_tree_plot_dat))
  mouse_mapping_occ = dplyr::left_join(mouse_mapping_occ,celltype_tree_plot_dat[,c("label","x","y")],by="label")
  # for empty ones -- species cor is na --> below threshold e.g 0.6
  mouse_mapping_occ$from_occ[is.na(mouse_mapping_occ$mouse_cluster)] =0
  
  mouse_mapping_occ$x = mouse_mapping_occ$x + current_offset
  #print(head(mouse_mapping_occ))
  # add geom_label_to plot
  celltype_tree_plot_5 = celltype_tree_plot_4 + 
    ggnewscale::new_scale_fill()+
    geom_label(mapping = aes(x=x,y=y,group=from_occ,label =from_occ,fill = species_cor),
               data = mouse_mapping_occ,
               size = label_textanno_size,label.padding= unit(0.05, "lines"),
               hjust = 0,
               #  fill="white",
               color="grey20",
               inherit.aes=F)+
    scale_fill_gradient2(low = "orange",mid = "yellow",high = "darkgreen",na.value = "darkorange",midpoint = 0.85,limits=c(0.6,1),name = "Species Cor")
  
  if(show_intermediate){celltype_tree_plot_5}
  
  current_offset = current_offset + 1.5
  
  ##########
  ### Plot mouse clusters names
  ##########
  
  # get data from existing plot --> need for a and y !
  celltype_tree_plot_dat = celltype_tree_plot$data
  
  # dataframe to use as data in additional plot
  mouse_mapping = gene_pct_anno_filter[,c("human_cluster","cluster_id","mouse_cluster","species_cor")]
  mouse_mapping$label = mouse_mapping$cluster_id
  mouse_mapping = dplyr::left_join(mouse_mapping,celltype_tree_plot_dat[,c("label","x","y")],by="label")
  # for empty ones -- species cor is na --> below threshold e.g 0.6
  mouse_mapping$mouse_cluster[is.na(mouse_mapping$species_cor)] = "                                       "
  
  mouse_mapping$x = mouse_mapping$x + current_offset
  
  # add geom_label_to plot
  celltype_tree_plot_6 = celltype_tree_plot_5 +
    # ggnewscale::new_scale_fill() +
    geom_text(mapping = aes(x=x,y=y,group=label,label =mouse_cluster), # ,  fill=species_cor
              data = mouse_mapping,
              size = label_textanno_size,#label.padding= unit(0.05, "lines"),
              hjust = 0,#,fill="white"
              inherit.aes=F)#+
  #scale_fill_gradient2(low = "orange",mid = "yellow",high = "darkgreen",na.value = "darkorange",midpoint = 0.75,limits=c(0.5,1),name = "Species Cor")
  #scale_fill_gradient(low = "grey90",high = "darkgreen",na.value = "grey90",limits=c(0.7,1))
  
  if(show_intermediate){celltype_tree_plot_6}
  
  
  current_offset = current_offset + mouse_label_offset
  
  ##########
  ### empty text for right margin
  ##########
  
  mouse_mapping_dummy = mouse_mapping
  #mouse_mapping_dummy$x
  mouse_mapping_dummy$x = mouse_mapping_dummy$x + mouse_label_offset
  # add geom_label_to plot
  celltype_tree_plot_7 = celltype_tree_plot_6 + 
    geom_blank(mapping = aes(x=x,y=y),data = mouse_mapping_dummy)
  
  # show
  return(celltype_tree_plot_7)
  
  
}
