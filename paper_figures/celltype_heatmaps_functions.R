##########
### Function to generate all plots
##########

save_cluster_minitree_helper = function(cluster_minitree,node_to_use,outpath,remove_name_part="",suffix=""){
  dir.create(outpath,showWarnings = F)
  if(suffix!= ""){suffix = paste0(suffix,"_")}
  filename = paste0(gsub("-","_",node_to_use),"_",gsub("-","_",remove_name_part),suffix,"plot_")
  
  # save plots
  ggsave(filename = paste0(outpath,filename,"tree_heatmap",".pdf"),plot = cluster_minitree$tree_heatmap, "pdf",dpi=300,width=200+10*nrow(cluster_minitree$tree_heatmap$layers[[10]]$data),height = 200,units="mm")#
  ggsave(filename = paste0(outpath,filename,"cluster_umap",".pdf"),plot = cluster_minitree$cluster_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
  ggsave(filename = paste0(outpath,filename,"overview_umap",".pdf"),plot = cluster_minitree$overview_umap, "pdf",dpi=300,width=200,height = 200,units="mm")#
}

make_cluster_minitree_helper = function(node_to_use,node_to_ignore=character(0),human_hypo_combined,combined_edgelist_mrtree,keep_quantile = 0.975,hjust_colnames=0,
                                        additional_nodes=character(0),additional_genes=character(0),label_size=6,label_size_tip=6,vjust_label=-0.25,hjust_label=0,
                                        tree_nodesize=5,tree_linesize = 1.5,heatmap_text_size=4,expr_type = "pct.exp",use_custom_scale =FALSE,gene_col="black",hilight_col = "red",
                                        remove_name_part="", heatmap_offset = 1, heatmap_offset_base = 4, matrix_width_per_gene = 0.12,avg_exp_cap = NULL,leaf_level=6,final_name_column="C4_named",seurat_ident_column="C4"){
  
  clusterlevel_to_use = stringr::str_extract(node_to_use[1],pattern = "C[0-9]+")
  relevantNodes = c(node_to_use,scUtils::find_children(nodes = node_to_use,edges = combined_edgelist_mrtree))
  relevantNodes = c(relevantNodes,additional_nodes) # manually add rootnode
  ignore_nodes = c(node_to_ignore,scUtils::find_children(nodes = node_to_ignore,edges = combined_edgelist_mrtree))
  relevantNodes = relevantNodes[!relevantNodes %in% ignore_nodes]
  ## subset data
  message("Subset data")
  cellsh = human_hypo_combined@meta.data$Cell_ID[ human_hypo_combined@meta.data[,clusterlevel_to_use] %in% c(node_to_use)]
  subset_seurat_plot = hypoMapUtils::subset_zoom_seurat(human_hypo_combined,cells = cellsh,keep_quantile = keep_quantile,reduction = "umap")
  subset_seurat = subset(human_hypo_combined,cells=cellsh) # to calc stuff
  
  ## make a highlight UMAP
  ## scattermore does not work with highlighting anymore
 # human_hypo_combined2 = human_hypo_combined
  human_hypo_combined@meta.data$hilight = NA
  human_hypo_combined@meta.data$hilight[rownames(human_hypo_combined@meta.data) %in% colnames(subset_seurat_plot)] = "hilight"
  message("overview_umap")
  overview_umap = DimPlot(human_hypo_combined,group.by = "hilight",reduction = "umap",label = F,raster = TRUE,raster.dpi = c(1024,1024),pt.size = 1.5,order = TRUE)+
    NoLegend()+NoAxes()+ggtitle(NULL)+scale_color_manual(values=c("hilight" = hilight_col),na.value = "grey80")
  #
  # old
  #overview_umap = DimPlot(human_hypo_combined,cells.highlight = colnames(subset_seurat_plot),sizes.highlight = 0.3,reduction = "umap",label = F,raster = F,raster.dpi = c(512,512),pt.size = 0.1)+NoLegend()+NoAxes()
  #DimPlot(subset_seurat_plot,group.by = "C4",label = TRUE,cols = color_order$color,reduction = "umap",raster.dpi = c(2048,2048),pt.size = 1,shuffle = TRUE,repel = TRUE,label.size = 5)+NoAxes()+NoLegend()+ggtitle(NULL)
  
  ### make a UMAP plot
  nodes_for_umap_annotation = data.frame(cluster_id = relevantNodes)
  nodes_for_umap_annotation = dplyr::left_join(nodes_for_umap_annotation,anno_df[,c("cluster_id","cluster_name","first_cluster_name","clusterlevel")])
  
  ## prepare colors and order
  color_order = nodes_for_umap_annotation[,c("cluster_id","first_cluster_name","cluster_name")]
  color_order = color_order[grepl(seurat_ident_column,color_order$cluster_id),]
  color_order$color = getOkabeItoPalette(nrow(color_order))
  #print(color_order)
  color_order$cluster_name = factor(color_order$cluster_name,levels = color_order$cluster_name )# make.unique
  color_order$final_name = gsub(remove_name_part,"",color_order$cluster_name)
  #print(color_order)
  #plot
  subset_seurat_plot@meta.data$final_name = gsub(remove_name_part,"",subset_seurat_plot@meta.data[,final_name_column])
  subset_seurat_plot@meta.data$final_name = factor(subset_seurat_plot@meta.data$final_name,levels = color_order$final_name)
  message("cluster_umap")
  cluster_umap = DimPlot(subset_seurat_plot,group.by = "final_name",label = TRUE,cols = color_order$color,reduction = "umap",raster.dpi = c(2048,2048),pt.size = 1,shuffle = TRUE,repel = TRUE,label.size = 5)+NoAxes()+NoLegend()+ggtitle(NULL)
  #FeaturePlot(subset_seurat_plot,reduction = "umap",order = TRUE,features = "QRFPR")
  
  ### make tree
  
  tree_color="grey40"
  
  message("tree_plot")
  message("relevantNodes: ",paste0(relevantNodes,collapse = " | "))
  subset_edge_list = combined_edgelist_mrtree[combined_edgelist_mrtree$from %in% c(relevantNodes),]
  mini_tree = plot_mini_tree_ov(edgelist = subset_edge_list,
                                leaf_level=leaf_level,
                                #layout="dendrogram",
                                anno_df = anno_df ,
                                metadata=human_hypo_combined@meta.data,
                                label_size = label_size, 
                                label_size_tip = label_size_tip,
                                show_genes = TRUE,
                                gene_col = gene_col,
                                vjust_label = vjust_label,
                                hjust_label = hjust_label,
                                node_size = tree_nodesize,
                                linesize = tree_linesize,
                                edge_color = tree_color, 
                                node_color = tree_color,
                                only_relevant_nodepoints =TRUE)
  
  #mini_tree
  gene_order = mini_tree$data %>% dplyr::arrange(branch.length,y) %>%
    dplyr::filter(first_cluster_name %in% c(unique(as.character(nodes_for_umap_annotation$first_cluster_name)))) %>%  
    dplyr::distinct(first_cluster_name,.keep_all = TRUE)
  #print(gene_order$first_cluster_name)
  ### cluster coloring 
  mini_treedata = mini_tree$data
  data_bg_fill <- color_order[color_order$cluster_id %in% mini_treedata$label,]
  data_bg_fill$node =mini_treedata$node[match(data_bg_fill$cluster_id,mini_treedata$label)]
  data_bg_fill$fill = data_bg_fill$color
  mini_tree2 <- mini_tree +
    ggtree::geom_hilight(
      data = data_bg_fill,
      mapping = aes(
        node = node,
        fill = I(fill)
      ),
      align = "none",
      alpha = 0.4,
      extendto = 2,
    )
  #mini_tree2
  # move rectangles as first layer so that tree is plotted on top
  mini_tree2$layers = c(mini_tree2$layers[[length(mini_tree2$layers)]],mini_tree2$layers[1:(length(mini_tree2$layers)-1)])
  
  #print(head(mini_tree2$layers[[1]]$data))
  mini_tree2$layers[[1]]$data$xmin = 5
  # save = mini_tree2$layers[[1]]$data$xmax
  # mini_tree2$layers[[1]]$data$xmax = mini_tree2$layers[[1]]$data$xmin
  # mini_tree2$layers[[1]]$data$xmin = save
  
  
  ## add gene expression as heatmap
  message("tree_plot heat map")
  # use seurat to get data
  Idents(subset_seurat) = seurat_ident_column # "C4"
  features_query = unique(c(unique(as.character(gene_order$first_cluster_name)),additional_genes))
  features_query = features_query[features_query %in% rownames(subset_seurat)]
  #dotplot= Seurat::DotPlot(subset_seurat,features = c(unique(as.character(nodes_for_umap_annotation$first_cluster_name)),"SST"),scale = FALSE)
  dotplot= Seurat::DotPlot(subset_seurat,features = features_query,scale = FALSE)
  dotplot_data = dotplot$data
  if(expr_type == "pct.exp"){
    heatmap_data = dotplot_data %>% dplyr::select(-avg.exp,-avg.exp.scaled) %>% tidyr::spread(key=features.plot,value=pct.exp)
    legend_title = "Pct"
    scale_lims = c(0,100)
  }else if(expr_type == "avg.exp"){ # with scale FALS avg.exp.scaled contains the normalized counts !!
    heatmap_data = dotplot_data %>% dplyr::select(-avg.exp,-pct.exp) %>% tidyr::spread(key=features.plot,value=avg.exp.scaled)
    legend_title = "Expr"
    if(is.null(avg_exp_cap)){
      scale_lims = c(0,max(dotplot_data$avg.exp.scaled))
    }else{
      scale_lims = c(0,avg_exp_cap)
    }
  }else if(expr_type == "avg.exp.scaled"){
    dotplot= Seurat::DotPlot(subset_seurat,features = features_query,scale = TRUE)
    dotplot_data = dotplot$data
    heatmap_data = dotplot_data %>% dplyr::select(-pct.exp,-avg.exp) %>% tidyr::spread(key=features.plot,value=avg.exp.scaled)
    legend_title = "Scaled Expr"
    limits_scaled = max(c(abs(min(dotplot_data$avg.exp.scaled)),abs(max(dotplot_data$avg.exp.scaled))))
    scale_lims =c(-1*limits_scaled,limits_scaled)
  }else if(expr_type == "specificity"){
    # calc specificty
    heatmap_data = calculate_specificity(nodes_for_umap_annotation,human_hypo_combined,downsample = 10000)
    legend_title = "Specificity"
    scale_lims =c(max(min(heatmap_data[,2:ncol(heatmap_data)]),-10),min(max(heatmap_data[,2:ncol(heatmap_data)]),10))
  }else{
    stop("Please provide a valid expr_type.") 
  }
  heatmap_matrix = heatmap_data[,2:ncol(heatmap_data)]
  rownames(heatmap_matrix) = heatmap_data[,1]
  
  # split matrix by level
  temp_annotation = gene_order %>% dplyr::distinct(first_cluster_name,.keep_all = TRUE) # use distinct to remove to one occurence of each gene
  genes_per_level = split(temp_annotation$first_cluster_name,temp_annotation$clusterlevel) # split into list of genes per level
  genes_per_level[["additional"]] = additional_genes
  
  ## loop over levels
  heatmap_matrix_list =list()
  for(i in 1:length(genes_per_level)){
    genes = genes_per_level[[i]]
    genes = genes[genes %in% rownames(subset_seurat_plot)] 
    if(length(genes) == 0){
      heatmap_matrix_list[[i]]=NULL
    }else{
      heatmap_matrix_list[[i]]= heatmap_matrix[,genes,drop=F]
    }
  }
  
  heatmap_matrix_list[sapply(heatmap_matrix_list, is.null)] <- NULL
  # loop over matrices
  mini_tree2_heat = mini_tree2
  library(scales)
  matrix_offset = heatmap_offset_base
  for(current_matrix in heatmap_matrix_list){
    # width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
    mini_tree2_heat = add_heatmap(circular_tree=mini_tree2_heat,
                                  heatmap_matrix = current_matrix,
                                  heatmap_colors= c("grey90","darkred"), # uses scale_fill_viridis_c anyway !!
                                  scale_limits = scale_lims,
                                  heatmap_colnames =TRUE, 
                                  legend_title = legend_title,
                                  matrix_offset = matrix_offset,
                                  matrix_width = ncol(current_matrix)*matrix_width_per_gene,
                                  colnames_angle=90,
                                  legend_text_size = 7,
                                  hjust_colnames = hjust_colnames,
                                  na_color = "white",
                                  heatmap_text_size=heatmap_text_size)+
      #  scale_fill_gradient(limits=c(0,4),oob=squish)
      scale_fill_viridis_c(limits=scale_lims,oob=squish,na.value = "white",name=legend_title)
    if(expr_type %in% c("specificity","avg.exp.scaled")){
      mini_tree2_heat = mini_tree2_heat + scale_fill_gradient2(low = "#4575B4",mid="#FFFFBF",high="#D73027",na.value = "white",midpoint = 0,limits=scale_lims,oob=squish) # similar to RdYlBu
    }else{
      # stick with viridis 
    }
    if(use_custom_scale){
      mini_tree2_heat = mini_tree2_heat + scale_fill_gradient2(low = "white",mid="#FFFFBF",high="#D73027",na.value = "white",midpoint = 0.001,limits=scale_lims,oob=squish) # similar to RdYlBu
    }
    # update offset
    matrix_offset = matrix_offset +heatmap_offset + (ncol(current_matrix)*(mini_tree2_heat$data$x %>% range(na.rm = TRUE) %>% diff)*matrix_width_per_gene ) #- (((mini_tree2_heat$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(current_matrix))*matrix_width_per_gene)
  }
  #rotate_tree(mini_tree2_heat,angle = 90)
  #mini_tree2_heat
  
  # return
  return(list(tree_heatmap = mini_tree2_heat, cluster_umap = cluster_umap,overview_umap = overview_umap, heatmap_data = heatmap_data))
  
}

##########
### Calculate specificty for heatmap
##########

calculate_specificity = function(nodes_for_umap_annotation,human_hypo_combined,cluster_level = "C4",offset = 0.01,downsample = 10000,seed =1234){
  
  all_clusters = unique(nodes_for_umap_annotation$cluster_id)[grepl(cluster_level,unique(nodes_for_umap_annotation$cluster_id))]
  all_specificities = list()
  message("Calculating specificity using ",downsample," cells for non-cluster reference")
  for(i in 1:length(all_clusters)){
    message("For cluster ",i," of ",length(all_clusters))
    current_cluster = all_clusters[i]
    current_cells = human_hypo_combined@meta.data$Cell_ID[human_hypo_combined@meta.data[,cluster_level] == current_cluster]
    other_cells = human_hypo_combined@meta.data$Cell_ID[!human_hypo_combined@meta.data[,cluster_level] == current_cluster & human_hypo_combined@meta.data$C1 %in% c("C1-1","C1-2","C1-3","C1-4","C1-5")]
    if(length(other_cells) > downsample){
      set.seed(seed) 
      other_cells =sample(other_cells,downsample)
    }
    fc_stats = FoldChange_seurat(human_hypo_combined@assays$RNA@data,current_cells,other_cells,features = unique(c(unique(as.character(nodes_for_umap_annotation$first_cluster_name)),additional_genes)))
    fc_stats$specificity = fc_stats$avg_log2FC * ((fc_stats$pct.1+offset) / (fc_stats$pct.2+offset))
    fc_stats$specificity[fc_stats$avg_log2FC < 0] = fc_stats$avg_log2FC[fc_stats$avg_log2FC < 0] * ((fc_stats$pct.2[fc_stats$avg_log2FC < 0]+offset) / (fc_stats$pct.1[fc_stats$avg_log2FC < 0]+offset))
    specificity_df = cbind(current_cluster,t(fc_stats$specificity))
    colnames(specificity_df)[1] = "id"
    colnames(specificity_df)[2:ncol(specificity_df)] = fc_stats$gene
    all_specificities[[current_cluster]] = specificity_df
  }
  heatmap_data = do.call(rbind,all_specificities,quote = F) %>% as.data.frame()
  heatmap_data[,2:ncol(heatmap_data)] = apply(heatmap_data[,2:ncol(heatmap_data)],2,as.numeric)
  
  return(heatmap_data)
}

##########
### Calculate specificty for heatmap
##########

calculate_specificity_general = function(cluster_ids=NULL,cluster_col,features,human_hypo_combined,filter_col="C1",cluster_filter = c("C1-1","C1-2","C1-3","C1-4","C1-5"),offset = 0.01,downsample = 10000,seed =1234){
  
  if(is.null(cluster_ids)){cluster_ids = unique(as.character(human_hypo_combined@meta.data[,cluster_col]))}
  if(is.null(cluster_filter)){cluster_filter = unique(as.character(human_hypo_combined@meta.data[,filter_col]))}
  all_specificities = list()
  message("Calculating specificity using ",downsample," cells for non-cluster reference")
  for(i in 1:length(cluster_ids)){
    message("For cluster ",i," of ",length(cluster_ids))
    current_cluster = cluster_ids[i]
    current_cells = human_hypo_combined@meta.data$Cell_ID[human_hypo_combined@meta.data[,cluster_col] == current_cluster]
    other_cells = human_hypo_combined@meta.data$Cell_ID[!human_hypo_combined@meta.data[,cluster_col] == current_cluster & human_hypo_combined@meta.data[,filter_col] %in% cluster_filter]
    if(length(other_cells) > downsample){
      set.seed(seed) 
      other_cells =sample(other_cells,downsample)
    }
    fc_stats = FoldChange_seurat(human_hypo_combined@assays$RNA@data,current_cells,other_cells,features = features)
    fc_stats$specificity = fc_stats$avg_log2FC * ((fc_stats$pct.1+offset) / (fc_stats$pct.2+offset))
    fc_stats$specificity[fc_stats$avg_log2FC < 0] = fc_stats$avg_log2FC[fc_stats$avg_log2FC < 0] * ((fc_stats$pct.2[fc_stats$avg_log2FC < 0]+offset) / (fc_stats$pct.1[fc_stats$avg_log2FC < 0]+offset))
    specificity_df = cbind(current_cluster,t(fc_stats$specificity))
    colnames(specificity_df)[1] = "cluster"
    colnames(specificity_df)[2:ncol(specificity_df)] = fc_stats$gene
    all_specificities[[current_cluster]] = specificity_df
  }
  heatmap_data = do.call(rbind,all_specificities,quote = F) %>% as.data.frame()
  heatmap_data[,2:ncol(heatmap_data)] = apply(heatmap_data[,2:ncol(heatmap_data)],2,as.numeric)
  
  return(heatmap_data)
}


##########
### foldchange from seurat
##########

#' Calculate general FC and pcts between two cell groups
#'
#' https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
#'
#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#'
#' @importFrom Matrix rowSums

FoldChange_seurat <- function(
    object,
    cells.1,
    cells.2,
    features = NULL
) {
  if(is.null(features)){
    features <- rownames(x = object)
  }else{
    features <- features[features %in% rownames(x = object)]
  }
  library(Matrix)
  # mean.fxn -> always use log2
  fc.name = "avg_log2FC"
  mean.fxn <- function(x) {
    base = 2
    pseudocount.use = 1
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
  }
  
  # Calculate percent expressed
  thresh.min <- 0
  pct.1 <- round(
    x = Matrix::rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = Matrix::rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
      length(x = cells.2),
    digits = 3
  )
  # Calculate fold change
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- data.frame(gene = features,fc = fc, pct.1 = pct.1, pct.2 = pct.2)
  colnames(fc.results)[2] <- fc.name
  return(fc.results)
}
