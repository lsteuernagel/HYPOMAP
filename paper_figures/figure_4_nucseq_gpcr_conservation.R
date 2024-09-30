
##########
###  Load data
##########

library(tidyverse)
library(ggplot2)
library(ggh4x)
source("utility_functions.R")
source("integration_pipeline/harmonization_functions.R")
source("merge_human_mouse_neurons/cluster_matching_functions.R")
source("paper_figures/tree_plotting_functions.R")

human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"

# result
matched_clusters_final = data.table::fread("merge_human_mouse_neurons/matched_clusters_scviCorPruned.tsv",data.table = F)


##########
###  load scseq objects
##########

# load seurat object
if(!exists("hypothalamus_neurons_cross_species", envir = .GlobalEnv )){
  hypothalamus_neurons_cross_species = readRDS(paste0(human_hypo_path,"hypothalamus_neurons_cross_species.rds")) 
}

## load markers
human_marker_genes = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_harmonization_3/human_hypo_combined/human_hypo_combined_comparisons_all_annotated.txt",data.table = F)
#human_marker_genes = data.table::fread("paper_figures/revision_figures/suppl_tables/marker_genes_all_top100.txt",data.table = F)

## mouse
hypoMap = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap.rds") # requires mouse HypoMap from: https://doi.org/10.17863/CAM.87955

hypoMap_annos = hypoMap@misc$annotation_result
hypoMap_edgelist =  hypoMap@misc$clustering_edgelist %>% 
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,from_named = clean_names_withID),by=c("from"="cluster_id")) %>%
  dplyr::left_join(hypoMap_annos %>% dplyr::select(cluster_id,to_named = clean_names_withID),by=c("to"="cluster_id"))


### Add clustering results from human neurons using the final object
if(!exists("human_hypo_combined", envir = .GlobalEnv )){
  human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP.rds"))
}
combined_edgelist_mrtree = human_hypo_combined@misc$edgelist

##########
###  GPCR comparison
##########

## load gpcrs and get human to mouse gene matching

gpcrs_human = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28G-protein%20coupled%20receptor%29%20AND%20(cc_scl_term%3ASL-0039)%20AND%20(organism_id%3A9606)%20NOT%20(cc_scl_term%3ASL-0191)%20AND%20(ft_transmem:*)" %>%
  url() %>%
  gzcon() %>%
  readr::read_tsv() %>%
  dplyr::filter(Reviewed == "reviewed")

gpcrs_mouse = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28G-protein%20coupled%20receptor%29%20AND%20(cc_scl_term%3ASL-0039)%20AND%20(organism_id%3A10090)%20NOT%20(cc_scl_term%3ASL-0191)%20AND%20(ft_transmem:*)" %>%
  url() %>%
  gzcon() %>%
  readr::read_tsv() %>%
  dplyr::filter(Reviewed == "reviewed")

gpcrs_human_vector = unlist(sapply(gpcrs_human$`Gene Names`,function(x){strsplit(x,split=" ")[[1]]}))
gpcrs_mouse_vector = stringr::str_to_upper(unlist(sapply(gpcrs_human$`Gene Names`,function(x){strsplit(x,split=" ")[[1]]})))
more_neuroactive_receptors = stringr::str_to_upper(hypoMapUtils::neuroactive_receptor_genes)
gpcrs_human_genes = unique(c(as.character(gpcrs_human_vector),more_neuroactive_receptors,gpcrs_mouse_vector))
gpcrs_human_genes = gpcrs_human_genes[!grepl("ADCY",gpcrs_human_genes)]
gpcrs_human_genes = gpcrs_human_genes[gpcrs_human_genes %in% rownames(human_hypo_combined@assays$RNA@counts)]

gpcrs_human_genes = unique(c(gpcrs_human_genes,c("HTR3A", "HTR3B","HTR3C", "HTR3D", "HTR3E")))


# get conserved genes
conversion_table = scUtils::get_mouse_human_genes_conversion()
gpcr_table = data.frame(human_gene = gpcrs_human_genes) %>% dplyr::left_join(conversion_table %>% dplyr::select(-ensembl_gene_id,-hsapiens_homolog_ensembl_gene),by=c("human_gene"="human_gene"),multiple = "all")

gpcr_table = gpcr_table[gpcr_table$mouse_gene != "Gm21973",]


## get all neruon cell ids
human_neurons = hypothalamus_neurons_cross_species@meta.data$Cell_ID[hypothalamus_neurons_cross_species@meta.data$species == "human"]
hypoMap_neurons = hypothalamus_neurons_cross_species@meta.data$Cell_ID[hypothalamus_neurons_cross_species@meta.data$species == "mouse"]
##

##########
###  GPCR conservation function
##########

# at the moment cannot handle 1:many gene relationships well !!!!

#' Plot conservation of marker / receptor genes of matched human and mouse clusters
#'
#' Retrieves the matching mosue or human clusters and queries the seurat object for gene expression in the cells of the selected clusters
#' genes are selected based on marker gene tbales that have to be provided
#' Then computes pcts and compares them to make a heatmap of up to 20 (+ additional manual genes) marker genes to compare across species.
#'
#' @param target_cluster_human a human cluster or NULL (if target_cluster_mouse is provided)
#' @param target_cluster_mouse a mouse cluster or NULL (if target_cluster_human is provided)
#' @param edgelist the human to muse mtached clusters edgelist
#' @param gene_table a conversion gene_table with human and mouse gene (can be subsetted to GPCRs, but also all genes)
#' @param additional_genes vector of manually added genes. have to be in gene_table
#' @param top_n_genes how many genes. defaults to 10
#' @param human_hypo_combined the human only seurat object
#' @param hypoMap the mouse only seurat object
#' @param human_neurons ids of neurons to subset to when calculating glbal pct. defaults to NULL which will use all cells. make sure to provide valid is when not NULL
#' @param hypoMap_neurons ids of neurons to subset to when calculating glbal pct. defaults to NULL which will use all cells. make sure to provide valid is when not NULL
#' @param human_hypo_combined_clustercol column to access in human_hypo_combined
#' @param hypoMap_clustercol column to access in hypoMap
#' @param human_marker_genes
#' @param mouse_marker_genes
#' @param min_spec specificity threshold when filtering markers. defaults to 1
#'
#' @return a list of ggplot and dataframe objects with comparison of top_n_genes between mouse and human clusters
#'

# human_marker_genes = comparisons_all_updated 
# mouse_marker_genes = hypoMap@misc$marker_genes_all

conserved_gene_expression <- function(target_cluster_human = NULL, target_cluster_mouse = NULL, edgelist,gene_table,additional_genes=character(0),top_n_genes = 10, human_hypo_combined , hypoMap ,human_neurons = NULL, hypoMap_neurons = NULL, human_hypo_combined_clustercol = "C6", hypoMap_clustercol = "C465_named", human_marker_genes, mouse_marker_genes , min_spec = 1,expr_scale_max=2){
  
  library(scales)
  # checks
  if(is.null(target_cluster_human) & is.null(target_cluster_mouse)){stop("Please specifiy cluster name")}
  if(! all(c("human_gene","mouse_gene") %in% colnames(gene_table))){stop("gene_table must have columns human_gene and mouse_gene")}
  if(! all(c("from","to","similarity") %in% colnames(edgelist))){stop("edgelist must have columns from and to")}
  
  edgelist = edgelist %>% dplyr::arrange(desc(similarity))  ## due to sorting the higher partner cluster is always the first
  
  # retrieve matched mouse clusters
  if(!is.null(target_cluster_human) & is.null(target_cluster_mouse)){
    if(! any(target_cluster_human %in% edgelist$from) ){stop("Provide valid human (from) cluster name.")}
    target_cluster_mouse = edgelist$to[edgelist$from %in% target_cluster_human]
    if(length(target_cluster_mouse) == 0){stop("Cannot find correlated mouse cluster")}
    message("Used provided human cluster ",target_cluster_human," and found ",length(target_cluster_mouse)," matched mouse cluster(s).")
  }else if(is.null(target_cluster_human) & !is.null(target_cluster_mouse)){
    if(! any(target_cluster_mouse %in% edgelist$to )){stop("Provide valid mouse (to) cluster name.")}
    target_cluster_human = edgelist$from[edgelist$to %in% target_cluster_mouse]
    if(length(target_cluster_human) == 0){stop("Cannot find correlated human cluster")}
    message("Used provided mouse cluster ",target_cluster_mouse," and found ",length(target_cluster_human)," matched human cluster(s).")
  }else{
    if(! any(target_cluster_human %in% unique(human_hypo_combined@meta.data[,human_hypo_combined_clustercol]))){stop("Provide valid human (from) cluster name or set not NULL.")}
    if(! any(target_cluster_mouse %in% unique(hypoMap@meta.data[,hypoMap_clustercol]))){stop("Provide valid mouse (to) cluster name or set not NULL.")}
    message("Using provided mouse and human cluster for comparison. The function does not check for correlation, so ensure that this makes sense!")
  }
  
  # choose receptors to prioritize
  human_gpcrs_cluster = human_marker_genes[human_marker_genes$gene %in% gene_table$human_gene & human_marker_genes$specificity > min_spec & human_marker_genes$p_val_adj < 0.001 & human_marker_genes$name %in% target_cluster_human ,] %>%
    dplyr::group_by(gene) %>% dplyr::slice_max(order_by = specificity,n = 1,with_ties = F) %>%
    dplyr::ungroup() %>%  dplyr::slice_max(order_by = specificity,n = top_n_genes,with_ties = F) %>%
    dplyr::arrange(desc(specificity))
  
  # choose mouse receptors to prioritize
  mouse_gpcrs_cluster = mouse_marker_genes[mouse_marker_genes$gene %in% gene_table$mouse_gene & mouse_marker_genes$specificity > min_spec & mouse_marker_genes$p_val_adj < 0.001 & mouse_marker_genes$cluster_name %in% target_cluster_mouse ,] %>%
    dplyr::group_by(gene) %>% dplyr::slice_max(order_by = specificity,n = 1,with_ties = F)  %>%
    dplyr::ungroup() %>%  dplyr::slice_max(order_by = specificity,n = top_n_genes,with_ties = F) %>%
    dplyr::arrange(desc(specificity))
  
  gene_table_subset = gene_table[gene_table$human_gene %in% human_gpcrs_cluster$gene | gene_table$mouse_gene %in% mouse_gpcrs_cluster$gene | gene_table$human_gene %in% additional_genes  | gene_table$mouse_gene %in% additional_genes,]
  
  ## human expression
  Idents(human_hypo_combined) = human_hypo_combined_clustercol
  dotplot_human = Seurat::DotPlot(human_hypo_combined,features = unique(gene_table_subset$human_gene),idents = target_cluster_human,scale = FALSE) # with scale FALS avg.exp.scaled contains the normalized counts !!
  dotplot_data_human = dotplot_human$data
  expr_genes_human_wide = dotplot_data_human %>% dplyr::select(cluster= id, gene= features.plot, avg.exp.scaled) %>% tidyr::spread(key = cluster,value = avg.exp.scaled)
  
  ## mouse expression
  Idents(hypoMap) = hypoMap_clustercol
  dotplot_mouse= Seurat::DotPlot(hypoMap,features = unique(gene_table_subset$mouse_gene),idents = target_cluster_mouse,scale = FALSE) # with scale FALS avg.exp.scaled contains the normalized counts !!
  dotplot_mouse_data = dotplot_mouse$data
  expr_genes_mouse_wide = dotplot_mouse_data %>% dplyr::select(cluster= id, gene= features.plot, avg.exp.scaled) %>% tidyr::spread(key = cluster,value = avg.exp.scaled)
  
  # put side-by-side
  genes_cross_species = expr_genes_human_wide %>%
    dplyr::left_join(gene_table_subset[,c("human_gene","mouse_gene")],by = c("gene"="human_gene")) %>%
    dplyr::left_join(expr_genes_mouse_wide,by = c("mouse_gene"="gene"))
  
  genes_cross_species_clean = genes_cross_species[!is.na(genes_cross_species$gene) & !is.na(genes_cross_species$mouse_gene),]
  genes_cross_species_clean[is.na(genes_cross_species_clean)] = 0
  
  # genes_cross_species_clean$pct_ratio = log2((genes_cross_species_clean[,2]+0.1) / (genes_cross_species_clean[,4]+0.1))
  # genes_cross_species_clean$pct_sum = genes_cross_species_clean[,2] + genes_cross_species_clean[,4]
  # 
  # get overall gene expression human
  human_hypo_combined@meta.data$neurons = "non"
  human_hypo_combined@meta.data$neurons[rownames(human_hypo_combined@meta.data) %in% human_neurons] = "neuron"
  Idents(human_hypo_combined) = "neurons"
  expr_genes_human_all = AverageExpression(object = human_hypo_combined,slot = "data",features = unique(gene_table_subset$human_gene))[["RNA"]] %>%  as.data.frame() 
  expr_genes_human_all$gene = rownames(expr_genes_human_all)
  expr_genes_human_all = expr_genes_human_all %>% dplyr::select(gene,expr = neuron) %>% dplyr::mutate(expr = log1p(expr))
  
  # get overall gene expression mouse
  hypoMap@meta.data$neurons = "non"
  hypoMap@meta.data$neurons[rownames(hypoMap@meta.data) %in% hypoMap_neurons] = "neuron"
  Idents(hypoMap) = "neurons"
  expr_genes_mouse_all = AverageExpression(object = hypoMap,slot = "data",features = unique(gene_table_subset$mouse_gene))[["RNA"]] %>%  as.data.frame() 
  expr_genes_mouse_all$gene = rownames(expr_genes_mouse_all)
  expr_genes_mouse_all = expr_genes_mouse_all %>% dplyr::select(gene,expr = neuron) %>% dplyr::mutate(expr = log1p(expr))
  
  ## join
  genes_cross_species_overall = expr_genes_human_all %>% dplyr::select(gene,expr_human_all = expr) %>%
    dplyr::left_join(gene_table_subset[,c("human_gene","mouse_gene")],by = c("gene"="human_gene")) %>%
    dplyr::left_join(expr_genes_mouse_all %>% dplyr::select(gene,expr_mouse_all = expr)  ,by = c("mouse_gene"="gene"))
  
  # subtract global pct to focus on more interesting receptors ?
  genes_cross_species_clean = dplyr::left_join(genes_cross_species_clean,genes_cross_species_overall,by=c("gene"="gene","mouse_gene"="mouse_gene"))
  # genes_cross_species_clean$pct_ratio_all = log2((genes_cross_species_clean[,"pct_human_all"]+0.1) / (genes_cross_species_clean[,"pct_mouse_all"]+0.1))
  #offset = 0.1
  #genes_cross_species_clean$pct_ratio_adj = log2((genes_cross_species_clean[,target_cluster_human[1]]+offset) /(genes_cross_species_clean[,"pct_human_all"]+offset)) - log2((genes_cross_species_clean[,target_cluster_mouse[1]]+offset) / (genes_cross_species_clean[,"pct_mouse_all"]+offset))
  
  #genes_cross_species_clean$pct_ratio_adj = log2((genes_cross_species_clean[,target_cluster_human[1]] - genes_cross_species_clean[,"expr_human_all"]+1) / (genes_cross_species_clean[,target_cluster_mouse[1]] - genes_cross_species_clean[,"expr_mouse_all"]+1))
  #genes_cross_species_clean$pct_ratio_adj = log2(pmax(0,genes_cross_species_clean[,target_cluster_human[1]] - genes_cross_species_clean[,"expr_human_all"]+ 0.1) / pmax(0,genes_cross_species_clean[,target_cluster_mouse[1]] - genes_cross_species_clean[,"expr_mouse_all"] + 0.1))
  genes_cross_species_clean$pct_ratio_adj = log2(genes_cross_species_clean[,target_cluster_human[1]] / genes_cross_species_clean[,target_cluster_mouse[1]])
  
  genes_cross_species_clean = genes_cross_species_clean %>% dplyr::arrange(desc(pct_ratio_adj))
  # genes_cross_species_clean$conservation = "similar"
  # genes_cross_species_clean$conservation[genes_cross_species_clean$pct_ratio_adj > 1] = "human"
  # genes_cross_species_clean$conservation[genes_cross_species_clean$pct_ratio_adj < -1] = "mouse"
  
  # heatmap
  genes_cross_species_clean_plot = genes_cross_species_clean %>%
    tidyr::gather(key = "celltype",value="expression",-gene,-mouse_gene,-pct_ratio_adj) # ,-conservation
  genes_cross_species_clean_plot$gene =factor(genes_cross_species_clean_plot$gene,levels = unique(genes_cross_species_clean$gene))
  genes_cross_species_clean_plot$species = "mouse"
  genes_cross_species_clean_plot$species[genes_cross_species_clean_plot$celltype %in% c(target_cluster_human,"expr_human_all")] = "human"
  genes_cross_species_clean_plot$celltype[genes_cross_species_clean_plot$celltype == "expr_human_all"] = "Human - all neurons"
  genes_cross_species_clean_plot$celltype[genes_cross_species_clean_plot$celltype == "expr_mouse_all"] = "Mouse - all neurons"
  
  ## 
  human_cluster_names = edgelist[edgelist$from %in% target_cluster_human,] %>%dplyr::ungroup() %>% dplyr::distinct(from) %>% as.data.frame()
  #human_cluster_names$full_name = paste0(human_cluster_names$from,": ",human_cluster_names$cluster_name)
  genes_cross_species_clean_plot = genes_cross_species_clean_plot %>% dplyr::left_join(human_cluster_names,by=c("celltype"="from"))
  genes_cross_species_clean_plot$celltype[!is.na(genes_cross_species_clean_plot$full_name)] =  genes_cross_species_clean_plot$full_name[!is.na(genes_cross_species_clean_plot$full_name)]
  
  # if auman gene has two matching mouse genes change the label to the mouse gene to properly accommodate that
  genes_cross_species_clean_plot$gene_label = genes_cross_species_clean_plot$gene
  genes_cross_species_clean_plot$gene_label = as.character(genes_cross_species_clean_plot$gene_label)
  #print(genes_cross_species_clean_plot)
  duplicated_genes = unique(as.character(genes_cross_species_clean_plot$gene_label[as.character(genes_cross_species_clean_plot$celltype) == as.character(genes_cross_species_clean_plot$celltype)[1]])[duplicated(as.character(genes_cross_species_clean_plot$gene_label[as.character(genes_cross_species_clean_plot$celltype) == as.character(genes_cross_species_clean_plot$celltype)[1]]))] )
  genes_cross_species_clean_plot$gene_label[genes_cross_species_clean_plot$gene_label %in% duplicated_genes] = genes_cross_species_clean_plot$mouse_gene[genes_cross_species_clean_plot$gene_label %in% duplicated_genes]
  genes_cross_species_clean_plot$gene_label =factor(genes_cross_species_clean_plot$gene_label,levels = unique(genes_cross_species_clean_plot$gene_label))
  
  # helper
  addline_format <- function(x,...){gsub('\\:\\s',':\n',x)}
  # genes_cross_species_clean_plot$Percentage = round(genes_cross_species_clean_plot$pct*100,3)
  genes_cross_species_clean_plot$expression = round(genes_cross_species_clean_plot$expression,2)
  # make heatmap
  heatmap = ggplot(data = genes_cross_species_clean_plot,aes(x = gene_label, y = addline_format(celltype), fill= expression))+
    geom_tile() + geom_text(aes(label = expression),size=5,color="white")+
    #  facet_wrap(~ conservation,scales = "free",ncol = 3)+
    ggh4x::facet_wrap2( ~ species,scales = "free_y",ncol = 1,
                        strip = ggh4x::strip_themed(
                          background_x =  ggh4x::elem_list_rect(fill = c("#E69F00","#009E73")),
                          by_layer_y = TRUE
                        ))+
    ggh4x::force_panelsizes(rows = c(length(unique(genes_cross_species_clean_plot$celltype[genes_cross_species_clean_plot$species=="human"])), 
                                     length(unique(genes_cross_species_clean_plot$celltype[genes_cross_species_clean_plot$species=="mouse"])))) +
    scale_fill_viridis_c(na.value = "grey50",limits=c(0,expr_scale_max),oob=squish,name="Expr")+
    cowplot::theme_cowplot()+ylab(NULL)+xlab(NULL)+
    theme(text=element_text(size=25),axis.text.x = element_text(angle = 90,size = 18),axis.text.y = element_text(size = 12))
  
  reslist = list(
    heatmap_plot = heatmap,
    heatmap_data = genes_cross_species_clean_plot,
    gene_data = genes_cross_species_clean
  )
  
  return(reslist)
}

save_conserved_gene_expression = function(conserved_gene_expression_output,output_path,name,suffix="gpcr"){
  filename = paste0(output_path,name,"_",suffix)
  # save table
  data.table::fwrite(conserved_gene_expression_output$gene_data,paste0(filename,"_table.txt"))
  
  # heatmap extent
  rownumber = length(unique(conserved_gene_expression_output$heatmap_data$celltype))
  colnumber = length(unique(conserved_gene_expression_output$heatmap_data$gene))
  width = 100 + colnumber*20
  height = 25 + rownumber*28
  
  # save heatmap
  ggsave(filename = paste0(paste0(filename,"_heatmap.pdf")),
         plot = conserved_gene_expression_output$heatmap_plot, "pdf",dpi=300,width=width,height = height,units="mm")
}

##########
###  Run function on target cluster
##########

save_gpcr_path = "paper_figures/revision_figures/gpcr_plots/"
dir.create(save_gpcr_path,showWarnings = TRUE)

# matched_clusters_pomc = matched_clusters_final[grepl("POMC",matched_clusters_final$from),]
# matched_clusters_pomc

gpcr_scenarios = list(
  pomc_lepr = "C4-373 Mid-2 GABA-GLU-3 POMC PRDM12",
  pomc_glp1r = "C4-374 Mid-2 GABA-GLU-3 POMC CALCR" ,
  pomc_glipr1 = "C4-375 Mid-2 GABA-GLU-3 POMC ANKRD30A" ,
  agrp= "C4-355 Mid-2 GABA-GLU-1 RGS22 AGRP"
  # htr3b_crabp1
  # sst_npy_res
  # nkx24_trh
)

addgenes_pomc = c("Glp1r","Lepr")

for(i in 1:length(gpcr_scenarios)){
  
  scenario_name = names(gpcr_scenarios)[i]
  message(scenario_name)
  if(grepl("pomc",scenario_name)){
    addgenes = addgenes_pomc 
  }else{
    addgenes = c() 
  }
  message("Using addgenes: ",addgenes )
  gpcr_res = conserved_gene_expression(target_cluster_human = gpcr_scenarios[[scenario_name]],
                                       edgelist = matched_clusters_final,
                                       gene_table = gpcr_table,
                                       additional_genes = addgenes,
                                       top_n_genes = 10,
                                       human_hypo_combined = human_hypo_combined,
                                       hypoMap = hypoMap,
                                       human_neurons = human_neurons,
                                       hypoMap_neurons = hypoMap_neurons,
                                       human_marker_genes = human_marker_genes,
                                       mouse_marker_genes = hypoMap@misc$marker_genes_all,
                                       human_hypo_combined_clustercol = "C4_named",
                                       hypoMap_clustercol = "C465_named",
                                       min_spec = 1)
  
  #pomc_lepr_res$heatmap_plot
  save_conserved_gene_expression(gpcr_res,output_path = save_gpcr_path,name = scenario_name)
  
}
