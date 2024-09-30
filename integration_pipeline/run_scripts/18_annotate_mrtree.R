##########
### Load parameters and packages
##########

message(Sys.time(),": Starting combined tree annotation .." )

message(" Load parameters and packages ")

require(tidyverse)
require(Seurat)
require(Matrix)

source("integration_pipeline/harmonization_functions.R")
source("integration_pipeline/mrtree_functions.R")

command_args<-commandArgs(TRUE)
param_file = command_args[1]
#param_file = "integration_pipeline/parameters/parameters_human_hypo_v1.json"
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### Load combined
##########

message(Sys.time(),": Load data .." )

subname = parameter_list$subname

# load subset specifc
parameter_list$harmonization_folder_path = paste0(parameter_list$harmonization_folder_path,"/",subname,"/")
parameter_list$new_name_suffix = subname

plot_path = paste0(parameter_list$harmonization_folder_path,"annotation_plots/")
dir.create(plot_path)

# read seurat:
human_hypo_combined = readRDS(paste0(parameter_list$harmonization_folder_path,subname,".rds"))
# other data
human_hypo_combined_edgelist = data.table::fread(paste0(parameter_list$harmonization_folder_path,subname,"_","edgelist_mrtree",".txt"),data.table = F)
#human_hypo_combined_edgelist = rbind(data.frame(from = "all",to=unique(human_hypo_combined_edgelist$from[grepl("C0-",human_hypo_combined_edgelist$from)]),isLeaf=FALSE,level),human_hypo_combined_edgelist)

human_hypo_combined_markers_all = data.table::fread(paste0(parameter_list$harmonization_folder_path,subname,"_","comparisons_all_updated",".txt"),data.table = F)
human_hypo_combined_markers_sibling = data.table::fread(paste0(parameter_list$harmonization_folder_path,subname,"_","comparisons_siblings_updated",".txt"),data.table = F)

mrtree_result_plotting_labelmat = human_hypo_combined@meta.data[,grepl("C[0-9]+",colnames(human_hypo_combined@meta.data))]

##########
### Define manual annotation : C1
##########

# top level non neuronal common names and basic neurons by main AP region
pre_annotated = c( "C0-1" = "AstroEpendymal",
                   "C0-2" = "Oligodendrocytes",
                   "C0-3" = "Neuron", 
                   "C0-4" = "Immune/Vascular", 
                   "C1-1" = "Astrocytes",
                   "C1-2" = "Ependymal",
                   "C1-3" = "Oligo-Precursor",
                   "C1-4" = "Oligo-Mature",
                   "C1-5" = "Pre-1",
                   "C1-6" = "Mid-1",
                   "C1-7" = "Pre-2",
                   "C1-8" = "Post-1",
                   "C1-9" = "Mid-3",
                   "C1-10" = "Mid-2",
                   "C1-11" = "Post-2",
                   "C1-12" = "Vascular",
                   "C1-13" = "Immune")

# Lower level non neuronal with common names
pre_annotated = c(pre_annotated,
                  "C3-141" = "Endothelial",
                  "C3-142" = "Fibroblasts" ,
                  #"C3-143" = #"Endothelial", C-3-143 could be astrocytes ?? or neutrophils ? 
                  "C3-144" = "Pericytes",
                  "C3-145" = "SMCs",
                  "C3-12" = "Ependymocytes",
                  "C3-13" = "Tanycytes",
                  "C3-14" = "Choroid", # check for Choroid plexus via Folr1, PRLR
                  "C2-52" = "Tcells",
                  "C2-51" = "Myocytes", #??
                  "C2-50" = "Macrophages",
                  "C2-49" = "Microglia"
)

# define a few C3 and C4 levels with common neurtransmitter genes, which are not first in the ranking but slightly below and are more helpful 
# --> see the C3 and C4 marker results
manual_C34_names = c(
  #  "C3-140" = "CD36", - use CD36 instead of Sim1 ? not relevant
  "C3-135" = "HDC", # overwrite HDC neurons, 
  "C4-355" = "AGRP", # these are the canonical AGRP/NPY neurons which have WIF1 with a slightly higher score than AGRP due to relatively low overall AGRP expression
  "C4-370"= "PMCH", # these are MCH neurons with a slightly higher ISL1 score than PMCH
  "C4-147" = "NPW" # these are the LH NPW neurons which are alternatively labeled using EBF1 (also in mouse data usually as a EBF1+ cluster --> sister-cluster of GLP1R/LEPR/TRH+ DMH cells 
) 


# get results / tranfrom to same df format as other levels
c1_level_annotation = data.frame(cluster = names(pre_annotated), anno_gene = pre_annotated)

##########
### Neurotransmitter status: C2
##########

message(Sys.time(),": Run annotation C2 .." )

neuron_celltypes = scUtils::find_children(nodes = c("C0-3"),human_hypo_combined_edgelist)

# construct annotation with neurortransmitters
thresh=2
neurotransmitter_expression_C2 = scUtils::gene_pct_cluster(human_hypo_combined,col_name = "C2",genes = c("SLC32A1","SLC17A6","SLC18A2","CHAT"),return_long = F) 
neurotransmitter_expression_C2$cluster = rownames(neurotransmitter_expression_C2)
neurotransmitter_expression_C2$ratio = log2((neurotransmitter_expression_C2$SLC32A1+0.01) / (neurotransmitter_expression_C2$SLC17A6+0.01))
#neurotransmitter_expression_C2 = neurotransmitter_expression_C2 %>% dplyr::left_join(human_hypo_combined_edgelist[,c("to","from")],by=c("cluster"="to"))
neurotransmitter_expression_C2$type = "NonNeuron"
neurotransmitter_expression_C2$type[neurotransmitter_expression_C2$cluster %in% neuron_celltypes] = "Neuron"
neurotransmitter_expression_C2$neuro_type = "GABA-GLU"
neurotransmitter_expression_C2$neuro_type[neurotransmitter_expression_C2$type == "NonNeuron"] = ""
neurotransmitter_expression_C2$neuro_type[neurotransmitter_expression_C2$ratio > thresh] = "GABA"
neurotransmitter_expression_C2$neuro_type[neurotransmitter_expression_C2$ratio < -1*thresh] = "GLU"
neurotransmitter_expression_C2$neuro_type[neurotransmitter_expression_C2$CHAT > neurotransmitter_expression_C2$SLC32A1 & neurotransmitter_expression_C2$CHAT > neurotransmitter_expression_C2$SLC17A6] = "CHOL"

## add full ano
neurotransmitter_anno_C2 = neurotransmitter_expression_C2 %>%
  dplyr::left_join(human_hypo_combined_edgelist[,c("to","from")],by=c("cluster"="to")) %>%
  dplyr::group_by(from,neuro_type) %>%
  dplyr::mutate(c2_anno = paste0(neuro_type,"-",seq_len(n())))

# ggplot(neurotransmitter_expression_C2,aes(x=ratio,fill=neuro_type,color=neuro_type))+geom_density(alpha=0.4)+
#   geom_point(y=0.1)+
#   facet_wrap(~type,nrow = 2)+
#   theme_light()+xlab("log2 SLC232A1 / SLC17A6 percentage")

c2_level_annotation = data.frame(cluster = neurotransmitter_anno_C2$cluster, anno_gene = neurotransmitter_anno_C2$c2_anno)

##########
### Run annotation
##########

manual_exclude_genes =c("ZNF385B","PLCG2","ELAVL2","SLC44A5","GALNTL6, MGAT4C")

message(Sys.time(),": Run annotation C3 .." )

annotation_result_C3 = annotate_level(edgelist = human_hypo_combined_edgelist, #[!grepl("C5|C6",human_hypo_combined_edgelist$to),],
                                      labelmat = mrtree_result_plotting_labelmat,
                                      level = "C3",
                                      markers_comparisons_all = human_hypo_combined_markers_all,
                                      markers_comparisons_siblings = human_hypo_combined_markers_sibling,
                                      manual_names=manual_C34_names ,
                                      overwrite_with_manual= TRUE, #FALSE,
                                      manual_exclude_genes=c(manual_exclude_genes),
                                      max_pval_adj=0.0001, 
                                      min_combined_score = 1.5,
                                      lower_pct_min=0.01,
                                      pct.2_max = 0.3,
                                      min_cells=20,
                                      max_score_siblings_children = 5,
                                      reverse_order = FALSE)
c3_level_allMarkers = do.call(rbind,annotation_result_C3$all_combined_markers_full) %>% as.data.frame()
c3_level_annotation = data.frame(cluster = names(annotation_result_C3$annotated_nodes_gene), anno_gene = unlist(annotation_result_C3$annotated_nodes_gene))

## manually add missing annos
message("Found ",length(unique(c3_level_annotation$cluster[c3_level_annotation$anno_gene=="NA"]))," un-annotated C3 cluster.")
for(cl in c3_level_annotation$cluster[c3_level_annotation$anno_gene=="NA"]){
  all_markers = c3_level_allMarkers[c3_level_allMarkers$cluster %in% cl & !c3_level_allMarkers$gene %in% manual_exclude_genes & 
                                      c3_level_allMarkers$pct.2 <= 0.6 & c3_level_allMarkers$combined_score > 1,] %>%
    dplyr::filter(! grepl("LINC|AL[0-9]{2,}|AP[0-9]{2,}|AC[0-9]{2,}|MIR|-AS|-PS|orf[0-9]{2,}|[0-9]{2,}\\.[0-9]",gene)) %>%
    dplyr::slice_max(order_by = combined_score,with_ties = F,n = 1)
  message("Using ",all_markers$gene[1]," for ",cl)
  c3_level_annotation$anno_gene[c3_level_annotation$cluster == cl] = all_markers$gene[1]
}

##########
### Run annotation
##########

message(Sys.time(),": Run annotation C4 .." )

annotation_result_C4 = annotate_level(edgelist = human_hypo_combined_edgelist, #[!grepl("C5|C6",human_hypo_combined_edgelist$to),],
                                      labelmat = mrtree_result_plotting_labelmat,
                                      level = "C4",
                                      markers_comparisons_all = human_hypo_combined_markers_all,
                                      markers_comparisons_siblings = human_hypo_combined_markers_sibling,
                                      manual_names= manual_C34_names,
                                      parent_names = annotation_result_C3$annotated_nodes_gene,
                                      overwrite_with_manual= TRUE, #FALSE,
                                      manual_exclude_genes=c(manual_exclude_genes),
                                      max_pval_adj=0.0001, 
                                      min_combined_score = 1.5,
                                      lower_pct_min=0.01,
                                      pct.2_max = 0.3,
                                      min_cells=20,
                                      max_score_siblings_children = 5,
                                      max_score_sibling = 20,
                                      max_score_self = 50,
                                      reverse_order = FALSE)

# get results
c4_level_allMarkers = do.call(rbind,annotation_result_C4$all_combined_markers_full) %>% as.data.frame()
c4_level_annotation = data.frame(cluster = names(annotation_result_C4$annotated_nodes_gene), anno_gene = unlist(annotation_result_C4$annotated_nodes_gene))

## manually add missing annos
message("Found ",length(unique(c4_level_annotation$cluster[c4_level_annotation$anno_gene=="NA"]))," un-annotated c4 cluster.")
for(cl in c4_level_annotation$cluster[c4_level_annotation$anno_gene=="NA"]){
  all_markers = c4_level_allMarkers[c4_level_allMarkers$cluster %in% cl & !c4_level_allMarkers$gene %in% manual_exclude_genes & 
                                      c4_level_allMarkers$pct.2 <= 0.6 & c4_level_allMarkers$combined_score > 1,] %>%
    dplyr::filter(! grepl("LINC|AL[0-9]{3,}|AP[0-9]{3,}|AC[0-9]{3,}|MIR|-AS|-PS|orf[0-9]{2,}|[0-9]{2,}\\.[0-9]",gene)) %>%
    dplyr::slice_max(order_by = combined_score,with_ties = F,n = 1)
  message("Using ",all_markers$gene[1]," for ",cl)
  c4_level_annotation$anno_gene[c4_level_annotation$cluster == cl] = all_markers$gene[1]
}

##########
### Build final annotation per cluster
##########

message(Sys.time(),": Combine and save .." )

## conactenate to names:

annotated_labelmat = mrtree_result_plotting_labelmat  %>%
  dplyr::left_join(c1_level_annotation %>% dplyr::rename(C0_name = anno_gene),by=c("C0"="cluster")) %>%
  dplyr::left_join(c1_level_annotation %>% dplyr::rename(C1_name = anno_gene),by=c("C1"="cluster")) %>%
  dplyr::left_join(c2_level_annotation %>% dplyr::rename(C2_name = anno_gene),by=c("C2"="cluster")) %>%
  dplyr::left_join(c3_level_annotation %>% dplyr::rename(C3_name = anno_gene),by=c("C3"="cluster")) %>%
  dplyr::left_join(c4_level_annotation %>% dplyr::rename(C4_name = anno_gene),by=c("C4"="cluster")) %>%
  dplyr::distinct()

# substitute pre annotated
for(i in 1:length(pre_annotated)){
  current_node = names(pre_annotated[i])
  current_level = stringr::str_extract(names(pre_annotated[i]),pattern = "C[0-9]")
  new_anno = pre_annotated[i]
  if(current_level != "C0"){
    annotated_labelmat[ annotated_labelmat[,current_level] == current_node,paste0(current_level,"_name")] = new_anno
  }
}
# get all nodes with no single
single_child_nodes = human_hypo_combined_edgelist$to[human_hypo_combined_edgelist$count <= 1 & !human_hypo_combined_edgelist$isLeaf]
children_of_single_child_nodes = human_hypo_combined_edgelist$to[human_hypo_combined_edgelist$from %in% single_child_nodes]
# delete annos if no siblings
for(col in 1:4){
  colname = paste0("C",col)
  annotated_labelmat[annotated_labelmat[,colname] %in% children_of_single_child_nodes,paste0("C",col,"_name")] =NA
}

## add full names:
for(col in 0:4){
  colname = paste0("C",col)
  if(col==0){
    annotated_labelmat[,paste0(colname,"_named")] = paste0(annotated_labelmat[,colname]," ",annotated_labelmat[,paste0(colname,"_name")])
  }else{
    annotated_labelmat[,paste0(colname,"_named")] = apply(annotated_labelmat,1,function(row,col){paste0(row[colname]," ",paste0(row[paste0("C",1:col,"_name")],collapse = " "))},col=col)
  }
}
# put C2 Non neurons together
for(col in 2:4){
  output_string <- gsub("([a-z]+)\\s+-([0-9]+)", "\\1-\\2", annotated_labelmat[grepl("[a-z]+ \\-[0-9]+",annotated_labelmat[,paste0("C",col,"_named")]),paste0("C",col,"_named")])
  annotated_labelmat[grepl("[a-z]+ -[0-9]+",annotated_labelmat[,paste0("C",col,"_named")]),paste0("C",col,"_named")]  = output_string
}
# delete NAs
for(col in 1:4){
  colname = paste0("C",col)
  # make sure to not gsub genes that start with NA
  annotated_labelmat[,paste0(colname,"_named")] = gsub(" NA(?!\\p{L})", "", annotated_labelmat[,paste0(colname,"_named")], perl = TRUE)
  # annotated_labelmat[,paste0(colname,"_named")] = gsub(" NA","",annotated_labelmat[,paste0(colname,"_named")])
}

##########
### Save results
##########

# save all results
all_results_list = list(
  annotated_labelmat = annotated_labelmat,
  c1_level_annotation = c1_level_annotation,
  c2_level_annotation = c2_level_annotation,
  c3_level_annotation = c3_level_annotation,
  c4_level_annotation = c4_level_annotation,
  c3_level_allMarkers = c3_level_allMarkers,
  c4_level_allMarkers = c4_level_allMarkers
)
saveRDS(all_results_list,paste0(parameter_list$harmonization_folder_path,subname,"_","cluster_annotation_allData",".rds"))

# save final annotated labelmat
annotated_labelmat_final = annotated_labelmat %>% dplyr::select(C0,C1,C2,C3,C4,C0_named,C1_named,C2_named,C3_named,C4_named) %>% as.data.frame()
data.table::fwrite(annotated_labelmat_final,file=paste0(parameter_list$harmonization_folder_path,subname,"_","cluster_annotation",".txt"),sep = "\t")

message(Sys.time(),":Complete." )

##########
### Combine and do some exploration
##########


#### some stats
# 
# neuron_celltypes = scUtils::find_children(nodes = c("C0-3"),human_hypo_combined_edgelist)
# cellsh=human_hypo_combined@meta.data$Cell_ID[human_hypo_combined@meta.data$C4 %in% neuron_celltypes]
# neuron_celltypes = unique((neuron_celltypes[grepl("C4",neuron_celltypes)]))
# Idents(human_hypo_combined) = "C4"
# genes_to_plot = unique(a2$anno_gene[a2$anno_gene %in% rownames(human_hypo_combined) & a2$cluster %in% neuron_celltypes])
# dp= Seurat::DotPlot(object = human_hypo_combined,features = genes_to_plot,group.by = "C3",idents = neuron_celltypes,scale = FALSE)
# dp+theme(axis.text.x = element_text(angle=90,size=7),axis.text.y = element_text(size=6))
# 
# dp_data = dp$data
# dp_data_gene = dp_data %>% dplyr::group_by(features.plot) %>% dplyr::summarise(mean_pct = mean(pct.exp)) #%>% dplyr::count()#
# gene_count = a2 %>% dplyr::group_by(anno_gene) %>% dplyr::count()
# 
# # combined_expr_anno
# 
# c4_genes = scUtils::gene_pct_cluster(human_hypo_combined,genes = annotated_labelmat$C4_name,col_name = "C4",return_long = TRUE)
# c4_genes1 = c4_genes %>% dplyr::rename(pct_c4 = pct, gene_c4=gene)
# c3_genes = scUtils::gene_pct_cluster(human_hypo_combined,genes = annotated_labelmat$C3_name,col_name = "C4",return_long = TRUE)
# c3_genes1 = c3_genes %>% dplyr::rename(pct_c3 = pct, gene_c3 = gene)
# 
# all_valid = paste0(annotated_labelmat$C3_name,"_",annotated_labelmat$C4_name)
# genes_both = c3_genes1 %>%
#   dplyr::left_join(c4_genes1,by=c("group"="group")) %>%
#   dplyr::mutate(concat = paste0(gene_c3,"_",gene_c4)) %>%
#   dplyr::filter(concat %in% all_valid)  %>%
#   dplyr::mutate(pct_product = pct_c3 * pct_c4)
# 
# 
# ggplot(genes_both,aes(x=concat,y=group,fill=pct_product))+
#   geom_tile()+scale_fill_gradient2(low="white",mid="orange",high="red",midpoint = 0.333)

##########
### Load parameters and packages
##########

# 
# # load test anno
# annores_list = annotation_result
# # load desc markers
# all_combined_markers_list = annores_list$all_combined_markers
# for(i in 1:length(all_combined_markers_list)){
#   if(nrow(all_combined_markers_list[[i]])>0){
#     all_combined_markers_list[[i]]$cluster = names(all_combined_markers_list)[i]
#   }
# }
# all_combined_markers_df = as.data.frame(do.call(rbind,all_combined_markers_list))
# 
# # load desc markers -full
# all_combined_markers_full = annores_list$all_combined_markers_full
# for(i in 1:length(all_combined_markers_full)){
#   if(nrow(all_combined_markers_full[[i]])>0){
#     all_combined_markers_full[[i]]$cluster = names(all_combined_markers_full)[i]
#   }
# }
# all_combined_markers_full_df = as.data.frame(do.call(rbind,all_combined_markers_full))
# 
# # amke a df with all genes
# annotated_nodes_gene = annores_list$annotated_nodes_gene
# annotated_nodes_gene_df = data.frame(cluster = names(annotated_nodes_gene),genes_conc = sapply(annotated_nodes_gene,function(x){paste0(x,collapse = "_")}))
# 
# annotated_nodes_gene_df
# 
# FeaturePlot(human_hypo_combined,"HS3ST4",order = TRUE,reduction = "umap_scvi_hypo")

##########
### Run annotation
##########

# theta <- acos( sum(A*B) / ( sqrt(sum(A ^ 2)) * sqrt(sum(B ^ 2)) ) )
# rad2deg <- function(rad) {(rad * 180) / (pi)}
# deg2rad <- function(deg) {(deg * pi) / (180)}

# # expects both numbers to be positive (so that max agle can be 90 degree (1.570796 rad)
# angular_similarity = function(A){
#   x1 = acos(sum(A*c(0,1)) / sqrt(sum(A^2)*sum(c(0,1)^2)))
#   x2 = acos(sum(A*c(1,0)) / sqrt(sum(A^2)*sum(c(1,0)^2)))
#   angular_sim = 1 - round((max(x1,x2)-0.7853982) / (1.570796-0.7853982),4) # -90 degree - 45 degree 
#   if(angular_sim < 0){angular_sim=0}
#   return(angular_sim)
# }
# 
# #a1 = do.call(rbind,annotation_result$all_combined_markers) %>% as.data.frame() %>% dplyr::filter(scor 
# a1$score2 = (apply(a1[,c("score_all_self","score_siblings_self")],1,angular_similarity)*(a1$score_all_self+a1$score_siblings_self)) - abs(a1$score_siblings_children)

