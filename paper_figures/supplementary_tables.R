
##########
### Load data
##########

table_dir = "paper_figures/revision_figures/suppl_tables/"
dir.create(table_dir)

# load files
human_hypo_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_HYPOMAP_publication/"
library(Seurat)
library(tidyverse)
library(ggtree)

human_hypo_combined = readRDS(paste0(human_hypo_path,"human_HYPOMAP_publication"))
combined_edgelist_mrtree = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_edgelist_mrtree_annotated.txt"),data.table = F)
comparisons_all_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_all_annotated.txt"),data.table = F)
comparisons_siblings_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_siblings_annotated.txt"),data.table = F)

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

stderror <- function(x) sd(x)/sqrt(length(x))

##########
### Prepare anno_df similar to Figure 1
##########

# make anno_df
anno_df =human_hypo_combined@misc$edgelist %>% dplyr::distinct(to,to_named) %>% dplyr::select(cluster_id = to,cluster_name = to_named)
anno_df$clusterlevel = stringr::str_extract(anno_df$cluster_id,pattern = "C[0-9]+")
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){tail(strsplit(x," ")[[1]],n=1)})
# manually shorten CO level
anno_df$first_cluster_name[anno_df$first_cluster_name == "AstroEpendymal"] = "AstroEpen"
anno_df$first_cluster_name[anno_df$first_cluster_name == "Oligodendrocytes"] = "Oligo"
anno_df$first_cluster_name[anno_df$first_cluster_name == "Immune/Vascular"] = "NN"
anno_df$first_cluster_name = stringr::str_remove(anno_df$first_cluster_name,pattern = "Oligo-")

##########
### donor overview
##########

donor_overview = data.table::fread("data/donor_overview.txt",data.table = F)

# and save in final dir
data.table::fwrite(donor_overview,paste0(table_dir,"donor_overview.txt"),sep="\t")

##########
### Per sample statistics
##########

sample_stats = human_hypo_combined@meta.data %>% dplyr::group_by(Sample_ID) %>%
  dplyr::add_count(name = "ncells") %>%
  dplyr::mutate(mean_nCount_RNA = mean(nCount_RNA),
                mean_nFeature_RNA = mean(nFeature_RNA),
                mean_percent.mt = mean(percent.mt),
                sem_nCount_RNA = stderror(nCount_RNA),
                sem_nFeature_RNA = stderror(nFeature_RNA),
                sem_percent.mt = stderror(percent.mt)) %>%
  dplyr::distinct(Sample_ID,Donor_ID,Dataset,mean_nCount_RNA,mean_nFeature_RNA,mean_percent.mt,sem_nCount_RNA,sem_nFeature_RNA,sem_percent.mt) %>%
  dplyr::mutate(facs = NA)

## add stats on probmeatic clusters using stats derived from pre-subsetting
status_meta_data_stats =  data.table::fread( "data/status_meta_data_statistics.txt",data.table = F)
# nee to map to new 

# add to sample stats
non_hypo_pct = status_meta_data_stats %>%
  dplyr::filter(celltype_status == "NonHypoNeuron") %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample_ID,pct_nonHypoNeuron = pct_status)
sample_stats = dplyr::left_join(sample_stats,non_hypo_pct,by="Sample_ID")

## save
data.table::fwrite(sample_stats,file = paste0(table_dir,"sup_sample_statistics.txt"),sep="\t")

WriteXLS::WriteXLS(x = sample_stats,ExcelFileName = paste0(table_dir,"sup_sample_statistics.xlsx"),col.names=TRUE)

# human_hypo_combined@meta.data %>% dplyr::group_by(Dataset) %>%
#   dplyr::add_count(name = "ncells") %>%
#   dplyr::summarise(mean_nCount_RNA = mean(nCount_RNA),
#                 mean_nFeature_RNA = mean(nFeature_RNA),
#                 mean_percent.mt = mean(percent.mt),
#                 sem_nCount_RNA = stderror(nCount_RNA),
#                 sem_nFeature_RNA = stderror(nFeature_RNA),
#                 sem_percent.mt = stderror(percent.mt))

##########
### human eval signatures
##########

human_celltype_signatures_2 = sapply(jsonlite::read_json("data/human_celltype_signatures_2.json"),unlist)
human_celltype_signatures_2_df = data.frame(name= names(human_celltype_signatures_2),signature_length= sapply(human_celltype_signatures_2,length),genes = sapply(human_celltype_signatures_2,paste0,collapse = ";"))
## save
data.table::fwrite(human_celltype_signatures_2_df,paste0(table_dir,"celltypes_evaluation_signatures.txt"),sep="\t")

##########
### scvi evaluation results human
##########

## load
scvi_eval_res = data.table::fread("data/complete_evaluation_results.txt",data.table = F)

## save
data.table::fwrite(scvi_eval_res,paste0(table_dir,"integration_evaluation_results.txt"),sep="\t")

##########
### Marker genes overviews
##########

## all
comparisons_all_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_all_annotated.txt"),data.table = F)
comparisons_all_updated_print = comparisons_all_updated %>% dplyr::filter(p_val_adj < 0.0001 & specificity > 1 )
comparisons_all_updated_print_topN = comparisons_all_updated_print %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(n = 100,order_by = specificity)

data.table::fwrite(comparisons_all_updated_print_topN,paste0(table_dir,"marker_genes_all_top100.txt"),sep="\t")
data.table::fwrite(comparisons_all_updated_print,paste0(table_dir,"marker_genes_all.txt"),sep="\t")

## siblings
comparisons_siblings_updated = data.table::fread(paste0(human_hypo_path,"human_hypo_combined_comparisons_siblings_annotated.txt"),data.table = F)
comparisons_siblings_updated_print = comparisons_siblings_updated %>% dplyr::filter(p_val_adj < 0.0001 & specificity > 1 )
comparisons_siblings_updated_print_topN = comparisons_siblings_updated_print %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(n = 100,order_by = specificity)

data.table::fwrite(comparisons_siblings_updated_print_topN,paste0(table_dir,"marker_genes_siblings_top100.txt"),sep="\t")
data.table::fwrite(comparisons_siblings_updated_print,paste0(table_dir,"marker_genes_siblings.txt"),sep="\t")

##########
### Table with human edgelist (of tree)
##########

## edgelist
combined_edgelist_mrtree = human_hypo_combined@misc$edgelist

data.table::fwrite(combined_edgelist_mrtree,paste0(table_dir,"edgelist_tree.txt"),sep="\t")

##########
### Table with human cluster annotations 
##########

### anno_df from above
annotation_table = anno_df

## top markers 
markers_per_cluster = comparisons_all_updated %>% 
  dplyr::filter(p_val_adj < 0.0001 &  specificity > 1 & !grepl("AC|AL|MIR|-AS|LINC",gene)) %>% 
  dplyr::select(gene,cluster,specificity) %>%
  #  dplyr::filter(cluster %in% glp1r_pct_per_cluster_top$group ) %>%
  # dplyr::filter(grepl("C6",cluster) ) %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = specificity,n = 5) %>%
  mutate(top_markers = paste0(gene, collapse = "|")) %>%
  dplyr::distinct(cluster,top_markers)

annotation_table = dplyr::left_join(annotation_table,markers_per_cluster,by=c("cluster_id"="cluster"))
annotation_table$top_markers[is.na(annotation_table$top_markers)]=""

entropy_fun = function(x,logfun ="log2"){
  log_vec = do.call(logfun,list(x))
  log_vec[is.infinite(log_vec)] = 0
  log_vec[is.nan(log_vec)] = 0
  return(-sum(x * log_vec))
}

# ncells
# make overall per cluster stats
all_cluster_levels = unique(annotation_table$clusterlevel)
general_cluster_stats_list = list()
for(current_level in all_cluster_levels){
  message(current_level)
  general_cluster_stats = human_hypo_combined@meta.data %>%
    dplyr::group_by(!!sym(current_level)) %>% 
    dplyr::mutate(mean_nCount_RNA = round(mean(nCount_RNA),2)) %>%
    dplyr::mutate(mean_nFeature_RNA = round(mean(nFeature_RNA),2)) %>%
    dplyr::mutate(mean_percent_mt = round(mean(percent.mt),4)) %>%
    dplyr::add_count(name="ncells") %>%
    dplyr::distinct(!!sym(current_level),ncells,mean_nCount_RNA,mean_nFeature_RNA,mean_percent_mt)# %>%
  target_vars = c("sex","Dataset","Donor_ID","Sample_ID")
  for(target_var in target_vars){
    #target_var = "sex"
    target_cluster_stats = human_hypo_combined@meta.data %>%
      dplyr::group_by(!!sym(current_level)) %>% 
      dplyr::add_count(name="ncells") %>%
      dplyr::group_by(!!sym(current_level),!!sym(target_var)) %>% 
      dplyr::add_count(name="cells_cluster_target") %>%
      dplyr::mutate(pct_cluster_target = round(cells_cluster_target / ncells *100,2)) #
    target_cluster_entropy = target_cluster_stats %>%
      dplyr::group_by(!!sym(current_level)) %>%
      dplyr::distinct(!!sym(current_level),pct_cluster_target) %>%
      dplyr::summarise(entropy_target = entropy_fun(pct_cluster_target / 100) / log2(length(unique(human_hypo_combined@meta.data[,target_var])))) #%>%
    colnames(target_cluster_entropy)[colnames(target_cluster_entropy) == "entropy_target"] = paste0("entropy_",target_var)
    target_cluster_wide = target_cluster_stats %>% 
      dplyr::distinct(!!sym(current_level),!!sym(target_var),pct_cluster_target) %>%
      tidyr::spread(key = !!sym(target_var),value = pct_cluster_target)
    target_cluster_final = target_cluster_wide#dplyr::left_join(target_cluster_entropy,target_cluster_wide)
    general_cluster_stats = dplyr::left_join(general_cluster_stats,target_cluster_final)
  }
  colnames(general_cluster_stats)[1] = "cluster_id"
  general_cluster_stats_list[[current_level]] = general_cluster_stats
}

general_cluster_stats_all = do.call(rbind,general_cluster_stats_list) %>% as.data.frame()
general_cluster_stats_all[is.na(general_cluster_stats_all)] = 0
annotation_table = dplyr::left_join(annotation_table,general_cluster_stats_all,by="cluster_id")

data.table::fwrite(annotation_table,paste0(table_dir,"cluster_annotations_and_statistics.txt"),sep="\t")

##########
### Table with mapping to herb
##########

#
herb_fetal_adult_hypo_meta = data.table::fread("siletti_processing/herb_fetal_adult_hypo_meta_table.txt",data.table = F)
#head(herb_fetal_adult_hypo_meta[herb_fetal_adult_hypo_meta$H108==31,] %>% na.omit())
neuron_metadata = human_hypo_combined@meta.data
colnames(neuron_metadata)
herb_fetal_adult_hypo_meta$V1 = stringr::str_replace(herb_fetal_adult_hypo_meta$V1,":","_")
neuron_metadata_wHerb = dplyr::left_join(neuron_metadata %>% dplyr::select(Cell_ID,Dataset,Sample_ID,Donor_ID,"C0","C1","C2","C3","C4","C0_named" ,"C1_named" ,"C2_named"  ,"C3_named", "C4_named"),
                                         herb_fetal_adult_hypo_meta %>% 
                                           dplyr::select(V1,H6,H53,H108,H369,H108_Class,H108_Nuclei,H108_MarkerGenes,H369_Class,H369_Nuclei,H369_MarkerGenes),
                                         by=c("Cell_ID"="V1"))
neuron_metadata_wHerb$H369 = as.character(neuron_metadata_wHerb$H369)
neuron_metadata_wHerb$H108 = as.character(neuron_metadata_wHerb$H108)
neuron_metadata_wHerb$H53 = as.character(neuron_metadata_wHerb$H53)
neuron_metadata_wHerb$H6 = as.character(neuron_metadata_wHerb$H6)

## make stats
c4_h369_stats = neuron_metadata_wHerb %>% 
  dplyr::filter(!grepl("C3-",C4)) %>%
  dplyr::group_by(C4) %>%
  dplyr::add_count(name="cells_C4_total") %>%
  dplyr::filter(Dataset=="Siletti" & !is.na(H369)) %>%
  dplyr::add_count(name="cells_C4") %>%
  dplyr::group_by(C4,H369) %>%
  dplyr::add_count(name="cells_C4_H369") %>% 
  dplyr::distinct(C4,H369,cells_C4_total,cells_C4,cells_C4_H369) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pct_C4_H369 = round(cells_C4_H369 / cells_C4 *100,2), pct_siletti = round(cells_C4 / cells_C4_total *100,2)) %>%
  dplyr::filter(pct_siletti >= 5) %>%
  dplyr::group_by(H369) %>%
  dplyr::arrange(H369,desc(pct_C4_H369)) %>% dplyr::mutate(cluster_order_H369 = row_number()) %>%
  dplyr::group_by(C4) %>%
  dplyr::arrange(C4,desc(pct_C4_H369)) %>% dplyr::mutate(cluster_order = row_number()) %>%
  dplyr::filter((cluster_order == 1 | pct_C4_H369 >= 25 | cluster_order_H369 == 1) & pct_siletti >= 5) %>%
  dplyr::select(Tadross_cluster = C4,Herb_cluster = H369,ncells_C_total = cells_C4_total, ncells_C_siletti = cells_C4, Percentage_H_Cluster = pct_C4_H369 )
#  dplyr::slice_max(order_by = pct_C4_H369,n = 3) #%>%

c3_H108_stats = neuron_metadata_wHerb %>% 
 # dplyr::filter(!grepl("C3-",C3)) %>%
  dplyr::group_by(C3) %>%
  dplyr::add_count(name="cells_C3_total") %>%
  dplyr::filter(Dataset=="Siletti" & !is.na(H108)) %>%
  dplyr::add_count(name="cells_C3") %>%
  dplyr::group_by(C3,H108) %>%
  dplyr::add_count(name="cells_C3_H108") %>% 
  dplyr::distinct(C3,H108,cells_C3_total,cells_C3,cells_C3_H108) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pct_C3_H108 = round(cells_C3_H108 / cells_C3 *100,2), pct_siletti = round(cells_C3 / cells_C3_total *100,2)) %>%
  dplyr::filter(pct_siletti >= 5) %>%
  dplyr::group_by(H108) %>%
  dplyr::arrange(H108,desc(pct_C3_H108)) %>% dplyr::mutate(cluster_order_h108 = row_number()) %>%
  dplyr::group_by(C3) %>%
  dplyr::arrange(C3,desc(pct_C3_H108)) %>% dplyr::mutate(cluster_order = row_number()) %>%
  dplyr::filter((cluster_order == 1 | cluster_order_h108 == 1 | pct_C3_H108 >= 25) & pct_siletti >= 5) %>%
  dplyr::select(Tadross_cluster = C3,Herb_cluster = H108,ncells_C_total = cells_C3_total, ncells_C_siletti = cells_C3, Percentage_H_Cluster = pct_C3_H108 )

C2_H53_stats = neuron_metadata_wHerb %>% 
  # dplyr::filter(!grepl("C2-",C2)) %>%
  dplyr::group_by(C2) %>%
  dplyr::add_count(name="cells_C2_total") %>%
  dplyr::filter(Dataset=="Siletti" & !is.na(H53)) %>%
  dplyr::add_count(name="cells_C2") %>%
  dplyr::group_by(C2,H53) %>%
  dplyr::add_count(name="cells_C2_H53") %>% 
  dplyr::distinct(C2,H53,cells_C2_total,cells_C2,cells_C2_H53) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pct_C2_H53 = round(cells_C2_H53 / cells_C2 *100,2), pct_siletti = round(cells_C2 / cells_C2_total *100,2)) %>%
  dplyr::filter(pct_siletti >= 5) %>%
  dplyr::group_by(H53) %>%
  dplyr::arrange(H53,desc(pct_C2_H53)) %>% dplyr::mutate(cluster_order_H53 = row_number()) %>%
  dplyr::group_by(C2) %>%
  dplyr::arrange(C2,desc(pct_C2_H53)) %>% dplyr::mutate(cluster_order = row_number()) %>%
  dplyr::filter((cluster_order == 1 | cluster_order_H53 == 1 | pct_C2_H53 >= 25) & pct_siletti >= 5) %>%
  dplyr::select(Tadross_cluster = C2,Herb_cluster = H53,ncells_C_total = cells_C2_total, ncells_C_siletti = cells_C2, Percentage_H_Cluster = pct_C2_H53 )


# double check:
length(unique(c3_H108_stats$Herb_cluster))
setdiff(neuron_metadata_wHerb$H108,c3_H108_stats$Herb_cluster)

# C4 & C3 & C2 cluster comparison:
herb_cluster_comparison = rbind(C2_H53_stats,rbind(c3_H108_stats,c4_h369_stats))

data.table::fwrite(herb_cluster_comparison,paste0(table_dir,"herb_cluster_comparison.txt"),sep="\t")

##########
### Table with human - mouse mapping
##########

human_mouse_edgelist = data.table::fread("merge_human_mouse_neurons/matched_clusters_scviCorPruned.tsv",data.table = F)
human_mouse_edgelist_print = human_mouse_edgelist %>% dplyr::select(human_cluster = from, 
                                                                    mouse_cluster = to, 
                                                                    human_parent = from_parent, 
                                                                    mouse_parent = to_parent, 
                                                                    similarity, 
                                                                    human_count = from_occ,
                                                                    mouse_count = to_occ  )

data.table::fwrite(human_mouse_edgelist_print,paste0(table_dir,"human_mouse_corresponding_clusters.txt"),sep="\t")

##########
### obesity cell type tables
##########

overview_MCR = data.table::fread("paper_figures/revision_figures/melanocortin_overview_table_filtered.txt",data.table = F)
data.table::fwrite(overview_MCR,paste0(table_dir,"melanocortin_celltypes.txt"),sep="\t")

overview_incretin = data.table::fread("paper_figures/revision_figures/incretin_overview_table_both_filtered.txt",data.table = F)
# overview_incretin2 = data.table::fread("paper_figures/revision_figures/incretin_overview_table_filtered.txt",data.table = F)
# overview_incretin2 = rbind(overview_incretin,overview_incretin2)
data.table::fwrite(overview_incretin,paste0(table_dir,"incretin_celltypes.txt"),sep="\t")

##########
### Load genetics
##########

## files included from genetics

# magma gwas
genetics_magma = data.table::fread(paste0(table_dir,"genetics_magma.txt"),data.table = F)

# effector genes
genetics_effector = data.table::fread(paste0(table_dir,"genetics_effector.txt"),data.table = F)

# effector genes
genetics_exome = data.table::fread(paste0(table_dir,"genetics_exome.txt"),data.table = F)

##########
### combined result
##########

supp_table_list = list(
  "1_donor_overview" = donor_overview %>% data.frame(), #donor_overview %>% as.data.frame(),
  "2_sample_overview" = sample_stats %>% as.data.frame(),
  "3_evaluation_celltypes" = human_celltype_signatures_2_df %>% as.data.frame(),
  "4_evaluation_results" = scvi_eval_res %>% as.data.frame(),
  "5_cluster_tree_edgelist" = combined_edgelist_mrtree %>% as.data.frame(),
  "6_cluster_tree_annotations" = annotation_table %>% as.data.frame(),
  "7_cluster_markers_global" = comparisons_all_updated_print_topN %>% as.data.frame(),
  "8_cluster_markers_siblings" = comparisons_siblings_updated_print_topN %>% as.data.frame(),
  "9_comparison_clusters_Herb" = herb_cluster_comparison %>% as.data.frame(),
  "10_human_mouse_corresponding" = human_mouse_edgelist_print %>% as.data.frame(),
  "11_region_assignment" = data.frame(),
  "12_region_abundances" = data.frame(),
  "13_melanocortin_celltypes" = overview_MCR %>% as.data.frame(),
  "14_incretin_celltypes" = overview_incretin %>% as.data.frame(),
  "15_ST_Magma" = genetics_magma %>% as.data.frame(),
  "16_ST_EffectorGenes" = genetics_effector %>% as.data.frame(),
  "17_ST_Exome_SignalLookUp" = genetics_exome %>% as.data.frame()
  # TODO: add all genetics tables
)

# make dictionary 
dictionary_list = lapply(names(supp_table_list),function(x,tlist){data.frame(table = c(x,colnames(tlist[[x]])))},tlist=supp_table_list)
#dictionary_list = dictionary_list[base::order(as.numeric(stringr::str_extract(names(dictionary_list),pattern = "[0-9]+")))]
dictionary_df = do.call("rbind",dictionary_list)
supp_table_list[["0_dictionary"]] = dictionary_df
# order numerically
supp_table_list = supp_table_list[base::order(as.numeric(stringr::str_extract(names(supp_table_list),pattern = "[0-9]+")))]
names(supp_table_list)

# shorten colnames
names(supp_table_list)[sapply(names(supp_table_list),nchar) > 31] = sapply(names(supp_table_list)[sapply(names(supp_table_list),nchar) > 31],substr,start=1,stop=31)
#write
WriteXLS::WriteXLS(x = supp_table_list,ExcelFileName = paste0(table_dir,"supplementary_tables_preliminary_revision.xlsx"),col.names=TRUE)




