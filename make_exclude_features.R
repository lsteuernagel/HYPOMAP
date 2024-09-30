##########
### Aim
##########

# Make list of genes to exclude from HVG selection

# lib
library(tidyverse)
library(Seurat)

##########
### Mouse exclude genes
##########

## Merged file: Take version from older scripts ?

## Feature exclude file:
# load short mouse version -> conbfert to human genes
mouse_exclude_genes = jsonlite::read_json("data/features_exclude_list_all2.json")
mouse_exclude_genes = sapply(mouse_exclude_genes$hvgs_exclude_long,unlist)

# convert
mouse_human = scUtils::get_mouse_human_genes_conversion()
mouse_exclude_genes_to_human = mouse_human$human_gene[mouse_human$mouse_gene %in% mouse_exclude_genes & mouse_human$hsapiens_homolog_perc_id > 0.5]

# additionally add seurat cell cycle markers
mouse_exclude_genes_to_human = unique(c(mouse_exclude_genes_to_human,unlist(Seurat::cc.genes.updated.2019)))

##########
### based on batch differences in gene expression
##########

# load data --> see merge script for creation of this!
gene_pcts_combined = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_merge/dataset_expression_genes.txt",data.table = F)

###### calculate ratios and filter baed on this
gene_pcts_combined$ratio1 = log2((gene_pcts_combined$pct_tadross+0.0001) / (gene_pcts_combined$pct_siletti+0.0001))

# I filter for genes with at least 10 occurences in one of the datasets and fold increase of 32 between tadross and siletti, or 129 between herb and any of the other two.
gene_pcts_combined_filter = gene_pcts_combined[abs(gene_pcts_combined$ratio1) > 5 & (gene_pcts_combined$occ_tadross > 10 | gene_pcts_combined$occ_siletti > 10),]

gene_dataset_differences = unique(gene_pcts_combined_filter$gene)

##########
### make complete list and save
##########

# make list with all genes to remove
all_human_exclude_genes = unique(c(mouse_exclude_genes_to_human,gene_dataset_differences))

all_human_exclude_genes = unique(c(all_human_exclude_genes,c("SELENOM","SELENOW","SELENOP")))



# and save
#scUtils::writeList_to_JSON(all_human_exclude_genes,filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_features_exclude_list_2.json")
scUtils::writeList_to_JSON(all_human_exclude_genes,filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scIntegration/data/human_features_exclude_list_2.json")
scUtils::writeList_to_JSON(all_human_exclude_genes,filename = "data/human_features_exclude_list_2.json")

##########
### also save key genes as jsons
##########

source("utility_functions.R")

# make list with all genes to remove
human_key_genes = list(keyGenes_neuron = KeyGenes, keyGenes_nonneuron = KeyGenes_nonNeuron)
# and save
scUtils::writeList_to_JSON(human_key_genes,filename = "data/human_key_genes.json")



