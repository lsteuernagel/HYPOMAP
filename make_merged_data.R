
# sbatch run_slurm.sh ~/Documents/r_scvi_v3_42.simg make_merged_data.R

##########
### Aim
##########

# Make merged dataset

# load list of genes to update:

##########
### data path
##########

# tadross
# see here: opts$data_path
tadross_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/"

# siletti
siletti_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/linnarson_human_data/"

# output path
output_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/human_hypothalamus_merge/"
system(paste0("mkdir -p ",output_path))

# lib
library(tidyverse)
library(Seurat)

##########
### Load pre-processed Tadross data
##########

message("-----",Sys.time(),": Load Tadross data ")

tadross_processed = readRDS(paste0(tadross_data_path ,"HumanNucSeq_processed_filtered.rds"))

##########
### Curate Tadross meta data
##########

colnames(tadross_processed@meta.data)
tadross_meta_data = tadross_processed@meta.data %>% dplyr::select(Cell_ID,nCount_RNA,nFeature_RNA,Sample_ID = sample, Donor_ID = donor, Batch_ID = batch, sex , age_years = age,
                                                                  percent.mt, InitialCluster = SNN_scvi_res.4, InitalClass = updated_annotation, DoubletClass = scDblFinder.class)
tadross_meta_data$Technology = "10x v3"
tadross_meta_data$Donor_ID = paste0("tadross_",tadross_meta_data$Donor_ID)
tadross_meta_data$Sample_ID = paste0("tadross_",tadross_meta_data$Sample_ID)
tadross_meta_data$Dataset = "Tadross"
tadross_meta_data$Batch_ID = paste0("tadross_",tadross_meta_data$Batch_ID)

rownames(tadross_meta_data) = tadross_meta_data$Cell_ID
tadross_processed@meta.data = tadross_meta_data

##########
### Load pre-processed Siletti data
##########

message("-----",Sys.time(),": Load Siletti data ")

# "raw" but already processed by the authors!
siletti_processed = readRDS(paste0(siletti_data_path,"linnarson_hypothalamus_raw.rds"))

##########
### Curate Siletti meta data
##########

siletti_meta_data = siletti_processed@meta.data %>% dplyr::select(Cell_ID = CellID,nCount_RNA,nFeature_RNA,Sample_ID = SampleID, Donor_ID  = Donor, age_years = Age,
                                                                  percent.mt = MT_ratio, InitialCluster = Clusters, DoubletClass = DoubletFinderFlag,Technology = Chemistry,unspliced_ratio, Roi, Tissue )
#siletti_meta_data$Technology = "10x v3"

siletti_samples_other_information = readxl::read_xlsx("siletti_processing/Stiletti_samples.xlsx")
siletti_meta_data = dplyr::left_join(siletti_meta_data,siletti_samples_other_information %>% dplyr::select(sex = Sex,info_sample_id = `10X Sample ID`,Sample_Tube = `Sample Tube`),by=c("Sample_ID"="info_sample_id"))

siletti_meta_data$Donor_ID = paste0("siletti_",siletti_meta_data$Donor_ID)
siletti_meta_data$Sample_ID = paste0("siletti_",siletti_meta_data$Sample_ID)
siletti_meta_data$Dataset = "Siletti"
siletti_meta_data$Technology = paste0("10x ",siletti_meta_data$Technology)
siletti_meta_data$Batch_ID = paste0("siletti_",gsub("M[0-9]+TX_","",siletti_meta_data$Sample_Tube))

siletti_meta_data$sex[siletti_meta_data$sex == "M"] = "Male"
siletti_meta_data$sex[siletti_meta_data$sex == "F"] = "Female"

rownames(siletti_meta_data) = siletti_meta_data$Cell_ID
siletti_processed@meta.data = siletti_meta_data

##########
### Find & rename changed gene symbols
##########

message("-----",Sys.time(),": Check gene reference differences between datasets ")

# get all Tadross genes via biomart --> standard 10x:
human_mart_ens98 = biomaRt::useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl",host="https://sep2019.archive.ensembl.org/")
human_genes_98 = biomaRt::getBM(mart = human_mart_ens98,attributes = c("ensembl_gene_id","external_gene_name"))
ensembl_98_genes = unique(human_genes_98$external_gene_name)
tadross_genes = ensembl_98_genes

# siletti:
# The reference genome and transcript annotations were based on the human GRCh38.p13 gencode V35 primary sequence assembly. 
# However, we filtered the reference. Because our pipeline only counted reads that uniquely aligned to one gene, reads that aligned to more than one gene were lost.
# Altogether we filtered 387 fusion genes, 1140 overlapping transcripts, 414 non-coding transcripts, 1127 coding paralogs, and 350 non-coding paralogs.
# -- > I just use the rownames of the final object:
siletti_genes = rownames(siletti_processed@assays$RNA@counts)

## make comparison:
intersection_genes = intersect(ensembl_98_genes,siletti_genes)
ensembl98_but_not_siletti = setdiff(ensembl_98_genes,siletti_genes)
siletti_but_not_ensembl98 = setdiff(siletti_genes,ensembl_98_genes)
all_genes = unique(c(ensembl_98_genes,siletti_genes))

# upset
UpSetR::upset(UpSetR::fromList(list(siletti=siletti_genes,ensembl98 = ensembl_98_genes)),text.scale = 2)

# compare:
# make a list of mathcing gene names
all_gene_match_list = list()
max_dist = 1
for(genex in ensembl98_but_not_siletti){
  str_dists = stringdist::stringdist(siletti_but_not_ensembl98,genex,method='lv')
  min_str_idx = which(str_dists == min(str_dists))
  min_str = str_dists[min_str_idx]
  names(min_str) = siletti_but_not_ensembl98[min_str_idx]
  all_gene_match_list[[genex]] = min_str
  if(min(str_dists)[1] > max_dist){
    all_gene_match_list[[genex]] = NA
  }
}
gene_match_overview = do.call(rbind,all_gene_match_list[!is.na(all_gene_match_list)]) %>% as.data.frame()

relevant_genes = all_gene_match_list[!is.na(all_gene_match_list)]
relevant_genes_df = data.frame(ensembl_98 = names(relevant_genes),siletti = unlist(sapply(relevant_genes,function(x){return(names(x[1]))})))

# save this file and manually curate
data.table::fwrite(relevant_genes_df,paste0(output_path,"ambigous_genes.txt"),sep="\t")

# load manually curated version for further use --> cannot do this automatically to many genes which you have to decide manually
#curated_ambigous_genes_df = data.table::fread(paste0(output_path,"ambigous_genes_curated.txt"),data.table = F)
curated_ambigous_genes_df = data.table::fread(paste0("data/ambigous_genes_curated.txt"),data.table = F)


##########
### Update Tadross gene names in count matrix (new object!)
##########

message("-----",Sys.time(),": Update gene reference differences in datasets ")

# get gene names, compare to curated table and set new rownames on matrix
tadross_count_matrix_updated = tadross_processed@assays$RNA@counts
gene_names_tadross = data.frame(old_gene = as.character(rownames(tadross_count_matrix_updated)))
gene_names_tadross = dplyr::left_join(gene_names_tadross,curated_ambigous_genes_df,by=c("old_gene"="ensembl_98"))
gene_names_tadross$siletti[is.na(gene_names_tadross$siletti)] = gene_names_tadross$old_gene[is.na(gene_names_tadross$siletti)]
rownames(tadross_count_matrix_updated) = gene_names_tadross$siletti

# make new seurat object
tadross_processed_updated = SeuratObject::CreateSeuratObject(counts = tadross_count_matrix_updated,project = "Tadross",meta.data = tadross_processed@meta.data)
tadross_processed_updated@meta.data$Cell_ID = rownames(tadross_processed_updated@meta.data) # make sure cellids are valid
tadross_processed_updated = Seurat::NormalizeData(tadross_processed_updated) # basic normalization

##########
### Make overview of dataset-specific genes
##########

message("-----",Sys.time(),": Find gene expression differences between datasets ")

# find genes that occur only in one or two of the datasets or show very large batch differences and export them --> can be excluded from HVG detection leater on

# helper function
gene_pcts = function(seurat_object){
  # get pct per gene
  message("Preparing gene subsets: ")
  all_genes = rownames(seurat_object@assays$RNA@counts)
  cut_breaks=ceiling(ncol(seurat_object@assays$RNA@counts) / 20000) # cut to avoid crashes
  cut_genes = cut(1:length(all_genes),breaks = cut_breaks)
  cut_genes_levels = as.character(levels(cut_genes))
  all_pct_res = list()
  message("Iterate over ",cut_breaks," gene subsets: ")
  for(cut_gene_level in cut_genes_levels){ # iterate over subsets
    message(cut_gene_level)
    local_genes = all_genes[which(cut_genes == cut_gene_level)]
    local_mat_expr = seurat_object@assays$RNA@counts[local_genes,]
    local_mat_expr[local_mat_expr != 0] <- 1
    local_gene_pcts = data.frame(gene = rownames(local_mat_expr), occ_local = Matrix::rowSums(local_mat_expr) ) # / ncol(local_mat_expr)
    local_gene_pcts$pct_local = local_gene_pcts$occ_local  / ncol(local_mat_expr)
    all_pct_res[[cut_gene_level]] = local_gene_pcts
  }
  all_gene_pcts = do.call(rbind,all_pct_res) %>% as.data.frame()
  return(all_gene_pcts)
}

# apply
gene_pcts_tadross = gene_pcts(tadross_processed_updated)
colnames(gene_pcts_tadross)[2:3] =c("occ_tadross","pct_tadross")
gene_pcts_siletti = gene_pcts(siletti_processed)
colnames(gene_pcts_siletti)[2:3] =c("occ_siletti","pct_siletti")

# make combined df
gene_pcts_combined = dplyr::full_join(gene_pcts_tadross,gene_pcts_siletti,by="gene")
gene_pcts_combined$pct_tadross[is.na(gene_pcts_combined$pct_tadross)] = 0
gene_pcts_combined$pct_siletti[is.na(gene_pcts_combined$pct_siletti)] = 0

# save
data.table::fwrite(gene_pcts_combined,paste0(output_path,"dataset_expression_genes.txt"),sep="\t")

# TODO: make some stat

# TODO: sav json for HVGS

##########
### Merge objects
##########

message("-----",Sys.time(),": Merge datasets ")

merged_seurat_object = merge(x = tadross_processed_updated,y=siletti_processed)

#merged_seurat_object = readRDS(paste0(file_name_prefix,".rds"))

# curate metadata
temp_meta =merged_seurat_object@meta.data %>% dplyr::select(Cell_ID,Sample_ID,Donor_ID, Dataset,Batch_ID,sex ,age_years, percent.mt,nCount_RNA ,nFeature_RNA ,InitialCluster, InitalClass, 
                                                            Technology,unspliced_ratio, Roi, Tissue)
rownames(temp_meta) = temp_meta$Cell_ID
merged_seurat_object@meta.data = temp_meta

##########
### Add basic seurat recipe
##########

message("-----",Sys.time(),": seurat recipe ")

all_human_exclude_genes = unlist(jsonlite::read_json(path = "data/human_features_exclude_list_2.json"))

merged_seurat_object = scUtils::seurat_recipe(merged_seurat_object,
                                              nfeatures_vst = 2000,
                                              genes_to_remove = all_human_exclude_genes,
                                              npcs_PCA = 80,
                                              calcUMAP = TRUE,
                                              findClusters = TRUE,
                                              clusterRes = 6)

# clear object for export
dummy=matrix(data = as.numeric())
merged_seurat_object@assays[["RNA"]]@var.features = character()
merged_seurat_object@assays[["RNA"]]@scale.data <- dummy[,-1] 

##########
### save to rds and h5ad
##########

message("-----",Sys.time(),": Save merged data ")

file_name_prefix = paste0(output_path,"human_hypothalamus_merged")

# save data to rds
saveRDS(merged_seurat_object,paste0(file_name_prefix,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = merged_seurat_object,filename = paste0(file_name_prefix,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(file_name_prefix,".h5seurat"), dest =  paste0(file_name_prefix,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)
system(paste0("rm ",paste0(file_name_prefix,".h5seurat")))

message("-----",Sys.time(),": Complete ")


