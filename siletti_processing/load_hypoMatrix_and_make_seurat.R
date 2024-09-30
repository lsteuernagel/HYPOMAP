# libs
library(hdf5r)
library(Matrix)
library(loomR)
library(Seurat)
library(scUtils)

# requires to have our matrix export hypomatrix.hdf5 and the full loom file for the netadata in the path
siletti_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/linnarson_human_data/"

## now re-open it
file.h5 <- H5File$new(paste0(siletti_path,"hypomatrix.hdf5"), mode="r+")
## lets look at the content
file.h5$ls(recursive=TRUE)
# get sub_matrix
align_full_matrix_h5 <- file.h5[["sub_matrix"]]
# get into memory
align_full_matrix <- align_full_matrix_h5[,]
# make dgc matrix in seurat format
align_full_matrix_dgc = as(align_full_matrix, "dgCMatrix") # align_full_matrix
align_full_matrix_dgc = t(align_full_matrix_dgc)
# rm dense matrix
rm(align_full_matrix)
gc()
# get gene and cell names
# ---> not included
# close h5 file properly
file.h5$close_all()

# query metadata from loom file
# Connect to the loom file in read/write mode
loom_level5_linnarson <- connect(filename = paste0(siletti_path,"adult_human_20221007.loom"), mode = "r+", skip.validate = TRUE)
# Viewing a dataset in the 'row_attrs' group with S3 $ chaining
all_col_attributes = names(loom_level5_linnarson$col.attrs)
# since get.attribute.df is broken (https://github.com/mojaveazure/loomR/issues/34) I iterate over all columns and manually build the metadata
columndata_list = list()
for(i in 1:length(all_col_attributes)){
  columndata_list[[i]]=loom_level5_linnarson[[paste0("col_attrs/",all_col_attributes[i])]][]
}
metadata_linnarson_level5 = as.data.frame(do.call(cbind,columndata_list),stringsAsFactors=F)
colnames(metadata_linnarson_level5) = all_col_attributes
# get gene names:
genes = loom_level5_linnarson[[paste0("row_attrs/","Gene")]][]
#close loom 
loom_level5_linnarson$close_all()
# IMPORTANT: we assume that the order is the same for the hypothalamus cells! and also for genes
metadata_hypothalamus = metadata_linnarson_level5[metadata_linnarson_level5$ROIGroupCoarse == "Hypothalamus",]
rownames(metadata_hypothalamus) = metadata_hypothalamus$CellID
# set row and col names of matrix
rownames(align_full_matrix_dgc) = genes
colnames(align_full_matrix_dgc) = metadata_hypothalamus$CellID

# create seura object
linnarson_hypo_seurat = Seurat::CreateSeuratObject(counts = align_full_matrix_dgc,meta.data = metadata_hypothalamus,project = "Linnarson10x_Hs",min.cells=0,min.features=0)
linnarson_hypo_seurat = Seurat::NormalizeData(linnarson_hypo_seurat)

# simple processing:
# read features to excludes
genes_to_exclude_file = "features_exclude_list2.json"
features_exclude_list= jsonlite::read_json(genes_to_exclude_file)
features_exclude_list = unlist(lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}}))
features_exclude_list = toupper(features_exclude_list)
# seurat recipe:
linnarson_hypo_seurat = scUtils::seurat_recipe(linnarson_hypo_seurat,nfeatures_vst = 3000,sample_column = "SampleID",npcs_PCA = 60,
                                               remove_hvgs = TRUE,genes_to_remove = features_exclude_list,calcUMAP = TRUE,clusterRes = 3)

# save as rds
saveRDS(linnarson_hypo_seurat,paste0(siletti_path,"/linnarson_hypothalamus_raw.rds"))





