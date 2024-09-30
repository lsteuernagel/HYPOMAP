#Read the data and save the list file------------
message("-----",Sys.time(),"01.Read.Files.R")

# from terminal (when in project dir):
# sbatch run_slurm.sh ~/Documents/r_scvi_v3.simg 01_read_data.R

##########
### Load libs
##########

library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(scDblFinder)
library(BiocParallel)
library(scater)
library(scran)
library(qs)
#library(outliers)
library(scuttle)
source("utility_functions.R")

library(scUtils)

##########
### set options
##########

#opts <- list(project = "HuHy", out_path = "test/",seed = 1234567)
opts <- workflow_options(project = "HuHy", out_path = "test/")
opts$cores <- 56

# raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/"
raw_data_path = opts$data_path
message("-----",Sys.time(),": Using path: ",raw_data_path)

print(opts)

# read features to excludes
genes_to_exclude_file = "features_exclude_list2.json"
features_exclude_list= jsonlite::read_json(genes_to_exclude_file)
features_exclude_list = unlist(lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}}))
features_exclude_list = toupper(features_exclude_list)

##########
### load metadata
##########

message("-----",Sys.time(),": load metadata")

metacolumns <-
  read_csv(paste0(raw_data_path,"data/metadata_edit.csv"),
           comment = "#",
           col_types = cols_only(
             paths = col_character(),
             sample = col_character(),
             donor = col_factor(),
             batch = col_factor()
           )
  )

donor_metacolumns <- 
  read_csv(paste0(raw_data_path,"data/donor_metadata_edit.csv"),
           comment = "#",
           col_types = cols(
             age = col_integer(),
             donor = col_factor(),
             sex = col_factor()
           )
  )

metacolumns <- left_join(metacolumns, donor_metacolumns, by = "donor")


print(metacolumns)
#metacolumns = metacolumns[1:2,]

# skip this for now -- but should be fine
# stopifnot(any(
#   metacolumns$paths %in% list.files(
#     path = raw_data_path,
#     pattern = "*.h5",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# ))

##########
### load counts
##########

message("-----",Sys.time(),": load counts")

data.10x <- vector("list", nrow(metacolumns))
QCcounts <- vector("list", nrow(metacolumns))
scrna.list <- vector("list", nrow(metacolumns))

for (i in seq_along(metacolumns$paths)) {
  data.10x[[i]] <- Seurat::Read10X_h5(filename = paste0(raw_data_path,metacolumns$paths[[i]]), use.names = TRUE)
  cells_per_gene <- Matrix::rowSums(data.10x[[i]] > 0)
  counts_per_gene <- Matrix::rowSums(data.10x[[i]])
  QCcounts[[i]] <-
    list("Counts" = counts_per_gene[unique(KeyGenes)],
         "Cells" = cells_per_gene[unique(KeyGenes)])
  QCcounts[[i]][["ratio"]] <-
    signif(c(QCcounts[[i]][[1]] / QCcounts[[i]][[2]]), digits = 2)
}

names(QCcounts) <- metacolumns$sample
# QCcountDF <- cbind.data.frame(QCcounts)
# QCcountDF[is.na(QCcountDF)] <- 0

#To pull values try QCcountDF["POMC", grepl(".ratio", names(QCcountDF))]

##########
### Create seurat object list
##########

message("-----",Sys.time(),": Create seurat object list")

for (i in seq_along(data.10x)) {
  print(i)
  scrna.list[[i]] = CreateSeuratObject(
    counts = data.10x[[i]],
    min.cells = 3,
    min.features = 250,
    project = metacolumns$sample[[i]]
  )
  scrna.list[[i]][["sample"]] = metacolumns$sample[[i]]
  scrna.list[[i]][["donor"]] = metacolumns$donor[[i]]
  scrna.list[[i]][["batch"]] = metacolumns$batch[[i]]
  scrna.list[[i]][["sex"]] = metacolumns$sex[[i]]
  scrna.list[[i]][["age"]] = metacolumns$age[[i]]
}

# rm(data.10x)

names(scrna.list) <- metacolumns$sample

scrna.list <- lapply(scrna.list, function(x) {
  x[['percent.mt']] <- PercentageFeatureSet(x, pattern = "^MT-|^mt-", assay = "RNA"); x
  #  x[['percent.XY']] <- PercentageFeatureSet(x, features = intersect(expressed_genes(x), c(.Y_genes_raw$hgnc_symbol, "XIST")), assay = RNA); x
})

# filter mt- genes
scrna.list <- lapply(scrna.list, myMTfilt)

##########
### merge seurat object list
##########

message("-----",Sys.time(),": merge seurat object list")

scrna <- merge(
  x = scrna.list[[1]],
  y = scrna.list[2:length(scrna.list)],
  add.cell.ids = metacolumns$sample,
  project = opts$project
)

print("merged file")
scrna
rm(scrna.list)
gc()

#write.csv(t(QCcountDF), paste0(opts$out_path,"results/01.SampleQCcountsDF.csv"))

##########
### scran qclust normalization
##########

message("-----",Sys.time(),": scran qclust normalization")

# spligt again
scrna.list <- SplitObject(scrna, split.by = "batch")
# rm(scrna, QCcountDF, QCcounts, cells_per_gene, counts_per_gene); gc()

# write_dims(scrnalist = scrna.list,out_filename = "02.batchdims.preQC")
# BatchQCcountDF <- counts_keygenes(seurat_list = scrna.list, gene_list = KeyGenes)
# write.csv(BatchQCcountDF, paste0(opts$out_path,"results/02.BatchQCcountDF.csv"))

# convert to sce
#print("converting to SCE")
scrna.list <- lapply(scrna.list, as.SingleCellExperiment)

# run scran::quickCluster and computeSumFactors
# print("min cluster sizes without floor value")
if (length(scrna.list) == 1 && length(unique(scrna.list[[1]]$sample)) == 1) {
  scrna.list <- lapply(scrna.list, function(x) {
    print(round(0.05*table(x[['sample']])))
    x[['qclust_batch']] <- scran::quickCluster(x,min.size = max(30, min(round(0.05*table(x[['sample']])))))
    x <- scran::computeSumFactors(x,cluster = x[['qclust_batch']],  min.mean = 0.1)
  })
} else {
  scrna.list <- lapply(scrna.list, function(x) {
    print(round(0.05*table(x[['sample']])))
    x[['qclust_batch']] <- scran::quickCluster(x,
                                               block = x[['sample']],
                                               min.size = max(30, min(round(0.05*table(x[['sample']])))),
                                               block.BPPARAM = MulticoreParam(opts$cores, RNGseed = opts$seed),
                                               BPPARAM = MulticoreParam(opts$cores, RNGseed = opts$seed))
    x <- scran::computeSumFactors(x,
                                  cluster = x[['qclust_batch']],
                                  min.mean = 0.1,
                                  BPPARAM = MulticoreParam(opts$cores, RNGseed = opts$seed))
  })
  
}

##########
### log normalization
##########

scrna.list <- lapply(scrna.list, scuttle::logNormCounts)

##########
### Doublet detection
##########

message("-----",Sys.time(),": Doublet detection")

#TODO must compare clusters vs no-clustering as there's some bad clustering 
scrna.list <- lapply(scrna.list, function(x) {
  print(unique(x[["sample"]]))
  x <- scDblFinder::scDblFinder(x,
                                clusters = x[['qclust_batch']],
                                samples = "sample",
                                multiSampleMode = "split"#,
                                #BPPARAM = MulticoreParam(opts$cores, RNGseed = opts$seed) 
                                #this worked with the default test 924S5, "C064.A.01" "C064.A.06"
                                #Doesn't work on the HPC, don't know why....
  )
})

##########
### revert back to seurat
##########

message("-----",Sys.time(),": revert back to seurat")

scrna.list <- lapply(scrna.list, function(x) {
  x <- as.Seurat(x)
  #  x <- RenameAssays(object = x, originalexp = 'RNA') #Fixed in bioc 3.15
})

##########
### doublet detection with DoubletFinder
##########
# 
# message("-----",Sys.time(),":  Doublet detection 2 ")
# doublet_info_list =list()
# for(i in 1:length(scrna.list)){
#   message(i, " of ",length(scrna.list))
#   current_seurat = scrna.list[[i]]
#   # remove features that should not be HVGs
#   keep_genes = rownames(current_seurat)[!rownames(current_seurat) %in% features_exclude_list]
#   current_seurat_new = Seurat::CreateSeuratObject(current_seurat@assays$RNA@counts[keep_genes,],meta.data = current_seurat@meta.data)
#   current_seurat_new@meta.data$Cell_ID = rownames(current_seurat_new@meta.data)
#   # run preprocessing for the current batch
#   # we remove unwanted genes beforehand so that doublet fidner won't use them --> also exclude below
#   current_seurat_new = scUtils::seurat_recipe(current_seurat_new,
#                                               nfeatures_vst = opts$nfeatures,
#                                               normalize_data = TRUE,
#                                               remove_hvgs = FALSE,
#                                               calcUMAP = FALSE,
#                                               findClusters = TRUE,
#                                               npcs_PCA = 70,
#                                               clusterRes = 1,
#                                               k.param = 30,
#                                               seed = opts$seed)
# 
#   # run doublet detection
#   current_seurat_new = scUtils::apply_DoubletFinder( # scUtils::
#     seurat_object = current_seurat_new,
#     npcs_PCA = 70,
#     pN_fixed = 0.25,
#     pK_max = 0.1,
#     doublet_formation_rate = opts$doublet_form_rate,
#     adjust_nExp = FALSE,
#     doublet_cluster_tresh = NULL,
#     return_seurat = TRUE
#   ) %>% suppressMessages()
# 
#   doublet_info_list[[i]] = current_seurat_new@meta.data[,c("Cell_ID","Doublet")]
# }

##########
### and merge again
##########

message("-----",Sys.time(),": merge again")

if (length(scrna.list) > 1) {
  scrna <- merge(x = scrna.list[[1]],
                 y = scrna.list[2:length(scrna.list)],
                 project = opts$project)
} else if (length(scrna.list) == 1) {
  scrna <- scrna.list[[1]]
}

rm(scrna.list)


### add doublet information back into processed seurat via join
# scrna@meta.data$Cell_ID = rownames(scrna@meta.data)
# message("-----",Sys.time(),":  Add Doublet annotation into merged")
# all_doublets = do.call(rbind,doublet_info_list) %>% as.data.frame()# %>% dplyr::rename(DoubletFinder = Doublet)
# print(colnames(all_doublets))
# print(head(all_doublets))
# colnames(all_doublets) = c("Cell_ID","DoubletFinder")
# meta_temp = dplyr::left_join(scrna@meta.data, all_doublets ,by="Cell_ID")
# rownames(meta_temp) = meta_temp[,"Cell_ID"]
# scrna@meta.data = meta_temp

#remove second list

# qsave(scrna, file = paste0(opts$out_path,"results/02.scrna.doublet"), nthreads = opts$cores)
# print("doublet file saved")

##########
### QC Tables and plots
##########

message("-----",Sys.time(),":  QC Tables and plots")

# scrna <- qread(paste0(opts$out_path, "results/02.scrna.doublet"), nthreads = opts$cores)

print(scrna)

colsOI <- c("nCount_RNA", "nFeature_RNA", "sample", "donor",
            "batch", "percent.mt", "miQC.keep",
            "qclust_batch", "scDblFinder.class")#, "DoubletFinder"   #"sizeFactor", "miQC.probability"

QCtable <- as_tibble(scrna@meta.data[colsOI],
                     rownames = "Cell.Barcode") %>%
  dplyr::mutate(
    percent.mt = signif(percent.mt, 4),
    nFeat.log = signif(log2(nFeature_RNA), 4),
    nCount.log = signif(log2(nCount_RNA), 4),
    MTlog1p = signif(log1p(percent.mt), 4),
    FCR = signif(nCount.log/nFeat.log, 4)
  )

#TODO feature QCs ----
# sort(Matrix::rowSums(scrna))
# table(cut(QCtable$FCR, 20, right = FALSE))

message("DCs <-")
print(QCtable)

DCs <- QCtable %>%
  dplyr::group_by(sample, qclust_batch, scDblFinder.class) %>%
  dplyr::count() %>%
  dplyr::group_by(sample, qclust_batch) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::filter(scDblFinder.class == "doublet") %>%
  dplyr::mutate(clust.dbl = if_else(percent >= 75, TRUE, FALSE),percent = signif(percent, 3)) %>%  #check this, 50% seemed to cut off a lot
  dplyr::select(-scDblFinder.class)

QCtable <- dplyr::left_join(QCtable, DCs, by = c("sample" = "sample", "qclust_batch" = "qclust_batch"))

QCtable$clust.dbl <- QCtable$clust.dbl %>%
  replace_na(FALSE)

QCtable <- QCtable %>%
  dplyr::mutate(doublet = if_else(scDblFinder.class == "doublet", TRUE, FALSE),
                doublet = doublet | clust.dbl,
                mtHigh = if_else(miQC.keep == "discard", TRUE, FALSE),
                #         FCRhigh = if_else(FCR >= 1.2, TRUE, FALSE),
                outlier.any = doublet | mtHigh)

##########
###plots stage 1----
##########
# 
# pdf(paste0(opts$out_path,"plots/03a.QC.PreFiltering.pdf"), width = opts$width, height = opts$height)
# 
# db1 <- QCtable %>%
#   mutate(sample = fct_reorder(sample, as.numeric(batch))) %>% 
#   ggplot(aes(x = sample, y = nFeature_RNA, fill = doublet)) +
#   geom_boxplot(notch = FALSE) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_fill_manual(values = c("light blue", "red")) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank())
# 
# db2 <- QCtable %>%
#   mutate(sample = fct_reorder(sample, as.numeric(batch))) %>% 
#   ggplot(aes(x = sample, y = nCount_RNA, fill = doublet)) +
#   geom_boxplot(notch = FALSE) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_fill_manual(values = c("light blue", "red")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
# 
# db1 / db2 + plot_layout(guides = "collect"); rm(db1, db2)
# 
# box_count_featureplot(group_x = "sample", filtered = FALSE)
# 
# box_count_featureplot(group_x = "batch", filtered = FALSE)
# 
# # GGally::ggpairs(
# #   QCtable,
# #   columns = c("nFeat.log",
# #               "nCount.log",
# #               "MTlog1p",
# #               "batch"),
# #   lower = list(continuous = GGally::wrap(
# #     "smooth", alpha = 0.3, size = 0.1
# #   )),
# #   aes(colour = !QCtable$outlier.any)
# # )
# 
# QCtable %>%
#   #  filter(outlier.any == FALSE) %>%
#   ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.mt)) +
#   geom_point(alpha = 1 / 5, size = 0.3) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   facet_wrap(~ batch) +
#   scale_color_viridis_c(trans = "log1p", breaks = c(0:5, 7,10, 15, 25, 50))
# 
# QCtable %>%
#   #  filter(outlier.any == FALSE) %>%
#   ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.mt)) +
#   geom_point(alpha = 1 / 5, size = 0.3) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   facet_wrap(~ batch + sample, labeller = labeller(.multi_line = FALSE)) +
#   scale_color_viridis_c(trans = "log1p", breaks = c(0:5, 7,10, 15, 25, 50))
# 
# QCtable %>%
#   ggplot(aes(x = nFeature_RNA, y = percent.mt, colour = miQC.keep)) +
#   geom_point(size = 0.3) +
#   facet_wrap(~ batch + sample, labeller = labeller(.multi_line = FALSE), scales = "free") +
#   scale_fill_manual(values = c("light blue", "red"))

##########
### Sample thresholds------
##########

message("-----",Sys.time(),":  Sample thresholds")

thresholdsSample <- thresholds(split.by = "sample", datum = QCtable)

sample.cuts <- thresholdsSample %>%
  dplyr::transmute(
    sample = sample,
    nCountHigh = nCount_RNA_q995, 
    nCountLow = if_else(nCount_RNA_q5 >= opts$lowUMI, nCount_RNA_q5, opts$lowUMI),
    nFeatureHigh = nFeature_RNA_q995,
    nFeatureLow = if_else(nFeature_RNA_q5 >= opts$lowGENES, nFeature_RNA_q5, opts$lowGENES) #?change
  )

##########
### SMore stage 1 plots
##########
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(nFeature_RNA)) +
#   geom_density(trim = TRUE) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   geom_vline(aes(xintercept = nFeatureHigh), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nFeatureLow), colour = "blue", sample.cuts) +
#   facet_wrap(~ sample)
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(nCount_RNA)) +
#   geom_density(trim = TRUE) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   geom_vline(aes(xintercept = nCountHigh), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nCountLow), colour = "blue", sample.cuts) +
#   facet_wrap(~ sample)
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.mt)) +
#   geom_point(alpha = 1 / 5, size = 0.3) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   facet_wrap(~ sample) +
#   scale_color_viridis_c(trans = "log1p", breaks = c(0:5, 7,10, 15, 25, 50)) + 
#   geom_hline(aes(yintercept = nCountHigh), colour = "blue", sample.cuts) +
#   geom_hline(aes(yintercept = nCountLow), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nFeatureHigh), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nFeatureLow), colour = "blue", sample.cuts) 
# 
# dev.off()

##########
### filters
##########

message("-----",Sys.time(),":  filter data")

#thresholdsSample %>% select(contains(c( "sample", "nCount")) & !contains(c("_MAD", "_MED", "_MIN"))) %>% View()

out <- vector("list", length = length(unique(QCtable$sample)))
names(out) <- unique(QCtable$sample)

#can't use the log calcs here as you've transformed them back to percent
for (samp in unique(QCtable$sample)) {
  i <- sample.cuts %>%
    filter(sample == samp)
  
  out[[samp]] <- QCtable %>% filter(sample == samp) %>%
    mutate(
      nCountHigh = if_else(nCount_RNA >= i$nCountHigh, TRUE, FALSE),
      nCountLow = if_else(nCount_RNA <= i$nCountLow, TRUE, FALSE),
      nFeatureHigh = if_else(nFeature_RNA >= i$nFeatureHigh, TRUE, FALSE),
      nFeatureLow = if_else(nFeature_RNA <= i$nFeatureLow, TRUE, FALSE),
    )
}

sampleQC <- dplyr::bind_rows(out)

#Fail if not identical-------
stopifnot(identical(sampleQC[1:ncol(QCtable)], QCtable[1:ncol(QCtable)]))

QCtable <- sampleQC %>% dplyr::mutate(
  outlierSample = mtHigh|nCountHigh|nCountLow|nFeatureHigh|nFeatureLow|doublet,
  outlier.any = outlierSample,
  nFeat.log = log2(nFeature_RNA),
  nCount.log = log2(nCount_RNA),
  MTlog1p = log1p(percent.mt)
)

sampleQC.dims <- table(QCtable$sample, QCtable$outlierSample)
#round(prop.table(sampleQC.dims, 1), 2)

##########
### group filters next
##########

message("-----",Sys.time(),":  filter data")

# TODO: threshodl function ?
thresholdsGroup <- thresholds(split.by = "batch", datum = QCtable)

#thresholdsGroup %>% select(contains(c("batch", "MT")) & !contains(c("_MIN"))) %>% View()

batch.cuts <- thresholdsGroup %>%
  dplyr::transmute(
    batch = batch,
    mtHigh = percent.mt_q995
  )

out <- vector("list", length = length(unique(QCtable$batch)))
names(out) <- unique(QCtable$batch)

#can't use the log calcs here as you've transformed them back to percent
for (bat in unique(QCtable$batch)) {
  i <- batch.cuts %>% dplyr::filter(batch == bat)
  
  out[[bat]] <- QCtable %>% filter(batch == bat) %>%
    dplyr::mutate(bat.mtHigh = if_else(percent.mt >= i$mtHigh, TRUE, FALSE))
}

batchQC <- dplyr::bind_rows(out)

QCtable <- batchQC %>% dplyr::mutate(outlier.any = outlierSample|bat.mtHigh)

batchQC.dims <- table(QCtable$batch, QCtable$outlier.any)

##########
### finish the QC steps PLOTS
##########

# pdf(paste0(opts$out_path,"plots/03b.QC.PostFiltering.pdf"), width = opts$width, height = opts$height)
# 
# box_count_featureplot(group_x = "sample", filtered = TRUE)
# 
# box_count_featureplot(group_x = "batch", filtered = TRUE)

# GGally::ggpairs(
#   QCtable,
#   columns = c("nFeat.log",
#               "nCount.log",
#               "MTlog1p",
#               "batch"),
#   lower = list(continuous = GGally::wrap(
#     "smooth", alpha = 0.3, size = 0.1
#   )),
#   aes(colour = !QCtable$outlier.any)
# )
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.mt)) +
#   geom_point(alpha = 1 / 5, size = 0.3) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   facet_wrap(~ batch) +
#   scale_color_viridis_c(trans = "log1p", breaks = c(0:5, 7,10, 15, 25, 50)) + 
#   theme_bw()
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.mt)) +
#   geom_point(alpha = 1 / 5, size = 0.3) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   facet_wrap(~ batch + sample, labeller = labeller(.multi_line = FALSE)) +
#   scale_color_viridis_c(trans = "log1p", breaks = c(0:5, 7,10, 15, 25, 50)) + 
#   theme_bw()
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(nFeature_RNA)) +
#   geom_density(trim = TRUE) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   geom_vline(aes(xintercept = nFeatureHigh), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nFeatureLow), colour = "blue", sample.cuts) +
#   facet_wrap(~ sample)
# 
# QCtable %>%
#   filter(outlier.any == FALSE) %>%
#   ggplot(aes(nCount_RNA)) +
#   geom_density(trim = TRUE) +
#   scale_x_continuous(trans = "log2", breaks = scales::breaks_log(n = 8, base = 2)) +
#   geom_vline(aes(xintercept = nCountHigh), colour = "blue", sample.cuts) +
#   geom_vline(aes(xintercept = nCountLow), colour = "blue", sample.cuts) +
#   facet_wrap(~ sample)
# 
# dev.off()

##########
### finish the QC steps
##########

#Filter the outliers - this could be delicate if the order changes, careful
#Try to switch to a vectorised approach e.g. ID %in% cellstoKeep

message("-----",Sys.time(),":  apply QC")

# old version: 
# stopifnot(identical(Cells(scrna), QCtable$Cell.Barcode))
#scrna@meta.data$outlier <- QCtable$outlier.any

# rewrite this part with a join to ensure it's ok:
temp_meta = scrna@meta.data
temp_meta$Cell_ID = rownames(temp_meta)
temp_meta = dplyr::left_join(temp_meta,QCtable %>% dplyr::select(Cell.Barcode,outlier = outlier.any),by=c("Cell_ID"="Cell.Barcode"))
# set rownames again:
rownames(temp_meta) = temp_meta$Cell_ID
# check that join worked:
stopifnot(nrow(temp_meta) == nrow(scrna@meta.data))
# add back to meatdata
scrna@meta.data = temp_meta

# subset
# scrna.outliers <- subset(scrna, outlier == TRUE)
# scrna <- subset(scrna, outlier == FALSE)


##########
### add simple non-integrated seurat
##########

scrna = scUtils::seurat_recipe(scrna,
                               nfeatures_vst = opts$nfeatures,
                               normalize_data = TRUE,
                               remove_hvgs = TRUE,
                               genes_to_remove = features_exclude_list,
                               sample_column = "sample",
                               calcUMAP = TRUE,
                               findClusters = TRUE,
                               npcs_PCA = 100,
                               clusterRes = 4,
                               k.param = 30,
                               seed = opts$seed)

##########
### save data - just a seurat rds
##########

message("-----",Sys.time(),":  save data")

# save version without subset:
saveRDS(scrna,paste0(raw_data_path,"human_hypo_raw.rds"))

# saveRDS(scrna,paste0(raw_data_path,"human_hypo_raw_filtered.rds"))
# saveRDS(scrna.outliers,paste0(raw_data_path,"human_hypo_raw_outliers.rds"))

#write.csv(QCtable, paste0(raw_data_path,"QCtable.csv"))
#write.csv(outlier_reason, paste0(raw_data_path,"outlierReason.csv"))
write.csv(QCtable, paste0(raw_data_path,"QCtable.csv"))
write.csv(batch.cuts, paste0(raw_data_path,"batch.cuts.csv"))
write.csv(sample.cuts, paste0(raw_data_path,"sample.cuts.csv"))
#write.csv(thresholdsSampleDF, paste0(raw_data_path,"thresholdsSampleDF.csv"))
write.csv(sampleQC.dims, paste0(raw_data_path,"sampleQC.dims.csv"))
#write.csv(thresholdsGroupDF, paste0(raw_data_path,"thresholdsGroupDF.csv"))
write.csv(batchQC.dims, paste0(raw_data_path,"batchQC.dims.csv"))

message("-----",Sys.time(),":  finished")



##########
### save data John's version
##########
# 
# qsave(scrna, paste0(opts$out_path,"results/03.scrna.list.filtered"), 
#       preset = "fast", 
#       nthreads = opts$cores)
# 
# qsave(scrna.outliers, paste0(opts$out_path,"results/03.scrna.outliers"), 
#       preset = "fast", 
#       nthreads = opts$cores)
# 
# print(scrna)
# 
# #save twice, one copy for records
# save(opts, file = paste0(opts$out_path,"opts.RData"))
# save(opts, file = "opts.RData")
# 
# thresholdsSampleDF <- data.frame(t(thresholdsSample))
# names(thresholdsSampleDF) <- thresholdsSample$sample
# 
# thresholdsGroupDF <- data.frame(t(thresholdsGroup))
# names(thresholdsGroupDF) <- thresholdsGroup$batch
# 
# outlier_reason <- QCtable %>%
#   group_by(sample) %>%
#   summarise_if(is.logical, mean) %>% #sum for count
#   mutate_if(is.numeric, ~ . * 100) 
# outlier_reason[,2:ncol(outlier_reason)] <- round(outlier_reason[,2:ncol(outlier_reason)], 2)
# 
# write.csv(outlier_reason, paste0(opts$out_path,"results/03.outlierReason.csv"))
# write.csv(QCtable, paste0(opts$out_path,"results/03.QCtable.csv"))
# write.csv(batch.cuts, paste0(opts$out_path,"results/03.batch.cuts.csv"))
# write.csv(sample.cuts, paste0(opts$out_path,"results/03.sample.cuts.csv"))
# write.csv(thresholdsSampleDF, paste0(opts$out_path,"results/03.thresholdsSampleDF.csv"))
# write.csv(sampleQC.dims, paste0(opts$out_path,"results/03.sampleQC.dims.csv"))
# write.csv(thresholdsGroupDF, paste0(opts$out_path,"results/03.thresholdsGroupDF.csv"))
# write.csv(batchQC.dims, paste0(opts$out_path,"results/03.batchQC.dims.csv"))
# 
#print("01.QC.R finished") ; Sys.time()
#sessionInfo()