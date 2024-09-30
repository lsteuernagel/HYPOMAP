# test_run <- function(metadata = metacolumns, 
#                      subset = TRUE, 
#                      datasets = c("924S5", "C064.A.01", "C064.A.06")) {
#   if (pacman::p_detectOS() == "Darwin" && subset == TRUE) {
#     x <- dplyr::filter(metadata, sample %in% datasets)
#     return(x)
#   } else {
#     return(metadata)
#   }
# }

workflow_options <- function(width = 18,
                             height = 18,
                             lowUMI = 500,
                             lowGENES = 450, 
                             batch_nPCS = 1:20, 
                             batch_RES = c(0.4, 1.2),
                             doublet_form_rate = 0.05,
                             seed = 1960, 
                             project = "HuHy", 
                             integration_nPCS = 1:50, 
                             nfeatures = 3000,
                             SIZE = 0.6, 
                             dpi = 320, 
                             out_path = NULL,
                             #^data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq_test/"
                             data_path= "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/human_nucseq/",
                             merged_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_human_data/humanHypoMap2/"
) {
  out_path <- paste0("results/", out_path)
  lapply(c(paste0(out_path, "plots"), paste0(out_path, "results")), function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  set.seed(seed)
  opts <- list(
    width = width,
    height = height,
    lowUMI = lowUMI,
    lowGENES = lowGENES, 
    batch_nPCS = batch_nPCS,
    seed = seed, 
    project = project,
    batch_RES = batch_RES, 
    doublet_form_rate = doublet_form_rate,
    integ_nPCS = integration_nPCS, 
    nfeatures = nfeatures,
    SIZE = SIZE, 
    dpi = dpi,
    out_path = out_path,
    base_path = out_path,
    data_path = data_path,
    merged_path = merged_path
  )
  return(opts)
}

KeyGenes_nonNeuron  <-
  unique(c("RBFOX3",
           "AGT",
           "OPALIN",
           "PECAM1",
           "LUM", # fibro
           "TAGLN", # mural
           "CD3G", # immune
           "F13A1",
           "TMEM119",
           "PDGFRA",
           "PLN", # ependymal
           "COL23A1"
           
  ))

KeyGenes <-
  unique(c("AGRP",
           "POMC",
           #    "GAL",
           "NPY",
           "CARTPT",
           "PMCH",
           "GNRH1",
           "HCRT",
           "AVP",
           "VIP",
           "TRH",
           "CCK",
           "MC3R",
           "MC4R",
           "GLP1R",
           "GLP2R",
           "GIPR",
           "TAC3",
           "PDYN",
           "CRH",
           "SIM1",
           "OTP",
           "OXT",
           "KISS1",
           "GHRH",
           "NR5A1",
           "GHSR",
           "SST",
           #    "LEPR", "INSR",
           #    "GPX3" #human hypothalamus
           "PNOC"
           #    "OPRL1"
  ))

# myMTfilt <- function(object, probab = 0.98, backup.pct = 15, ...){ 
#   if (sum(object[["percent.mt"]]) <= 0) {
#     print("percent.mt is 0 for all cells")
#     object[["miQC.keep"]] <- "keep"
#     return(object)
#   }
#   y <- outliers::scores(object[["percent.mt"]], type = "t", prob = probab)
#   y1 <- object[["percent.mt"]]
#   backup.pct = min(min(y1$percent.mt[y$percent.mt]), backup.pct)
#   try(print(paste("Threshold", unique(object[["sample"]]), "=", signif(backup.pct, digits = 2 ))), silent = TRUE)
#   object[["miQC.keep"]] <- dplyr::if_else(object[["percent.mt"]] > backup.pct, "discard", "keep")
#   return(object)
# }

#Calculate QC thresholds
thresholds <- function(split.by = c("batch", "sample"),datum = QCtable) {
  split.by <- match.arg(split.by)
  QCcal <- list(
    L2SD = ~ mean(.x, trim = 0.1) - 2 * sd(.x),
    #      L3SD = ~ mean(.x, trim = 0.1) - 3 * sd(.x),
    L2MAD = ~ median(.x) - 2 * mad(.x),
    #      L3MAD = ~ median(.x) - 3 * mad(.x),
    q5 =  ~ quantile(.x, probs = c(.05), na.rm = T),
    q1 =  ~ quantile(.x, probs = c(.01), na.rm = T),
    U3MAD = ~ median(.x) + 3 * mad(.x),
    U4MAD = ~ median(.x) + 4 * mad(.x),
    U5MAD = ~ median(.x) + 5 * mad(.x),
    U2SD = ~ mean(.x, trim = 0.1) + 2 * sd(.x),
    q995 = ~ quantile(.x, probs = c(.995), na.rm = T),
    #      MAD = ~ mad(.x),
    MED = ~ median(.x),
    MIN = ~ min(.x),
    MAX = ~ max(.x),
    MEAN = ~ mean(.x, trim = 0.1)
  )
  
  if (split.by == "sample") {
    thresholdsSample <- datum %>%
      #        filter(outlier.any == FALSE) %>%
      group_by(sample) %>%
      select(sample, 
             nCount_RNA,
             nCount.log,
             nFeature_RNA,
             nFeat.log) %>%
      #               percent.mt,
      #               MTlog1p) %>%
      summarise(across(where(is.numeric), QCcal)) %>% # , .names = "{.fn}.{.col}"
      mutate(across(contains(c(
        "nCount.log", "nFeat.log"
      )), ~ 2 ^ .x)) %>%
      mutate(across(contains("log1p"), ~ expm1(.x))) %>%
      select(
        c("sample",
          "nCount_RNA_MED", 
          "nCount_RNA_MEAN", 
          "nCount_RNA_MAX",
          "nCount_RNA_q995",
          "nCount_RNA_U4MAD", 
          "nCount_RNA_U5MAD", 
          "nCount.log_MEAN",
          "nCount.log_U2SD",
          "nCount.log_MED",
          "nCount_RNA_MIN",
          "nCount_RNA_q5",
          #            "nCount_RNA_L2MAD", 
          "nCount.log_L2SD",
          "nFeature_RNA_MED", 
          "nFeature_RNA_MEAN", 
          "nFeature_RNA_MAX",
          "nFeature_RNA_q995",
          "nFeature_RNA_U4MAD", 
          "nFeature_RNA_U5MAD", 
          "nFeat.log_MEAN",
          "nFeat.log_U2SD",
          "nFeat.log_MED",
          "nFeature_RNA_MIN",
          "nFeature_RNA_q5",
          #            "nFeature_RNA_L2MAD", 
          "nFeat.log_L2SD"
        )
      ) %>%
      #        select(order(colnames(.))) %>%
      mutate(across(where(is.numeric), round, 1))
    # thresholdsSampletemp <- data.frame(t(thresholdsSample[,-1]))
    # names(thresholdsSampletemp) <- thresholdsSample$sample
    return(thresholdsSample)
  }
  if (split.by == "batch") {
    thresholdsGroup <- datum %>%
      filter(outlier.any == FALSE) %>%
      group_by(batch) %>%
      select(batch, 
             nCount_RNA,
             nCount.log,
             nFeature_RNA,
             nFeat.log,
             percent.mt,
             MTlog1p
      ) %>%
      summarise(across(where(is.numeric), QCcal)) %>% # , .names = "{.fn}.{.col}"
      mutate(across(contains(c(
        "nCount.log", "nFeat.log"
      )), ~ 2 ^ .x)) %>%
      mutate(across(contains("log1p"), ~ expm1(.x))) %>%
      select(
        c("batch",
          "percent.mt_MEAN",
          "MTlog1p_MEAN",
          "percent.mt_MED",
          "MTlog1p_MED",
          "percent.mt_MAX",
          "percent.mt_q995", 
          "percent.mt_U3MAD",
          "percent.mt_U4MAD",
          "MTlog1p_U2SD"
        )
      ) %>%
      #        select(order(colnames(.))) %>%
      mutate(across(where(is.numeric), round, 1))
    # thresholdsGrouptemp <- data.frame(t(thresholdsGroup[,-1]))
    # names(thresholdsGrouptemp) <- thresholdsGroup$batch
    return(thresholdsGroup)
  }
}

##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){
  
  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames
  if(is.character(current_embedding[,1])){
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      message("Found ",length(rnames)," rows in embedding and ",length(rownames(metadata))," rows in metadata.")
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
}


##########
### clear_clustering
##########

#' Eliminate small clusters from a vector of labels substituting with the labels of NN
#' @param x vector of labels
#' @param min_cells minimum number of cells to keep cluster
#' @param nn_idx matrix of cells x k NN idx --> e.g. output of annoy or rann
#' @return updated vector of labels

clear_clustering = function(x,min_cells,nn_idx){
  x = as.character(x)
  new_x=x
  # which clusters are too small ?
  small_clusters = names(table(x)[table(x) < min_cells])
  # go through small clusters and move cells to laregr clusters based on neighbors
  if(length(small_clusters)>0){
    #message("Removing ",length(small_clusters)," clusters")
    for(i in 1:length(small_clusters)){
      current_cluster = small_clusters[i]
      which_idx = which(x == current_cluster)
      # get idx for k NN
      neighbor_idx = nn_idx[which_idx,]
      # substitute with cluster labels
      if(length(which_idx)>1){
        neighbor_clusters = apply(neighbor_idx,1,function(z,cluster_vector){return(cluster_vector[z])},cluster_vector=x)
        # extract that most common label in neighbors
        clusters_vote = apply(neighbor_clusters,2,function(z,exclude_clusters){
          return(names(sort(table(z),decreasing = TRUE))[! names(sort(table(z),decreasing = TRUE)) %in% exclude_clusters][1])
        },exclude_clusters=small_clusters)
      }else{
        neighbor_clusters = x[neighbor_idx]
        clusters_vote = names(sort(table(neighbor_clusters),decreasing = TRUE))[! names(sort(table(neighbor_clusters),decreasing = TRUE)) %in% small_clusters][1]
      }
      # overwrite cluster label
      new_x[which_idx]=clusters_vote
    }
  }
  
  return(new_x)
}

# set_parallel_cores <- function(num_cores = NULL) {
#   options(future.rng.onMisuse = "ignore")
#   if (!is.null(num_cores) && pacman::p_detectOS() == "Linux") {
#     options(future.globals.maxSize = 180000 * 1024 ^ 2)
#     cores_to_use <- num_cores
#     return(cores_to_use)
#   } else if (pacman::p_detectOS() == "Darwin") {
#     options(future.globals.maxSize = 12288 * 1024 ^ 2)
#     x <- parallel::detectCores(logical = FALSE) - 1
#     cores_to_use <- max(x, 1)
#     return(cores_to_use)
#   } else if (pacman::p_detectOS() == "Linux") {
#     options(future.globals.maxSize = 162000 * 1024 ^ 2)
#     x <- parallel::detectCores(logical = FALSE) - 2 #switch to parallelly::availableCores(logical = FALSE, methods = "Slurm")
#     cores_to_use <- min(max(x, 1), 30)
#     return(cores_to_use)
#   } else {
#     stop("OS not clear, set cores manually")
#   }
# }