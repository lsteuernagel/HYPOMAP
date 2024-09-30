
# from: https://github.com/rbpatt2019/chooseR/blob/master/R/helper_functions.R
#' Find clusters in Seurat object
#'
#' This function combines the necessary and sequential calls
#'  to Seurat::FindNeighbors and Seurat::FindClusters
#'
#' @param obj Seurat object
#' @param reduction String naming the reduction for \code{Seurat::FindNeighbors}
#' @param npcs Integer of number of dimensions for \code{Seurat::FindNeighbors}
#' @param assay String naming the assay for \code{Seurat::FindNeighbors}
#'  This assay is also used to define features in the default settings
#' @param features Character vector containing features for \code{Seurat::FindNeighbors}
#' @param resolution Numeric of resolution to use for \code{Seurat::FindClusters}
#' @param verbose Boolean controlling verbosity
#'  If \code{TRUE}, print updates to terminal
#'
#' @return A Seurat object with the new clusters stored in Idents and the
#'  old clusters stored in \code{seurat_clusters}
#'
#' @section Warning:
#' Each call to this function overwrites \code{seurat_clusters}
#'
#' @seealso \url{https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FindClusters}
#'  and\url{https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FindNeighbors}
#'  for further documentation on the component functions.
#'
#' @examples
#' data <- find_clusters(obj)
#' data <- find_clusters(obj, npcs = 30, verbose = TRUE)
#' \dontrun{
#' data <- find_clusters()
#' }
#' @export
find_clusters <- function(
    obj,
    reduction = "pca",
    npcs = 100,
    assay = "integrated",
    features = NULL,
    resolution = 0.8,
    verbose = FALSE) {
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = reduction,
    dims = 1:npcs,
    assay = integtrated,
    features = features,
    verbose = verbose,
    graph.name = paste(reduction, assay, sep = ".")
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = resolution,
    graph.name = paste(reduction, assay, sep = "."),
    verbose = verbose
  )
  return(obj)
}

#' Generate n sub-samples
#'
#' Wrapper function around rplicate and split to generate
#'  \code{n} subsamples of size \code{size} from \code{input}
#'
#' @param n Integer number of replicates
#' @param input Vector to be subset
#' @param size Numeric fraction of input to take for subsample
#' @param replace Boolean whether to sample with or without replacement
#' @param simplify Whether to simplify output to array or matrix, if possible
#'
#' @return A list of subsamples of \code{input}.
#'  If \code{simplify}, then in a simplfied array or matrix
#'
#' @seealso \code{\link[base]{replicate}} and \code{\link[base]{sample}}
#'  for documentation on the component functions.
#
#' @examples
#' n_samples(3, c(1, 3, 5, 2))
#' n_samples(5, c('a', 'b', 'c'), size = 0.5, replace = TRUE)
#'
#' \dontrun{
#' n_samples(c(1, 2, 3))
#' }
#' @export
n_samples <- function(
    n,
    input,
    size = 0.8,
    replace = FALSE,
    simplify = FALSE) {
  splits <- replicate(
    n,
    sample(
      input,
      as.integer(length(input) * size),
      replace = replace
    ),
    simplify = simplify
  )
}

# from: https://github.com/rbpatt2019/chooseR/blob/master/R/pipeline.R
# Functions used directly in the analysis and example script

`%>%` <- magrittr::`%>%`

#' Run multiple clusters on the same Seurat Object
#'
#' Given a Seurat object, cluster a random subset of the object n times
#'  This assumes that the object has already been apropriately
#'  normalised, scaled, and reduced.
#'
#' @param obj Seurat object containing apropriately preprocessed cells
#' @param n Integer number of clusterings to perform
#' @param size Numeric fraction of obj to use for random subset
#' @param npcs Integer number of dimensions from the reduction to use
#' @param res Numeric resolution to use for clustering
#' @param reduction String reduction slot of obj to use for clustering
#' @param assay String assay slot of obj to use for clustering
#'
#' @return clusters dplyr::tibble with cells in first column and clustering
#'          results in consecutive columns
#'
#' @examples
#' results <- multiple_clusters(obj)
#' results <- multiple_clusters(obj, n = 10, res = 4, reduction = "my_pca")
#' \dontrun{
#' multiple_cluster()
#' }
#' @export
multiple_cluster <- function(
    obj,
    n = 100,
    size = 0.8,
    npcs = 100,
    res = 1.2,
    reduction = "pca",
    assay = "SCT") {
  
  # Initialise tibble for data
  clusters <- dplyr::as_tibble(Seurat::Cells(obj))
  clusters <- dplyr::rename(clusters, "cell" = value)
  
  # Get samples
  samples <- n_samples(n, Seurat::Cells(obj), size = size)
  
  # Repeated clusters
  j <- 1
  for (idx in samples) {
    message(paste0("\tClustering ", j, "..."))
    small_obj <- obj[, idx]
    small_obj <- find_clusters(
      small_obj,
      reduction = reduction,
      npcs = npcs,
      resolution = res,
      assay = assay
    )
    clusters <- dplyr::left_join(
      clusters,
      dplyr::as_tibble(Seurat::Idents(small_obj), rownames = "cell"),
      by = "cell"
    )
    j <- j + 1
  }
  return(clusters)
}

#' Find matches for a given clustering resolution
#'
#' Given a dataframe where clustering results are stored in columns,
#'  find which cells were included in the same cluster.
#'  This function scores matches as 1 (TRUE), mis-matches as 0 (FALSE)
#'  and drops as 1i (imaginary). This allows further calcualtions to
#'  distinguish between mismatches and drops.
#'
#' @param col String column name containing results
#' @param df Dataframe where each row is a cell and each column is the
#'  results of a clustering
#'
#' @return dgCMatrix of dimension \code{n_cells x n_cells} encoding matches
#'
#' @examples
#' find_matches("first_cluster", data)
#' \dontrun{
#' find_matches(first_cluster, data)
#' }
#' @export
find_matches <- function(col, df) {
  mtchs <- outer(df[[col]], df[[col]], "==")
  # Records drops as imaginary, mtchs as 1, not mtchs as 0
  mtchs[is.na(mtchs)] <- 1i
  return(mtchs)
}

# sparse version
find_matches_sparse <- function(col, df) {
  mtchs <- outer(df[[col]], df[[col]], "==")
  # Records drops as imaginary, mtchs as 1, not mtchs as 0
  # mtchs[is.na(mtchs)] <- 1i
  mtchs = as(mtchs, "dgCMatrix")
  return(mtchs)
}


#' Score the number of matches
#'
#' Given a complex number x, calculate the fraction of time that it was a match
#'  Assuming the reall part represents matches, the imaginary part represents
#'  drops, and that there were \code{n} repeats
#'
#' @param x Complex number encoding the number of matches and drops
#' @param n Integer number of repeats
#'
#' @return Float fraction of times that a match occurred
#'
#' @examples
#' percent_match(5 + 3i, 10)
#' \dontrun{
#' percent_match("abc")
#' }
#' @export
percent_match <- function(x, n = 100) {
  return(Re(x) / (n - Im(x)))
}

#' Compute group average frequencies
#'
#' Given a tibble of pairwise colcustering frequencies, compute the average
#'  pairwise co-clustering frequency for each cluster.
#'
#' @param tbl Tibble with pairwise co-clustering frequency per cell
#' @param clusters Vector of cluster labels. This should be for each
#'  cell. That is, if there are 100 cells, this should have
#'  \code{length(clusters) == 100}
#'
#' @return Tibble with clusters in both column and rows.
#'  This format makes it ideal for plotting heatmaps
#'
#' @examples
#' group_scores(data, Seurat::Idents(data))
#' \dontrun{
#' group_scores(data, unique(Seurat::Idents(data)))
#' }
#' @export
group_scores <- function(tbl, clusters) {
  colnames(tbl) <- clusters
  data <- tbl %>%
    tibble::add_column("cell_1" = clusters) %>%
    tidyr::pivot_longer(-cell_1, names_to = "cell_2", values_to = "percent") %>%
    dplyr::group_by(cell_1, cell_2) %>%
    dplyr::summarise("avg_percent" = mean(percent)) %>%
    dplyr::ungroup()
  return(data)
}

#' Compute group average silhouette scores
#'
#' Given a path to an RDS object containing silhouette scores, calculate the
#'  cluster-wise average silhouette score
#'
#' @param sil Silhouette score output
#' @param res Numeric resolution at which the silhouette score was calculated
#'
#' @return Tibble with 3 columns: cluster, avg_sil, and res
#'
#' @examples
#' scores <- group_sil("/path/to/my/silhouette.rds", 0.8)
#' \dontrun{
#' group_scores(sil_scores, 1.0)
#' }
#' @export
group_sil <- function(sil, res) {
  sil <- tibble::as_tibble(sil[, ]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise("avg_sil" = mean(sil_width)) %>%
    tibble::add_column("res" = res)
  return(sil)
}

#' Compute confidence intervals on the median
#'
#' Given a numeric vector, calculate confidence intervals on the median
#'  The user may specify the interval, number of replicates used for
#'  bootstrapping, and the type of bootstrap to use
#'
#' @param x Numeric vector containing distribution of values to analysed
#' @param interval Numeric confidence interval to calculate. 0 < interval < 1
#' @param R Integer number of bootstrap replicates to use
#' @param type String type of confidence interval to calculate
#'
#' @return Tibble with clusters in both column and rows.
#'  This format makes it ideal for plotting heatmaps
#'
#' @seealso \code{\link[boot]{boot}} and \code{\link[boot]{boot.ci}}
#'  for more information on their parameters and meanings
#'
#' @examples
#' ci <- boot_median(c(1, 2, 3, 4, 5), interval = 0.9)
#' \dontrun{
#' group_scores(1, 2, 3, 4, 5)
#' }
#' @export
boot_median <- function(x, interval = 0.95, R = 25000, type = "bca") {
  # Define median to take data and indices for use with boot::
  med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }
  
  # Calculate intervals
  boot_data <- boot::boot(data = x, statistic = med, R = R)
  boot_ci <- boot::boot.ci(boot_data, conf = interval, type = type)
  
  # Extract desired statistics
  ci <- list(
    low_med = boot_ci$bca[4],
    med = boot_ci$t0,
    high_med = boot_ci$bca[5]
  )
  return(ci)
}