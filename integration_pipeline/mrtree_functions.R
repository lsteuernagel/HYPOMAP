##### These are functions copied from: https://github.com/pengminshi/MRtree
#### Installation of the full repo required too many dependencies, so I am only using the code of the core functions
#### Includes some modifications

#' MRtree main function
#'
#' The multi-resolution clusters is input into \code{mrtree}, and the output is the constructed hierarchical cluster tree.
#' modfified!!!???
#'
#' @param x object containing multi-resolution clustering results
#' **data source**
#' building a reconciled hierarchical clustering tree requires information about which cluster each sample has been assinged
#'  to at different resolutions. This information can be supplied in various forms, as a matrix, data.frame or more specialised
#'   object. In all cases the object provided must contain numeric columns with the naming structure `PXS` where `P` is
#'   a prefix indicating that the column contains clustering information, `X` is a numeric value indicating the clustering
#'   resolution and `S` is any additional suffix to be removed. For `SingleCellExperiment` objects this information must
#'   be in the `colData` slot and for `Seurat` objects it must be in the `meta.data` slot. For all objects except matrices
#'   any additional columns can be used as aesthetics, for matrices an additional metadata data.frame can be supplied if
#'   required.
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param max.k the maximum resolution (number of clusters) to consider in building the tree
#' @param consensus boolean, whether to perform consensus clustering within the clusterings with the same number of clusters.
#' @param verbose boolean, whether to show logmsgs
#' @param n.cores integer, number of cores for parallel computation
#'
#' @return A list containing  \describe{
#'   \item{labelmat.mrtree}{The Reconciled tree saved as a label matrix, with duplicated layers omited.}
#'   \item{labelmat.recon}{The full reconciled tree save in a label matrix}
#'   \item{labelmat.flat}{The initial flat clustering (cluster tree) as input for MRtree algorithm}
#'   \item{resolutions}{The corresponding clustering resolution of the initial cluster tree}
#'   \item{paths}{The unique path in the resulted reconciled hierarchical cluster tree}
#' }
#'


mrtree <- function(x, ...) {
  UseMethod("mrtree", x)
}

# require(tidytree)
# require(bnstruct)

#' \code{mrtree} with label saved in a matrix as input
#' @rdname mrtree
#' @import checkmate parallel
#' @importFrom bnstruct knn.impute

mrtree.matrix <- function(labelmat, prefix = NULL, suffix = NULL,  max.k = Inf,
                          consensus = F, sample.weighted = F, augment.path = F,
                          verbose = F, n.cores = parallel::detectCores()-1) {
  if (n.cores > parallel::detectCores()-1){
    warnings('Use ', parallel::detectCores()-1, 'cores instead!')
    n.cores = parallel::detectCores()-1
  }

  Ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)]))) # number of clusters per resolution

  if (is.unsorted(Ks)){
    ord = order(Ks, decreasing = F)
    labelmat = labelmat[, ord]
    Ks = Ks[ord]
  }

  if (max.k != Inf){
    is.effective = Ks < max.k
    labelmat = labelmat[,is.effective, drop=F]
    Ks = Ks[is.effective]
  }

  labelmat.flat = labelmat        # save the initial flat clustering results
  nk = length(Ks)                 # number of initial clusterings

  if (consensus){
    # within resolution consensus clustering
    ccs.out = consensus_clustering_within_resolution(labelmat,
                                                     ks = Ks,
                                                     sample.weighted=sample.weighted,
                                                     n.cores = n.cores)
    labelmat = ccs.out$labelmat.ccs     # after consensus cluster, each layer has unique number of clusters
    labelmat.ccs = labelmat
  } else {
    labelmat.ccs = NULL
  }

  if (sample.weighted){
    labelmat.encoded = do.call(cbind, lapply(1:ncol(labelmat), function(l){
      label_to_membership(labelmat[,l])
    }))

    colsum = colSums(labelmat.encoded)
    sample_weights = 1 / sqrt(apply(labelmat.encoded, 1, function(x) sum(colsum[x]))) # W^{-1/2}

  } else {
    sample_weights = rep(1, nrow(labelmat))         # weights for samples
  }


  if (any(is.na(labelmat))){
    # if there are missing labels, input via nearest neighbor
    labelmat = bnstruct::knn.impute(labelmat, k=5)
    labelmat.imputed = labelmat
  } else {
    labelmat.imputed = NULL
  }

  ###############
  # initialize
  if (verbose)
    logmsg("initilize the tree ...")

  tree = construct_tree_from_labelmat(labelmat)  # initialize tree
  labelmat.in.paths = labelmat
  paths = unique(labelmat.in.paths)  # start from all existing paths
  bad.node = get_bad_nodeset(tree)
  if (verbose) {
    logmsg("initial size of bad nodes:", length(bad.node))
  }

  candidate.ind = which(tree$end %in% bad.node)  # include all edges to the nodes not visited
  candidate = tree[candidate.ind, ]

  # progress bar
  message('Run MrTree:')
  total_bad = length(bad.node)
  pb <- txtProgressBar(min = -total_bad-1, max =0, style = 3)
  lowest.layer.last =  NULL


  while (length(bad.node) > 0) {
    # progress bar
    setTxtProgressBar(pb, -length(bad.node))

    # top down-collection, lowest three layers
    layer = as.numeric(sapply(candidate$start, function(x) strsplit(x, split = ";")[[1]][1]))
    layer.unique.sorted = sort(unique(layer), decreasing = F)
    lowest.layer = layer.unique.sorted[1: min(2, length(layer.unique.sorted))] # among three lowest layers
    ind.lowest.layer = which(layer %in% lowest.layer)
    candidate = candidate[ind.lowest.layer, ]

    if (augment.path){
      paths = augment_path(paths, layers = setdiff(lowest.layer, lowest.layer.last))
      lowest.layer.last = unique(c(lowest.layer.last, lowest.layer))
    }

    # calculate the cost for candidates
    candidate$cost = unlist(parallel::mclapply(1:nrow(candidate), function(i) {
      cost(node.start = candidate$start[i], node.end = candidate$end[i], paths = paths,
           labelmat = labelmat, labelmat.in.paths = labelmat.in.paths, sample_weights = sample_weights)
    }, mc.cores = n.cores))  #, mc.cores =  parallel::detectCores()-1)

    # choose the edge with minimum cost
    ind.min = which.min(candidate$cost)
    if (length(ind.min) > 1) {
      # choose the one with max count
      ind.min = ind.min[which.max(candidate$count[ind.min])]
    }
    node.start = candidate$start[ind.min]
    node.end = candidate$en[ind.min]

    if (verbose) {
      logmsg("Select edge-> start:", node.start, ", end:", node.end, ", cost:",
             candidate$cost[ind.min])
    }

    # remove paths that has the same end node but different start node
    if (verbose)
      logmsg("prune path ...")
    paths = prune_paths(paths = paths, node.start = node.start, node.end = node.end)

    # only update the nodes that are affected
    if (verbose)
      logmsg("assign sample to the path ...")
    node.end.decoded = decode(node.end)
    node.start.decoded = decode(node.start)
    node.ind.affected = which(labelmat.in.paths[, node.end.decoded$layer] ==
                                node.end.decoded$label & labelmat.in.paths[, node.start.decoded$layer] !=
                                node.start.decoded$label)  # only the node on the eliminated paths are affected

    if (verbose) {
      logmsg("length(node.ind.affected)=", length(node.ind.affected))
    }

    if (length(node.ind.affected) > 0) {
      checkmate::assert_true(nrow(labelmat.in.paths) == nrow(labelmat))

      labelmat.in.paths[node.ind.affected, ] = paths[assign_samples_to_paths(labelmat = labelmat[node.ind.affected, ,drop=F],
                                                                             paths = paths), , drop=F]
    }

    # update tree
    if (verbose)
      logmsg("update the tree ...")

    tree = construct_tree_from_labelmat(labelmat.in.paths)
    bad.node = get_bad_nodeset(tree)

    if (verbose) {
      logmsg("number of bad.node = ", length(bad.node))
    }

    candidate.ind = which((tree$end %in% bad.node))  # include all edges to the nodes not visited
    candidate = tree[candidate.ind, ]

    if (augment.path){
      paths = rbind(unique(labelmat.in.paths),  # remove the path with 0 data points assigned
                    unique(paths[apply(paths[,lowest.layer, drop=F], 1, function(x) any(x == -1)),, drop=F]) # the path that has -1 in lowest layers
      )
      paths = unique(paths)
    } else {
      paths = unique(labelmat.in.paths)
    }

    if (verbose) {
      logmsg("number of remaining paths is ", nrow(paths))
    }
  }
  setTxtProgressBar(pb, -length(bad.node))
  close(pb)

  if (consensus){
    # map the reconciled layers with the original layers
    labelmat.recon = labelmat.in.paths[, ccs.out$k.idx]
    colnames(labelmat.recon) = colnames(labelmat.flat)
    Ks.recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
  } else {
    labelmat.recon =  labelmat.in.paths
    Ks.recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
  }

  unique.idx = which(!duplicated(Ks.recon[-length(Ks.recon)], MARGIN = 2)) # remove the last column in labelmat.recon
  labelmat.tree = labelmat.recon[, unique.idx, drop=F]
  colnames(labelmat.tree) = paste0("K", Ks.recon[unique.idx])

  resolutions = colnames(labelmat.flat)

  return(list(labelmat.mrtree = labelmat.tree,
              labelmat.recon = labelmat.recon,
              labelmat.flat = labelmat.flat,
              labelmat.ccs = labelmat.ccs,
              labelmat.imputed = labelmat.imputed,
              resolutions = resolutions,
              paths = paths,
              params = list(prefix=prefix, suffix=suffix,
                            consensus=consensus, sample.weighted=sample.weighted, max.k=max.k)))
}


#' \code{mrtree} with labelmatrix save as a data frame as input
#' @rdname mrtree
#' @import checkmate
#' @export
mrtree.data.frame <- function(x, prefix, suffix = NULL, ...) {
  checkmate::assert_character(prefix, any.missing = FALSE, len = 1)

  clust_cols <- grepl(prefix, colnames(x), fixed = TRUE)
  if (!is.null(suffix)) {
    clust_match_suffix <- grepl(suffix, colnames(x), fixed = TRUE)
    clust_cols = clust_cols & clust_match_suffix
  }

  if (sum(clust_cols) < 2) {
    stop(paste("Less than two column names matched the prefix: ", prefix, "and suffix: ",
               suffix), call. = FALSE)
  }

  clusterings <- as.matrix(x[, clust_cols])
  mode(clusterings) <- "numeric"

  mrtree(clusterings, ...)
}


#' \code{mrtree} with SingleCellExperiment object as input
#' @rdname mrtree
#' @import checkmate
#' @export
mrtree.SingleCellExperiment <- function(x, prefix='sc3_', suffix = "_clusters",...) {

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("The SingleCellExperiment package is missing, this must be",
         "installed for clustree to use SingleCellExperiment objects",
         call. = FALSE)
  }

  checkmate::assert_class(x, "SingleCellExperiment")
  # checkmate::assert_character(exprs, any.missing = FALSE, len = 1)

  mrtree(data.frame(x@colData), prefix, suffix, ...)
}


#' \code{mrtree} with Seurat object as input
#' @rdname mrtree
#' @import checkmate
#' @export
mrtree.Seurat <- function(x, prefix = 'RNA_snn_res.', suffix = NULL,...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is missing, this must be installed for ",
         "clustree to use seurat objects",
         call. = FALSE)
  }

  checkmate::assert_class(x, "Seurat")
  # checkmate::assert_character(exprs, any.missing = FALSE)

  mrtree(x@meta.data, prefix, suffix, ...)
}


#' Construct the cluster tree from multi-resolution clustering save as a label matrix
#'
#' @param labelmat sample by m label matrix, each columns represents a clustering for a certain resolution
#' @return tree saved as an edge list containing \describe{
#'     \item{start}{label of start node}
#'     \item{end}{label of end node}
#' }
construct_tree_from_labelmat <- function(labelmat) {

  nb.layer = ncol(labelmat)

  if (nb.layer <= 1) {
    stop("Error! Can not construct tree with less than two layers")
  }

  tree = NULL
  for (l in 1:(nb.layer - 1)) {

    edgelist = table_to_edgelist(labels.start = labelmat[, l], labels.end = labelmat[,
                                                                                     l + 1])

    edgelist$start = paste(l, edgelist$start, sep = ";")
    edgelist$end = paste(l + 1, edgelist$end, sep = ";")

    tree = rbind(tree, edgelist)
  }

  return(tree)
}

#' Augment the path to include on additional node as alternative to each layer
augment_path <- function(paths, layers, aug.name=-1){

  if (length(layers)==0){
    return (paths)
  }

  nb.layers = ncol(paths)
  nb.paths = nrow(paths)

  paths.augmented = NULL
  for (layer in layers){
    if (layer != nb.layers & !(aug.name %in% paths[,layer])){
      if (layer == 1){
        aug_ind2 = apply(paths[,2:nb.layers, drop=F], 1, function(x) any(x==-1))
        paths2 = paths[!aug_ind2, 2:nb.layers, drop=F]

        paths.augmented = rbind(path_dot_product(paths1=matrix(aug.name, nrow=1, ncol=1),
                                                 paths2=paths2),
                                paths.augmented)
      } else{
        aug_ind1 = apply(paths[,1:(layer-1), drop=F], 1,  function(x) any(x==-1))
        aug_ind2 = apply(paths[,(layer+1):nb.layers,drop=F], 1, function(x) any(x==-1))

        paths1 = cbind(paths[!aug_ind1,1:(layer-1), drop=F],
                       rep(aug.name, sum(!aug_ind1)))
        paths2 = paths[!aug_ind2, (layer+1):nb.layers,drop=F]

        paths.augmented = rbind(path_dot_product(paths1=paths1, paths2=paths2),
                                paths.augmented)
      }
    }
  }
  paths = unique(rbind(paths, paths.augmented))

  return(paths)
}


#' Calculate the cost generate by adding the edge = (start, end) and removing all the other conflicting edge
#'
#' @param node.start start node of the edge (out-vertex)
#' @param node.end end node of the edge (in-vertex)
#' @param paths viable path
#' @param labelmat.in.paths assigned labels of data points to the optimimum paths
#' @param sample_weights weights per sample (clustering) used in loss function (all ones by default)
#'
#' @return a scalar representing the cost
cost <- function(node.start, node.end, paths, labelmat, labelmat.in.paths, sample_weights) {

  # prune path
  paths.new = prune_paths(paths = paths, node.start = node.start, node.end = node.end)

  # only update the nodes that are affected
  node.end.decoded = decode(node.end)
  node.ind.affected = which(labelmat.in.paths[, node.end.decoded$layer] == node.end.decoded$label)

  labelmat.in.paths.affected.old = labelmat.in.paths[node.ind.affected,,drop=F]
  labelmat.in.paths.affected.new = paths.new[assign_samples_to_paths(labelmat = labelmat[node.ind.affected,,drop=F],
                                                                     paths = paths.new), ]

  cost = sum(sample_weights[node.ind.affected] * (labelmat.in.paths.affected.old != labelmat.in.paths.affected.new))

  return(cost)
}


#' Convert the start end table to an edgelist
#'
#' @param labels.start label of the start nodes
#' @param labels.end labels of the end nodes
#'
#' @return and edge list \describe{
#'     \item{start}{label of start node}
#'     \item{end}{label of end node}
#'     \item{count}{number of data points in the edge}
#' }
table_to_edgelist <- function(labels.start, labels.end) {
  tab = table(labels.start, labels.end)
  rowsum = rowSums(tab)
  colsum = colSums(tab)
  norm = matrix(rep(rowsum, ncol(tab)), ncol = ncol(tab)) + matrix(rep(colsum, nrow(tab)), nrow = nrow(tab), byrow = T)
  norm[norm == 0] = Inf
  tab = tab/norm
  nonzero.ind = which(tab > 0)
  edgelist = data.frame(start = rownames(tab)[row(tab)[nonzero.ind]],
                        end = colnames(tab)[col(tab)[nonzero.ind]],
                        count = tab[nonzero.ind], cost = Inf)

  return(edgelist)
}


#' parse the node name as the layer and cluster
#'
#' @param node.name name of the node
#'
#' @return A list \describe{
#'     \item{layer}{the layer in the tree the node belongs to}
#'     \item{label}{the label of the cluster in the layer the node belongs to}
#' }
decode <- function(node.name) {

  layer_label = as.numeric(unlist(strsplit(node.name, ";")))

  layer = layer_label[1]
  label = layer_label[2]

  return(list(layer = layer, label = label))
}


#' Remove the paths that include the same node.end but different node.start
#'
#' @param paths Initial paths to be pruned
#' @param node.start the label of the starting node
#' @param node.end the label of the ending node
#'
#' @return A matrix whose rows contains the remaining paths after prunning
prune_paths <- function(paths, node.start, node.end) {

  node.start.decoded = decode(node.start)
  node.end.decoded = decode(node.end)

  is.conflict = (paths[, node.start.decoded$layer] != node.start.decoded$label) &
    (paths[, node.end.decoded$layer] == node.end.decoded$label)

  if (sum(is.conflict) > 0 & node.end.decoded$layer < ncol(paths)){
    is.selected = (paths[, node.start.decoded$layer] == node.start.decoded$label) &
      (paths[, node.end.decoded$layer] == node.end.decoded$label)
    paths.to.remove.modified = path_dot_product(paths1=paths[is.selected, 1:node.end.decoded$layer, drop=F],
                                                paths2=paths[is.conflict, (1+node.end.decoded$layer):ncol(paths), drop=F])
  } else {
    paths.to.remove.modified = NULL
  }

  paths = unique(rbind(paths[!is.conflict, ], paths.to.remove.modified))

  return(paths)
}

#' outer product of path prefix and suffix
path_dot_product <- function(paths1, paths2){
  paths1 = unique(paths1)
  paths2 = unique(paths2)

  nb.paths1 = nrow(paths1)
  nb.paths2 = nrow(paths2)

  if (nb.paths1 ==0 | nb.paths2==0){
    return(NULL)
  }
  path.product = cbind(paths1[rep(1:nb.paths1, nb.paths2), , drop=F],
                       paths2[rep(1:nb.paths2, each=nb.paths1), , drop=F])
  return(path.product)
}


#' Assign data points to the viable path that are closest to the original clustering assignment
#'
#' @param labelmat Initial multi-resolution cluster assignment, n-by-m, corresponding to m resolutions
#' @param paths The viable paths
#'
#' @return A vector of length n containing the path index to which the data points are assigned to
#'
#' @import checkmate parallel
#' @export
assign_samples_to_paths <- function(labelmat, paths, n.cores=parallel::detectCores() - 1) {

  n.layers = ncol(paths)
  paths = t(paths)

  path.labels = apply(labelmat, 1, function(label){
    min.ind = which.min(colSums((abs(paths - label) > 0)))
    if (length(min.ind) > 1) {
      min.ind = sample(min.ind, 1)
    }
    return(min.ind)
  })

  return(path.labels)
}


#' Output the set of bad node in the given tree
#'
#' @param tree a list containg start & and pairs representing a cluster tree
#'
#' @return A vector containing the name of the bad vertex
get_bad_nodeset <- function(tree) {
  nb_in_edge.out = aggregate(tree$start, by = list(tree$end),
                             FUN = function(x) length(unique(x)))
  bad.node = nb_in_edge.out$Group.1[nb_in_edge.out$x > 1]

  return(bad.node)
}


#' Perform consensus clustering among clusterings of same number of clusters
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param ks vector of length m, number of clustering for each clustering
#' @param n.cores maximum number of cores used for parallel computation
#'
#' @return a list containing \describe{
#' \item{labelmat.ccs}{consensus clustering results saved in a label matrix}
#' \item{k.idx}{idx mapping the original columns to the new columns}
#' }
#' @import parallel
#' @export
consensus_clustering_within_resolution <- function(labelmat, ks=NULL, sample.weighted=F, n.cores = NULL){
  if (is.null(ks)){
    ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
  }

  unique.out = unique_ind(ks)
  k_unique = unique.out$C
  k_idx = unique.out$ic
  nk = length(k_unique)

  if (is.null(n.cores)){
    n.cores = min(nk, parallel::detectCores() - 4)
  }

  labelmat.ccs = do.call(cbind, lapply(1:nk, function(i){
    k = k_unique[i]
    idx = which(k_idx==i)
    if (length(idx) > 1){
      if (sample.weighted){
        consensus_clustering_weighted(labelmat[,idx], k=k)
      } else {
        consensus_clustering(labelmat[,idx], k=k)
      }
    } else {
      labelmat[,idx, drop=F]
    }
  }))
  return(list(labelmat.ccs=labelmat.ccs, k.idx = k_idx))
}


#' Perform consensus clustering
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param k number of clusters for the final consensus clustering
#'
#' @return consensus clustering result saved as a label vector (NA is labels missing for all clusterings)
#'
#' @import RSpectra
#' @export
consensus_clustering <- function(labelmat, k=NULL){
  m = ncol(labelmat)
  n = nrow(labelmat)
  if (m==1){
    return (matrix(labelmat, ncol=1))
  }

  if (is.null(k)){ # if k not provided
    ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
    k = unique(ks)
    if (length(k)>1){
      warning('Error, clustering applied to non-k clustering')
      k = max(ks)
    }
  }

  if (k==1){
    return (matrix(1, nrow=nrow(labelmat), ncol=1))
  }

  # one-hot encode the label matrix
  labelmat.encoded = do.call(cbind, lapply(1:m, function(l){
    label_to_membership(labelmat[,l])
  }))

  # if there are missing values, weight rows by 1/\sqrt{#assigned clustering}
  if(any(is.na(labelmat))){
    w = sqrt(rowSums(labelmat.encoded))
    allmissing = w==0

    labels = rep(NA, n)
    svd.out = RSpectra::svds(labelmat.encoded[!allmissing,] / w[!allmissing], k=k, nv=0)
    u.weigted = svd.out$u * rep(svd.out$d, each=nrow(svd.out$u))
    kmeans.out = kmeans(u.weigted, centers=k, nstart = 5, iter.max = 30)
    labels[!allmissing] = kmeans.out$cluster
  } else {
    svd.out = RSpectra::svds(labelmat.encoded, k=k, nv=0)
    u.weigted = svd.out$u * rep(svd.out$d, each=nrow(svd.out$u))
    kmeans.out = kmeans(u.weigted, centers=k, nstart = 5, iter.max = 30)
    labels = kmeans.out$cluster
  }

  return(labels)
}

#' weighted consensus clustering
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param k number of clusters for the final consensus clustering
#'
#' @return consensus clustering result saved as a label vector (NA is labels missing for all clusterings)
#' @import RSpectra
#' @export
consensus_clustering_weighted <- function(labelmat, k=NULL){
  m = ncol(labelmat)
  n = nrow(labelmat)
  if (m==1){
    return (matrix(labelmat, ncol=1))
  }

  if (is.null(k)){ # if k not provided
    ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
    k = unique(ks)
    if (length(k)>1){
      warning('Error, clustering applied to non-k clustering')
      k = max(ks)
    }
  }

  if (k==1){
    return (matrix(1, nrow=nrow(labelmat), ncol=1))
  }

  # one-hot encode the label matrix
  labelmat.encoded = do.call(cbind, lapply(1:m, function(l){
    label_to_membership(labelmat[,l])
  }))

  colsum = colSums(labelmat.encoded)
  weights_sqrt = 1 / sqrt(apply(labelmat.encoded, 1, function(x) sum(colsum[x]))) # W^{-1/2}

  # if there are missing values, weight rows by 1/\sqrt{#assigned clustering}
  if(any(is.na(labelmat))){
    allmissing = rowSums(labelmat.encoded) == 0
    labels = rep(NA, n)
    svd.out = RSpectra::svds(weights_sqrt[!allmissing] * labelmat.encoded[!allmissing,], k=k, nv=0)
    u.weigted = svd.out$u * rep(svd.out$d, each=nrow(svd.out$u))
    kmeans.out = kmeans(u.weigted, centers=k, nstart = 5, iter.max = 30)
    labels[!allmissing] = kmeans.out$cluster
  } else {
    svd.out = RSpectra::svds(weights_sqrt * labelmat.encoded, k=k, nv=0)
    u.weigted = svd.out$u * rep(svd.out$d, each=nrow(svd.out$u))
    kmeans.out = kmeans(u.weigted, centers=k, nstart = 5, iter.max = 30)
    labels = kmeans.out$cluster
  }

  return(labels)
}


#' Impute the missing labels by K nearest neighbor
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param k number of nearest neighbor to consider
#'
#' @return Imputed labelmatrix
#' @importFrom FNN get.knnx
#' @export
knn_impute <- function(labelmat){
  labelmat.encoded = do.call(cbind, lapply(1:ncol(labelmat),
                                           function(j) label_to_membership(labelmat[,j]))) # convert chategorical matirx to numeric
  missing.row = apply(labelmat, 1, function(x) any(is.na(x)))

  labelmat.encoded.nonmissing = labelmat.encoded[!missing.row,]
  labelmat.nonmissing = labelmat[!missing.row,]
  for (missing.row.idx in which(missing.row)){
    knn.idx = FNN::get.knnx(labelmat.encoded.nonmissing,
                            labelmat.encoded[missing.row.idx, ,drop=FALSE],
                            k=1)$nn.index
    missing.col = is.na(labelmat[missing.row.idx, ]) # missing columns for this missing sample
    labelmat[missing.row.idx, missing.col] = labelmat.nonmissing[knn.idx, missing.col, drop=F]
  }

  return (labelmat)
}


#' Convert label vector to membership matrix
#'
#' @param labels a vector of labels from K distint classes
#' @param labels.names (optional) alternative label names, used for naming columns of membership matrix
#'
#' @return an n-by-K binary membership matrix
#' @import checkmate
#' @export
label_to_membership <- function(labels, label.names = NULL) {
  if (is.null(label.names)) {
    label.names = sort(unique(labels[!is.na(labels)]))
  } else {
    checkmate::assert_true(all(labels[!is.na(labels)] %in% label.names))
  }

  K = length(label.names)
  memb = t(sapply(labels, function(lab) {
    if (is.na(lab)){
      rep(0, K) # if label missing then return 0 row
    } else {
      as.numeric(label.names == lab)
    }
  } ))
  return(matrix(memb, ncol=K))
}


#### From plot

#' Plot MRtree results as a dendrogram. If reference labels are provided, a pie chart is
#' shown at each tree node, detailing the label proprotions.
#'
#' @param labelmat a n by m label matrix, in ouput of \code{mrtree} function, \code{labelmat.mrtree}
#' @param ref.labels a factor or characteristic vector specifying the reference labels of n data points
#' @param show.ref.labels boolean, whether to show the labels of major type at tree nodes and tips
#' @param label.order a vector specifying the order of labels for default colors
#' @param node.size scalar, the size of the pie chart / node
#' @param cols a vector of colors, one per each label
#' @param plot.piechart boolean, whether to draw the pie chart in the tree plot
#' @param tip.labels a vector of strings specifying the labels of tree leafs
#' @param tip.label.dist distance of the tip labels to the tree tips
#' @param show.branch.labels boolean, whether to show the branch labels for convenience of flipping branches
#' @param fip.branch a list of vectors each of size 2, indicating the branch labels to flip. Each time two branches are flipped.
#' @param legend.title string as legend title, default is an empty string
#' @param bottom.margin size of the bottom margin, need to be adjusted to show the full labels
#'
#' @importFrom data.tree as.Node as.phylo.Node
#' @importFrom ape Ntip Nnode
#' @importFrom tibble as_tibble
#' @importFrom tidytree full_join as.treedata
#' @importFrom ggtree ggtree nodepie layout_dendrogram
#' @importFrom scales hue_pal
#' @export
plot_tree <- function(labelmat, ref.labels = NULL, show.ref.labels = T,
                      label.order = NULL, node.size = 0.2, cols = NULL,
                      plot.piechart = T, tip.labels = NULL, tip.label.dist = 2,
                      show.branch.labels = F, flip.branch=NULL,
                      legend.title = "",  bottom.margin = 25){
  if(is.null(ref.labels)){
    ref.labels = rep('', nrow(labelmat))
    show.ref.labels=F
    draw.pie.chart=F
  } else {
    ref.labels = as.character(ref.labels)
    ref.labels = gsub('-','_', ref.labels)
  }

  if(is.null(label.order)){
    label.order = sort(unique(ref.labels))
  } else {
    label.order = gsub('-','_', label.order)
    if (! all(label.order %in% ref.labels)){
      warnings(sum(!label.order %in% ref.labels),'label name not if the reference labels!')
    }
  }

  if (plot.piechart){
    pointsize = 0.01
  } else{
    pointsize = 5
  }

  n = nrow(labelmat)
  p = ncol(labelmat)

  # save in data.tree format
  labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
  df = as.data.frame(unique(labelmat), stringsAsFactors = F)
  df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
  tree.datatree = data.tree::as.Node(df)

  # phylo tree for visualization
  tree.phylo = data.tree::as.phylo.Node(tree.datatree)

  # reference type per node
  ord = data.frame(node=1:(ape::Ntip(tree.phylo)+ape::Nnode(tree.phylo)),
                   row.names =c(tree.phylo$tip.label, tree.phylo$node.label))
  df = data.frame(labelmat = c(labelmat),
                  ref.labels = rep(ref.labels, p))
  df = rbind(df, data.frame(labelmat = 'all', ref.labels))

  # calculate per type percentage
  pct = aggregate(df$ref.labels, by = list(node = df$labelmat), FUN = function(x) {
    t = table(x)
    t/sum(t)
  })
  pct = data.frame(pct$x, row.names = pct$node, stringsAsFactors = F)
  pct = transform(merge(pct,ord,by="row.names",all=TRUE), row.names=Row.names, Row.names=NULL) # use transform to remove the rownames

  # set the node size
  nodesize = aggregate(df$labelmat, by = list(node = df$labelmat), FUN = function(x) length(x))
  nodesize = data.frame(nodesize = nodesize$x/max(nodesize$x),
                        node = ord[as.character(nodesize$node),],
                        row.names = ord[as.character(nodesize$node),])
  nodesize$nodesize = nodesize$nodesize^(1/8) * node.size  # rescale to reduce the difference

  # set the major label of the node and tips
  major.labels = data.frame(major.labels = colnames(pct[,colnames(pct)!='node'])[apply(pct[,1:(ncol(pct)-1)], 1, which.max)],
                            node = pct$node,
                            row.names = pct$node)


  # only plot the splits and leaf
  tab = table(tibble::as_tibble(tree.phylo)$parent)
  issplit = setdiff(names(tab[tab>1]),ord['all',1])
  isleaf = 1:ape::Ntip(tree.phylo)
  nodesize = nodesize[c(issplit, isleaf),]
  major.labels = major.labels[c(issplit, isleaf),]
  #  major.labels$major.labels = factor(major.labels$major.labels, levels = label.order)

  # tree with label and size
  tree.plot = tidytree::full_join(tidytree::as.treedata(tree.phylo),
                                  merge(major.labels, nodesize, by="node"),
                                  by = 'node')

  # set the order of labels
  if (!is.null(cols)){
    if (length(cols)!= length(label.order)){
      warnings('Number of color does not match the number of labels!')
    }
  } else {
    cols = scales::hue_pal()(length(label.order)) # approximate the default color
  }

  # plot the tree
  gg = ggtree::ggtree(tree.plot, size = 1)+ ggtree::layout_dendrogram() + xlim(bottom.margin,-110)

  # if (!is.null(flip.branch)){
  #   for (i in 1:length(flip.branch)){
  #     gg =  ggtree::flip(tree_view = gg, node1=which(gg$data$label==flip.branch[[i]][1]),
  #                        node2=which(gg$data$label==flip.branch[[i]][2]))
  #   }
  # }

  if (show.ref.labels){
    gg = gg +
      # geom_point(aes(color = major.labels, size=nodesize), stroke = 0)+
      ggtree::geom_tippoint(aes(color = major.labels, size=nodesize), stroke = 0)+
      ggtree::geom_nodepoint(aes(color = major.labels, size=nodesize), stroke = 0)+
      scale_color_manual(values = cols, labels = label.order, drop = FALSE)
    if (!is.null(tip.labels)){
      if (length(tip.labels)!= sum(gg$data$isTip)){
        stop('Error: leaf labels of different size with number of leaf: ', ape::Ntip(tree.phylo),'!')
      }
      gg = gg+ ggtree::geom_tiplab(aes(x = x+tip.label.dist,
                                       label = c(tip.labels[rank(gg$data$y[gg$data$isTip])],
                                                 rep(NA, sum(!gg$data$isTip)))),
                                   angle = 270, color='black')
    } else {
      gg = gg + ggtree::geom_tiplab(aes(x = x+tip.label.dist, label = major.labels),  angle = 270, color='black')
    }

    if (show.branch.labels){
      gg = gg + ggtree::geom_nodelab(aes(x = x-10, label=label),  angle = 0, color='black')+
        ggtree::geom_tiplab(aes(x = x-10, label=label),  angle = 0, color='black')
    }
    gg = gg + guides(colour = guide_legend(override.aes = list(size=5)), size=FALSE)+
      labs( color = legend.title)
  }

  # if (plot.piechart){
  #   pies = ggtree::nodepie(pct, cols = 1:(ncol(pct)-1), color= cols[order(label.order)])
  #   pies = pies[c(issplit, isleaf)]
  #   piesize = nodesize$nodesize
  #   gg = gg + ggtree::geom_inset(pies, reverse_x = T, height = piesize, width = piesize)
  # }

  gg
}

#' Log message
#'
#' @param level log streshold
#' @param l... log messages
#'
#' @import futile.logger
#' @export
logmsg <- function(..., level = "info") {
  msg = paste0(...)
  if (level == "debug") {
    futile.logger::flog.debug(msg)
  } else if (level == "info") {
    futile.logger::flog.info(msg)
  } else if (level == "warn") {
    futile.logger::flog.warn(msg)
  } else if (level == "error") {
    futile.logger::flog.error(msg)
  } else {
    stop("Log level not found!")
  }
}


# statistical mode, for ordered categorical vectors. 
# returns the first mode in case of multiple modes
stat.mode.ord <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# from: https://github.com/sambofra/bnstruct/blob/master/R/imputation.R
#' @param data a numerical matrix.
#' @param k number of neighbours to be used; for categorical variables 
#'        the mode of the neighbours is used, for continuous variables 
#'        the median value is used instead. Default: 10.
#' @param cat.var vector containing the indices of the variables to be 
#'        considered as categorical. Default: all variables.
#' @param to.impute vector indicating which rows of the dataset are to be imputed. 
#'        Default: impute all rows.
#' @param using vector indicating which rows of the dataset are to be used 
#'        to search for neighbours. Default: use all rows.
#'      
#' @return imputed matrix.
#' 
knn.impute <- function( data, k = 10, cat.var = 1:ncol(data), 
                        to.impute = 1:nrow(data), using = 1:nrow(data) )
{
  n.var <- dim(data)[2]
  num.var <- setdiff(1:n.var,cat.var)
  
  use.data <- data[using,,drop=FALSE] # retain dimensions even for one row
  imp.data <- data[to.impute,,drop=FALSE] 
  
  storage.mode(use.data) <- "double"
  storage.mode(imp.data) <- "double"
  
  use.num.var.max <- apply(use.data[,num.var],2,max,na.rm=TRUE)
  use.num.var.min <- apply(use.data[,num.var],2,min,na.rm=TRUE)
  
  na.cases <- (rowSums(is.na(imp.data)) > 0)
  
  neigh <- rep(0,k)
  storage.mode(num.var) <- "integer"
  
  for( i in 1:nrow(imp.data) )
  {
    if( na.cases[i] )
    {
      num.var.max <- pmax(use.num.var.max, imp.data[i,num.var], na.rm=TRUE)
      num.var.min <- pmin(use.num.var.min, imp.data[i,num.var], na.rm=TRUE)
      num.var.range <- num.var.max - num.var.min
      
      d <- .Call( "bnstruct_heom_dist", imp.data[i,], use.data, 
                  num.var, num.var.range, PACKAGE = "bnstruct" )
      s <- sort(d, index.return=TRUE)$ix
      
      for( j in which(is.na(imp.data[i,])) )
      {
        # find the k closest neighbours with nonmissing value for j
        ind.neigh <- 1
        # if to.impute[i] is in _using_, it is the first of the list
        if (to.impute[i] %in% using)
          ind.s <- 2
        else
          ind.s <- 1
        while( ind.neigh < k + 1)
        {
          if( !is.na(use.data[s[ind.s],j]) )
          {
            neigh[ind.neigh] <- use.data[s[ind.s],j]
            ind.neigh <- ind.neigh + 1
          }
          ind.s <- ind.s + 1
        }
        
        # impute from the neighbours
        if( j %in% cat.var )
          imp.data[i,j] <- stat.mode.ord(neigh)
        else
          imp.data[i,j] <- median(neigh)      
      }
    }
  }
  
  return(imp.data)
}

