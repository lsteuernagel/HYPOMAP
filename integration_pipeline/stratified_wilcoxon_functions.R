# batch-aware wilcox / vanelteren:

# From: https://github.com/KChen-lab/stratified-tests-for-seurat/blob/master/R/stratified_test_seurat_rename.R

## Paper: https://www.biorxiv.org/content/10.1101/2021.06.08.447617v1.full
## git: https://github.com/KChen-lab/stratified-tests-for-seurat

## all credit to authors

# ve_test_core modified from r-source/src/library/stats/R/wilcox.test.R
# Other functions modified from Seurat/R/differential_expression.R
# Per GNU GPL v3.0, modification and distribution is allowed.

library(pbapply)
library(future.apply)
library(Matrix)

#`%||%` <- Seurat:::`%||%`
WilcoxDETest <- Seurat:::WilcoxDETest
CheckDots <- Seurat:::CheckDots
GLMDETest <- Seurat:::GLMDETest
DiffExpTest <- Seurat:::DiffExpTest
MarkerTest <- Seurat:::MarkerTest
LRDETest <- Seurat:::LRDETest
MASTDETest <- Seurat:::MASTDETest
DESeq2DETest <- Seurat:::DESeq2DETest

# Compared with the original Wilcoxon U-test,
# we drop the paired test, the mu parameter,
# and only retains the rank sum test (i.e. drop signed rank test)
ve_test_core <- function (x, y, correct = TRUE)
{
  alternative = "two.sided"
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (!is.numeric(y))
    stop("'y' must be numeric")

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 1L)
    return(list(U = 0, V = 0))

  if (length(y) < 1L)
    return(list(U = 0, V = 0))

  r <- rank(c(x, y))
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))

  CORRECTION <- 0
  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)

  TIES <- (length(r) != length(unique(r)))
  NTIES <- table(r)

  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 1))))
  if (correct) CORRECTION <- sign(z) * 0.5

  # instead of continue and calculate the wilcoxon statistics, we return the U and V needed by Van Eteren test.
  list(U = unname(z - CORRECTION), V = SIGMA, U_raw = unname(STATISTIC))
}

ve_test <- function (x, label, strata, ...)
{
  param <- list(...)
  if ('correct' %in% names(param)) correct <- param$correct
  else correct <- TRUE
  if ('genre' %in% names(param)) genre <- param$genre
  else genre <- 'locally-best'

  stopifnot(genre %in% c('locally-best', 'design-free'))

  Us = c()
  U_raws = c()
  Vs = c()
  Ms = c()
  Ns = c()
  Ss = c()
  for (s in unique(strata))
  {
    mask1 <- label$group == "Group1" & strata == s
    mask2 <- label$group == "Group2" & strata == s
    temp <- ve_test_core(x=x[mask1],
                         y=x[mask2],
                         correct=correct)
    Us = c(Us, temp$U)
    U_raws = c(U_raws, temp$U_raw)
    Vs = c(Vs, temp$V)
    Ms = c(Ms, sum(mask1))
    Ns = c(Ns, sum(mask2))
    Ss = c(Ss, s)
  }

  if (genre == 'locally-best'){
    res <- list(p.value = pchisq(sum(Us / (Ms + Ns + 1)) ^ 2 / sum((Vs / (Ms + Ns + 1)) ^ 2),
                                 df = 1, lower.tail = F),
                effect.size = sum(U_raws / (Ms + Ns + 1)) / sum(Ms * Ns / (Ms + Ns + 1)))
  } else if (genre == 'design-free'){
    res <- list(p.value = pchisq(sum(Us / (Ms * Ns)) ^ 2 / sum((Vs / (Ms * Ns)) ^ 2),
                                 df = 1, lower.tail = F),
                effect.size = sum(U_raws / (Ms * Ns)) / length(Ms))
  }
  return(res)
}

VEDETest <- function (data.use, cells.1, cells.2, latent.vars, verbose = TRUE, ...)
{
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(x = group.info), drop = FALSE]

  latent.vars <- latent.vars[rownames(x = group.info), , drop = FALSE]
  latent.vars <- do.call(paste, latent.vars)

  my.sapply <- ifelse(test = verbose && nbrOfWorkers() == 1,
                      yes = pbsapply, no = future_sapply)

  test_res <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
    return(ve_test(data.use[x, ], group.info, latent.vars, ...))
  })

  p.val <- unlist(test_res['p.value', ])
  effect.size <- unlist(test_res['effect.size', ])
  #print(test_res)
  #print(p.val)
  #print(effect.size)
  #print(rownames(x = data.use))

  return(data.frame(p_val=p.val,
                    effect_size=effect.size,
                    row.names = rownames(x = data.use)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param counts Count matrix if using scale.data for DE tests. This is used for
#' computing pct.1 and pct.2 and for filtering features based on fraction
#' expressing
#' @param features Genes to test. Default is to use all genes
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.25
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  \item{"negbinom"} : Identifies differentially expressed genes between two
#'   groups of cells using a negative binomial generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"poisson"} : Identifies differentially expressed genes between two
#'   groups of cells using a poisson generalized linear model.
#'   Use only for UMI-based datasets
#'  \item{"LR"} : Uses a logistic regression framework to determine differentially
#'  expressed genes. Constructs a logistic regression model predicting group
#'  membership based on each feature individually and compares this to a null
#'  model with a likelihood ratio test.
#'  \item{"MAST"} : Identifies differentially expressed genes between two groups
#'  of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST
#'  package to run the DE testing.
#'  \item{"DESeq2"} : Identifies differentially expressed genes between two groups
#'  of cells based on a model using DESeq2 which uses a negative binomial
#'  distribution (Love et al, Genome Biology, 2014).This test does not support
#'  pre-filtering of genes based on average difference (or percent detection rate)
#'  between cell groups. However, genes may be pre-filtered based on their
#'  minimum detection rate (min.pct) across both cell groups. To use this method,
#'  please install DESeq2, using the instructions at
#'  https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#'  \item{"VE"} : Stratified Wilcoxon Test
#' }
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.1
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos Only return positive markers (FALSE by default)
#' @param verbose Print a progress bar once expression testing begins
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test, used only when \code{test.use} is one of
#' 'LR', 'negbinom', 'poisson', 'MAST', or 'VE'
#' @param min.cells.feature Minimum number of cells expressing the feature in at least one
#' of the two groups, currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#'
#' @importFrom Matrix rowSums
#' @importFrom stats p.adjust
#'
#' @rdname FindMarkers2
#' @export
#' @method FindMarkers2 default
#'
FindMarkers2.default <- function(
  object,
  slot = "data",
  counts = numeric(),
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  reduction = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  base = 2,
  ...
) {
  features <- features %||% rownames(x = object)
  methods.noprefiliter <- c("DESeq2")
  if (test.use %in% methods.noprefiliter) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  # feature selection (based on percentages)
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  #message("   FindMarkers2.default: ","Calculate pcts")
  # subset only once to speed up!
  data_for_cells1 = data[features, cells.1, drop = FALSE]
  data_for_cells2 = data[features, cells.2, drop = FALSE]

  if (is.null(x = reduction)) {
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data_for_cells1 > thresh.min) /
        length(x = cells.1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data_for_cells2 > thresh.min) /
        length(x = cells.2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > min.pct))
    if (length(x = features) == 0) {
      stop("No features pass min.pct threshold")
    }
    alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
    features <- names(
      x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
    )
    if (length(x = features) == 0) {
      stop("No features pass min.diff.pct threshold")
    }
  } else {
    data.alpha <- data.frame(
      pct.1 = rep(x = NA, times = length(x = features)),
      pct.2 = rep(x = NA, times = length(x = features))
    )
  }
  #message("   FindMarkers2.default: ","Calculate mean diff")
  # feature selection (based on average difference)
  # edited this to resemble the latest seurat versions (skip apply to avoid implicit cast to dense ! )
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use,base = base))
      },
      function(x) {
        return(log(x = rowMeans(x = x) + pseudocount.use,base = base))
      }
    )
  } else {
    rowMeans
  }
  # also updated mean calc
  data.1 = mean.fxn(data_for_cells1)
  data.2 = mean.fxn(data_for_cells2)
  # get total diff
  total.diff <- (data.1 - data.2)
  if (is.null(x = reduction) && slot != "scale.data") {
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff > logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) > logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      stop("No features pass logfc.threshold threshold")
    }
  }
  #message("   FindMarkers2.default: ","Subset to max.cells.per.ident")
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  # perform DE
  #message("   FindMarkers2.default: ","perform DE")
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST', "LR", 'VE')) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for 'negbinom', 'poisson', 'LR', 'VE', and 'MAST' tests",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!test.use %in% c('wilcox', 'VE', 'MAST', 'DESeq2')) {
    CheckDots(...)
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    'bimod' = DiffExpTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'roc' = MarkerTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    't' = DiffTTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    'negbinom' = GLMDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'poisson' = GLMDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    'MAST' = MASTDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = DESeq2DETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = LRDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    "VE" = VEDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    stop("Unknown test: ", test.use)
  )
  if (is.null(x = reduction)) {
    diff.col <- ifelse(
      test = slot == "scale.data" || test.use == 'roc',
      yes = "avg_diff",
      no = "avg_logFC"
    )
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
    de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  } else {
    diff.col <- "avg_diff"
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
  }
  if (only.pos) {
    de.results <- de.results[de.results[, diff.col] > 0, , drop = FALSE]
  }
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, diff.col]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, diff.col]), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }
  return(de.results)
}

#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#'
#' @importFrom methods is
#'
#' @rdname FindMarkers2
#' @export
#' @method FindMarkers2 Seurat
#'
FindMarkers2.Seurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  base = 2,
  ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = 'counts',
    no = slot
  )
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)
  } else {
    if (data.slot == "counts") {
      stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
    }
    data.use <- t(x = Embeddings(object = object, reduction = reduction))
  }
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
  }
  counts <- switch(
    EXPR = data.slot,
    'scale.data' = GetAssayData(object = object[[assay]], slot = "counts"),
    numeric()
  )
  de.results <- FindMarkers2.default(
    object = data.use,
    slot = data.slot,
    counts = counts,
    cells.1 = ident.1,
    cells.2 = ident.2,
    features = features,
    reduction = reduction,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    verbose = verbose,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    latent.vars = latent.vars,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
  return(de.results)
}
