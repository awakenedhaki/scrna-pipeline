# Dependencies =================================================================

# Constants ====================================================================
SEED <- 100

# Helpers ======================================================================

#' @title Find the elbow point in the curve of variance explained by each successive PC. This can be used to determine the number of PCs to retain.
#'
#' @param variance Numeric vector containing the variance explained by each PC.
#'   Should be monotonic decreasing.
#'
#' @details Find the elbow point in the curve of variance explained by each successive PC. This can be used to determine the number of PCs to retain.
#'
#' @return An integer scalar specifying the number of PCs at the elbow point.
#'
#' @author Aaron Lun
.findElbowPoint <- function(variance) {
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }

  # Finding distance from each point on the curve to the diagonal.
  dy <- -diff(range(variance))
  dx <- length(variance) - 1
  l2 <- sqrt(dx^2 + dy^2)
  dx <- dx/l2
  dy <- dy/l2

  dy0 <- variance - variance[1]
  dx0 <- seq_along(variance) - 1

  parallel.l2 <- sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
  normal.x <- dx0 - dx * parallel.l2
  normal.y <- dy0 - dy * parallel.l2
  normal.l2 <- sqrt(normal.x^2 + normal.y^2)

  # Picking the maximum normal that lies below the line.
  # If the entire curve is above the line, we just pick the last point.
  below.line <- normal.x < 0 & normal.y < 0
  if (!any(below.line)) {
    length(variance)
  } else {
    which(below.line)[which.max(normal.l2[below.line])]
  }
}

#' @title Load variance model generated during feature selection.
#'
#' TODO: Select the correct model file if many are present.
.load.variance.model <- function() {
  diagnostic.files <- list.files(here(getwd(), 'data', 'diagnostic'))

  variance.model.idx <- which(grepl('featureSelect-.*\\.rds',
                                     diagnostic.files))
  if (length(variance.model.idx) == 0) {
    stop('Feature selection must be run prior to dimensionality reduction. Ensure that a variance model is present in the diagnostic directory.')
  }
  model <- diagnostic.files[variance.model.idx]

  return(model)
}

#' @title Percent of explained variance
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#'
#' @importFrom SingleCellExperiment reducedDim
#'
#' @return Numeric vector
.explained.variance <- function(sce) {
  attr(SingleCellExperiment::reducedDim(sce), 'percentVar')
}

#' @title Subset highly variably genes
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param subset.name String, name of subset
#'
#' @importFrom SingleCellExperiment rowSubset
#'
#' @return Logical vector
.hvg.subset <- function(sce, subset.name) {
  SingleCellExperiment::rowSubset(sce, field = subset.name)
}

# Functions ====================================================================

#' @title Perform principal component analysis on cells
#'
#' TODO: Refactor nested if statement
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param hvariable.genes String, location of HVGs
#' @param ncomponents Integer, number of components to retain
#' @param infer Logical, if the number of components should be inferred
#' @param variance.model DataFrame, Generated variance model from feature selection
#'
#' @importFrom scater runPCA
#' @importFrom scran denoisePCA
#' @importFrom SingleCellExperiment reducedDim
#'
#' @return \code{\link{SingleCellExperiment}} object
.pca <- function(sce,
		             hvariable.genes,
		             ncomponents,
		             infer,
		             variance.model) {
  logger('Running PCA.')
  # TODO: Add `ifelse` in case `altExp` was used insted of `rowSubset`
  sce <- scater::runPCA(sce,
                        exprs_values = 'logcounts',
                        subset_row = hvariable.genes)

  # Infer the number of components using the elbow method
  if (!is.null(infer)) {
    if (infer == 'elbow') {
      logger('Inferring number of components by "elbow" method.')
      ncomponents <- ncomponents <- .findElbowPoint(.explained.variance(sce))
    } else if (infer == 'denoise') {
      logger('Inferring number of components by "denoise" method.')
      ncomponents <- ncol(scran::denoisePCA(sce,
                                            technical = variance.model,
                                            subset.row = hvariable.genes))
    }
  }

  logger('Saving PCA data in diagnostics directory.')
  save.diagnostic(SingleCellExperiment::reducedDim(sce, 'PCA'), 'reducedDim-pca.rds')

  # Subset PCA using specified or infered number of components
  SingleCellExperiment::reducedDim(sce, 'PCA') <- SingleCellExperiment::reducedDim(sce, 'PCA')[ , 1:ncomponents]

  return(sce)
}

#' @title Conventional dimensionality reduction (PCA, tSNE, and UMAP)
#'
#' @description Performs PCA on read count matrix and selects user-defined, or
#'   inferred, top number of components. The selected PCA subspace is then used
#'   to initialize tSNE and UMAP algorithms.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param subset.name String, name of subset to be used
#' @param pca.ncomponents Integer, number of components to retain
#' @param pca.infer String, method of inference for ncomponents
#' @param ... Additional arguments for .pca
#' @param tsne.perplexity Numeric, perplexity (size of search radius)
#' @param tsne.ncomponents Integer, number of components
#' @param umap.min.dist Numeric, minimum distance between two nodes
#' @param umap.nneighbours Integer, number of neighbours
#' @param umap.metric String, distance metric
#' @param umap.ncomponents Integer, number of components
#'
#' @importFrom scater runTSNE runUMAP
#'
#' @return \code{\link{SingleCellExperiment}} object
.reducedVanilla <- function(sce,
                            subset.name,
                            pca.ncomponents = 30,
                            pca.infer = NULL,
                            ...,
                            tsne.perplexity = 20,
                            tsne.ncomponents = 3,
                            umap.min.dist = 0.5,
                            umap.nneighbours = 15,
                            umap.metric = 'euclidean',
                            umap.ncomponents = 3) {
  variance.model <- .load.variance.model()
  hvariable.genes <- .hvg.subset(sce, subset.name = subset.name)

  sce <- .pca(sce,
  	    		  ncomponents = pca.ncomponents,
  	    		  infer = pca.infer,
  			      hvariable.genes = hvariable.genes,
  			      variance.model = variance.model,
  			      ...)
  logger('Completed PCA.')

  logger('Running TSNE.')
  sce <- scater::runTSNE(sce,
                         dimred = 'PCA',
                         exprs_values = 'logcounts',
                         ncomponents = umap.ncomponents)
  logger('Completed TSNE.')
  logger('Running UMAP.')
  sce <- scater::runUMAP(sce,
                         dimred = 'PCA',
                         exprs_values = 'logcounts',
                         ncomponents = umap.ncomponents,
                         min_dist = umap.min.dist,
                         n_neighbors = umap.nneighbours,
                         metric = umap.metric)
  return(sce)
}

#' @title Dimensionality reduction (dispatcher)
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param seed Integer, random number generator initial state
#' @param method String, method of dimensionality reduction
#' @param subset.name String, name of subset to be used
#' @param ... Additional parameters for dispatched functions
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
dimReduced <- function(sce,
                       seed = SEED,
                       method = 'vanilla',
                       subset.name = 'subset',
                       ...) {
  set.seed(seed)
  logger('Initializing dimensionality reduction...')

  if (method == 'vanilla') {
    sce <- .reducedVanilla(sce, subset.name = subset.name, ...)
  }

  logger('Saving dimensionality reduction data in "processed" directory.')
  save.processed(sce, 'dimReduced.rds')

  logger('Dimensionality reduction has ended.')
  return(sce)
}
