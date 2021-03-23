# Dependencies =================================================================

# Constants ====================================================================
# . General
ASSAY.TYPE <- 'counts'

# . scran::quickCluster
SEED <- 100

# Helpers ======================================================================

#' @title Maximum window size.
#'
#' @description The number of cells to be pooled together.
#'
#' @param ncolumns Integer, number of columns
#' @param min.cells Integer, minimum number of cells
#' @param min.window.size Numeric, minimum upper limit of window size
#' @param proportion Numeric, proportion of cells to be pooled given ncolumns
#'
#' @return Integer
.max.window.size <- function(ncolumns,
                             min.cells,
                             min.window.size,
                             proportion) {
  # Minimum size for the max size of the window
  min.size <- min(min.cells, floor(ncolumns * proportion))

  # Select the max size of the window
  max.window <- min(min.window.size, min.size + 1)

  return(max.window)
}

# Functions ====================================================================

#' @title Scaling normalization by library size.
#'
#' @details May not yield accurate normalized expression values for downstream analysis, but
#'   sufficient for exploratory analysis.
#'
#'   Each cells is a library, and the size of the library is the total sum of counts
#'   across all genes in that cell. The library size factor is directly proportional to
#'   library size, relative to the mean size factor across all cells equaling 1.
#'
#'   Normalized values will be on the same scale as the original counts
#'
#'   Any upregulation is met with an eqaul magnitude downregulation, since this method
#'   assumes that the dataset is perfectly normal. This prevents composition effects.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param ... \code{\link[scater]{librarySizeFactors}} arguments
#'
#' @importFrom scater librarySizeFactors
#'
#' @return \code{\link{SingleCellExperiment}} object
.normVanilla <- function(sce, ...) {
  # Calculate size factors
  sce <- scater::librarySizeFactors(sce, ...)
  logger('Size factors have been calculated.')

  # Calculate log2-transformed normalized counts
  # . Divide gene counts by their respective size factors
  # . Log2-transform normalized counts
  # . New counts are stored as `logcounts`
  sce <- scater::logNormCounts(sce)
  logger('Counts have been normalized and log2-transformed.')

  return(sce)
}

#' @title Pooled-counts based normalization.
#'
#' @description Counts from multiple cells are pooled to increase the size of the counts for
#'   accurate size factor estimation. Pooling is done by clustering of cells, using the
#'   \code{igraph} method. Each cluster is normalized separately and size factors are rescaled
#'   to make clusters comparable. Clustering is performed on a PCA subspace, and must #' be seeded.
#'
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param seed Random number generator initial state
#' @param method String, clustering method
#' @param min.size Integer, minimum number of clusters
#' @param use.ranks Logical, use rank matrix for clustering
#' @param bsparam Algorithm to be used for PCA
#' @param min.cells Integer, TODO
#' @param min.proportion Numeric, TODO
#' @param min.window.size Integer, TODO
#' @param min.mean Numeric, TODO
#' @param steps Integer, TODO
#' @param initial.pool Integer, TODO
#'
#' @importFrom BiocSingular IrlbaParam
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom scater logNormCounts
#' @importFrom SingleCellExperiment sizeFactors
#'
#' @return \code{\link{SingleCellExperiment}} object
.deconvolution <- function(sce,
					  	   seed = SEED,
					   	   method = 'igraph',
					       min.size = 10,
		          		   use.ranks = FALSE,
		          		   bsparam = BiocSingular::IrlbaParam(),
		          		   min.cells = 150,
					       min.proportion = 0.3,
		          	       min.window.size = 101,
		          		   min.mean = 0.1,
		          		   steps = 5,
		          		   initial.pool = 21) {
  # Seeding random number generator
  set.seed(seed)

  # Clustering cells to pool counts
  clustered.sce <-  scran::quickCluster(sce,
                                        assay.type = ASSAY.TYPE,
                                        method = method,
                                        min.size = min.size,
                                        use.ranks = use.ranks,
                                        BSPARAM = bsparam)
  logger('Clusters have been generated.')

  logger('Saving cluster table in "diagnostic" directory.')
  save.diagnostic(clustered.sce, 'normQC-clustered.rds')

  # Max window size determining the number of cells per pool
  max.window <- .max.window.size(ncol(sce),
                                 min.cells = min.cells,
                                 min.window.size = min.window.size,
                                 proportion = min.proportion)

  # Deconvoluted size factors from pooled cells
  sce <- scran::computeSumFactors(sce,
                                  cluster = clustered.sce,
                                  assay.type = ASSAY.TYPE,
                                  sizes = seq(initial.pool, max.window, by = steps),
                                  min.mean = min.mean)
  logger('Sum factors have been computed.')

  if (summary(sizeFactors(sce))[['Min.']] <= 0) {
    logger('Negative size factors encountered. Check number of cells, and revisit QC.')
    stop('Cannot have negative size factors.')
  }

  # Calculate log2-transformed normalized counts
  # . Divide gene counts by their respective size factors
  # . Log2-transform normalized counts
  # . New counts are stored as `logcounts`
  sce <- scater::logNormCounts(sce)
  logger('Normalized counts have been log2-transformed.')

  return(sce)
}

#' @title Normalization of read counts (dispatcher)
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param method String, normalization method
#' @param ... Additional arguments specific to normalization method
#'
#' @return \code{\link{SingleCellExperiment}}
#' @export
normCounts <- function(sce, method = 'deconvolution', ...) {
  logger('Initializing normalization of counts...')
  if (method == 'vanilla') {
    logger('Starting library size normalization.')
    sce <- .normVanilla(sce, ...)
  } else if (method == 'deconvolution') {
    logger('Starting library size normalization by deconvolution.')
    sce <- .deconvolution(sce, ...)
  }

  logger('Saving normalized data in "processed" directory.')
  save.processed(sce, 'normQC.rds')

  logger('Normalization has ended.')
  return(sce)
}
