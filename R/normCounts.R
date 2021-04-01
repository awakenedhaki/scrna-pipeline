#' @title Dispatcher function for normalization of read counts.
normCounts <- function(sce, protocol, seed = 100, ...) {
  set.seed(seed)
	if (protocol == 'vanilla') {
		norm.sce <- .norm.vanilla(sce, ...)
	} else if (protocol == 'deconvolution') {
		norm.sce <- .deconvolution(sce, ...)
	}

  save.processed(norm.sce, paste0('norm-', protocol))

	return(norm.sce)
}

# Vanilla ======================================================================

#' @importFrom scater librarySizeFactors logNormCounts
.norm.vanilla <- function(sce, ...) {
	norm.sce <- scater::librarySizeFactors(sce, ...) %>%
		scater::logNormCounts()

	return(norm.sce)
}

# Deconvolution ================================================================

#' @title Deconvolution-based read count normalization based on pooling of cells
#'
#' @description TODO
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param clustering A named list of parameters for \code{\link[scran]{quickCluster}}.
#'   * \code{method}: The method of clustering to be performed (ie. \code{'igraph'}, \code{'hclust'})
#'   * \code{min.size}: The minimum number of clusters.
#'   * \code{use.ranks}: Whether or not to use the rank matrix for clustering.
#'   * \code{BSPARAM}: The PCA algorithm to be used, selected from \code{link{BiocSingular}}.
#' @param deconv A named list of parameters for \code{\link[scran]{computeSumFactors}}.
#'   * \code{assay}: The counts matrix to be used for normalization.
#'   * \code{min.window.size}: The minimum number of cells in one pool.
#'   * \code{max.window.size}: A list of parameters determining the maximum number of cells in a pool.
#'     * \code{ncells}: Number of cells in pool.
#'     * \code{size}: User-specified maximum window size.
#'     * \code{proportion}: Determining max size by proportion of cells
#'   * \code{steps}: Sliding window.
#'   * \code{min.mean}: The minimum mean per pool of cells.
#'
#' @importFrom BiocSingular IrlbaParam
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom scater logNormCounts
#'
#' @return \code{\link{SingleCellExperiment}} object.
.deconvolution <- function(sce,
                           clustering = list(method = 'igraph',
                                             min.size = 10,
                                             use.ranks = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam()),
                           deconv = list(assay = 'counts',
                                         min.window.size = 21,
                                         max.window.size = list(ncells = 150,
                                                                size = 101,
                                                                proportion = 0.3),
                                         steps = 5,
                                         min.mean = 0.1)) {
  clust.sce <- .kwargs(
    scran::quickCluster,
    sce,
    clustering
  )

  save.diagnostic(table(clust.sce), paste0('norm-deconv-clust'))

	max.window.size <- .kwargs(
	  .get.max.window.size,
	  ncol(sce),
	  deconv$max.window.size
	)
	deconv.sce <- scran::computeSumFactors(sce,
	                                       cluster = clust.sce,
	                                       assay.type = deconv$assay,
	                                       sizes = seq(deconv$min.window.size,
	                                                   max.window.size,
	                                                   deconv$steps),
	                                       min.mean = deconv$min.mean)
  deconv.sce <- scater::logNormCounts(deconv.sce)

	return(deconv.sce)
}

.get.max.window.size <- function(ncolumns,
                                 ncells,
                                 size,
                                 proportion) {
  # Minimum size for the max size of the window
  infer.size <- min(ncells, floor(ncolumns * proportion))

  # Select the max size of the window
  max.window <- min(size, infer.size + 1)

  return(max.window)
}

# . Contants ===================================================================

CLUSTERING <- list(method = 'igraph',
                   min.size = 10,
                   use.ranks = FALSE,
                   BSPARAM = BiocSingular::IrlbaParam())

DECONV <- list(assay = 'counts',
               min.window.size = 21,
               max.window.size = list(ncells = 150,
                                      size = 101,
                                      proportion = 0.3),
               steps = 5,
               min.mean = 0.1)
