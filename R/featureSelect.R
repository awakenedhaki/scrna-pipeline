# Dependencies =================================================================

# Constants ====================================================================
SEED <- 100
ASSAY.COUNTS <- 'logcounts'

# Helpers ======================================================================

# Functions ====================================================================

#' @title Model gene variance by log2-transformed normalized counts.
#'
#' @details Gene-base trend decomposition of technical and biological
#'   components via fitting of mean-variance trend to model the variance of log2
#'   expression profiles of each genes. This model assumes that most genes do not
#'   vary due to biological heterogeneity.
#'
#'   Default arguments for \code{\link[scran]{modelGeneVar}} may lead to overfitting of the data,
#'   when the high-abundance genes are also highly variable.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param density.weight Logical, use of inverse density weights
#' @param ... Additional \code{\link[scran]{modelGeneVar}} arguments
#'
#' @importFrom scran modelGeneVar
#'
#' @return \code{\link{SingleCellExperiment}}
.vanilla <- function(sce,
	          				 density.weight = TRUE,
					           ...) {
    model <- scran::modelGeneVar(sce,
                                 assay.type = ASSAY.COUNTS,
                                 density.weight = density.weight,
                                 ...)
    logger('Saving mean-variance fitTrendVar data in "diagnostic" directory.')
    save.diagnostic(model, 'featureSelect-vanilla.rds')

    return(model)
  }

#' @title Model gene variance by fitting into poisson distribution.
#'
#' @details Without a standard (ie. spike-in) the technical noise can be modeled by our
#'   understanding of the data-generation process. We assume that technical noise
#'   manifests as a Poisson distribution.
#'
#'   This model assumes that most genes do not vary due to biological heterogeneity.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param seed Random number generator initial state
#' @param ... Additional \code{\link[scran]{modelGeneVarByPoisson}} arguments
#'
#' @importFrom scran modelGeneVarByPoisson
#'
#' @return \code{\link{SingleCellExperiment}}
.poisson <- function(sce, seed = SEED, ...) {
  set.seed(seed)

  model <- scran::modelGeneVarByPoisson(sce, assay.type = ASSAY.COUNTS, ...)

  logger('Saving mean-variance fitTrendPoisson data in "diagnostic" directory.')
  save.diagnostic(model, 'featureSelect-poisson.rds')

  return(model)
}

#' @title Select highly variable genes based on per-gene variation model (dispatcher).
#'
#' @details Noise to signal ratio is a function of the number of highly variable genes
#'   retained, given that the larger \code{n} is, the more biologically disinteresting
#'   genes there will be in the final dataset.
#'
#'   Set \code{n} as the number of genes that you expect to be differentially expressed
#'   in your system. Ballpark a number between 500 to 5,000 to start, if you do not
#'   know how many genes you would expect to be differentially expressed.
#'
#'   Otherwise, select a false discovery rate (FDR) threshold (\code{fdr.threshold}),
#'   where the proportion of false positive genes are expected to have an FDR
#'   above the specified threshold. Number of genes will vary depending on the
#'   false positive error rate and severity of multiple test correction. Genes
#'   could also be ranked by their p-value, and obtain the top hits.
#'
#'   Selecting for genes above a certain variance will maximize noise, but minimize
#'   bias. This is default setting in \code{\link[scran]{getTopHVGs}}.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param fit.trend String, the trend tha data will be fit into
#' @param subset.name String, name highly variable gene subset
#' @param n Integer, number of genes to select
#' @param prop Numeric, proportion of genes to select
#' @param fdr.threshold Numeric, false discovery rate threshold
#'
#' @importFrom scran getTopHVGs
#' @importFrom SingleCellExperiment rowSubset
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
featureSelect <- function(sce,
                          n = NULL,
                          prop = NULL,
                          fdr.threshold = NULL,
                          fit.trend = 'poisson',
                          subset.name = 'subset') {
  logger('Initializing feature selection...')

  if (fit.trend == 'vanilla') {
    logger('Performing mean-variance modeling with fitTrendVar.')
    gene.variance <- .vanilla(sce)
  } else if (fit.trend == 'poisson') {
    logger('Performing mean-variance modeling with fitTrendPoisson.')
    gene.variance <- .poisson(sce)
  }

  # Selecting the top highly variable genes, depending on parameterization
  selected <- scran::getTopHVGs(gene.variance,
                                n = n,
                                prop = prop,
                                fdr.threshold = fdr.threshold)
  logger('Top highly variable genes have been selected.')

  # Stores logical vector to a `rowSubset` of the SingleCellExperiment object
  SingleCellExperiment::rowSubset(sce, field = subset.name) <- selected
  logger('Data has been added as a subset.')

  logger('Saving feature selected data in "processed" directory.')
  save.processed(sce, 'featureSelect.rds')

  return(sce)
}
