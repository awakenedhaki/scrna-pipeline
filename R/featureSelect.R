#' @title Dispatcher function for feature selection protocols.
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
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param seed Initial state of random number generator.
#' @param protocol Gene variance model to be fitted.
#'   Current protocols include:
#'   * 'vanilla': \code{\link[scran]{modelGeneVar}}.
#'   * 'poisson': \code{\link[scran]{modelGeneVarByPoisson}}
#' @param model List of parameters for the gene variance models.
#' @param ... Parameters for \code{\link[scran]{getTopHVGs}}
#'
#' @importFrom scran getTopHVGs
#' @importFrom SingleCellExperiment rowSubset
#'
#' @return \code{\link{SingleCellExperiment}} object.
#'   Highly variable genes are added to SingleCellExperiment::rowData.
#' @export
featureSelect <- function(sce, seed = 100, protocol, model, ...) {
  if (protocol == 'vanilla') {
    gene.variance <- .mode.vanilla(sce, model)
  } else if (protocol == 'poisson') {
    gene.variance <- .model.poisson(sce, model)
  }

  hvgs <- scran::getTopHVGs(gene.variance, ...)
  SingleCellExperiment::rowSubset(sce, field = protocol) <- hvgs

  return(sce)
}

# Vanilla ======================================================================

#' @importFrom scran modelGeneVar
.model.vanilla <- function(sce, model) {
  .kwargs(
    scran::modelGeneVar,
    sce,
    model
  )
}

# Poisson ======================================================================

#' @importFrom scran modelGeneVarByPoisson
.model.poisson <- function(sce, model) {
  .kwargs(
    scran::modelGeneVarByPoisson,
    sce,
    model
  )
}
