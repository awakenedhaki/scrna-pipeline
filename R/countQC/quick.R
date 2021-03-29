# scater::quikckPerCellQCMetrics ===============================================

#' @title Quick quality control using \code{\link[scater]quickPerCellQCMetrics}.
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param ... Additional parameters for \code{\link[scater]{quickPerCellQCMetrics}}.
#'
#' @importFrom scater quikckPerCellQCMetrics
.scater <- function(sce, ...) {
  scater::quickPerCellQCMetrics(sce, ...)
}
