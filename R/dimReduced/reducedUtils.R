#' @importFrom SingleCellExperiment rowSubset
.get.hv.genes <- function(sce, subset.name) {
  SingleCellExperiment::rowSubset(sce, field = subset.name)
}

.read.variance.model <- function() {
  diagnostic.files <- list.files(here(getwd(), 'data', 'diagnostic'))
  variance.model.idx <- .where('featureSelect-.*\\.rds', diagnostic.files, matcher = grepl)

  if (length(variance.model.idx) == 0) {
    stop('Fitted variance model data must be within diagnostic directory.')
  }

  model <- readRDS(diagnostic.files[variance.model.idx])
  return(model)
}

#' @importFrom SingleCellExperiment reducedDim
.explained.variance <- function(sce) {
  attr(SingleCellExperiment::reducedDim(sce), 'percentVar')
}
