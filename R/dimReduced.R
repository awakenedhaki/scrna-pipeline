dimReduced <- function(sce, protocol, seed = 100, params) {
  set.seed(seed)
  if (protocol == 'vanilla') {
    sce <- .reduced.vanilla(sce, params)
  }

  return(sce)
}

# Vanilla ======================================================================

#' @importFrom scater runTSNE runUMAP
.reduced.vanilla <- function(sce, params) {
  sce <- .pca(sce, hv.genes = hv.genes, params$pca)
  sce <- .kwargs(scater::runTSNE, sce, params$tsne)
  sce <- .kwargs(scater::runUMAP, sce, params$umap)
}

#' @importFrom scater runPCA
#' @importFrom PCAtools findElbowPoint
#' @importFrom scran denoisePCA
#' @importFrom SingleCellExperiment reducedDim
.pca <- function(sce, hv.genes, params) {
  params$scater$subset_row <- .get.hvg.subset(sce, params$subset.name)
  sce <- .kwargs(scater::runPCA, sce, params$scater)

  if ('infer' %in% names(params) & !(ncomponents %in% names(params))) {
    if (params$infer == 'elbow') {
      ncomponents <- findElbowPoint(.explained.variance(sce))
    } else if (params$infer == 'denoise') {
      variance.model <- .read.variance.model(...)
      ncomponents <- ncol(scran::denoisePCA(sce,
                                            technical = variance.model,
                                            subset.row = hv.genes))
    }
  }

  SingleCellExperiment::reducedDim(sce, 'PCA') <-
    SingleCellExperiment::reducedDim(sce, 'PCA')[, 1:params$ncomponents]
  return(sce)
}

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
