#' TODO: documentation
dimReduced <- function(sce,
                       protocol,
                       seed = 100,
                       params) {
  set.seed(seed)
  if (protocol == 'vanilla') {
    sce <- .reduced.vanilla(sce, params)
  }

  return(sce)
}

# Vanilla ======================================================================

#' TODO: documentation
#' @importFrom scater runTSNE runUMAP
.reduced.vanilla <- function(sce, params) {
  sce <- .pca(sce, params$pca)
  sce <- .kwargs(scater::runTSNE, sce, params$tsne)
  sce <- .kwargs(scater::runUMAP, sce, params$umap)
}

#' TODO: documentation
#' @importFrom scater runPCA
#' @importFrom PCAtools findElbowPoint
#' @importFrom scran denoisePCA
#' @importFrom SingleCellExperiment reducedDim
.pca <- function(sce, params) {
  params$scater$subset_row <- .get.hvg.subset(sce, params$subset.name)

  pca.sce <- .kwargs(scater::runPCA, sce, params$scater)

  ncomponents <- .get.ncomponents(sce, params)
  SingleCellExperiment::reducedDim(pca.sce, 'PCA') <-
    SingleCellExperiment::reducedDim(pca.sce, 'PCA')[, 1:ncomponents]

  return(pca.sce)
}

#' TODO: documentation
.get.ncomponents <- function(sce, params) {
  if ('infer' %in% names(params) & !('ncomponents' %in% names(params))) {
    ncomponents <- .infer.ncomponents(sce, params)
  } else {
    ncomponents <- pca$ncomponents
  }

  return(ncomponents)
}

#' TODO: documentation
.infer.ncomponents <- function(sce, params) {
  if (params$infer == 'elbow') {
    ncomponents <- findElbowPoint(.explained.variance(sce))
  } else if (params$infer == 'denoise') {
    variance.model <- .read.variance.model()
    # TODO: Finding better way of passing hv.genes logical vector to denoisePCA
    denoised <- scran::denoisePCA(sce,
                                  technical = variance.model,
                                  subset.row = params$scater$subset_row)
    ncomponents <- ncol(reducedDim(denoised))
  }
  return(ncomponents)
}

#' TODO: documentation
#' @importFrom SingleCellExperiment rowSubset
.get.hvg.subset <- function(sce, subset.name) {
  SingleCellExperiment::rowSubset(sce, field = subset.name)
}

#' TODO: documentation
.read.variance.model <- function() {
  # TODO: Handle multiple model files.
  diagnostic.path <- here(getwd(), 'data', 'diagnostic')

  diagnostic.files <- list.files(diagnostic.path)
  variance.model.idx <- .where('featureSelect-.*\\.rds', diagnostic.files, matcher = grepl)
  if (length(variance.model.idx) == 0) {
    stop('Fitted variance model data must be within diagnostic directory.')
  }

  model.path <- here(diagnostic.path, diagnostic.files[variance.model.idx])
  model <- readRDS(model.path)
  return(model)
}

#' TODO: documentation
#' @importFrom SingleCellExperiment reducedDim
.explained.variance <- function(sce) {
  . <- attr(SingleCellExperiment::reducedDim(sce), 'percentVar')
  save.diagnostic(., 'dimReduced-explained-variance')
  return(.)
}
