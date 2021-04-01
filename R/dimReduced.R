#' TODO: documentation
dimReduced <- function(sce,
                       protocol,
                       seed = 100,
                       params) {
  IDENTIFIER <<- .get.identifier(sce)

  set.seed(seed)
  if (protocol == 'vanilla') {
    sce <- .reduced.vanilla(sce, params)
  }

  save.processed(sce, paste('dimReduced', protocol, IDENTIFIER, sep = '-'))
  rm(IDENTIFIER, pos = ".GlobalEnv")
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
    variance.model <- .read.variance.model(IDENTIFIER)
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
#' @importFrom SingleCellExperiment reducedDim
.explained.variance <- function(sce) {
  . <- attr(SingleCellExperiment::reducedDim(sce), 'percentVar')
  save.diagnostic(., paste('dimReduced-explained-variance', IDENTIFIER, sep = '-'))
  return(.)
}
