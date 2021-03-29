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
