#' @title Dispatcher function for clustering protocols.
clustering <- function(sce, protocol, params) {
  if (protocol == 'graph') {
    sce <- .graph(sce, params)
  } else if (protocol == 'kmeans') {
    sce <- .kmeans(sce, params)
  }

  return(sce)
}

# Graph-based clustering =======================================================

#' @import igraph
#' @importFrom scran buildSNNGraph
#' @importFrom SingleCellExperiment reducedDim colLabels
.graph <- function(sce, params) {
  # TODO: Change naming scheme allowing for different graph-based approaches.
  graph <- .kwargs(scran::buildSNNGraph, sce, params$graphing)

  if (!missing(layout)) {
    sce <- .generate.layout(sce, graph, params$layout)
  }

  f.clust <- .get.cluster.method(params$method)
  clusters <- f.cluster(graph)$membership
  save.diagnostic(clusters, paste0('clustering-graph-', params$grpah$k))

  SingleCellExperiment::colLabels(sce) <- factor(clusters)
  return(sce)
}

#' @import igraph
#' @importFrom SingleCellExperiment reducedDim
.generate.layout <- function(sce, graph, layout) {
  f <- base::get(paste0('with_', layout),
                 envir = as.environment('package:igraph'))

  field.name <- paste('force', layout, sep = '-')
  SingleCellExperiment::reducedDim(sce, field.name) <- igraph::layout_(graph, layout = f())

  return(sce)
}

#' @import igraph
.get.cluster.method <- function(method) {
  function.name <- paste0('cluster_', method)
  f <- base::get(function.name,
                 envir = as.environment('package:igraph'))
  return(f)
}

# K-means clustering ===========================================================

#' @importFrom SingleCellExperiment reducedDim colLabels
.kmeans <- function(sce, params) {
  clust.kmeans <- stats::kmeans(SingleCellExperiment::reducedDim(sce, params$field),
                                clusters = params$k)

  save.diagnostic(clust.kmeans, paste0('clustering-kmeans-', params$k))
  SingleCellExperiment::colLabels(sce) <- factor(clust.kmeans$cluster)
  return(sce)
}
