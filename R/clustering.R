#' @title Dispatcher function for clustering protocols.
clustering <- function(sce, protocol, params) {
  IDENTIFIER <<- .get.identifier(sce)

  if (protocol == 'graph') {
    sce <- .graph(sce, params)
  } else if (protocol == 'kmeans') {
    sce <- .kmeans(sce, params)
  }

  save.processed(sce, paste('clustering', 'protocol', IDENTIFIER, sep = '-'))

  rm(IDENTIFIER, pos = ".GlobalEnv")
  return(sce)
}

# Graph-based clustering =======================================================

#' @import igraph
#' @importFrom scran buildSNNGraph
#' @importFrom SingleCellExperiment reducedDim colLabels
.graph <- function(sce, params) {
  # TODO: Change naming scheme allowing for different graph-based approaches.
  graph <- .kwargs(scran::buildSNNGraph, sce, params$graphing)

  if (('layout' %in% names(params))) {
    sce <- .generate.layout(sce, graph, params$layout)
  }

  f.clust <- .get.cluster.method(params$method)
  clusters <- f.clust(graph)$membership
  save.diagnostic(clusters, paste('clustering-graph', params$graph$k, IDENTIFIER, sep = '-'))

  SingleCellExperiment::colLabels(sce) <- factor(clusters)
  return(sce)
}

#' @import igraph
#' @importFrom SingleCellExperiment reducedDim
.generate.layout <- function(sce, graph, layout) {
  function.name <- paste0('with_', layout)
  f <- .blackbox.namespace(function.name)

  field.name <- paste('force', layout, sep = '-')
  SingleCellExperiment::reducedDim(sce, field.name) <- igraph::layout_(graph, layout = f())

  return(sce)
}

#' @import igraph
.get.cluster.method <- function(method) {
  function.name <- paste0('cluster_', method)
  f <- .blackbox.namespace(function.name)
  return(f)
}

# K-means clustering ===========================================================

#' @importFrom SingleCellExperiment reducedDim colLabels
.kmeans <- function(sce, params) {
  clust.kmeans <- stats::kmeans(SingleCellExperiment::reducedDim(sce, params$field),
                                clusters = params$k)

  save.diagnostic(clust.kmeans, paste('clustering-kmeans', params$k, IDENTIFIER, sep = '-'))
  SingleCellExperiment::colLabels(sce) <- factor(clust.kmeans$cluster)
  return(sce)
}
