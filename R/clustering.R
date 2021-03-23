# Dependencies =================================================================

# Constants ====================================================================
SEED <- 100

# Helpers ======================================================================

#' @title Generates layout for force-directed visualization.
#'
#' TODO: Bad practice, call function from string...
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param layout String, layout method
#' @param graph Graph object (\code{\link[igraph]{make_graph}})
#'
#' @import igraph
#' @importFrom SingleCellExperiment reducedDim
#'
#' @return \code{\link{SingleCellExperiment}} object
.generate.layout <- function(sce, graph, layout) {
  f <- get(paste0('with_', layout))

  field.name <- paste('force', layout, sep = '-')
  SingleCellExperiment::reducedDim(sce, field.name) <- igraph::layout_(graph, layout = f())

  return(sce)
}

#' @title Get cluster method from igraph
#'
#' @param cluster String, community detection method
#'
#' @import igraph
#'
#' @return closure
.get.cluster.method <- function(cluster) {
  function.name <- paste0('cluster_', cluster)
  return(get(function.name))
}

# Functions ====================================================================

#' @title Graph-based clustering.
#'
#' TODO: Separate graph generation and community detection.
#' TODO: Assess cluster separation
#'
#' @description Generate a k-nearest neighbours graph for downstream community detection.
#'   \code{\link[igraph]{cluster_walktrap}} detects densely connected subgraphs via a
#'   random walk. The intuition for this algorithm is that random walks tend to get
#'   trapped in regions of high density.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param k Integer, number of nearest neighbours (resolution)
#' @param type String, edge weighting scheme
#' @param cluster String, community detection method
#' @param layout String, method for force-directed graph layout
#' @param dimred String, dimensionality reduced subspace
#'
#' @import igraph
#' @importFrom scran buildSNNGraph
#' @importFrom SingleCellExperiment reducedDim colLabels
#'
#' @return \code{\link{SingleCellExperiment}} object
.graph <- function(sce,
                   k,
                   type = 'number',
                   cluster = 'walktrap',
                   layout = 'fr',
                   dimred = 'PCA') {
  graph <- scran::buildSNNGraph(sce,
                                k = k,
                                type = type,
                                use.dimred = dimred)

  sce <- .generate.layout(sce, graph, layout = layout)

  # Cluster (community) detection
  f.cluster <- .get.cluster.method(cluster)
  clusters <- f.cluster(graph)$membership

  logger('Saving cluster data in diagnostic directory.')
  save.diagnostic(clusters, paste0('clustering-graph-', k, '.rds'))

  SingleCellExperiment::colLabels(sce) <- factor(clusters)
  return(sce)
}

#' @title k-means clustering
#'
#' TODO: Assess cluster gap statistic
#' TODO: Assess cluster separation
#'
#' @details Clustering cells based on centroids, minimizing for within-cluster sum of squares.
#'   In other words, it groups cells that are closests to a given centroid with could be
#'   any given cell or a "proto"-cell.
#'
#'   This method is biased towards normally distributed, spherical, clusters.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param k Integer, number of clusters (centroids)
#' @param field String, dimensionality-reduced subspace to use
#'
#' @importFrom SingleCellExperiment reducedDim colLabels
#'
#' @return \code{\link{SingleCellExperiment}} object
.kmeans <- function(sce, k, field) {
  cluster.kmeans <- stats::kmeans(SingleCellExperiment::reducedDim(sce, field),
                                  clusters = k)
  save.diagnostic(clusters, paste0('clustering-kmeans', k, '.rds'))

  SingleCellExperiment::colLabels(sce) <- factor(cluster.kmeans$cluster)
  return(sce)
}

#' @title Subsclustering
#'
#' TODO: Implement
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param cluster Integer, cluster of interest
.subcluster <- function(sce, cluster) {

}

#' @title Cell clustering (dispatcher)
#'
#' TODO: Add hierarchical clustering
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param method String, method of clustering
#' @param ... Additional parameters for specified clustering method
#'
#' @return \code{\link{SingleCellExperiment}} object
clustering <- function(sce, method = 'graph', ...) {
  if (method == 'graph') {
    sce <- .graph(sce, ...)
  } else if (method == 'kmeans') {
    sce <- .kmeans(sce, ...)
  }

  return(sce)
}
