# Constants ====================================================================
BASE.CRITERIA <- c('^sum$', '^detected$')

# Helpers ======================================================================
# . General ====================================================================
# . . Getters ==================================================================

#' @title Get criteria patterns.
#'
#' @description Get the first elements of the `criteria` list of lists, which corresponds to
#'   a gene symbol pattern (ex. '^MT-').
#'
#' @param criteria A list of lists.
#'
#' @return List of gene symbol patterns.
.get.gene.symbol.patterns <- function(criteria) {
  return(lapply(criteria, FUN = `[[`, 1))
}

#' @title Get criteria ranges.
#'
#' @description Get the second elements of the `criteria` list of lists, which corresponds to
#'   a criteria range. The range is the lower and upper bounds a criteria value in a cell.
#'
#' @param criteria A list of lists.
#'
#' @return List of numeric vectors.
.get.criteria.ranges <- function(criteria) {
  return(lapply(criteria, FUN = `[[`, 2))
}

#' @title Get the columns for the criteria that will be used for quality control.
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param criteria.names A string vector with the user-defined criteria (ex. mito).
#'
#' @importFrom SummarizedExperiment colData
#'
#' @return S4Vectors object, each column representing one criteria.
.get.criteria.data <- function(sce, criteria.names) {
  coldata <- SummarizedExperiment::colData(sce)

  patterns <- .build.criteria.patterns(criteria.names)
  idxs <- vapply(patterns, FUN = function(pattern) {
    .where(pattern, colnames(coldata), matcher = grepl)
  }, numeric(1))

  return(coldata[idxs])
}

# . . Builders =================================================================

#' @title Builds logical vector defining the a subset (ie. criteria).
#'
#' @description Builds a logical vector that corresponds to matches to the gene symbol pattern
#'   (ex. '^MT-') specified by the user. These logical vectors will define which genes (rows)
#'   will be used to compute the criteria values (ie. mitochondrial proportion).
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param gene.symbol.patterns A list of gene symbol patterns.
#'
#' @return A list of gene symbol patterns.
.build.subsets <- function(sce, gene.symbol.patterns) {
  subsets <- base::Map(f = function(pattern) {grepl(pattern, rownames(sce))},
                       gene.symbol.patterns)

  return(subsets)
}

#' @title Build criteria name patterns corresponding to the percent column.
#'
#' @description Build a string pattern for criteria that corresponds to the calculated
#'   percentage that criteria occupy within a cell (ie. mitochondrial count proportion).
#'
#' @param criteria.names String vector of criteria names.
#'
#' @return List of modified criteria names. The list is appended with `BASE.CRITERIA`,
#'   which include general criteria (ex. `sum`).
.build.criteria.patterns <- function(criteria.names) {
  user.specified <- vapply(criteria.names, FUN = function(name) {
    paste0(name, '_percent')
  }, FUN.VALUE = character(1))

  return(append(BASE.CRITERIA, user.specified))
}

# . . Miscellaneous ============================================================

#' @title Collapses a logical matrix, row-wise, into a logical vector.
#'
#' @description Performs a \code{\link{|}} operation using row-wise sums. If the value is greater than
#'   \code{0}, then the result is \code{TRUE}. All \code{TRUE} elements correspond to cells that are to be
#'   discarded.
#'
#' @param mtx A logical matrix.
#'
#' @seealso [Reduce] for a conceptual understanding.
#'
#' @return A logical vector.
.reduce.flag.matrix <- function(mtx) {
  return(!!base::rowSums(mtx))
}
