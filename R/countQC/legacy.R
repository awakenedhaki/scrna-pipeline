# Legacy Thresholds Protocol ===================================================

#' @title The go to quality control protocol in our lab.
#'
#' @details The legacy protocol is a mixture of fixed and adaptive thresholds,
#'   specific to a given criteria.
#'
#'   Mitochondrial and ribosomal gene percentages are subject to user-defined boundaries,
#'   which are specified within the list of lists discussed in \code{countQC}. Cells whose
#'   percentages are within a given range are kept.
#'
#'   Library size (\code{sum}) and number of detected genes (\code{detected}) are
#'   flagged by the \code{\link[scuttle]isOutlier} function. The method relies on median
#'   absolute deviation, which is the median distance any given cell has to the median
#'   value. This method should not be influenced by extreme outliers. The outlier detection
#'   is parameterized to look at the lower tail (ex. small library sizes) in a log2 scale.
#'
#' @param criteria.data An S4Vector with quality control criteria of interest.
#' @param ... Additional parameters for \code{.legacy.flag.cells}.
#'
#' @importFrom scater isOutlier
#' @importFrom stringr str_detect
#' @importFrom purrr partial
#'
#' @return Logical vector, flagging cells to be removed.
.legacy <- function(criteria.data, ...) {
  # TODO: Add nmads into `criteria` param list
  flag.matrix <- sapply(names(criteria.data), FUN = function(name) {
    .legacy.flag.cells(criteria.data[[name]], name, ...)
  })

  flags <- .reduce.flag.matrix(flag.matrix)
  return(flags)
}

#' @title Flag cells to be removes from \code{\link{SingleCellExperiment}} object.
#'
#' @details Individual QC criteria columns are being processed in this function. For example,
#'   the `criteria.column` vector would correspond to the number of detected genes column.
#'
#' @param criteria.column Numeric vector for one QC criteria.
#' @param criteria.name Name of the criteria currnetly being operated on.
#' @param nmads The median absolute deviation for \code{\link[scuttle]{isOutlier}}.
#'
#' @importFrom scater isOutlier
#' @importFrom stringr str_detect
#'
#' @return Logical matrix of cell to be removed.
.legacy.flag.cells <- function(criteria.column, criteria.name, .ranges, nmads) {
  if (.is.base.criteria(criteria.name)) {
    flags <- scater::isOutlier(criteria.column,
                               nmads,
                               type = 'lower',
                               log = TRUE)
  } else {
    idx <- .where(criteria.name, names(.ranges),matcher = stringr::str_detect)
    .range <- .ranges[[idx]]
    flags <- .between(criteria.column, .range)
  }

  return(flags)
}

#' @title Determine if a criteria name is from the base criteria.
#'
#' @param criteria.name A string for the name of a criteria
#'
#' @return logical
.is.base.criteria <- function(criteria.name) {
  . <- grepl(criteria.name, BASE.CRITERIA)
  return(base::Reduce(`|`, .))
}

#' @title Determine if a value is within a given range.
#'
#' @param criteria.data A numeric vector
#' @param .range A numeric vector of length 2.
#'
#' @return A logical vector of cells to be removed.
.between <- function(criteria.data, .range) {
  .range <- unlist(.range)
  return(min(.range) > criteria.data |
           max(.range) < criteria.data
  )
}
