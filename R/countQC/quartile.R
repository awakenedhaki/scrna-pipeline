# Quartile Thresholds Protocol =================================================

#' @title Flag cells that are in the first quartile of given criteria.
#'
#' @details This method does not rely on user-defined boundaries for criteria values.
#'   Therefore, the user does *not* have to specify a range.
#'
#'   The quartile protocol is adaptive. It will infer low quality cells from the data itself.
#'   This could be beneficial since it removes human bias. However, the adaptive protocols
#'   may remove more cells that necessary/desired. It is best to cross-reference this protocol
#'   with the \code{legacy} implementation.
#'
#' @param criteria.data An S4Vector with quality control criteria of interest.
#'
#' @return Logical vector, flagging cells to be removed.
.quartile <- function(criteria.data) {
  flag.matrix <- sapply(names(criteria.data), FUN = function(name) {
    .quartile.flag.cells(criteria.data[[name]])
  })

  flags <- .reduce.flag.matrix(flag.matrix)
  return(flags)
}

#' @title Determine if cells are in the first quartile range.
#'
#' @param criteria.column Numeric vector for one QC criteria.
#'
#' @return Logical matrix, flagging cells to be removed.
.quartile.flag.cells <- function(criteria.column) {
  quartile <- .build.quartile(criteria.column)
  first.quartile.range <- quartile[['25%']]

  flags <- criteria.column < first.quartile.range
  return(flags)
}

#' @title Builds the quartile ranges for a given criteria.
#'
#' @param criteria.column Numeric vector for one QC criteria.
#'
#' @return Named vector corresponding the quartile ranges.
.build.quartile <- function(criteria.column) {
  quartile <- stats::quantile(criteria.column, probs = seq(0, 1, by = 0.25))
  return(quartile)
}
