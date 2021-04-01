#' @title Dispatcher function for quality control protocols.
#'
#' @description Perfoms user-specified quality control (QC) protocols on
#'   \code{\link{SingleCellExperiment}} objects.
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param criteria A list of a list containig details for each user-defined criteria.
#'
#'   The first element of a sublist must contain the gene symbol pattern to be searched.
#'   The second element is only required if `legacy` protocol is being used, and it must contain
#'   a vector of two numeric elements specifying the lower and upper bounds for the criteria value.
#'
#'   Common criteria include:
#'   * Mitochondrial gene percentage (ex. list(mito = list('^MT-', c(1, 60))))
#'   * Ribosomal gene percentage (ex. list(ribo = list('^RP[SL][[:digit:]], c(0, 60))))
#' @param protocol String specifying QC protocol.
#'   Current protocols include:
#'   * 'legacy'
#'   * 'quartile'
#'   * 'scater'
#' @param rm.genes Determine if rows corresponding to QC criteria should be removed.
#'
#' @details Implementated quality control methods are:
#' * Fixed tresholds: User-defined threshold for specified QC metrics.
#' * Quartile ranges: Removal of cells that fit within the first quartile range.
#' * \code{\link[scater]{quickPerCellQC}}: Relies on QC protocol implemented within \code{\link{scater}} .
#'
#' @importFrom scater addPerCellQC
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
countQC <- function(sce,
                    criteria,
                    protocol = 'legacy',
                    rm.genes = TRUE,
                    ...) {
	gene.symbol.patterns <- .get.gene.symbol.patterns(criteria)
	subsets <- .build.subsets(sce, gene.symbol.patterns)
	qc.sce <- scater::addPerCellQC(sce, subsets = subsets)

	criteria.data <- .get.criteria.data(qc.sce, names(criteria))
	if (protocol == 'legacy') {
		ranges <- .get.criteria.ranges(criteria)
		flagged.cells <- .legacy(criteria.data, ranges, ...)
	} else if (protocol == 'quartile') {
		flagged.cells <- .quartile(criteria.data, ...)
	} else if (protocol == 'scater') {
		flagged.cells <- .scater(sce, ...)
	}

	if (rm.genes) {
	  qc.sce <- qc.sce[!Reduce(f = `|`, subsets), ]
	}

  qc.sce <- qc.sce[, !flagged.cells]
  save.processed(qc.sce, paste0('qc-', protocol))
	return(qc.sce)
}

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

  save.diagnostic((summary(flag.matrix)), 'qc-legacy')

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

  save.diagnostic((summary(flag.matrix)), 'qc-quartile')

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

# scater::quickPerCellQC ===============================================

#' @title Quick quality control using \code{\link[scater]quickPerCellQC}.
#'
#' @param sce \code{\link{SingleCellExperiment}} object.
#' @param ... Additional parameters for \code{\link[scater]{quickPerCellQC}}.
#'
#' @importFrom scater quickPerCellQC
.scater <- function(sce, ...) {
  scater::quickPerCellQC(sce, ...)
}
