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
#' * \code{\link[scater]{quikckPerCellQCMetrics}}: Relies on QC protocol implemented within \code{\link{scater}} .
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

	return(qc.sce[ , !flagged.cells])
}
