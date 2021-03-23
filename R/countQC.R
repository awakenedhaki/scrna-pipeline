# Dependencies =================================================================

# Constants ====================================================================
MITOCHONDRIAL.GENE <- '^MT-'
RIBOSOMAL.GENE <- '^RP[SL][[:digit:]]'

# Helpers ======================================================================

#' @title Determines if numeric is outside a given range, exclusively.
#'
#' @details Given a numeric vector, determine which values are outside a range
#'   specified by the 'range' vector. The min/max comparison are combined using a
#'   logical operator (ex. `|` or `&`).
#'
#' @param numerics Numeric vector
#' @param range Numeric vector representing the min/max range
#' @param operator Comparison (logical) operator function
#'
#' @return Logical vector
.outside.range <- function(numerics, range, operator) {
  operator(numerics < min(range), numerics > max(range))
}

#' @title Generate quartiles for given numeric vector.
#'
#' @param data Numeric vector
#'
#' @return Named vector of quartile ranges
.quartile <- function(data) {
  stats::quantile(x = data, probs = base::seq(0, 1, 0.25))
}

# Functions ====================================================================

#' @title Identify low QC cells with frequent QC metrics.
#'
#' @param sce \code{\link{SingleCellExperiment}}
#' @param ... TODO
#'
#' @importFrom scater quickPerCellQC
#'
#' @return DataFrame
.quickPerCellQCMetrics <- function(sce, ...) {
  qc <- scater::quickPerCellQC(sce, ...)
}

#' @title Filters cells based on fixed thresholds.
#'
#' @details Removes cells (columns) that do not meet the paramaterized criteria. These cells
#'   can be referred to as low-quality libraries that may arise due to cell damage, or
#'   failure during library preparation (ex. drop-out events)
#'
#'   The median absolute deviation is a metric for variability, As the name suggests,
#'   it is the deviation of the absolute deviation (sum of differences between point X
#'   and the median) of a vector of values. This method should not be influenced by
#'   extreme outlier values, and can be used instead of z-scores. Assumes a normal
#'   distribution, without consideration of outlier events.
#'
#'   Library size filtering performed by \code{\link[scater]{isOutlier}}, looking at the 'lower'
#'   end of the distribution, in the log scale. Number of detected genes is also
#'   filtered under the same criteria.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param nmads Median absolute deviation
#' @param mito.range Min/max vector for mitochondrial percentage
#' @param ribo.range Min/max vector for ribosomal percentage
#'
#' @importFrom scater isOutlier
#'
#' @return SingleCellExperiment object
.fixedFilter <- function(sce,
                         nmads = 3,
                         mito.range = c(1, 60),
                         ribo.range = c(0, 60)) {
  # Flagging cells (columns) with high or low mitochondrial reads
  out.mitochondrial <- .outside.range(sce$subsets_mitochondrial_percent,
                                      range = mito.range,
                                      operator = `|`)


  # Flagging cells (columns) with high or low ribosomal reads
  out.ribosomal <- .outside.range(sce$subsets_ribosomal_percent,
                                  range = ribo.range,
                                  operator = `|`)

  # Flagging cells (columns) by library size's deviation from the median
  out.libsize <- scater::isOutlier(sce$sum, nmads = nmads, type = 'lower', log = TRUE)
  # Flagging cells (columns) by number of detected genes deviation from the median
  out.gene <- scater::isOutlier(sce$detected, nmads = nmads, type = 'lower', log = TRUE)

  # Cells (columns) to be removed
  rm.cells <- (out.mitochondrial | out.ribosomal | out.libsize | out.gene)

  save.diagnostic(data.frame(Library.Size = sum(out.libsize),
                             Detected.Genes = sum(out.gene),
                             Mitochondrial = sum(out.mitochondrial),
                             Ribosomal = sum(out.ribosomal),
                             Total = sum(rm.cells)),
                  'countQC-fixed.rds')
  logger('Number of cells filtered per criteria were saved in "diagnostic" directory.')

  # Filtering low-quality cells
  return(sce[ , !rm.cells])
}

#' @title Filters cells based on quartiles.
#'
#' @description Generate quartiles for mitochondrial percent, library size, and number of
#'   genes detected. Cells in the first quartile of any feature is removed.
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#'
#' @return \code{\link{SingleCellExperiment}} object
.adaptiveFilter <- function(sce) {
  # Flagging cells (columns) in the first quartile of mitochondrial content
  out.mitochondrial <- sce$subsets_mitochondrial_percent < .quartile(sce$subsets_mitochondrial_percent)[['25%']]

  # Flagging cells (columns) in the first quartile of library size
  out.libsize <- sce$sum < .quartile(sce$sum)[['25%']]

  # Flagging cells (columns) in the first quartile of genes detected
  out.gene <- sce$detected < .quartile(sce$detected)[['25%']]

  # Cells (columns to be removed)
  rm.cells <- (out.mitochondrial | out.libsize | out.gene)

  logger('Saving number of cells flagged per criteria in "diagnostic" directory.')
  save.diagnostic(data.frame(Library.Size = sum(out.libsize),
                             Detected.Genes = sum(out.gene),
                             Mitochondrial = sum(out.mitochondrial),
                             Total = sum(rm.cells)),
                  'countQC-quartile.rds')

  return(sce[ , !rm.cells])
}

#' @title Main function for read count quality control (dispatcher).
#'
#' TODO: Look into inheritParams tag
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#' @param method Character specifying which quality control method to use
#' @param rm.genes Remove mitochondrial and ribosomal genes
#' @param ... Additional parameters for method function
#'
#' @importFrom scater addPerCellQC
#' @importFrom SingleCellExperiment colData
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
countQC <- function(sce,
                    method,
                    rm.genes = TRUE,
                    ...) {
  logger('Initializing quality control...')

  # Flag mitochondrial genes
  is.mito <- grepl(MITOCHONDRIAL.GENE, rownames(sce))
  # Flag ribosomal genes
  is.ribo <- grepl(RIBOSOMAL.GENE, rownames(sce))

  # Adding QC metrics
  # . Sum: Sum of counts per cells
  # . Detected: Number of genes detected per cell
  # . Mitochondrial: Proportion of mitochondrial genes
  # . Ribosomal: Proportion of ribosomal genes
  sce <- scater::addPerCellQC(sce, subsets = list(mitochondrial = is.mito, ribosomal = is.ribo))

  # Dispatcher
  if (method == 'fixed') {
    logger('Initializing cell filtering by fixed thresholds.')
    sce <- .fixedFilter(sce, ...)
  } else if (method == 'quartile') {
    logger('Initializing cell filtering by adaptive thresholds.')
    sce <- .adaptiveFilter(sce)
  } else if (method == 'scater') {
    logger('Initializing cell filtering by quickPerCellQCMetrics.')
    qc <- .quickPerCellQCMetrics(SingleCellExperiment::colData(sce), ...)
    sce <- sce[ , !qc$discard]
  }

  # Remove mitochondrial and ribosomal genes
  if (rm.genes) {
    logger('Removing mitochondrial and ribosomal genes from dataset.')
    sce <- sce[-c(which(is.mito), which(is.ribo)), ]
  }

  logger('Saving quality controlled data in "processed" directory.')
  save.processed(sce, 'qc.rds')

  logger('Quality control has ended.')
  return(sce)
}
