# Dependencies =================================================================

# Constants ====================================================================

# Helpers ======================================================================

# Functions ====================================================================
#' @title Prepare data for processing.
#'
#' @details Reading 10x sample directory containing `barcodes.tsv`, `genes.tsv` and `matrix.mtx`
#'   and encapsulating the data as a \code{\link{SingleCellExperiment}} object.
#'
#'   The data is prepared for downstream processing by:
#'   * Renaming column names to the `Barcode` IDs.
#'   * Renaming row names in to ENSEMBL ID-gene symbol pairs.
#'   * Adding a experiment identifier column for ease of dataset integration.
#'   * Replacing \code{NA} values with character representation \code{"NA"}.
#'
#'   \code{\link{SummarizedExperiment}} is imported since \code{\link[SummarizedExperiment]{rowData}}
#'   is not within the \code{\link{SingleCellExperiment}} namespace.
#'
#' @param path10x String, path to 10X sample data.
#' @param identifier String, experiment identifier
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom scater uniquifyFeatureNames
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
prep.data <- function(path10x, identifier) {
  # TODO:
  # . Include @seealso for functions used in implementation?
  sce <- DropletUtils::read10xCounts(samples = path10x)

  colnames(sce) <- sce$Barcode
  rownames(sce) <- scater::uniquifyFeatureNames(ID = rownames(sce),
											                          names = SummarizedExperiment::rowData(sce)$Symbol)

  sce$experiment.name <- identifier

  colnames(SummarizedExperiment::rowData(sce)) <- sce %>%
    SummarizedExperiment::rowData() %>%
    names() %>%
    tidyr::replace_na(replace = 'NA')

  # Save data in ./data/processed/ directory
  save.processed(sce, 'prep', identifier)
  return(sce)
}
