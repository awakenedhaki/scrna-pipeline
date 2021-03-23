# Dependencies =================================================================

# Constants ====================================================================

# Helpers ======================================================================

# Functions ====================================================================

#' @title Prepare data for processing
#'
#' @param path10x Path to CellRanger output directory
#' @param experiment.name String, experiment identifier
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom SummarizedExperiment rowData
#' @importFrom scater uniquifyFeatureNames
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#'
#' @return \code{\link{SingleCellExperiment}} object
#' @export
prep.data <- function(path10x, experiment.name) {
  sce <- DropletUtils::read10xCounts(samples = path10x)

  # Rename column names in colData to Barcode IDs
  colnames(sce) <- sce$Barcode

  # Rename row names in colData to ensembl ID-gene symbol pairs
  rownames(sce) <- scater::uniquifyFeatureNames(ID = rownames(sce),
											                          names = SummarizedExperiment::rowData(sce)$Symbol)

  # Create column in colData with experiment name for downstream group_by
  sce$experment.name <- experiment.name

  # Replace NA with 'NA' in rowData column names
  colnames(SummarizedExperiment::rowData(sce)) <- sce %>%
    SummarizedExperiment::rowData() %>%
    names() %>%
    tidyr::replace_na(replace = 'NA')

  save.processed(sce, 'prep.rds')

  return(sce)
}
