# Dependencies =================================================================

# Constants ====================================================================

# Helpers ======================================================================

# Functions ====================================================================

#' @title Find markers within a dataset
#'
#' @param sce \code{\link{SingleCellExperiment}} object
#'
#' @importFrom scran findMarkers
.findMarkers <- function(sce) {
  markers.sce <- scran::findMarkers(sce)

}

#' @title Find markers
diffExp <- function() {}
