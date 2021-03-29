# Vanilla ======================================================================

#' @importFrom scater librarySizeFactors logNormCounts
.norm.vanilla <- function(sce, ...) {
  norm.sce <- scater::librarySizeFactors(sce, ...) %>%
    scater::logNormCounts()

  return(norm.sce)
}
