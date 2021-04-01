#' TODO: Add documentation
.get.identifier <- function(sce) {
  identifier <- unique(sce$identifier)
  if (length(identifier) > 1) {
    identifier <- paste('integrated',
                        paste(identifier, collapse = '-'),
                        sep = '-')
  }

  return(identifier)
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.processed <- function(obj, name) {
  file.name <- paste0(name, '.rds')
  saveRDS(obj, file = here(getwd(), 'data', 'processed', file.name))

  console.logger(paste('Saved', name, 'within data/processed directory as RDS.'))
  file.logger(paste('Saved', name, 'within data/processed directory as RDS.'))
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.diagnostic <- function(obj, name) {
  file.name <- paste0(name, '.rds')
  saveRDS(obj, file = here(getwd(), 'data', 'diagnostic', file.name))

  console.logger(paste('Saved', name, 'within data/diagnostic directory as RDS.'))
  file.logger(paste('Saved', name, 'within data/diagnostic directory as RDS.'))
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.output <- function(obj, name) {
  saveRDS(obj, file = here(getwd(), 'output', name))

  console.logger(paste('Saved', name, 'within data/output directory as RDS.'))
  file.logger(paste('Saved', name, 'within data/output directory as RDS.'))
}
