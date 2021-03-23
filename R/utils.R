# Dependencies =================================================================

# Helpers ======================================================================

#' TODO: Add documentation
#'
#' @param name String, name for logger
#'
#' @importFrom here here
#' @importFrom futile.logger flog.logger appender.file
generate.log <- function(name) {
  filename <- paste0(name, '.log')

  if (!file.exists(here(filename))) {
    file.create(filename)
  }
  futile.logger::flog.logger(name, appender = futile.logger::appender.file(filename))
}

#' TODO: Add documentation
#'
#' @param msg String, output message
#'
#' @importFrom futile.logger flog.info
logger <- function(msg) {
  futile.logger::flog.info(msg)
  futile.logger::flog.info(msg, name = 'analysis')
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.processed <- function(obj, name) {
  saveRDS(obj, file = here(getwd(), 'data', 'processed', name))
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.diagnostic <- function(obj, name) {
  saveRDS(obj, file = here(getwd(), 'data', 'diagnostic', name))
}

#' TODO: Add documentation
#'
#' @param obj \code{\link{SingleCellExperiment}} object
#' @param name String, name of file
#'
#' @importFrom here here
save.output <- function(obj, name) {
  saveRDS(obj, file = here(getwd(), 'output', name))
}

