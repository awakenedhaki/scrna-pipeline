# Data processing ==============================================================

#' TODO: Documentation
.kwargs <- function(f, param, kwargs) {
  do.call(purrr::partial(f, param), kwargs)
}

#' @title Find the index of strings that match a given pattern.
#'
#' @param pattern A string pattern.
#' @param strings A vector of strings.
#' @param matcher A matching function (ex. \code{\link{grepl}})
#'
#' @return Vector of numerics corresponding to indexes where the pattern matched in the string vector.
.where <- function(pattern, strings, matcher) {
  return(which(matcher(pattern, strings)))
}

# Logging ======================================================================

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
console.logger <- function(msg) {
  futile.logger::flog.info(msg)
}

#' TODO: Add documentation
#'
#' @param msg String, output message
#'
#' @importFrom futile.logger flog.info
file.logger <- function(msg) {
  futile.logger::flog.info(msg, name = 'analysis')
}

# Saving =======================================================================

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
