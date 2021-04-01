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
