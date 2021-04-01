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

#' TODO: documentation
.read.variance.model <- function(identifier) {
  # TODO: Handle multiple model files.
  diagnostic.path <- here(getwd(), 'data', 'diagnostic')

  diagnostic.files <- list.files(diagnostic.path)
  pattern <- paste0(paste('featureSelect-.*', identifier, sep = '-'), '.rds')
  variance.model.idx <- .where('featureSelect-.*\\.rds', diagnostic.files, matcher = grepl)
  if (length(variance.model.idx) == 0) {
    stop('Fitted variance model data must be within diagnostic directory.')
  }

  model.path <- here(diagnostic.path, diagnostic.files[variance.model.idx])
  model <- readRDS(model.path)
  return(model)
}

#' TODO: Add documentation
.blackbox.namespace <- function(name) {
  base::get(name, envir = as.environment('package:blackbox'))
}
