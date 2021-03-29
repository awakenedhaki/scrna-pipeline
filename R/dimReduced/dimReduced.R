dimReduced <- function(sce, protocol, seed = 100, params) {
  set.seed(seed)
  if (protocol == 'vanilla') {
    sce <- .reduced.vanilla(sce, params)
  }

  return(sce)
}
