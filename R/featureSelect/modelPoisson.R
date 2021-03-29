# Poisson ======================================================================

#' @importFrom scran modelGeneVarByPoisson
.model.poisson <- function(sce, model) {
  .kwargs(
    scran::modelGeneVarByPoisson,
    sce,
    model
  )
}
