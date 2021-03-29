# Vanilla ======================================================================

#' @importFrom scran modelGeneVar
.model.vanilla <- function(sce, model) {
  .kwargs(
    scran::modelGeneVar,
    sce,
    model
  )
}
