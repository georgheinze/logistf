#' @exportS3Method terms logistf
terms.logistf <- function(x) {
  attr(x$model, which = 'terms', exact = TRUE)
}