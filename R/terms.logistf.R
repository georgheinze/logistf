#' @exportS3Method terms logistf
terms.logistf <- function(x, ...) {
  attr(x$model, which = 'terms', exact = TRUE)
}

#' @exportS3Method terms logistf
terms.flic <- function(x, ...) {
  attr(x$model, which = 'terms', exact = TRUE)
}

#' @exportS3Method terms logistf
terms.flac <- function(x, ...) {
  attr(x$model, which = 'terms', exact = TRUE)
}