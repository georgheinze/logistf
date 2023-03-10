#' Recover data method for logistf objects
#' 
#' @description This function provides a \code{recover_data} method for \code{logistf} objects. This is required for \code{emmeans} support.
#' 
#' @param object a \code{logistf} object.
#' @param frame the model frame.
#' @param ... additional arguments.
#'
#' @export recover_data.logistf

recover_data.logistf <- function(object, frame = object$model, ...){
    fcall = object$call
    emmeans::recover_data(fcall, delete.response(terms(object)),
                          object$na.action, frame = frame, ...)
}


