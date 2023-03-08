#' Recover data method for logistf objects
#' 
#' @description This function provides a \code{recover_data} method for \code{logistf} objects. This is required for \code{emmeans} support.
#' 
#' @export

recover_data.logistf <- function(object, frame = object$model, ...){
    fcall = object$call
    recover_data(fcall, delete.response(terms(object)), object$na.action,
                 frame = frame, ...)
}


