#' @exportS3Method predict logistf
predict.logistf <- function (lfobject, newdata, type = c("link", "response"), flic=FALSE) 
{
  type <- match.arg(type)
  if (missing(newdata)) {#no data - return linear.predictors or response according to type
    if (flic) {
      pred <- switch(type, link = lfobject$flic.linear.predictors, response = lfobject$flic.predict)
    }
    else {
      pred <- switch(type, link = lfobject$linear.predictors, response = lfobject$predict)
    }
  }
  else {
    newlin <- c(lfobject$coefficients[-1] %*% t(newdata))
    if (flic) {
    pred <- switch(type, link = lfobject$flic.coefficients[1]+newlin, 
                   response = 1/(1+exp(-(lfobject$flic.coefficients[1]+newlin))) )
    }
    else {
      pred <- switch(type, link = lfobject$coefficients[1]+newlin , 
                     response = 1/(1+exp(-(lfobject$coefficients[1]+newlin))))
    }
  }
  pred
}