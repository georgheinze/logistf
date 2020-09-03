#' @exportS3Method extractAIC logistf
extractAIC.logistf<-function(fit, scale, k=2, ...){
  dev<- -2*diff(fit$loglik)
  AIC<- dev+k*fit$df
  edf<- fit$df
  return(c(edf,AIC))
}

#' @exportS3Method extractAIC flic
extractAIC.flic<-function(fit, scale, k=2, ...){
  dev<- -2*diff(fit$loglik)
  AIC<- dev+k*fit$df
  edf<- fit$df
  return(c(edf,AIC))
}
#' @exportS3Method extractAIC flac
extractAIC.flac<-function(fit, scale, k=2, ...){
  dev<- -2*diff(fit$loglik)
  AIC<- dev+k*fit$df
  edf<- fit$df
  return(c(edf,AIC))
}

#' @exportS3Method nobs logistf
nobs.logistf<-function(object, ...){
  return(object$n)
}

#' @exportS3Method nobs flic
nobs.flic<-function(object, ...){
  return(object$n)
}
#' @exportS3Method nobs flac
nobs.flac<-function(object, ...){
  return(object$n)
}

# @exportS3Method terms logistf
#terms.logistf<-function(x, ...){
#  return(terms(formula(x)))
#}