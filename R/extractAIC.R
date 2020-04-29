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

#' @exportS3Method terms logistf
terms.logistf<-function(x, ...){
  object<-x
  options(warn=-1)
  fakeglm<-try(glm(formula=object$formula,data=object$data, family=binomial,maxit=1))
  options(warn=0)
  return(terms(fakeglm))
}
#' @exportS3Method terms flic
terms.flic<-function(x, ...){
  object<-x
  options(warn=-1)
  fakeglm<-try(glm(formula=object$formula,data=object$data, family=binomial,maxit=1))
  options(warn=0)
  return(terms(fakeglm))
}
#' @exportS3Method terms flic
terms.flac<-function(x, ...){
  object<-x
  options(warn=-1)
  fakeglm<-try(glm(formula=object$formula,data=object$data, family=binomial,maxit=1))
  options(warn=0)
  return(terms(fakeglm))
}