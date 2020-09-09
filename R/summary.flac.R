#' @exportS3Method summary flac
summary.flac <- function(object, ...)
{
  cat("Firth's logistic regression with added covariate\n\n")
  print(object$call)
  cat("\nModel fitted by", object$method)
  
  #check if flac was called with logistf object:
  if(!is.null(object$call$lfobject)){
    lfobject <- eval(object$call$lfobject, parent.frame())
    if(!is.null(lfobject$call$terms.fit)){
      terms.fit <- lfobject$call$terms.fit
    }
    else terms.fit <- NULL
  }
  else terms.fit <- NULL
  
  
  #consider for wald only covariance matrix with columns corresponding to variables in terms.fit
  call <- object$call
  if(!is.null(call$terms.fit) | !is.null(terms.fit)){
    terms.fit <- c(call$terms.fit,terms.fit)[!is.null(c(call$terms.fit,terms.fit))]
    terms.fit <- eval(terms.fit, parent.frame())
    loc <- match(terms.fit, object$terms)
    var.red <- object$var[loc,loc]
    coefs <- coef(object)[loc]
    chi2 <- vector(length=length(object$terms))
    chi2[loc] <- qchisq(1 - object$prob[loc], 1)
    chi2[-loc] <- 0
  }
  else {
    var.red <- object$var
    coefs <- coef(object)
    chi2 <- qchisq(1 - object$prob, 1)
  }
  
  cat("\nCoefficients:\n")
  out <- cbind(object$coefficients, diag(object$var)^0.5, object$ci.lower,object$ci.upper, qchisq(1 - object$prob, 1), object$prob, ifelse(object$method.ci=="Wald", 1, ifelse(object$method.ci=="-", 3, 2)))
  dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method"))
  print(out)
  
  cat("\nMethod: 1-Wald, 2-Profile penalized log-likelihood, 3-None\n")
  
  LL <- 2 * diff(object$loglik)
  cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 1 -pchisq(LL, object$df), ", n=",object$n, sep = "")
  if(object$terms[1]!="(Intercept)"){
    wald.z <- tryCatch({
      t(coefs) %*% solve(var.red) %*% coefs
    }, 
    error=function(cond){
      message("\n Variance-Covariance matrix is singular \n")
      return(NA)
    }
    )
  }
  else{
    wald.z <- tryCatch({
      t(coefs[2:nrow(var.red)]) %*%
        solve(var.red[2:nrow(var.red),2:nrow(var.red)]) %*%
        coefs[2:nrow(var.red)]
    }, 
    error=function(cond){
      message("\n Variance-Covariance matrix is singular \n")
      return(NA)
    }
    )
  }
  cat("\nWald test =", wald.z, "on", object$df, "df, p =", 1 - pchisq(wald.z, object$df))

  invisible(object)
}