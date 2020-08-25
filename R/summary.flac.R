#' @exportS3Method summary flac
summary.flac <- function(object, ...)
{
  # x ... object of class flac
  cat("Firth's logistic regression with added covariate\n\n")
  print(object$call)
  cat("\nModel fitted by", object$method)
  cat("Coefficients:\n")
  out <- cbind(object$coefficients, object$var, object$ci.lower,object$ci.upper, qchisq(1 - object$prob, 1), object$prob, ifelse(object$method.ci=="Wald", 1, 2))
  dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method"))
  print(out)
  
  cat("\n Method: 1-Wald, 2-Profile penalized log-likelihood\n")
  
  LL <- 2 * diff(object$loglik)
  cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 1 -pchisq(LL, object$df), ", n=",object$n, sep = "")
  if(object$terms[1]!="(Intercept)"){
    wald.z <- tryCatch({
      t(coef(object)) %*% solve(object$var) %*% coef(object)
    }, 
    error=function(cond){
      message("\n Variance-Covariance matrix is singular \n")
      return(NA)
    }
    )
  }
  else{
    wald.z <- tryCatch({
      t(coef(object)[2:(object$df+1)]) %*%
        solve(object$var[2:(object$df+1),2:(object$df+1)]) %*%
        coef(object)[2:(object$df+1)]
    }, 
    error=function(cond){
      message("\n Variance-Covariance matrix is singular \n")
      return(NA)
    }
    )
  }
  cat("\nWald test =", wald.z, "on", object$df, "df, p =", 1 - pchisq(wald.z, object$df))
  #cat("\n\nCovariance-Matrix:\n")
  #print(object$var)
  invisible(object)
}