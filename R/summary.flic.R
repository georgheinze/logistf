#' @exportS3Method summary flic
summary.flic <- function(object, ...)
{ 
  cat("Firth's logistic regression with intercept correction\n\n")
  print(object$call)
  cat("\nModel fitted by", object$method)
  cat("Coefficients:\n")
  out <- cbind(object$coefficients,object$var, object$ci.lower,object$ci.upper, qchisq(1 - object$prob, 1), object$prob, ifelse(object$method.ci=="Wald", 1, 2))
  dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method"))
  print(out)

  cat("\n Method: 1-Wald, 2-Profile penalized log-likelihood\n")
  
  LL <- 2 * diff(object$loglik)
  cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 1 -pchisq(LL, object$df), ", n=",object$n, sep = "")

  invisible(object) 
}