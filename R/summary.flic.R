#' @exportS3Method summary flic
summary.flic <- function(x, ...)
{
  # x ... object of class flic
  cat("Firth's logistic regression with intercept correction\n\n")
  cat("Call:\n")
  print(x$call)
  out <- cbind(x$coefficients,x$var, x$ci.lower, x$ci.upper, x$probabilities)
  dimnames(out) <- list(names(x$coefficients), c("coef","se(coef)",paste(c("lower", "upper"),1 - x$alpha), "p"))
  print(out)
}