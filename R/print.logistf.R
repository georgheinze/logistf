#' @exportS3Method print logistf
print.logistf <-
function(x, ...)
{
# x ... object of class logistf
 print(x$call)
 cat("Model fitted by", x$method)
 cat("\nConfidence intervals and p-values by", paste(unique(x$method.ci), sep="/"), "\n\n")
 cat("Coefficients:\n")
 out <- x$coefficients
 print(out)
 LL <- 2 * diff(x$loglik)
 cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", 1 -
pchisq(LL, x$df), ", n=",
  x$n, "\n\n", sep = "")
 invisible(x)
}

