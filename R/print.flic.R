#' @exportS3Method print flic
print.flic <-
  function(x, ...)
  {
    # x ... object of class flic
    cat("Firth's logistic regression with intercept correction\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\n\nCoefficients:\n")
    out <- x$coefficients
    print(out)
    LL <- 2 * diff(x$loglik)
    cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", 1 -
          pchisq(LL, x$df), ", n=",
        x$n, "\n\n", sep = "")
    invisible(x)
  }