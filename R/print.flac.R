#' @exportS3Method print flac
print.flac <-
  function(x, ...)
  {
    # x ... object of class flac
    cat("Firth's logistic regression with added covariate\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\n\nCoefficients:\n")
    out <- x$coefficients
    print(out)
    LL <- -2 * (x$loglik['null']-x$loglik['full'])
    cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", 1 -
          pchisq(LL, x$df), ", n=",
        x$n, "\n\n", sep = "")
    invisible(x)
  }