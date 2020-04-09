#' @exportS3Method print flac
print.flac <-
  function(x, ...)
  {
    # x ... object of class flac
    cat("Firth's logistic regression with added covariate\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\n\nCoefficients:\n")
    out <- cbind(x$coefficients)
    print(out)
  }