#' @exportS3Method print logistftest
print.logistftest <-
function(x, ...){
# x ... object of class logistftest
 print(x$call)
 cat("Model fitted by", x$method, "\n\nFactors fixed as follows:\n")
 print(x$testcov)
 LL <- -2 * (x$loglik['null']-x$loglik['full'])
 out <- c(x$loglik['null'], x$loglik['full'], LL/2)
 names(out) <- c("Restricted model", "Full model", "difference")
 cat("\nLikelihoods:\n")
 print(out)
 cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", x$prob, "\n", sep = "")
 invisible(x)
}

