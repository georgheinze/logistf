#' @exportS3Method summary logistf
summary.logistf <-function(object,...){
  # object ... object of class logistf
   print(object$call)
   cat("\nModel fitted by", object$method)
   cat("\nCoefficients:\n")
   
   #consider for wald only covariance matrix with columns corresponding to variables in terms.fit
   call <- object$call
   if(!is.null(call$terms.fit)){
      loc <- match(call$terms.fit, object$terms)
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

   out <- cbind(object$coefficients, diag(object$var)^0.5, object$ci.lower,object$ci.upper, chi2, object$prob, ifelse(object$method.ci=="Wald", 1, ifelse(object$method.ci=="-", 3, 2)))
   dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method"))
  # if(object$method.ci != "Wald")
  #  dimnames(out)[[2]][5] <- "Chisq"
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
   #cat("\n\nCovariance-Matrix:\n")
   #print(object$var)
   invisible(object) 
}