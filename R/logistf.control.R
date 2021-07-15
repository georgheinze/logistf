#' Control Parameters for \code{logistf}
#'
#' Sets parameters for iterations in Firth's penalized-likelihood logistic regression.
#' 
#' \code{logistf.control()} is used by \code{logistf} and \code{logistftest} to set control parameters to default values. 
#' Different values can be specified, e. g., by \code{logistf(..., control= logistf.control(maxstep=1))}.
#'
#' @param maxit The maximum number of iterations
#' @param maxhs The maximum number of step-halvings in one iteration. The increment of the 
#' beta vector within one iteration is divided by 2 if the new beta leads to a decrease 
#' in log likelihood.
#' @param maxstep Specifies the maximum step size in the beta vector within one iteration.
#' @param lconv Specifies the convergence criterion for the log likelihood.
#' @param gconv Specifies the convergence criterion for the first derivative of the log likelihood (the score vector).
#' @param xconv Specifies the convergence criterion for the parameter estimates.
#' @param collapse If \code{TRUE}, evaluates all unique combinations of x and y and collapses data set.
#' @param fit  Fitting method used. One of Newton-Raphson: "NR" or Iteratively reweighted least squares: "IRLS"
#'
#' @return
#'    \item{maxit}{The maximum number of iterations}
#'    \item{maxhs}{The maximum number of step-halvings in one iteration. The increment of the 
#' beta vector within one iteration is divided by 2 if the new beta leads to a decrease 
#' in log likelihood.}
#'    \item{maxstep}{Specifies the maximum step size in the beta vector within one iteration.}
#'    \item{lconv}{Specifies the convergence criterion for the log likelihood.}
#'    \item{gconv}{Specifies the convergence criterion for the first derivative of the log likelihood (the score vector).}
#'    \item{xconv}{Specifies the convergence criterion for the parameter estimates.}
#'    \item{collapse}{If \code{TRUE}, evaluates all unique combinations of x and y and collapses data set.}
#'    \item{fit}{Fitting method used. One of Newton-Raphson: "NR" or Iteratively reweighted least squares: "IRLS"}
#'    \item{call}{The function call.}
#' @export
#' 
#' @encoding UTF-8
#' @examples
#' data(sexagg)
#' fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, control=logistf.control(maxstep=1))
#' summary(fit2)
#' 
logistf.control <-
function(maxit=25, maxhs=0, maxstep=5, lconv=0.00001, gconv=0.00001, xconv=0.00001, collapse=TRUE, fit = "NR"){
  res<-list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, gconv=gconv, xconv=xconv, collapse=collapse, fit = fit, call=match.call())
  attr(res, "class")<-"logistf.control"
  return(res)
}

