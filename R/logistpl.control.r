#' Control Parameters for logistf Profile Likelihood Confidence Interval Estimation
#'
#' Sets parameters for modified Newton-Raphson iteration for finding 
#' profile likelihood confidence intervals in Firth's penalized likelihood logistic regression
#' 
#' \code{logistpl.control()} is used by \code{logistf} to set control parameters to default values 
#' when computing profile likelihood confidence intervals. 
#' Different values can be specified, e. g., by \code{logistf(..., control= logistf.control(maxstep=1))}.
#'
#' @param maxit The maximum number of iterations
#' @param maxhs The maximum number of step-halvings in one iteration. The increment of the 
#' beta vector within one iteration is divided by 2 if the new beta leads to a decrease 
#' in log likelihood.
#' @param maxstep Specifies the maximum step size in the beta vector within one iteration. Set to -1 for infinite stepsize.
#' @param lconv Specifies the convergence criterion for the log likelihood.
#' @param xconv Specifies the convergence criterion for the parameter estimates.
#' @param ortho Requests orthogonalization of variable for which confidence intervals are computed with respect to other covariates
#' @param pr Request rotation of the matrix spanned by the covariates
#'
#' @return
#'    \item{maxit}{The maximum number of iterations}
#'    \item{maxhs}{The maximum number of step-halvings in one iteration. The increment of the 
#' beta vector within one iteration is divided by 2 if the new beta leads to a decrease 
#' in log likelihood.}
#'    \item{maxstep}{Specifies the maximum step size in the beta vector within one iteration.}
#'    \item{lconv}{Specifies the convergence criterion for the log likelihood.}
#'    \item{xconv}{Specifies the convergence criterion for the parameter estimates.}
#'    \item{ortho}{specifies if orthogonalization is requested.}
#'    \item{pr}{specifies if rotation is requested}
#' @export
#' 
#' @author Georg Heinze
#'
#' @examples
#' data(sexagg)
#' fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, 
#'     plcontrol=logistpl.control(maxstep=1))
#' summary(fit2)
#'
logistpl.control<-function(maxit=100, maxhs=0, maxstep=5, lconv=0.00001, xconv=0.00001, ortho=FALSE, pr=FALSE){
  list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, xconv=xconv, ortho=ortho, pr=pr)
}
