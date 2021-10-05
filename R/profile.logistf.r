#' Compute Profile Penalized Likelihood
#' 
#' Evaluates the profile penalized likelihood of a variable based on a logistf model fit
#'
#' @param fitted An object fitted by \code{logistf}
#' @param which A righthand formula to specify the variable for which the profile should be evaluated, e.g., which=~X).
#' @param variable Alternatively to which, a variable name can be given, e.g., variable="X"
#' @param steps Number of steps in evaluating the profile likelihood
#' @param pitch Alternatively to steps, one may specify the step width in multiples of standard errors
#' @param limits Lower and upper limits of parameter values at which profile likelihood is to be evaluated
#' @param alpha The significance level (1-\eqn{\alpha} the confidence level, 0.05 as default).
#' @param firth Use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) 
#' or the standard maximum likelihood method (\code{firth=FALSE}) for the logistic regression.
#' @param legends legends to be included in the optional plot
#' @param control Controls Newton-Raphson iteration. Default is \code{control= logistf.control(maxstep, 
#' maxit, maxhs, lconv, gconv, xconv)}
#' @param plcontrol Controls Newton-Raphson iteration for the estimation of the profile likelihood 
#' confidence intervals. Default is \code{plcontrol= logistpl.control(maxstep, maxit, maxhs, lconv, xconv, ortho, pr)}
#' @param ... Further arguments to be passed.
#'
#' @return An object of class \code{logistf.profile} with the following items:
#'   \item{beta}{Parameter values at which likelihood was evaluated}
#'   \item{stdbeta}{Parameter values divided by standard error}
#'   \item{profile}{profile likelihood, standardized to 0 at maximum of likelihood. The values in 
#'   profile are given as minus \eqn{\chi^2}}
#'   \item{loglik}{Unstandardized profile likelihood}
#'   \item{signed.root}{signed root (z) of \eqn{\chi^2} values (negative for values below the maximum likelihood
#'    estimate, positive for values above the maximum likelihood estimate)}
#'    \item{cdf}{profile likelihood expressed as cumulative distribution function, obtained as
#'    \eqn{\Phi(z)}, where \eqn{\Phi} denotes the standard normal distribution function.}
#' 
#' @export
#' @references Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining 
#' profile likelihood information from logistic regressions. Statistics in Medicine, to appear.
#' 
#' @rdname profile.logistf
#' @encoding UTF-8
#'
#' @examples
#' data(sex2)
#' fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
#' plot(profile(fit,variable="dia"))
#' plot(profile(fit,variable="dia"), "cdf")
#' plot(profile(fit,variable="dia"), "density")
#' 
#' @importFrom utils data
#' 
#' @method profile logistf
#' @exportS3Method profile logistf
profile.logistf <-
function(fitted,  which, variable, steps=100, pitch = 0.05, limits,
                    alpha = 0.05,  firth = TRUE,
                    legends = TRUE,  control, plcontrol, plot=FALSE, ...){

  # by MP, 06.02.01
  # adapted and renamed by GH, 10.03.11 (incredible! 10 years!)
  # again adapted and renamed by GH, 13.05.13
  # which  ... righthand formula des zu plottenden Term (z.B. ~B oder ~A:D)
  # pitch  ... distances between points in std's
  # limits ... vector of MIN & MAX in std's, default=extremes of both CI's
  #            +- 0.5 std. of beta
  #
  
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("fitted","which"), names(mf), 0L)
  mf <- mf[c(1, m)]
  fitted <- eval(mf$fitted, parent.frame())
  variables <- fitted$terms[-1]
  
  formula<-fitted$formula

  # Next line added by Harry Southworth, 22/10/02.
  if (missing(which) & missing(variable)) stop("You must specify a variable: either by which (a one-sided formula) or by variable.")
  if (missing(control)) control<-fitted$control
  if (missing(plcontrol)) plcontrol<-logistpl.control()

  modcontrol <- fitted$modcontrol 
   
  call <- match.call()

  x <- model.matrix(fitted$formula, model.frame(fitted))
  n <- nrow(x)
  
  cov.name <- labels(x)[[2]]

  offset <- model.offset(mf)   
  weight <- model.weights(mf)
  if(is.null(offset)) {
    offset <- rep(0,n)
  }
  else {
    offset<-as.vector(offset)
  }
  if(is.null(weight)) {
    weight<-rep(1,n)
  }
      
  cov.name <- labels(x)[[2]]
  k <- ncol(x)
  
  
  
  if(!(identical(modcontrol$terms.fit, 1:k) | is.null(modcontrol$terms.fit))){
    stop("Please call profile on a logistf-object with all terms fitted.")
  }
  
  if(!missing(which)) {
    cov.name2 <- colnames(get_all_vars(which, data = model.frame(fitted))) 
  } else {
    cov.name2 <- variable
  }
  pos <- match(cov.name2, cov.name) 
  if(is.na(pos)) {
    stop(paste(variable,"is not a model term."))
  }
  
  
  std.pos <- diag(fitted$var)[pos]^0.5
   
  coefs <- fitted$coefficients 
  covs <- fitted$var 

  LL.0 <- fit$loglik['full'] - qchisq(1 - alpha, 1)/2
  if(missing(limits)) {
    lower.fit <- logistpl(x, y, init=fitted$coefficients, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=-1, i=pos, plcontrol=plcontrol, modcontrol = modcontrol)
    upper.fit <- logistpl(x, y, init=fit$coefficients, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=1, i=pos, plcontrol=plcontrol, modcontrol = modcontrol)
    limits <- c(lower.fit$beta, upper.fit$beta)
  }
  
  limits <- (limits - coefs[pos])/std.pos
  limits <- c(min(qnorm(alpha/2), limits[1]) - 0.5, max(qnorm(1 - alpha/2), limits[2]) + 0.5)
  limits <- c(floor(limits[1]/pitch) * pitch, ceiling(limits[2]/pitch) * pitch)
  
  knots <- seq(limits[1], limits[2], diff(limits)/steps)
  nn <- length(knots)
  res <- matrix(NA, nn, 3) 
  dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
  res[,1] <- knots
  res[,2] <- coefs[pos] + covs[pos, pos]^0.5 * knots
   
  #first iteration: 
  init<-fitted$coefficients
  init[pos]<-res[1,2]
  xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, init=init,
                     control=control, modcontrol = update(fitted$modcontrol, terms.fit = (1:k)[-pos])) 
  res[1, 3] <- xx$loglik
  for(i in 2:nn) {
    if(xx$iter == 0 && xx$loglik == 0){
      init<-numeric(length(xx$beta))
      } else {
        init<-xx$beta
      }
    init[pos]<-res[i,2]
    xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, init=init,
                     control=control, modcontrol = update(fitted$modcontrol, terms.fit = (1:k)[-pos])) 
    res[i, 3] <- xx$loglik
   }
  
   signed.root<-sqrt(2*(-res[,3]+max(res[,3])))*sign(res[,2]-fitted$coefficients[pos])
   cdf<-pnorm(signed.root)
   
   results<-list(beta=res[,2], 
                 stdbeta=res[,1], 
                 profile=2*(res[,3]-max(res[,3])), 
                 loglik=res[,3], 
                 signed.root=signed.root, 
                 cdf=cdf)
   attr(results,"class")<-"logistf.profile"
   results
}

