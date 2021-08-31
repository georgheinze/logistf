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
#' @param plot If \code{TRUE}, profile likelihood is plotted. This parameter becomes obsolete as a 
#' generic plot function is now provided.
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
#' @author Georg Heinze and Meinhard Ploner
#' @export
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
  #if(is.null(data)) stop("Call logistf with dataout=TRUE.\n")
  
  # Next line added by Harry Southworth, 22/10/02.
   if (missing(which) & missing(variable)) stop("You must specify a variable: either by which (a one-sided formula) or by variable.")
   if (missing(control)) control<-logistf.control()
   if (missing(plcontrol)) plcontrol<-logistpl.control()
   
     call <- match.call()

  y <- fitted$y
  n <- length(y)
  x <- model.matrix(fitted$formula, model.frame(fitted)) ## Model-Matrix 
  cov.name <- labels(x)[[2]]
  # weight <- as.vector(model.weights(mf)  )
  offset <- model.offset(mf)   
  weight<-model.weights(mf)
  if (is.null(offset)) offset<-rep(0,n)
  else offset<-as.vector(offset)
  if (is.null(weight)) weight<-rep(1,n)
  
      
  cov.name <- labels(x)[[2]]
  k <- ncol(x)
      if (dimnames(x)[[2]][1] == "(Intercept)")  {
          int <- 1
          coltotest <- 2:k
      }
  
      else {
          int <- 0
          coltotest <-1:k
      }
    if(!missing(which)) cov.name2 <- colnames(get_all_vars(which, data = model.frame(fitted))) ## Label des Test-Fakt.
    else cov.name2 <- variable
    pos <- match(cov.name2, cov.name) ## Position des Testfakors
    if(is.na(pos)) {
      stop(paste(variable,"is not a model term."))
    }
    fit<-logistf.fit(x, y, weight=weight, offset=offset, firth=firth, control=control) 
    std.pos <- diag(fit$var)[pos]^0.5
   
   coefs <- fit$beta ## "normale" Koeffizienten
   covs <- fit$var ## Varianzen
  # n <- nrow(data)
   n <- nrow(x)
   cov.name <- labels(x)[[2]]
   if(missing(limits)) {
    lim.pl<-numeric(0)
    LL.0 <- fit$loglik - qchisq(1 - alpha, 1)/2
    lower.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=-1, i=pos, plcontrol=plcontrol)
    lim.pl[1]<-lower.fit$beta
    upper.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=1, i=pos, plcontrol=plcontrol)
    lim.pl[2]<-upper.fit$beta
    lim.pl <- (lim.pl - coefs[pos])/std.pos
    limits <- c(min(qnorm(alpha/2), lim.pl[1]) - 0.5, max(qnorm(1 - alpha/2), lim.pl[2]) + 0.5)
   }
  
   limits <- c(floor(limits[1]/pitch) * pitch, ceiling(limits[2]/pitch) * pitch)
  
   knots <- seq(limits[1], limits[2], diff(limits)/steps)
   nn <- length(knots)
   res <- matrix(knots, nn, 3) #initialisiere Werte
   dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
   for(i in 1:nn) {
    res[i, 2] <- coefs[pos] + covs[pos, pos]^0.5 * knots[i]
    if(i == 1){
       init<-lower.fit$betahist[nrow(lower.fit$betahist),]
       init[pos]<-res[i,2]
       xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, init=init,
                     control=control, fitcontrol = update(fitted$fitcontrol, terms.fit = (1:k)[-pos])) 
    }     
    else {
       init<-xx$beta
       init[pos]<-res[i,2]
       xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, init=init,
                     control=control, fitcontrol = update(fitted$fitcontrol, terms.fit = (1:k)[-pos])) # use solution from last step
    }
    res[i, 3] <- xx$loglik
   }
  
   #### Graphischer Output:
  
   if(plot==TRUE){
     my.par <- act.par <- par()
     my.par$mai[3] <- 1.65 * act.par$mai[3]
  ## if(legends) my.par$mai[1] <- 2 * act.par$mai[1]
     par(mai = my.par$mai)
     ind <- (1:nn)[round(4 * res[, 1]) == round(4 * res[, 1], 10)]
     if(length(ind) == 0) ind <- 1:nn
     pp <- max(res[, 3]) - 0.5 * res[, 1]^2
  
     plot(res[, -1], type = "l", xlab=expression(beta)) ##Profile likelihood
  
   #lines(res[,2], pp, lty=4)  #<<<Wald approximative profile lik. >>>
  
     points(res[res[, 1] == 0, 2], max(res[, 3])) ##Maximum of likelihood
  
     segments(min(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1),
                  max(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1), lty = 3) ##refer.line
  
     yy <- par("usr")[4] - (par("usr")[4] - par("usr")[3]) * c(0.9, 0.95)
  
     segments(fit$beta[pos] - qnorm(alpha/2) * std.pos, yy[1], fit$beta[pos] - qnorm(1 - alpha/2) *
              std.pos, yy[1], lty = 6) ##Wald-CI
     segments(lower.fit$beta, yy[2], upper.fit$beta, yy[2], lty = 8) ##prof.pen.lik.-CI
  
     axis(side = 3, at = res[ind, 2], labels = res[ind, 1])
  
     mtext(expression(paste("distance from ", hat(beta)," in multiples of ", hat(sigma))), side = 3, line = 3)
    ## mtext(expression(paste(beta, " of ", cov.name2)), side=1, line = 3)
     par(mai = act.par$mai)
   
     if (legends)
      {
       legend(x=fit$beta[pos],
              y=min((min(res[,3])+max(res[,3]))/2,(max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1))),
          legend=c("Profile penalized likelihood",
                   paste(100 * (1 - alpha),"%-reference line"),
                   "Wald confidence interval",
                   "Profile likelihood confidence interval"),
          lty=c(1,3,6,8), 
          text.col=c("black","black","black","black"), ncol=1, bty="n", xjust=0.5)
      }
  
  
      title(paste("Profile of penalized likelihood for Variable",cov.name2))
   }
   signed.root<-sqrt(2*(-res[,3]+max(res[,3])))*sign(res[,1])
   cdf<-pnorm(signed.root)
   
   results<-list(beta=res[,2], stdbeta=res[,1], profile=2*(res[,3]-max(res[,3])), loglike=res[,3], signed.root=signed.root, cdf=cdf)
   attr(results,"class")<-"logistf.profile"
   results
}

