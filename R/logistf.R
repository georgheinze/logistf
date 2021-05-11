#' Firth's Bias-Reduced Logistic Regression
#'
#' Implements Firth's bias-Reduced penalized-likelihood logistic regression. 
#'
#' \code{logistf} is the main function of the package. It fits a logistic regression 
#' model applying Firth's correction to the likelihood. The following generic methods are available for logistf's output 
#' object: \code{print, summary, coef, vcov, confint, anova, extractAIC, add1, drop1, 
#' profile, terms, nobs, predict}. Furthermore, forward and backward functions perform convenient variable selection. Note 
#' that anova, extractAIC, add1, drop1, forward and backward are based on penalized likelihood 
#' ratios.
#'
#' @param formula A formula object, with the response on the left of the operator, 
#' and the model terms on the right. The response must be a vector with 0 and 1 or \code{FALSE} and 
#' \code{TRUE} for the outcome, where the higher value (1 or \code{TRUE}) is modeled. It is 
#' possible to include contrasts, interactions, nested effects, cubic or polynomial 
#' splines and all S features as well, e.g. Y ~ X1*X2 + ns(X3, df=4).
#' @param data A data.frame where the variables named in the formula can be found, 
#' i. e. the variables containing the binary response and the covariates.
#' @param pl Specifies if confidence intervals and tests should be based on the profile 
#' penalized log likelihood (\code{pl=TRUE}, the default) or on the Wald method (\code{pl=FALSE}).
#' @param alpha The significance level (1-\eqn{\alpha} the confidence level, 0.05 as default).
#' @param control Controls Newton-Raphson iteration. Default is \code{control= logistf.control(maxstep, 
#' maxit, maxhs, lconv, gconv, xconv)}
#' @param plcontrol Controls Newton-Raphson iteration for the estimation of the profile 
#' likelihood confidence intervals. Default is \code{plcontrol= logistpl.control(maxstep, maxit, 
#' maxhs, lconv, xconv, ortho, pr)}
#' @param firth Use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the 
#' standard maximum likelihood method (\code{firth=FALSE}) for the logistic regression. 
#' Note that by specifying \code{pl=TRUE} and \code{firth=FALSE} (and probably a lower number 
#' of iterations) one obtains profile likelihood confidence intervals for maximum likelihood 
#' logistic regression parameters.
#' @param init Specifies the initial values of the coefficients for the fitting algorithm
#' @param weights specifies case weights. Each line of the input data set is multiplied 
#' by the corresponding element of weights
#' @param na.action a function which indicates what should happen when the data contain NAs
#' @param plconf specifies the variables (as vector of their indices) for which profile likelihood 
#' confidence intervals should be computed. Default is to compute for all variables.
#' @param flic If \code{TRUE}, intercept is altered such that the predicted probabilities become unbiased while 
#' keeping all other coefficients constant
#' @param model  If TRUE the corresponding components of the fit are returned.
#' @param tau  Penalization parameter (default = 0.5)
#' @param ... Further arguments to be passed to \code{logistf}
#' 
#' @return The object returned is of the class \code{logistf} and has the following attributes:
#'    \item{coefficients}{the coefficients of the parameter in the fitted model.}
#'    \item{alpha}{the significance level (1- the confidence level) as specified in the input.}
#'    \item{terms}{the column names of the design matrix}
#'    \item{var}{the variance-covariance-matrix of the parameters.}
#'    \item{df}{the number of degrees of freedom in the model.}
#'    \item{loglik}{a vector of the (penalized) log-likelihood of the restricted and the full models.}
#'    \item{iter}{the number of iterations needed in the fitting process.}
#'    \item{n}{the number of observations.}
#'    \item{y}{the response-vector, i. e. 1 for successes (events) and 0 for failures.}
#'    \item{formula}{the formula object.}
#'    \item{call}{the call object.}
#'    \item{terms}{the model terms (column names of design matrix).}
#'    \item{linear.predictors}{ a vector with the linear predictor of each observation.}
#'    \item{predict}{a vector with the predicted probability of each observation.}
#'    \item{hat.diag}{a vector with the diagonal elements of the Hat Matrix.}
#'    \item{conv}{the convergence status at last iteration: a vector of length 3 with elements: last change in log likelihood, max(abs(score vector)), max change in beta at last iteration.}
#'    \item{method}{depending on the fitting method 'Penalized ML' or `Standard ML'.}
#'    \item{method.ci}{the method in calculating the confidence intervals, i.e. `profile likelihood' or `Wald', depending on the argument pl and plconf.}
#'    \item{ci.lower}{the lower confidence limits of the parameter.}
#'    \item{ci.upper}{the upper confidence limits of the parameter.}
#'    \item{prob}{ the p-values of the specific parameters.}
#'    \item{pl.iter}{only if pl==TRUE: the number of iterations needed for each confidence limit.}
#'    \item{betahist}{only if pl==TRUE: the complete history of beta estimates for each confidence limit.}
#'    \item{pl.conv}{only if pl==TRUE: the convergence status (deviation of log likelihood from target value, last maximum change in beta) for each confidence limit.}
#'    \item{control}{a copy of the control parameters.}
#'    \item{flic}{logical, is TRUE  if intercept was altered such that the predicted probabilities become unbiased while 
#' keeping all other coefficients constant. According to input of logistf.}  
#'    \item{model}{if requested (the default), the model frame used.}
#'    \item{fit}{The fitting algorithm used.}
#'    \item{na.action}{information returned by model.frame on the special handling of NAs}
#'     
#' @export
#'
#' @encoding UTF-8
#' @examples
#' data(sex2)
#' fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
#' summary(fit)
#' nobs(fit)
#' drop1(fit)
#' plot(profile(fit,variable="dia"))
#' extractAIC(fit)
#' 
#' fit1<-update(fit, case ~ age+oc+vic+vicl+vis)
#' extractAIC(fit1)
#' anova(fit,fit1)
#' 
#' data(sexagg)
#' fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT)
#' summary(fit2)
#' 
#' # simulated SNP example
#' set.seed(72341)
#' snpdata<-rbind(
#'   matrix(rbinom(2000,2,runif(2000)*0.3),100,20),
#'   matrix(rbinom(2000,2,runif(2000)*0.5),100,20))
#' colnames(snpdata)<-paste("SNP",1:20,"_",sep="")
#' snpdata<-as.data.frame(snpdata)
#' for(i in 1:20) snpdata[,i]<-as.factor(snpdata[,i])
#' snpdata$case<-c(rep(0,100),rep(1,100))
#' 
#' fitsnp<-logistf(data=snpdata, formula=case~1, pl=FALSE)
#' add1(fitsnp, scope=paste("SNP",1:20,"_",sep=""))
#' fitf<-forward(fitsnp, scope = paste("SNP",1:20,"_",sep=""))
#' fitf
#' 
#' @author Georg Heinze and Meinhard Ploner
#' @references Firth D (1993). Bias reduction of maximum likelihood estimates. Biometrika 80, 27-38. 
#' Heinze G, Schemper M (2002). A solution to the problem of separation in logistic regression. 
#' Statistics in Medicine 21: 2409-2419.
#' 
#' Heinze G, Ploner M (2003). Fixing the nonconvergence bug in logistic regression with SPLUS and 
#' SAS. Computer Methods and Programs in Biomedicine 71: 181-187.
#' 
#' Heinze G, Ploner M (2004). Technical Report 2/2004: A SAS-macro, S-PLUS library and R
#' package to perform logistic regression without convergence problems. Section of Clinical Biometrics, Department of Medical Computer Sciences, Medical University of Vienna, Vienna, Austria. 
#' http://www.meduniwien.ac.at/user/georg.heinze/techreps/tr2_2004.pdf
#' 
#' Heinze G (2006). A comparative investigation of methods for logistic regression with separated or 
#' nearly separated data. Statistics in Medicine 25: 4216-4226. 
#' 
#' Venzon DJ, Moolgavkar AH (1988). A method for computing profile-likelihood based confidence 
#' intervals. Applied Statistics 37:87-94.
#' @seealso [add1.logistf()], [drop1.logistf()], [anova.logistf()]
#' @rdname logistf
logistf <-
function(formula, data, pl = TRUE, alpha = 0.05, control, plcontrol, firth = TRUE, init, weights,na.action, plconf=NULL,flic=FALSE, model = TRUE,tau=0.5,  ...){
   call <- match.call()
   extras <- list(...)
   call_out <- match.call()
   if(missing(control)) control<-logistf.control()
   if(pl==TRUE & missing(plcontrol)) plcontrol<-logistpl.control()
   
    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula", "data","weights","na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, type="any")
    if(is.logical(y)){
      y <- as.numeric(y)
    }
    else if(is.factor(y)){
      if(length(levels(y))==2){
        warning(paste0("Ivalid response variable: Changing value ", levels(y)[1], " to 0 and value ",  levels(y)[2], " to 1." ))
        y <- as.numeric(y)-1
      }
      
    }else if(!is.numeric(y)){
      stop("Invalid response variable: must be logical or numeric")
    }
    
    n <- length(y)
    x <- model.matrix(mt, mf)
    
    k <- ncol(x)
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf)) # avoid any problems with 1D or nx1 arrays by as.vector
    offset <- as.vector(model.offset(mf))
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    if (missing(init)) init<-rep(0,k)
    if (is.null(plconf)) {  #if only intercept has to be fitted: calculate Wald CI for intercept
      if(isspecnum(cov.name,"(Intercept)")){
        plconf <- NULL
        rest_plconf <- 1
      }
      else{
      plconf<-1:k
      rest_plconf <- NULL
      }
    }
    else { #if plconf is passed to logistf write NA to variables not specified
      rest_plconf <- (1:k)[-plconf]
    }
    if (cov.name[1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }
    else {
        int <- 0
        coltotest <-1:k
    }
    
    #for backward function to update logistf object
    if (!is.null(extras$terms.fit)){
      colfit <- eval(extras$terms.fit)
      matched <- match(colfit, cov.name)
      if(any(is.na(matched))) stop(paste0("term(s): ", paste(colfit[is.na(matched)],collapse=", ")," not found in formula.\n"))
      nterms <- length(matched)
      colfit <- (1:k)[matched]
      call_out$terms.fit <- extras$terms.fit
      rest_plconf <- (1:k)[-matched] #compute confidence intervals for variables not specified in terms.fit with Wald
      plconf <- intersect(plconf, matched)
    }
    else {
      colfit <- 1:k
      nterms <- k
    }
    
    if(!firth & control$fit == "IRLS"){
        warning("Fitting method IRLS with firth = FALSE currently not implemented. Using Newton-Raphson.")
        control$fit <- "NR"
    }

    fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=colfit, init, control=control, tau=tau)
    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=int, init, control=control, tau=tau)

    
    if(fit.full$iter>=control$maxit){
      warning(paste("logistf.fit: Maximum number of iterations for full model exceeded. Try to increase the number of iterations or alter step size by passing 'logistf.control(maxit=..., maxstep=...)' to parameter control"))
    }
    if(fit.null$iter>=control$maxit){
      warning(paste("logistf.fit: Maximum number of iterations for null model exceeded. Try to increase the number of iterations or alter step size by passing 'logistf.control(maxit=..., maxstep=...)' to parameter control"))
    }
    
    fit <- list(coefficients = fit.full$beta, alpha = alpha, terms=colnames(x), var = fit.full$var, df = nterms-int, loglik =c(fit.null$loglik, fit.full$loglik),
        iter = fit.full$iter, n = sum(weight), y = y, formula = formula(formula), call=call_out, conv=fit.full$conv, tau = tau, fit = fit.full$fit)
    
    names(fit$conv)<-c("LL change","max abs score","beta change")
    beta<-fit.full$beta
    covs<-fit.full$var
    pi<-fit.full$pi
    fit$firth<-firth
    fit$linear.predictors <- as.vector(x %*% beta + offset)
    fit$predict <- fit.full$pi
    fit$hat.diag <- fit.full$Hdiag
    if(firth) fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    
    if (!is.null(extras$terms.fit)){ #consider for wald only covariance matrix with columns corresponding to variables in terms.fit
      loc <- match(extras$terms.fit, fit$terms)
      var.red <- fit$var[loc,loc]
      vars <- diag(as.matrix(var.red))
      waldprob <- wald_ci.lower <- wald_ci.upper <- vector(length = k)
      waldprob[loc] <- 1 - pchisq((beta[loc]^2/vars), 1)
      wald_ci.lower[loc] <- as.vector(beta[loc] + qnorm(alpha/2) * vars^0.5)
      wald_ci.upper[loc] <- as.vector(beta[loc] + qnorm(1 - alpha/2) * vars^0.5)
    }
    else {
      vars <- diag(covs)
      waldprob <- 1 - pchisq((beta^2/vars), 1)
      wald_ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
      wald_ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * vars^0.5)
    }

    fit$alpha<-alpha
    fit$conflev<-1-alpha
    
    if(pl) {
        #initialisation
        betahist.lo<-vector(length(plconf),mode="list")
        betahist.up<-vector(length(plconf),mode="list")
        pl.conv<-matrix(0,length(plconf),4)
        dimnames(pl.conv)[[1]]<-as.list(plconf)
        dimnames(pl.conv)[[2]]<-as.list(c("lower, loglik","lower, beta", "upper, loglik", "upper, beta"))
        LL.0 <- fit.full$loglik - qchisq(1 - alpha, 1)/2
        pl.iter<-matrix(0,k,3)
        icount<-0
        iters <- vector() #number of iterations of fit.i per variable 
        for(i in plconf) {
            icount<-icount+1
            inter<-logistpl(x, y, beta, i, LL.0, firth, -1, offset, weight, plcontrol)
            fit$ci.lower[i] <- inter$beta
            pl.iter[i,1]<-inter$iter
            betahist.lo[[icount]]<-inter$betahist
            pl.conv.lower<-t(inter$conv)
            inter<-logistpl(x, y, beta, i, LL.0, firth, 1, offset, weight, plcontrol)
            fit$ci.upper[i] <- inter$beta
            pl.iter[i,2]<-inter$iter
            betahist.up[[icount]]<-inter$betahist
            pl.conv.upper<-t(inter$conv)
            pl.conv[icount,]<-cbind(pl.conv.lower,pl.conv.upper)
            fit.i<-logistf.fit(x,y, weight=weight, offset=offset, firth, col.fit=setdiff(colfit,i), control=control, tau=tau)
            pl.iter[i,3]<-fit.i$iter
            fit$prob[i] <- 1-pchisq(2*(fit.full$loglik-fit.i$loglik),1)
            fit$method.ci[i] <- "Profile Likelihood"
        }
        fit$pl.iter <- pl.iter
        colnames(fit$pl.iter)<-c("Lower", "Upper", "Null model")
        if(sum(fit$pl.iter>=plcontrol$maxit)){
          notconv <- cov.name[apply(fit$pl.iter>=plcontrol$maxit, 1, sum)>0]
          warning(paste("Nonconverged PL confidence limits: maximum number of iterations for variables:",paste0(notconv,collapse=", ")), " exceeded. Try to increase the number of iterations by passing 'logistpl.control(maxit=...)' to parameter plcontrol")
        }
          
        fit$betahist<-list(lower=betahist.lo, upper=betahist.up)
        fit$pl.conv<-pl.conv
        if(sum(iters>=control$maxit)>0){ #check if algorithms for all models have converged
          notconv <- cov.name[iters>=control$maxit]
          warning(paste("Maximum number of iterations for PLR test for variables:",paste0(notconv,collapse=", ")), " exceeded. P-value may be incorrect. Try to increase the number of iterations by passing 'logistf.control(maxit=...)' to parameter control")
        }
    }
    else { 
      fit$prob <- waldprob
      fit$method.ci <- rep("Wald",k)
      fit$ci.lower <- wald_ci.lower
      fit$ci.upper <- wald_ci.upper
    }
    for (i in rest_plconf){
      fit$ci.lower[i] <- 0
      fit$ci.upper[i] <- 0
      fit$prob[i] <- 0
      fit$method.ci[i] <- "-"
    }
    names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$coefficients) <- dimnames(x)[[2]]
    #flic: 
    if (flic){
        fit$flic <- TRUE
        #calculate linear predictors ommiting the intercept
        lp_flic <-  fit$linear.predictors-fit$coef[1]
        #determine ML estimate of intercept 
        response <- formula.tools::lhs.vars(formula)
        fit_flic <- glm(as.formula(paste(response, paste("1"), sep=" ~ ")), family=binomial(link=logit), 
                   data=data, offset=lp_flic)
        W <- diag(fit_flic$fitted.values*(1-fit_flic$fitted.values))
        tmp.var <- solve(t(x)%*%W%*%x)
        beta0.se <- sqrt(tmp.var[1,1])
        fit$coefficients <- c(fit_flic$coef, fit$coef[-1])
        fit$ci.lower <- c(fit_flic$coef-beta0.se*1.96, fit$ci.lower[-1])
        fit$ci.upper <- c(fit_flic$coef+beta0.se*1.96, fit$ci.upper[-1])
        fit$linear.predictors <- fit_flic$linear
        fit$predict <-fit_flic$fitted
        fit$var <- tmp.var
        fit$method.ci[1] <- "Wald"
      }
    else fit$flic <- FALSE
    
    fit$control <- control
    if (model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    attr(fit, "class") <- c("logistf")
    fit
}


is.logistf<-function(object){
  if(attr(object,"class")=="logistf") return(TRUE)
  else return(FALSE)
  }
  
#' @method coef logistf  
coef.logistf<-function(object,...){
  return(object$coefficients)
}

#' @method confint logistf
#' @exportS3Method confint logistf
confint.logistf<-function(object,parm, level=0.95, exp=FALSE, ...){
  # in fact, level is already determined in the object and will be overwritten
  level<-object$conflev
  levstr<-paste(level*100,"%",sep="")
  cimat<-cbind(object$ci.lower,object$ci.upper)
  rownames(cimat)<-names(object$ci.lower)
  colnames(cimat)<-c(paste("Lower ",levstr,sep=""),paste("Upper ",levstr,sep=""))
  if(exp) cimat<-exp(cimat)
  return(cimat)
}

#' @method vcov logistf
#' @exportS3Method vcov logistf
vcov.logistf<-function(object,...){
  var <- object$var
  colnames(var) <- rownames(var) <- object$terms
  return(var)
}

#' @method family logistf
#' @exportS3Method family logistf
family.logistf<-function(object, ...){
  return(stats::binomial())
}
  