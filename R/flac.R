#' FLAC - Firth's logistic regression with added covariate
#'
#' \code{flac} implements Firth's bias-reduced penalized-likelihood logistic regression with added covariate.
#'
#' Flac is a simple modification of Firth's logistic regression which provides average predicted
#' probabilities equal to the observed proportion of events, while preserving the ability to deal
#' with seperation.
#'
#' The modified score equation to estimate coefficients for Firth's logistic regression can be 
#' interpreted as score equations for ML estimates for an augmented data set. This data set can be 
#' created by complementing each original observation i with two pseudo-observations weighted by 
#' \eqn{h_i/2} with unchanged covariate values and with response values set to \eqn{y=0} and \eqn{y=1}
#' respectively. The basic idea of Flac is to discriminate between original and pseudo-observations
#' in the alternative formulation of Firth's estimation as an iterative data augmentation procedure.
#' The following generic methods are available for flac's output object: \code{print, summary, coef, vcov, confint, anova, extractAIC, add1, drop1, 
#' profile, terms, nobs, predict}. Furthermore, forward and backward functions perform convenient variable selection. Note 
#' that anova, extractAIC, add1, drop1, forward and backward are based on penalized likelihood 
#' ratios.
#' @encoding UTF-8
#' 
#' @param x Either formula or \code{\link{logistf}} object
#' @param data If using with formula, a data framecontaining the variables in the model. 
#' @param ... Further arguments passed to the method or \code{\link{logistf}}-call.
#'
#' @return A \code{flac} object with components:
#'   \item{coefficients}{The coefficients of the parameter in the fitted model.}
#'   \item{predicted.probabilities}{A vector with the predicted probability of each observation}
#'   \item{linear.predictions}{A vector with the linear predictor of each observation.}
#'   \item{probabilities}{The p-values of the specific parameters}
#'   \item{ci.lower}{The lower confidence limits of the parameter.}
#'   \item{ci.upper}{The upper confidence limits of the parameter.}
#'   \item{call}{The call object.}
#'   \item{alpha}{The significance level: 0.95}
#'   \item{var}{The variance-covariance-matrix of the parameters.}
#'   \item{loglik}{A vector of the (penalized) log-likelihood of the restricted and the full models.}
#'   \item{n}{The number of observations.}
#'   \item{formula}{The formula object.}
#'   \item{augmented.data}{The augmented dataset used}
#'   \item{terms}{The model terms (column names of design matrix).}
#'   \item{df}{The number of degrees of freedom in the model.}
#'   \item{method}{depending on the fitting method 'Penalized ML' or `Standard ML'.}
#'   \item{method.ci}{the method in calculating the confidence intervals, i.e. `profile likelihood' or `Wald', depending on the argument pl and plconf.}
#'   
#' 
#' @export
#'
#' @examples 
#' #With formula and data:
#' data(sex2)
#' flac(case ~ age + oc + vic + vicl + vis + dia, sex2)
#' 
#' #With a logistf object:
#' lf <- logistf(formula = case ~ age + oc + vic + vicl + vis + dia, data = sex2)
#' flac(lf)
#' 
#' @references Puhr, R., Heinze, G., Nold, M., Lusa, L., and Geroldinger, A. (2017) Firth's logistic regression with rare events: accurate effect estimates and predictions?. Statist. Med., 36: 2302– 2317. doi: 10.1002/sim.7273.
#' 
#' @seealso [logistf()] for Firth's bias-Reduced penalized-likelihood logistic regression.
#' 
#' @rdname flac
#' @export flac
flac <- function(x,...){
  UseMethod("flac",x)
}
#' @method flac formula
#' @exportS3Method flac formula
#' @describeIn flac With formula and data
flac.formula <- function(x, data, ...){
  call <- match.call()
  extras <- list(...)
  call_out <- match.call()  
  
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("x", "data","weights","na.action","offset"), names(mf), 0L)
  
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  formula <- x
  y <- model.response(mf)
  n <- length(y)
  x <- model.matrix(mt, mf)
  
  response <- formula.tools::lhs.vars(formula)
  scope <- formula.tools::rhs.vars(formula)
  
  if (!is.null(extras$terms.fit)){
    termsfit <- eval(extras$terms.fit)
    n_termsfit <- length(termsfit)
    call_out$terms.fit <- extras$terms.fit
    temp.fit1 <- logistf(formula, data=mf, terms.fit=termsfit, pl=F,...)
  }
  else temp.fit1 <- logistf(formula, data=mf, pl=FALSE, ...)
  #apply firths logistic regression and calculate diagonal elements h_i of hat matrix
  #and construct augmented dataset and definition of indicator variable g
  temp.pseudo <- c(rep(0,length(y)), rep(1,2*length(y)))
  temp.neww <- c(rep(1,length(y)), temp.fit1$hat/2, temp.fit1$hat/2)
  newdat <- data.frame("newresp"=c(y, y, 1-y),
                       rbind(mf[-which(names(mf) %in% c(response))],
                             mf[-which(names(mf) %in% c(response))],
                             mf[-which(names(mf) %in% c(response))]), 
                       temp.pseudo=temp.pseudo, temp.neww=temp.neww)
  #ML estimation on augmented dataset
  rhs <- paste(paste(scope, collapse="+"),"temp.pseudo", sep="+")
  newform <- paste("newresp", "~", rhs)
  if (!is.null(extras$terms.fit)) {
    temp.fit2 <- logistf(newform,data=newdat, weights=temp.neww,terms.fit=termsfit, firth=FALSE, plconf=plconf, ...)
  }
  else {
    temp.fit2 <- logistf(newform, data=newdat, weights=temp.neww, firth=FALSE, ...)
  }
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients) & "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(temp.fit1$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(y)]
  prob <- temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))]
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  var <- diag(temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)])^0.5
  
  res <- list(coefficients=coefficients,
              fitted = fitted, 
              method = temp.fit2$method,
              linear.predictions=linear.predictors, 
              probabilities=prob,
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              call=match.call(), alpha = temp.fit1$alpha, 
              var=var, 
              loglik = temp.fit1$loglik, 
              n=temp.fit1$n, formula=formula(formula), augmented.data = newdat, 
              terms=colnames(x)[-1], df = (temp.fit2$df-1), #-1 because temp.fit2 has one variable more: weights 
              formula=formula, 
              method.ci = temp.fit2$method.ci[-length(temp.fit2$method.ci)], 
              control = temp.fit2$control
              )
  attr(res, "class") <- c("flac")
  res
}

#' @method flac logistf
#' @exportS3Method flac logistf
#' @describeIn flac With logistf object
flac.logistf <- function(x, ... ){
  
  mf <- match.call(expand.dots =FALSE)
  m <- match("x", names(mf), 0L)
  mf <- mf[c(1, m)]
  lfobject <- eval(mf$x, parent.frame())
  variables <- lfobject$terms[-1]
  data <- model.frame(lfobject)
  
  response <- formula.tools::lhs.vars(lfobject$formula)
  scope <- formula.tools::rhs.vars(lfobject$formula)
  
  temp.pseudo <- c(rep(0,length(lfobject$y)), rep(1,2*length(lfobject$y)))
  temp.neww <- c(rep(1,length(lfobject$y)), lfobject$hat/2, lfobject$hat/2)
  newdat <- data.frame("newresp"=c(lfobject$y, lfobject$y, 1-lfobject$y), 
                       rbind(data[-which(names(data) %in% c(response))],
                             data[-which(names(data) %in% c(response))],
                             data[-which(names(data) %in% c(response))]),
                       temp.pseudo=temp.pseudo, temp.neww=temp.neww)
  #ML estimation on augmented dataset
  rhs <- paste(paste(scope, collapse="+"),"temp.pseudo", sep="+")
  newform <- paste("newresp", "~", rhs)
  temp.fit2 <- logistf(newform,data=newdat, weights=temp.neww, firth=FALSE, ...)
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients)& "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(lfobject$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(lfobject$y)]
  prob <- lfobject$prob
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  var <- diag(temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)])^0.5
  
  
  res <- list(coefficients=coefficients,
              fitted = fitted, 
              method = temp.fit2$method,
              linear.predictions=linear.predictors, 
              probabilities=prob,
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              call=match.call(), alpha = lfobject$alpha, 
              var=var, 
              loglik = lfobject$loglik, 
              n=lfobject$n, formula=lfobject$formula, augmented.data = newdat, 
              terms=lfobject$terms, df = lfobject$df, 
              formula=lfobject$formula, 
              method.ci = temp.fit2$method.ci[-length(temp.fit2$method.ci)], #exclude estimation of CI of temp.pseudo
              control = temp.fit2$control
              )
  attr(res, "class") <- c("flac")
  res
}
