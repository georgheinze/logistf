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
#' The following generic methods are available for flac's output object: \code{print, summary, coef, confint, anova, extractAIC, add1, drop1, 
#' profile, terms, nobs, predict}. Furthermore, forward and backward functions perform convenient variable selection. Note 
#' that anova, extractAIC, add1, drop1, forward and backward are based on penalized likelihood 
#' ratios.
#' @encoding UTF-8
#' 
#' @param formula A formula object, with the response on the left of the operator, 
#' and the model terms on the right. The response must be a vector with 0 and 1 or \code{FALSE} and 
#' \code{TRUE} for the outcome, where the higher value (1 or \code{TRUE}) is modeled.
#' @param data If using with formula, a data frame containing the variables in the model. 
#' @param lfobject A fitted \code{\link{logistf}} object
#' @param ... Further arguments passed to the method or \code{\link{logistf}}-call.
#'
#' @return A \code{flac} object with components:
#'   \item{coefficients}{The coefficients of the parameter in the fitted model.}
#'   \item{predict}{A vector with the predicted probability of each observation}
#'   \item{linear.predictors}{A vector with the linear predictor of each observation.}
#'   \item{prob}{The p-values of the specific parameters}
#'   \item{ci.lower}{The lower confidence limits of the parameter.}
#'   \item{ci.upper}{The upper confidence limits of the parameter.}
#'   \item{call}{The call object.}
#'   \item{alpha}{The significance level: 0.95}
#'   \item{var}{The variance-covariance-matrix of the parameters.}
#'   \item{loglik}{A vector of the (penalized) log-likelihood of the restricted and the full models.}
#'   \item{n}{The number of observations.}
#'   \item{formula}{The formula object.}
#'   \item{augmented.data}{The augmented dataset used}
#'   \item{df}{The number of degrees of freedom in the model.}
#'   \item{method}{depending on the fitting method 'Penalized ML' or `Standard ML'.}
#'   \item{method.ci}{the method in calculating the confidence intervals, i.e. `profile likelihood' or `Wald', depending on the argument pl and plconf.}
#'   \item{control}{a copy of the control parameters.}  
#'   \item{terms}{the model terms (column names of design matrix).}
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
#' @importFrom stats formula
#' @rdname flac
#' @export flac
flac <- function(...){
  UseMethod("flac")
}
#' @method flac formula
#' @exportS3Method flac formula
#' @describeIn flac With formula and data
flac.formula <- function(formula, data, ...){
  extras <- list(...)
  call_out <- match.call()
  
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.response(mf)
  n <- length(y)
  x <- model.matrix(mt, mf)
  
  if (!is.null(extras$terms.fit)){
    colfit <- eval(extras$terms.fit)
    call_out$terms.fit <- extras$terms.fit
  }
  
  response <- formula.tools::lhs.vars(formula)
  scope <- formula.tools::rhs.vars(formula)
  
  temp.fit1 <- logistf(formula, data=mf, pl=F,...)
  
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
  temp.fit2 <- logistf(newform,data=newdat, weights=temp.neww, firth=FALSE, ...)
  temp.fit3 <- logistf(newresp ~ temp.pseudo,data=newdat, weights=temp.neww, firth=FALSE, ...)
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients) & "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(temp.fit1$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(y)]
  prob <- temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))]
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  
  res <- list(coefficients=coefficients,
              alpha = temp.fit1$alpha, 
              terms = colnames(x),
              var=temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)], 
              df = (temp.fit1$df),
              loglik = c(temp.fit2$loglik[2],temp.fit3$loglik[2]),
              n=temp.fit1$n,
              formula=formula(formula), 
              call=call_out,
              linear.predictors=linear.predictors, 
              predict = fitted, 
              prob=prob,
              method = temp.fit2$method,
              method.ci = temp.fit2$method.ci[-length(temp.fit2$method.ci)], 
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              control = temp.fit2$control, 
              augmented.data = newdat)
  attr(res, "class") <- c("flac")
  res
}

#' @method flac logistf
#' @exportS3Method flac logistf
#' @describeIn flac With logistf object
flac.logistf <- function(lfobject, ... ){
  extras <- list(...)
  call_out <- match.call()
  
  mf <- match.call(expand.dots =FALSE)
  m <- match("lfobject", names(mf), 0L)
  mf <- mf[c(1, m)]
  lfobject <- eval(mf$lfobject, parent.frame())
  variables <- lfobject$terms[-1]
  data <- model.frame(lfobject)
  
  if (!is.null(extras$terms.fit)){
    colfit <- eval(extras$terms.fit)
    call_out$terms.fit <- extras$terms.fit
  }
  
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
  temp.fit2 <- update(lfobject, formula. = newform, data=newdat, weights=temp.neww, firth=FALSE)
  temp.fit3 <- update(lfobject, formula. = newresp ~ temp.pseudo, data=newdat, weights=temp.neww, firth=FALSE)
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients)& "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(lfobject$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(lfobject$y)]
  prob <- temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))]
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  
  res <- list(coefficients=coefficients, 
              alpha = lfobject$alpha,
              terms=lfobject$terms, 
              var=temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)], 
              df = lfobject$df,  
              loglik = c(temp.fit2$loglik[2],temp.fit3$loglik[2]),
              n=lfobject$n,
              formula=formula(lfobject$formula),
              call=call_out,
              linear.predictors=linear.predictors, 
              predict = fitted, prob=prob,
              method = temp.fit2$method,
              method.ci = temp.fit2$method.ci[-length(temp.fit2$method.ci)],
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              augmented.data = newdat, 
              control = temp.fit2$control
              )
  attr(res, "class") <- c("flac")
  res
}
