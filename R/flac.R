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
#' @param data A data frame containing the variables in the model. 
#' @param lfobject A fitted \code{\link{logistf}} object
#' @param model If TRUE the corresponding components of the fit are returned.
#' @param control Controls iteration parameter. Taken from \code{logistf}-object when specified. Otherwise default is \code{control= logistf.control()}
#' @param fitcontrol  Controls additional parameter for fitting. Taken from \code{logistf}-object when specified. Otherwise default is \code{logistf.fit.control = logistf.fit.control()}
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
#'   \item{fitcontrol}{a copy of the fitcontrol parameters.}  
#'   \item{terms}{the model terms (column names of design matrix).}
#'   \item{model}{if requested (the default), the model frame used.}
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
#' flac(lf, data=sex2)
#' 
#' @references Puhr, R., Heinze, G., Nold, M., Lusa, L., and Geroldinger, A. (2017) Firth's logistic regression with rare events: accurate effect estimates and predictions?. Statist. Med., 36: 2302-2317. doi: 10.1002/sim.7273.
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
flac.formula <- function(formula, data, model=TRUE, control, fitcontrol, ...){
  extras <- list(...)

  if(missing(control)){
    control <- logistf.control()
  }
  if(missing(fitcontrol)){
    fitcontrol <- logistf.fit.control()
  }
  
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
  
  response <- formula.tools::lhs.vars(formula)
  scope <- formula.tools::rhs.vars(formula)
  
  temp.fit1 <- logistf(formula, data=data, pl=F, control = control, fitcontrol = fitcontrol, ...)
  
  #apply firths logistic regression and calculate diagonal elements h_i of hat matrix
  #and construct augmented dataset and definition of indicator variable g
  temp.pseudo <- c(rep(0,length(y)), rep(1,2*length(y)))
  temp.neww <- c(rep(1,length(y)), temp.fit1$hat/2, temp.fit1$hat/2)

  data.inverse <- data
  data.inverse[,match(response, names(data.inverse))] <- 1- data.inverse[,match(response, names(data.inverse))]
  newdat <- data.frame(rbind(data, data, data.inverse), temp.pseudo=temp.pseudo, temp.neww=temp.neww)
  names(newdat)[match(names(newdat), response)] <- "newresp"
  
  #ML estimation on augmented dataset
  newform <- update(formula, ~ .+temp.pseudo)
  newform <- update(newform, newresp ~ . )
  temp.fit2 <- logistf(newform,data=newdat, weights=temp.neww, firth=FALSE, control = control, fitcontrol = fitcontrol, ...)
  temp.fit3 <- logistf(newresp ~ temp.pseudo,data=newdat, weights=temp.neww, firth=FALSE, ...)
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients) & "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(temp.fit1$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(y)]
  prob <- temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))]
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  var <- temp.fit2$var
  rownames(var) <- colnames(var) <- names(temp.fit2$coefficients)
  var <- var[which("temp.pseudo"!=names(temp.fit2$coefficients)), which("temp.pseudo"!=names(temp.fit2$coefficients))]
  method.ci <- temp.fit2$method.ci[which("temp.pseudo"!=names(temp.fit2$coefficients))]
  
  res <- list(coefficients=coefficients,
              alpha = temp.fit1$alpha, 
              terms = colnames(x),
              var=var,
              df = (temp.fit1$df),
              loglik = c(temp.fit3$loglik[2],temp.fit2$loglik[2]),
              n=temp.fit1$n,
              formula=formula(formula), 
              call=match.call(),
              linear.predictors=linear.predictors, 
              predict = fitted, 
              prob=prob,
              method = temp.fit2$method,
              method.ci = method.ci, 
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              control = control, 
              fitcontrol = fitcontrol,
              augmented.data = newdat)
  if(model) res$model <- mf
  attr(res, "class") <- c("flac")
  res
}

#' @method flac logistf
#' @exportS3Method flac logistf
#' @describeIn flac With logistf object
flac.logistf <- function(lfobject, data, model=TRUE, ... ){
  extras <- list(...)

  mf <- match.call(expand.dots =FALSE)
  m <- match(c("lfobject", "data"), names(mf), 0L)
  mf <- mf[c(1, m)]
  lfobject <- eval(mf$lfobject, parent.frame())
  
  fitcontrol <- lfobject$fitcontrol
  
  if(!is.null(extras$formula)) lfobject <- update(lfobject, formula. = extras$formula) #to update flac.logistf objects at least with formula

  variables <- lfobject$terms[-1]
  
  if(lfobject$flic) stop("Please call flac() only on logistf-objects with flic=FALSE")
  
  response <- formula.tools::lhs.vars(lfobject$formula)
  scope <- formula.tools::rhs.vars(lfobject$formula)
  
  temp.pseudo <- c(rep(0,length(lfobject$y)), rep(1,2*length(lfobject$y)))
  temp.neww <- c(rep(1,length(lfobject$y)), lfobject$hat/2, lfobject$hat/2)
  
  data.inverse <- data
  data.inverse[,match(response, names(data.inverse))] <- 1- data.inverse[,match(response, names(data.inverse))]
  newdat <- data.frame(rbind(data, data, data.inverse), temp.pseudo=temp.pseudo, temp.neww=temp.neww)
  names(newdat)[match(names(newdat), response)] <- "newresp"
  
  #ML estimation on augmented dataset
  newform <- update(lfobject$formula, ~ .+temp.pseudo)
  newform <- update(newform, newresp ~ . )
  temp.fit2 <- update(lfobject, formula. = newform, data=newdat, weights=temp.neww, firth=FALSE)
  temp.fit3 <- update(lfobject, formula. = newresp ~ temp.pseudo, data=newdat, weights=temp.neww, firth=FALSE)
  
  #outputs
  coefficients <- temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients)& "`(weights)`"!=names(temp.fit2$coefficients))]
  fitted <- temp.fit2$predict[1:length(lfobject$y)]
  linear.predictors <- temp.fit2$linear.predictors[1:length(lfobject$y)]
  prob <- temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))]
  ci.lower <- temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))]
  ci.upper <- temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))]
  var <- temp.fit2$var
  rownames(var) <- colnames(var) <- names(temp.fit2$coefficients)
  var <- var[which("temp.pseudo"!=names(temp.fit2$coefficients)), which("temp.pseudo"!=names(temp.fit2$coefficients))]
  method.ci <- temp.fit2$method.ci[which("temp.pseudo"!=names(temp.fit2$coefficients))]
  
  res <- list(coefficients=coefficients, 
              alpha = lfobject$alpha,
              terms=lfobject$terms, 
              var=var, 
              df = lfobject$df,  
              loglik = c(temp.fit3$loglik[2],temp.fit2$loglik[2]),
              n=lfobject$n,
              formula=formula(lfobject$formula),
              call=match.call(),
              linear.predictors=linear.predictors, 
              predict = fitted, prob=prob,
              method = temp.fit2$method,
              method.ci = method.ci,
              ci.lower=ci.lower,
              ci.upper=ci.upper,
              augmented.data = newdat, 
              control = lfobject$control, 
              fitcontrol = fitcontrol
              )
  
  if(model) {
    res$model <-lfobject$model
  }
  attr(res, "class") <- c("flac")
  res
}
