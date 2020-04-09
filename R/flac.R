#' FLAC - Firth's logistic regression with added covariate
#'
#' \code{flac} implements Firth's bias-Reduced penalized-likelihood logistic regression with added covariate.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#' 
#' @param x Either formula and data or \code{\link{logistf}} object
#' @param ... Further arguments passed to the method.
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
#'   \item{method}{The fitting method: 'Penalized ML'}
#'   \item{method.ci}{The method calculating the confidence intervals of the parameter (without intercept.)}
#'   \item{var}{The variance-covariance-matrix of the parameters.}
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
flac.formula <- function(formula = attr(data, "formula"), data = sys.parent()){
  
  #apply firths logistic regression and calculate diagonal elements h_i of hat matrix
  #and construct augmented dataset and definition of indicator variable g
  temp.fit1 <- logistf(formula, data=data, pl=FALSE)
  temp.pseudo <- c(rep(0,length(temp.fit1$y)), rep(1,2*length(temp.fit1$y)))
  temp.neww <- c(rep(1,length(temp.fit1$y)), temp.fit1$hat/2, temp.fit1$hat/2)
  response <- formula.tools::lhs.vars(formula)
  newdat <- data.frame(y=c(temp.fit1$y, temp.fit1$y, 1-temp.fit1$y),
                       rbind(temp.fit1$data[-which(names(temp.fit1$data) %in% c(response))],
                             temp.fit1$data[-which(names(temp.fit1$data) %in% c(response))],
                             temp.fit1$data[-which(names(temp.fit1$data) %in% c(response))]), 
                       temp.pseudo=temp.pseudo)
  #ML estimation on augmented dataset
  temp.fit2 <- logistf(y ~.,data=newdat, weights=temp.neww, family=binomial(logit), firth=FALSE, pl=TRUE)
  
  res <- list(coefficients=temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients))],
              fitted = temp.fit2$predict[1:length(temp.fit1$y)], 
              linear.predictions=temp.fit2$linear.predictors[1:length(temp.fit1$y)], 
              probabilities=temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))],
              ci.lower=temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))],
              ci.upper=temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))],
              call=match.call(), alpha = temp.fit1$alpha, 
              var=diag(temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)])^0.5)
  attr(res, "class") <- c("flac")
  res
}

#' @method flac logistf
#' @exportS3Method flac logistf
#' @describeIn flac With logistf object
flac.logistf <- function(lfobject){
  response <- formula.tools::lhs.vars(lfobject$formula)
  temp.pseudo <- c(rep(0,length(lfobject$y)), rep(1,2*length(lfobject$y)))
  temp.neww <- c(rep(1,length(lfobject$y)), lfobject$hat/2, lfobject$hat/2)
  newdat <- data.frame(y=c(lfobject$y, lfobject$y, 1-lfobject$y), 
                       rbind(lfobject$data[-which(names(temp.fit1$data) %in% c(response))],
                             lfobject$data[-which(names(temp.fit1$data) %in% c(response))],
                             lfobject1$data[-which(names(temp.fit1$data) %in% c(response))]),
                       temp.pseudo=temp.pseudo)
  #ML estimation on augmented dataset
  temp.fit2 <- logistf(y ~.,data=newdat, weights=temp.neww, family=binomial(logit), firth=FALSE, pl=TRUE)
  
  res <- list(coefficients=temp.fit2$coefficients[which("temp.pseudo"!=names(temp.fit2$coefficients))],
              fitted = temp.fit2$predict[1:length(lfobject$y)], 
              linear.predictions=temp.fit2$linear.predictors[1:length(lfobject$y)], 
              probabilities=temp.fit2$prob[which("temp.pseudo"!=names(temp.fit2$prob))],
              ci.lower=temp.fit2$ci.lower[which("temp.pseudo"!=names(temp.fit2$ci.lower))],
              ci.upper=temp.fit2$ci.upper[which("temp.pseudo"!=names(temp.fit2$ci.upper))],
              call=match.call(), alpha = lfobject$alpha, 
              var=diag(temp.fit2$var[-nrow(temp.fit2$var), -ncol(temp.fit2$var)])^0.5)
  attr(res, "class") <- c("flac")
  res
}
