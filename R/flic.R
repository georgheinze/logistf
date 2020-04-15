#' FLIC - Firth's logistic regression with intercept correction
#' 
#' \code{flic} implements Firth's bias-Reduced penalized-likelihood logistic regression with intercept correction.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#' 
#' @param x Either formula and data or \code{\link{logistf}} object
#' @param ... Further arguments passed to the method.
#'
#' @return A \code{flic} object with components:
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
#' flic(case ~ age + oc + vic + vicl + vis + dia, sex2)
#' 
#' #With a logistf object:
#' lf <- logistf(formula = case ~ age + oc + vic + vicl + vis + dia, data = sex2)
#' flic(lf)
#' 
#' @seealso \code{\link{logistf}} for Firth's bias-Reduced penalized-likelihood logistic regression.
#' 
#' @rdname flic
#' @export flic
flic <- function(x,...){
  UseMethod("flic",x)
}
#'
#' @method flic formula
#' @exportS3Method flic formula
#' @describeIn flic With formula and data
flic.formula <- function(formula = attr(data, "formula"), data = sys.parent()){
  
  #Determine coefficient estimates by Firths penalization
  FL <- logistf(formula, data=data)
  designmat <- model.matrix(formula, data)
  response <- lhs.vars(formula) 
  #calculate linear predictors ommiting the intercept
  lp <- FL$linear.predictors-FL$coefficients[1]
  #determine ML estimate of intercept 
  fit <- glm(as.formula(paste(response, paste("1"), sep=" ~ ")), family=binomial(link=logit), 
             data=data, offset=lp)
  #se of intercept
  W <- diag(fit$fitted.values*(1-fit$fitted.values))
  tmp.var <- solve(t(designmat)%*%W%*%designmat)
  beta0.se <- sqrt(tmp.var[1,1])
  
  ic <- fit$coef
  res <- list(coefficients=c(ic, FL$coef[-1]), predicted.probabilities = fit$fitted, linear.predictions=fit$linear, 
              probabilities=c(summary(fit)$coef[, "Pr(>|z|)"], FL$prob[-1]),ci.lower=c(ic-beta0.se*1.96, FL$ci.lower[-1]),
              ci.upper=c(ic+beta0.se*1.96, FL$ci.upper[-1]),call=match.call(), alpha = FL$alpha, 
              method=FL$method, method.ci=FL$method.ci, var=c(beta0.se, diag(FL$var)[-1]^0.5))
  attr(res, "class") <- c("flic")
  res
}

#'
#'
#' @describeIn flic With logistf object
#' @method flic logistf
#' @exportS3Method flic logistf
flic.logistf <- function(lfobject){
  #calculate linear predictors ommiting the intercept
  lp <- lfobject$linear.predictors-lfobject$coefficients[1]
  #determine ML estimate of intercept 
  response <- all.vars(lfobject$formula)[1]
  lfformula <- as.formula(paste(response, paste("1"), sep=" ~ "))
  fit <- glm(lfformula, family=binomial(link=logit), data=lfobject$data, offset=lp)
  #se of intercept
  W <- diag(fit$fitted.values*(1-fit$fitted.values))
  designmat <- model.matrix(lfobject$formula, lfobject$data)
  if( det(t(designmat)%*%W%*%designmat) ) stop('Fisher Information matrix is singular')
  tmp.var <- solve(t(designmat)%*%W%*%designmat)
  beta0.se <- sqrt(tmp.var[1,1])
  ic <- fit$coef
  res <- list(coefficients=c(ic, lfobject$coef[-1]), fitted = fit$fitted, linear.predictions=fit$linear, 
              probabilities=c(summary(fit)$coef[, "Pr(>|z|)"], lfobject$prob[-1]),ci.lower=c(ic-beta0.se*1.96, lfobject$ci.lower[-1]),
              ci.upper=c(ic+beta0.se*1.96, lfobject$ci.upper[-1]),call=match.call(), alpha = lfobject$alpha, 
              method=lfobject$method, method.ci=lfobject$method.ci, var=c(beta0.se, diag(lfobject$var)[-1]^0.5))
  attr(res, "class") <- c("flic")
  res
}
