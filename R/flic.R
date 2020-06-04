#' FLIC - Firth's logistic regression with intercept correction
#' 
#' \code{flic} implements Firth's bias-Reduced penalized-likelihood logistic regression with intercept correction.
#'
#' Flic is a simple modification of Firth's logistic regression which provides average predicted
#' probabilities equal to the observed proportion of events, while preserving the ability to deal
#' with seperation.
#' 
#' In general the average predicted probability in FL regression is not equal to the observed 
#' proportion of events. Because the determinant of the Fisher-Information matrix is maximized 
#' for \eqn{\pi_i = \frac{1}{2}} it is concluded that Firth's penalization tends to push the 
#' predicted probabilities towards one-half compared with ML-estimation.
#' Flic fits a logistic regression model applying Firth’s correction to the likelihood with a 
#' correction of the intercept, such that the predicted probabilities become unbiased while
#' keeping all other coefficients constant.
#' The following generic methods are available for flic's output object: \code{print, summary, coef, vcov, confint, anova, extractAIC, add1, drop1, 
#' profile, terms, nobs, predict}. Furthermore, forward and backward functions perform convenient variable selection. Note 
#' that anova, extractAIC, add1, drop1, forward and backward are based on penalized likelihood 
#' ratios.
#' 
#' @param x Either formula and data or \code{\link{logistf}} object
#' @param ... Further arguments passed to the method.
#'
#' @return A \code{flic} object with components:
#'   \item{coefficients}{The coefficients of the parameter in the fitted model.}
#'   \item{terms}{The model terms (column names of design matrix).}
#'   \item{predicted.probabilities}{A vector with the predicted probability of each observation}
#'   \item{linear.predictors}{A vector with the linear predictor of each observation.}
#'   \item{probabilities}{The p-values of the specific parameters}
#'   \item{ci.lower}{The lower confidence limits of the parameter.}
#'   \item{ci.upper}{The upper confidence limits of the parameter.}
#'   \item{call}{The call object.}
#'   \item{alpha}{The significance level: 0.95}
#'   \item{method}{The fitting method: 'Penalized ML'}
#'   \item{method.ci}{The method calculating the confidence intervals of the parameter (without intercept.)}
#'   \item{var}{The variance-covariance-matrix of the parameters.}
#'   \item{df}{The number of degrees of freedom in the model.}
#'   \item{loglik}{A vector of the (penalized) log-likelihood of the restricted and the full models.}
#'   \item{n}{The number of observations.}
#'   \item{formula}{The formula object.}
#'   \item{data}{A copy of the input dataset.}
#'   
#' 
#' @export
#' 
#' @encoding UTF-8
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
#' @references Puhr, R., Heinze, G., Nold, M., Lusa, L., and Geroldinger, A. (2017) Firth's logistic regression with rare events: accurate effect estimates and predictions?. Statist. Med., 36: 2302– 2317. doi: 10.1002/sim.7273.
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
#' @export flic.formula
flic.formula <- function(formula = attr(data, "formula"), data = sys.parent(),...){
  #Determine coefficient estimates by Firths penalization
  extras <- list(...)
  call_out <- match.call()
  if (!is.null(extras$terms.fit)){
    termsfit <- eval(extras$terms.fit)
    call_out$terms.fit <- extras$terms.fit
    # estimate profile likelihood confidence intervals only for variables in terms.fit
    plconf <- match(termsfit, formula.tools::rhs.vars(formula))
    # logistf call
    FL <- logistf(formula, data=data, terms.fit=termsfit,plconf = plconf)
  }
  else FL <- logistf(formula, data=data)
  designmat <- model.matrix(formula, data)
  response <- lhs.vars(formula) 
  #calculate linear predictors ommiting the intercept
  lp <- FL$linear.predictors-FL$coefficients[1]
  #determine ML estimate of intercept 
  fit <- glm(as.formula(paste(response, paste("1"), sep=" ~ ")), family=binomial(link=logit), 
             data=data, offset=lp)
  #se of intercept
  W <- diag(fit$fitted.values*(1-fit$fitted.values))
  XWX <- t(designmat)%*%W%*%designmat
  tmp.var <- solve(XWX)
  beta0.se <- sqrt(tmp.var[1,1])
  
  #compute penalized likelihood: 
  loglik <- sum(log(fit$fitted^(FL$y)*(1-fit$fitted)^(1-FL$y)))
  I <- 0.5*log(det(XWX))
  full_loglik <- loglik+I
  
  ic <- fit$coef
  res <- list(coefficients=c(ic, FL$coef[-1]),terms=colnames(designmat), predicted.probabilities = fit$fitted, linear.predictions=fit$linear, 
              probabilities=c(summary(fit)$coef[, "Pr(>|z|)"], FL$prob[-1]),ci.lower=c(ic-beta0.se*1.96, FL$ci.lower[-1]),
              ci.upper=c(ic+beta0.se*1.96, FL$ci.upper[-1]),call=call_out, alpha = FL$alpha, 
              method=FL$method, method.ci=FL$method.ci, var=c(beta0.se, diag(FL$var)[-1]^0.5), df=FL$df, loglik=c(FL$loglik[1], full_loglik), n=FL$n, 
              formula=formula(formula), data = data)
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
  if( det(t(designmat)%*%W%*%designmat) == 0) stop('Fisher Information matrix is singular')
  XWX <- t(designmat)%*%W%*%designmat
  tmp.var <- solve(XWX)
  beta0.se <- sqrt(tmp.var[1,1])
  #compute penalized likelihood: 
  loglik <- sum(log(fit$fitted^(lfobject$y)*(1-fit$fitted)^(1-lfobject$y)))
  I <- 0.5*log(det(XWX))
  full_loglik <- loglik+I
  
  ic <- fit$coef
  res <- list(coefficients=c(ic, lfobject$coef[-1]), fitted = fit$fitted, linear.predictors=fit$linear, 
              probabilities=c(summary(fit)$coef[, "Pr(>|z|)"], lfobject$prob[-1]),ci.lower=c(ic-beta0.se*1.96, lfobject$ci.lower[-1]),
              ci.upper=c(ic+beta0.se*1.96, lfobject$ci.upper[-1]),call=match.call(), alpha = lfobject$alpha, 
              method=lfobject$method, method.ci=lfobject$method.ci, var=c(beta0.se, diag(lfobject$var)[-1]^0.5), df=lfobject$df-1, loglik=c(lfobject$loglik[1], full_loglik), n=lfobject$n, 
              formula=lfobject$formula, data = lfobject$data)
  attr(res, "class") <- c("flic")
  res
}
