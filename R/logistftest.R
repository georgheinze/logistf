#' Penalized likelihood ratio test
#' 
#' This function performs a penalized likelihood ratio test on some (or all) selected factors. 
#' The resulting object is of the class logistftest and includes the information printed by the
#' proper print method.
#' 
#' This function performs a penalized likelihood ratio test on some (or all) selected factors. The resulting object is of the class logistftest and includes the information printed by the proper print 
#' method. Further documentation can be found in Heinze & Ploner (2004). 
#' In most cases, the functionality of the logistftest function is replaced by anova.logistf, which 
#' is a more standard way to perform likelihood ratio tests. However, as shown in the example below, logistftest provides some specials such as testing against non-zero values. (By the way, 
#' anova.logistf calls logistftest. 
#' 
#' @param object A fitted \code{logistf} object
#' @param test righthand formula of parameters to test (e.g. ~ B + D - 1). As default 
#' all parameter apart from the intercept are tested. If the formula includes -1, the 
#' intercept is omitted from testing. As alternative to the formula one can give 
#' the indexes of the ordered effects to test (a vector of integers). To test only the 
#' intercept specify test = ~ - . or test = 1.
#' @param values Null hypothesis values, default values are 0. For testing the specific hypothesis 
#' B1=1, B4=2, B5=0 we specify test= ~B1+B4+B5-1 and values=c(1, 2,0).
#' @param firth Use of Firth's (1993) penalized maximum likelihood (firth=TRUE, default) or 
#' the standard maximum likelihood method (firth=FALSE) for the logistic regression. 
#' Note that by specifying pl=TRUE and firth=FALSE (and probably lower number of iterations) 
#' one obtains profile likelihood confidence intervals for maximum likelihood logistic 
#' regression parameters.
#' @param beta0 Specifies the initial values of the coefficients for the fitting algorithm
#' @param weights Case weights
#' @param control Controls parameters for iterative fitting
#' @param col.fit.object Numerical vector containing the positions of the variables to fit, if not
#' specified: all variables are taken
#' @param ... further arguments passed to logistf.fit
#'
#' @return The object returned is of the class logistf and has the following attributes:
#'   \item{testcov}{A vector of the fixed values of each covariate; NA stands for a parameter which is not tested.}
#'   \item{loglik}{A vector of the (penalized) log-likelihood of the full and the restricted models. If 
#'   the argument beta0 not missing, the full model isn't evaluated}
#'   \item{df}{The number of degrees of freedom in the model}
#'   \item{prob}{The p-value of the test}
#'   \item{call}{The call object}
#'   \item{method}{Depending on the fitting method 'Penalized ML' or 'Standard ML'}
#'   \item{beta}{The coefficients of the restricted solution}
#'
#' @author Georg Heinze 
#' @references 
#' Firth D (1993). Bias reduction of maximum likelihood estimates. Biometrika 80, 27-38.
#' 
#' 
#' Heinze G, Ploner M (2004). Technical Report 2/2004: A SAS-macro, S-PLUS library and R 
#' package to perform logistic regression without convergence problems. Section of Clinical Biometrics, Department of Medical Computer Sciences, Medical University of Vienna, Vienna, Austria. 
#' http://www.meduniwien.ac.at/user/georg.heinze/techreps/tr2_2004.pdf 
#' 
#' Heinze G (2006). A comparative investigation of methods for logistic regression with separated or 
#' nearly separated data. Statistics in Medicine 25: 4216-4226
#' 
#'  
#'   
#' @export
#'
#' @examples
#' data(sex2) 
#' fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
#' logistftest(fit, test = ~ vic + vicl - 1, values = c(2, 0))
#' 
#' 
logistftest <-
function(object, test, values, firth = TRUE, beta0, weights, control, col.fit.object = NULL, ...)
{
    call <- match.call()
    formula<-object$formula
    
    mf <- match.call(expand.dots =FALSE)
    m <- match(c("object"), names(mf), 0L)
    mf <- mf[c(1, m)]
    object <- eval(mf$object, parent.frame())
    
    y <- model.response(model.frame(object))
    n <- length(y)
    x <- model.matrix(object$formula, model.frame(object))

    if (missing(control)) control<-object$control

    cov.name <- labels(x)[[2]]
    if(missing(weights) & !is.null(object$weights)) weight <- object$weights
    else weight<-NULL
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    k <- ncol(x)
    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }

###    fit.full<-logistf.fit(    ) # unrestricted, define init and col.fit from values, beta0 and test
###    fit.null<-logistf.fit(    ) # restricted, define init and col.fit from values, beta0 and test
    if(!is.null(col.fit.object)){
        fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=col.fit.object, control=control)
    }
    else {
        fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, control=control)
    }
    
    if(fit.full$iter>=control$maxit){
        warning(paste("Maximum number of iterations exceeded. Try to increase the number of iterations by passing 'logistf.control(maxit=...)' to parameter control"))
    }

    pos<-coltotest
    if(missing(test)) {
        test <- coltotest
    }
    if(is.vector(test)) {
        cov.name2 <- cov.name[test]
    }
    else {
        test <- eval(test, parent.frame())
        cov.name2 <- labels(model.matrix(test, model.frame(object)))[[2]]
        #cov.name2 <- attr(terms(test), "term.labels")
    }
    pos <- match(cov.name2, cov.name)   ## Position der Testfakt.
    OK <- !is.na(pos)
    pos <- pos[OK]
    cov.name2 <- cov.name2[OK]
    k2 <- length(cov.name2) ## Anzahl Faktoren
    if(!missing(beta0)) {
        offset1 <- beta0
    }
    else {
        offset1 <- rep(0, k)    ## Vektor der fixierten Werte
    }
    if(!missing(values)) {
        offset1[pos] <- values
    }
    beta <- offset1  ########################################
    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=(1:k)[-pos], control=control, init=beta)
    loglik<-c(fit.null$loglik,fit.full$loglik)
    
    offset1[ - pos] <- NA
    names(offset1) <- cov.name
    fit <- list(testcov = offset1, loglik = loglik, df = k2, prob = 1 - pchisq(2 *
        diff(loglik), k2), call = match.call(), beta = beta)
    if(firth) {
        fit$method <- "Penalized ML"
    }
    else {
        fit$method <- "Standard ML"
    }
    attr(fit, "class") <- "logistftest"
    fit
}

