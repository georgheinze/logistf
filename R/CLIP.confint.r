#' Confidence Intervals after Multiple Imputation: Combination of Likelihood Profiles
#' 
#' This function implements the new combination of likelihood profiles (CLIP) method described in 
#' Heinze, Ploner and Beyea (2013). This method is useful for computing confidence intervals for 
#' parameters after multiple imputation of data sets, if the normality assumption on parameter estimates and consequently the validity of applying Rubin‘s rules (pooling of variances) is in doubt. It 
#' consists of combining the profile likelihoods into a posterior. The function CLIP.confint searches 
#' for those values of a regression coefficient, at which the cumulative distribution function of the 
#' posterior is equal to the values specified in the argument ci.level (usually 0.025 and 0.975). The 
#' search is performed using R's optimize function.
#'
#' For each confidence limit, this function performs a binary search to evaluate the combined posterior, 
#' which is obtained by first transforming the imputed-data likelihood profiles into cumulative distribution functions (CDFs), and then averaging the CDFs to obtain the CDF of the posterior. Usually, 
#' the binary search manages to find the confidence intervals very quickly. The number of iterations 
#' (mean and maximum) will be supplied in the output object. Further details on the method can be 
#' found in Heinze, Ploner and Beyea (2013).
#'
#' @param obj Either a list of logistf fits (on multiple imputed data sets), or the result of analysis of a \code{mice} (multiply imputed) object using \code{with.mids}
#' @param variable The variable of interest, for which confidence intervals should be computed. If missing, confidence intervals for all variables will be computed.
#' @param data A list of data set corresponding to the model fits. Can be left blank if obj was obtained with the \code{dataout=TRUE} option or if obj was obtained by mice
#' @param firth If \code{TRUE}, applies the Firth correction. Should correspond to the entry in obj.
#' @param weightvar An optional weighting variable for each observation. 
#' @param control Control parameters for \code{logistf}, usually obtained by \code{logistf.control()}
#' @param ci.level The two confidence levels for each tail of the posterior distribution.
#' @param pvalue If \code{TRUE}, will also compute a P-value from the posterior.
#' @param offset An optional offset variable
#' @param bound.lo Bounds (vector of length 2) for the lower limit. Can be left blank. Use only if problems are encountered.
#' @param bound.up Bounds (vector of length 2) for the upper limit. Can be left blank. Use only if problems are encountered.
#' @param legacy If \code{TRUE}, will use pure R code for all model fitting. Can be slow. Not recommended.
#'
#' @return An object of class \code{CLIP.confint}, with items:
#'    \item{variable}{The variable(s) which were analyzed}
#'    \item{estimate}{The pooled estimate (average over imputations)}
#'    \item{ci}{The confidence interval(s)}
#'    \item{pvalue}{The p-value(s)}
#'    \item{imputations}{The number of imputed datasets}
#'    \item{ci.level}{The confidence level (input)}
#'    \item{bound.lo}{The bounds used for finding the lower confidence limit; usually not of interest. May be useful for error-tracing.}
#'    \item{bound.up}{The bounds used for finding the upper confidence limit}
#'    \item{iter}{The number of iterations (for each variable and each tail)}
#'    \item{call}{The call object}
#' 
#' @export CLIP.confint
#'
#' @examples
#' #generate data set with NAs 
#' freq=c(5,2,2,7,5,4)
#' y<-c(rep(1,freq[1]+freq[2]), rep(0,freq[3]+freq[4]), rep(1,freq[5]), rep(0,freq[6]))
#' x<-c(rep(1,freq[1]), rep(0,freq[2]), rep(1,freq[3]), rep(0,freq[4]), rep(NA,freq[5]),
#' rep(NA,freq[6]))
#' toy<-data.frame(x=x,y=y)
#' 
#' # impute data set 5 times
#' set.seed(169)
#' toymi<-list(0)
#' for(i in 1:5){
#'   toymi[[i]]<-toy
#'   y1<-toymi[[i]]$y==1 & is.na(toymi[[i]]$x)
#'   y0<-toymi[[i]]$y==0 & is.na(toymi[[i]]$x) 
#'   xnew1<-rbinom(sum(y1),1,freq[1]/(freq[1]+freq[2]))
#'   xnew0<-rbinom(sum(y0),1,freq[3]/(freq[3]+freq[4]))
#'   toymi[[i]]$x[y1==TRUE]<-xnew1
#'   toymi[[i]]$x[y0==TRUE]<-xnew0
#'   }
#'   
#'   # logistf analyses of each imputed data set
#'   fit.list<-lapply(1:5, function(X) logistf(data=toymi[[X]], y~x, pl=TRUE, dataout=TRUE))
#'   
#'   # CLIP confidence limits
#'   CLIP.confint(obj=fit.list)
#' @author Georg Heinze and Meinhard Ploner
#' @references Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining 
#' profile likelihood information from logistic regressions. Statistics in Medicine, to appear.
#' 
#' @seealso [logistf()] for Firth's bias-Reduced penalized-likelihood logistic regression.
#' @encoding UTF-8
#' 
#' @rdname CLIP.confint
CLIP.confint <- function(obj=NULL, variable=NULL, data, firth=TRUE, weightvar=NULL, 
                         control=logistf.control(), ci.level=c(0.025,0.975), pvalue=TRUE, 
                         offset=NULL, bound.lo=NULL, bound.up=NULL, legacy=FALSE) {
  which <- NULL
  if(length(ci.level)==1) ci.level<-c((1-ci.level)/2, 1-(1-ci.level)/2)
  
  
  if (!is.null(obj)) {
    if (is.mira(obj)) {
      # assuming a mira object with fits on each imputed data set
      #          data.imp<-texteval(paste("",obj$call$data,"",sep=""))
      data.imp<-eval(obj$call$data)    #############GEHT NOCH NICHT!!
      nimp<-data.imp$m          
      data<-lapply(1:nimp, function(x) complete(data.imp, action=x))
      fits<-obj$analyses
      formula<-as.formula(fits[[1]]$call$formula)
    }
    else {
      # assuming as input a list of data sets (data) and a list of fits (obj)
      fits<-obj
      if(missing(data)) if(is.null(fits[[1]]$data)) stop("Please provide data either as list of imputed data sets or by calling logistf on the imputed data sets with dataout=TRUE.\n")
      else data<-lapply(1:length(fits), function(X) fits[[X]]$data)
      formula<-as.formula(fits[[1]]$call$formula)
      nimp<-length(data)
    } 
  }
  else  {
    nimp<-length(data)
    stop(paste("Please provide a list of fit objects or a mira object with fits on the ",nimp," imputed data sets.\n"))
  }     
  
  if(is.null(variable))       variable<-names(fits[[1]]$coefficients)
  nvar<-length(variable)    
  
  inner.CLIP<-function(myvar){
    variable<-myvar
    old<-legacy
    
    
    
    
    imputations<-nimp
    if (is.null(which) & is.null(variable))
      stop("You must specify which (a one-sided formula) or variable (the name of the variable).")
    res<-matrix(0,length(grid),3)
    
    nperimp<-unlist(lapply(1:imputations,function(x) nrow(data[[x]])))
    
    variable.names<-colnames(model.matrix(formula, data = data[[1]]))
    #    mat.data<-matrix(lapply(1:imputations,function(x) matrix(unlist(data[[x]]),nperimp[[x]],ncol(data[[x]]))),sum(nperimp),ncol(data[[1]]))
    
    
    imputation.indicator<-unlist(sapply(1:imputations, function(x) rep(x,nperimp[x]))[TRUE])       #[TRUE]makes vector out of matrix
    
    mat.data<-matrix(0,sum(nperimp),ncol(data[[1]]))   ### copy list of data set into a matrix
    for(i in 1:imputations) mat.data[imputation.indicator==i,1:ncol(data[[i]])]<-as.matrix(data[[i]])
    
    
    #  if(missing(weightvar)) {
    #     weightvar<-"weights"
    #     mat.data[,ncol(mat.data)]<-rep(1,nrow(mat.data))
    #   }
    big.data<-data.frame(mat.data)
    colnames(big.data)<-colnames(data[[1]])
    
    
    k<-length(variable.names)  
    xyw<-matrix(0,sum(nperimp),k+2)
    xyw[,1:k]<-  model.matrix(formula, data = big.data)
    xyw[,k+1]<- model.response(model.frame(as.formula(formula),data=big.data))
    if (!is.null(which)) {
      cov.name <- variable.names
      cov.name2 <- colnames(model.matrix(which, data = big.data))
      pos <- match(cov.name2, cov.name)
      variable <- cov.name2
    }
    else {
      cov.name <- variable.names 
      cov.name2 <- variable
      pos <- match(cov.name2, cov.name)
    }
    #    fit <- logistf.fit(x, y, weight = weight, offset = offset,
    #        firth = firth, control = control)
    #    std.pos <- diag(fit$var)[pos]^0.5
    #    coefs <- fit$beta
    #    covs <- fit$var
    #    n <- nrow(x)
    #    cov.name <- labels(x)[[2]]
    
    if (is.null(weightvar)) {
      # data<-lapply(1:imputations, function(z) {
      #  data[[z]]$weightvar<-rep(1,nrow(data[[z]]))
      #  data[[z]]
      #  })
      #  weightvar<-"weightvar"
      xyw[,k+2]<-rep(1,nrow(xyw))
    }
    else xyw[,k+2]<-big.data[,weightvar]
    posweight<-k+2
    
    
    lf<-function(index) logistf.fit(y=xyw[index,k+1], x=xyw[index,1:k], weight=xyw[index,k+2])
    
    if (is.null(fits))   fits<-lapply(1:imputations, function(z) lf(imputation.indicator==z))
    
    if(is.null(bound.lo)){
      lower.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.lower[pos])))
      lowerbound.lo<-min(lower.collect) ###von dem den index nehmen und davon das PL CI ausrechnen
      upperbound.lo<-max(lower.collect)
    }  else {
      lowerbound.lo<-min(bound.lo)
      upperbound.lo<-max(bound.lo)
    }
    
    if(is.null(bound.up)){
      upper.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.upper[pos])))
      lowerbound.up<-min(upper.collect) ###von dem den index nehmen und davon das PL CI ausrechnen
      upperbound.up<-max(upper.collect)
    }  else {
      lowerbound.up<-min(bound.up)
      upperbound.up<-max(bound.up)
    }
    
    if (lowerbound.lo == upperbound.lo){
      lowerbound.lo <- lowerbound.lo - 1/2   # very crude!
      upperbound.lo <- upperbound.lo + 1/2
    }
    
    if (lowerbound.up == upperbound.up){
      lowerbound.up <- lowerbound.up - 1/2   # very crude!
      upperbound.up <- upperbound.up + 1/2
    }
    
    
    estimate<-mean(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients[pos])))
    
    iter<-numeric(0)
    
    loglik<-unlist(lapply(1:imputations, function(x) fits[[x]]$loglik[2]))
    beta<-t(matrix(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients)),k,imputations))
    
    
    lpdf<-function(zz,z) logistf.pdf(x=xyw[imputation.indicator==zz,1:k], y=xyw[imputation.indicator==zz,k+1], 
                                     weight=xyw[imputation.indicator==zz,k+2], beta=beta[zz,],loglik=loglik[zz],
                                     pos=pos, firth=firth, offset=offset, control=control, b=z, old=old)$pdf
    
    f=function(z)  mean(unlist(lapply(1:imputations, function(zz) lpdf(zz,z))))
    
    f.lower<-f(lowerbound.lo)-ci.level[1]
    f.upper<-f(upperbound.lo)-ci.level[1]
    iter[1]<-2
    itwhile<-0
    while(f.lower > 0 & (upperbound.lo - lowerbound.lo) > 0 & itwhile<5) {
      itwhile<-itwhile+1
      lowerbound.lo<-lowerbound.lo - (upperbound.lo - lowerbound.lo)/2
      f.lower<-f(lowerbound.lo)-ci.level[1]
      iter[1]<-iter[1]+1
    }
    
    if (itwhile>=5 & f.upper<0) stop("pool.pl can not find a lower boundary for the lower confidence limit.\n Try to increase number of imputations or supply boundaries by bound.lo=c(x,xx).\n")
    
    itwhile<-0
    while(f.upper < 0 & (upperbound.lo - lowerbound.lo) > 0 & itwhile<5) {
      itwhile<-itwhile+1
      upperbound.lo<-upperbound.lo + (upperbound.lo - lowerbound.lo)/2
      f.upper<-f(upperbound.lo)-ci.level[1]
      iter[1]<-iter[1]+1
    }
    if (itwhile>=5 & f.upper<0) stop("pool.pl can not find an upper boundary for the lower confidence limit.\n Try to increase number of imputations or supply boundaries by bound.lo=c(x,xx).\n")
    
    ci<-numeric(0)
    res.ci<-uniroot(f=function(z) {f(z)-ci.level[1]}, 
                    lower=lowerbound.lo, upper=upperbound.lo, f.lower=f.lower, f.upper=f.upper)
    ci[1]<-res.ci$root
    iter[1]<-res.ci$iter+iter[1]
    
    
    
    f.lower<-f(lowerbound.up)-ci.level[2]
    f.upper<-f(upperbound.up)-ci.level[2]
    iter[2]<-2
    itwhile<-0
    while(f.lower > 0 & (upperbound.up - lowerbound.up) > 0 & itwhile<5) {
      itwhile<-itwhile+1
      lowerbound.up<-lowerbound.up - (upperbound.up - lowerbound.up)/2
      f.lower<-f(lowerbound.up)-ci.level[2]
      iter[2]<-iter[2]+1
    }
    
    if (itwhile>=5 & f.lower>0) stop("pool.pl can not find a lower boundary for the upper confidence limit.\n Try to increase number of imputations or supply boundaries by bound.up=c(x,xx).\n")
    
    
    itwhile<-0
    while(f.upper < 0 & (upperbound.up - lowerbound.up) > 0 & itwhile<5) {
      itwhile<-itwhile+1
      upperbound.up<-upperbound.up + (upperbound.up - lowerbound.up)/2
      f.upper<-f(upperbound.up)-ci.level[2]
      iter[2]<-iter[2]+1
    }
    
    if (itwhile>=5 & f.upper<0) stop("pool.pl can not find an upper boundary for the upper confidence limit.\n Try to increase number of imputations or supply boundaries by bound.up=c(x,xx).\n")
    
    res.ci<-uniroot(f=function(z) {f(z)-ci.level[2]}, 
                    lower=lowerbound.up, upper=upperbound.up, f.lower=f.lower, f.upper=f.upper)
    ci[2]<-res.ci$root
    iter[2]<-res.ci$iter+iter[2]
    
    if (pvalue=="TRUE") {
      pvalue1<-f(0)
      pvalue<-2*min(pvalue1, (1-pvalue1))   # two-sided = twice smaller tail
    }
    else pvalue<-NA
    #uniroot(f=function(x) pdf.imp(b=x, formula = formula, data = data,
    #    which=which, firth = firth, weights=weights, control=control, plcontrol=plcontrol)-ci.level[2], 
    #    lower=lowerbound.up, upper=upperbound.up)
    
    #ci[2]<-uniroot(f=function(x) pdf.combi(b=x, formula = formula, data = data,
    #    which=which, firth = firth, weights=weights, control=control, plcontrol=plcontrol)-ci.level[2], lower=lowerbound.up, upper=upperbound.up,
    #     formula = formula, data = data,
    #    which=which, firth = TRUE, weights=weights, control=control, plcontrol=plcontrol)
    
    res<-list(
      estimate=estimate,ci=ci,pvalue=pvalue, imputations=imputations, 
      ci.level=ci.level, myvar=myvar, call=match.call(), 
      bound.lo=c(lowerbound.lo, upperbound.lo),
      bound.up=c(lowerbound.up, upperbound.up), iter=iter)
    return(res)
  }
  
  estimate<-numeric(nvar)
  ci<-matrix(0,nvar,2)
  pvalue.out<-numeric(nvar)
  bound.lo<-matrix(0,nvar,2)
  bound.up<-matrix(0,nvar,2)
  iter<-matrix(0,nvar,2)
  for(i in 1:nvar)  {
    res.tmp<-inner.CLIP(myvar=variable[i])
    estimate[i]<-res.tmp$estimate
    ci[i,]<-res.tmp$ci
    pvalue.out[i]<-res.tmp$pvalue
    bound.lo[i,]<-res.tmp$bound.lo
    bound.up[i,]<-res.tmp$bound.up
    iter[i,]<-res.tmp$iter
  }
  
  res<-list(
    variable=variable, estimate=estimate, ci=ci, pvalue=pvalue.out,
    imputations=res.tmp$imputations, ci.level=ci.level, bound.lo=bound.lo, 
    bound.up=bound.up, iter=iter, call=match.call())
  
  attr(res,"class") <- "CLIP.confint"
  return(res)
}

#' @method print CLIP.confint
print.CLIP.confint <- function(
  x,
  exp=FALSE, 
  ...
) {
  object <- x
  mat <- cbind(object$estimate, object$ci[,1], object$ci[,2], object$pvalue)
  if(exp) 
    mat[,1:3]<-exp(mat[,1:3])
  colnames(mat)<-c("Estimate", "Lower", "Upper", "P-value")
  
  if(exp) 
    colnames(mat)[1]<-"Odds ratio"
  rownames(mat)<-object$variable
  
  print(object$call)
  cat("Number of imputations: ",object$imputations,"\n\n")
  cat("Iterations, mean: ", mean(object$iter), " max:", max(object$iter),"\n\n")
  cat("Confidence level, lower:", object$ci.level[1]*100,"%, upper:",object$ci.level[2]*100,"%\n")
  print(mat)
}

is.mira <- function(
  object
) {
  if(is.null(attr(object,"class"))) 
    return(FALSE)
  if(attr(object,"class")=="mira") 
    return(TRUE)
  return(FALSE)
}
