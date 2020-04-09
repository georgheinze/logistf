#' Pseudo Variance Modification of Rubin's Rule
#' 
#' The pseudo-variance modification proposed by Heinze, Ploner and Beyea (2013) provides a quick 
#' way to adapt Rubin‘s rules to situations of a non-normal distribution of a regression coefficient. 
#' However, the approxiation is less accurate than that of the CLIP method.
#' 
#' The pseudo-variance modification computes a lower and an upper pseudo-variance, which are based 
#' on the distance between profile likelihood limits and the parameter estimates. These are then 
#' plugged into the usual Rubin‘s rules method of variance combination
#'
#' @param obj A fitted \code{logisf} object
#' @param variable The variable(s) to compute the PVR confidence intervals, either provided as names or as numbers
#' @param skewbeta If \code{TRUE}, incorporates information on the skewness of the parameter estimates 
#' across the imputed data sets.
#'
#' @return An object of class \code{PVR.confint} with items:
#'   \item{estimate}{the pooled parameter estimate(s) (the average across completed-data estimates)} 
#'   \item{ci}{the confidence intervals based on the PVR method} 
#'   \item{lower.var}{the lower pseudo-variance(s)} 
#'   \item{upper.var}{the upper pseudo-variance(s)} 
#'   \item{conflev}{the confidence level: this is determined by the confidence level (1-alpha) used in the input fit objects} 
#'   \item{call}{the function call} 
#'   \item{variable}{the variable(s) for which confidence intervals were computed} 
#'   
#' @author Georg Heinze
#' @references Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining 
#' profile likelihood information from logistic regressions. Statistics in Medicine, to appear.   
#' @export
#'
#' @examples
#' #generate data set with NAs
#' freq=c(5,2,2,7,5,4)
#' y<-c(rep(1,freq[1]+freq[2]), rep(0,freq[3]+freq[4]), rep(1,freq[5]), rep(0,freq[6]))
#' x<-c(rep(1,freq[1]), rep(0,freq[2]), rep(1,freq[3]), rep(0,freq[4]), rep(NA,freq[5]),
#'    rep(NA,freq[6]))
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
#'# logistf analyses of each imputed data set
#'fit.list<-lapply(1:5, function(X) logistf(data=toymi[[X]], y~x, pl=TRUE, dataout=TRUE))
#'
#'# CLIP confidence limits
#'PVR.confint(obj=fit.list)
#'
#' @rdname PVR.confint
PVR.confint<-function(obj, variable, skewbeta=FALSE){
   if(missing(obj)) stop("Please provide an object with a list of logistf fits for analysis.\n")
   fit.list<-obj
   
   if(missing(variable)) variable=names(fit.list[[1]]$coefficients)
   nimp<-length(fit.list)
   nvar<-length(variable)
   zalpha<-qnorm(1-fit.list[[1]]$alpha/2)
   pseudo.var<-function(x,tail){
      meanx<-mean(x)
      if(tail==1) signx<--1
      else signx<-1
      diff<-x-meanx
      if(skewbeta) thistail<-(sign(diff)==signx) 
      else thistail<-(diff==diff)
      pseudovar<-sum(diff[thistail]^2)/(sum(thistail)-1)
      return(pseudovar)
    }
   if(nvar>1){
     lower.pse<-(t(sapply(1:nimp,function(X) (fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.lower[variable])/zalpha)))
    upper.pse<-(t(sapply(1:nimp,function(X) -(fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.upper[variable])/zalpha)))
    lower.var<-lower.pse^2
    upper.var<-upper.pse^2
    coef<-t(sapply(1:nimp,function(X) fit.list[[X]]$coefficient[variable]))
    var.coef.lo<-apply(coef,2,pseudo.var,1)
    var.coef.up<-apply(coef,2,pseudo.var,2)
    mean.coef<-apply(coef,2,mean)
    RR.lower.var<-apply(lower.var,2,mean)+(1+1/nimp)*var.coef.lo
    RR.lower.ci<-mean.coef-zalpha*sqrt(RR.lower.var)
    RR.upper.var<-apply(upper.var,2,mean)+(1+1/nimp)*var.coef.up
    RR.upper.ci<-mean.coef+zalpha*sqrt(RR.upper.var)
    ret<-cbind(RR.lower.ci, RR.upper.ci)
    rownames(ret)<-variable
    colnames(ret)<-c("lower","upper")
   } else {
    lower.pse<-((sapply(1:nimp,function(X) (fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.lower[variable])/zalpha)))
    upper.pse<-((sapply(1:nimp,function(X) -(fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.upper[variable])/zalpha)))
    lower.var<-lower.pse^2
    upper.var<-upper.pse^2
    coef<-(sapply(1:nimp,function(X) fit.list[[X]]$coefficient[variable]))
    var.coef.lo<-pseudo.var(coef,1)
    var.coef.up<-pseudo.var(coef,2)
    mean.coef<-mean(coef)
    RR.lower.var<-mean(lower.var)+(1+1/nimp)*var.coef.lo
    RR.lower.ci<-mean.coef-zalpha*sqrt(RR.lower.var)
    RR.upper.var<-mean(upper.var)+(1+1/nimp)*var.coef.up
    RR.upper.ci<-mean.coef+zalpha*sqrt(RR.upper.var)
    ret<-cbind(RR.lower.ci, RR.upper.ci)
    colnames(ret)<-c("lower","upper")
    rownames(ret)<-variable
   }
 res<-list(estimate=mean.coef, ci=ret,lower.var=RR.lower.var, upper.var=RR.upper.var,conflev=1-fit.list[[1]]$alpha, call=match.call(), variable=variable)
 attr(res,"class")="PVR.confint"
 return(res)
 }


print.PVR.confint<-function(x,exp=FALSE,...){
  object<-x
  print(object$call)
  cat("Pseudo-variance modification of Rubins Rules\n")
  cat("Confidence level: ", object$conflev*100, "%\n")
  mat<-cbind(object$estimate, object$ci[,1], object$ci[,2], object$lower.var, object$upper.var)
  colnames(mat)<-c("Estimate", "Lower", "Upper", "Lower pseudo variance", "Upper pseudo variance")
  if(exp) {
    mat[,1:3]<-exp(mat[,1:3])
    colnames(mat)[1]<-"Odds ratio"
    }
  rownames(mat)<-object$variable
  print(mat)
 }
 

