#' Combine Profile Likelihoods from Imputed-Data Model Fits
#' 
#' This function uses CLIP (combination of likelihood profiles) 
#' to compute the pooled profile of the posterior after multiple imputation.
#'
#' While CLIP.confint iterates to find those values at which the CDF of the 
#' pooled posterior equals the confidence levels, CLIP.profile will evaluate 
#' the whole profile, which enables plotting and evaluating the skewness of the combined and the completed-data profiles. The combined and completeddata profiles are available as cumulative distribution function (CDF) or in the scaling of relative 
#' profile likelihood (minus twice the likelihood ratio statistic compared to the maximum). Using a 
#' plot method, the pooled posterior can also be displayed as a density. 
#'
#' @param obj Either a list of logistf fits (on multiple imputed data sets), or the result 
#' of analysis of a \code{mice} (multiply imputed) object using \code{with.mids}.
#' @param variable The variable of interest, for which confidence intervals should be computed. 
#' If missing, confidence intervals for all variables will be computed.
#' @param data A list of data set corresponding to the model fits. Can be left blank if obj was 
#' obtained with the dataout=TRUE option or if obj was obtained by mice.
#' @param which Alternatively to variable, the argument which allows to specify the variable to 
#' compute the profile for as righthand formula, e.g. which=~X.
#' @param firth If \code{TRUE}, applies the Firth correction. Should correspond to the entry in obj.
#' @param weightvar An optional weighting variable for each observation 
#' @param control control parameters for \code{logistf}, usually obtained by \code{logistf.control()}
#' @param offset An optional offset variable
#' @param from Lowest value for the sequence of values for the regression coefficients for which the profile will be computed. Can be left blank.
#' @param to Highest value for the sequence of values for the regression coefficients for which the profile will be computed. Can be left blank
#' @param steps Number of steps for the sequence of values for the regression coefficients for which the profile will be computed
#' @param legacy If \code{TRUE}, only R code will be used. Should be avoided.
#' @param keep If \code{TRUE}, keeps the profiles for each imputed data sets in the output object.
#'
#' @return An object of class \code{CLIP.profile} with items:
#'    \item{beta}{The values of the regression coefficient}
#'    \item{cdf}{The cumulative distribution function of the posterior}
#'    \item{profile}{The profile of the posterior}
#'    \item{cdf.matrix}{An imputations x steps matrix with the values of the completed-data CDFs for each beta}
#'    \item{profile.matrix}{An imputations x steps matrix with the values of the completed-data profiles for each beta}
#'    \item{call}{The function call}
#' @export
#' @author Georg Heinze und Meinhard Plonar
#' @references Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining profile 
#' likelihood information from logistic regressions. Statistics in Medicine, to appear.
#' @examples
#' 
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
#' }
#' 
#' # logistf analyses of each imputed data set
#' fit.list<-lapply(1:5, function(X) logistf(data=toymi[[X]], y~x, pl=TRUE))
#' 
#' # CLIP profile
#' xprof<-CLIP.profile(obj=fit.list, variable="x",data =toymi, keep=TRUE)
#' plot(xprof)
#' 
#' #plot as CDF
#' plot(xprof, "cdf")
#' 
#' #plot as density
#' plot(xprof, "density")
#' 
#' @rdname CLIP.profile
CLIP.profile <- function(obj=NULL, variable, data, which, firth=TRUE, weightvar, 
  control=logistf.control(), offset=NULL, from=NULL, to=NULL, steps=101, legacy=FALSE, keep=FALSE) {
  old <- legacy
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
      if(missing(data)) stop("Please provide data either as list of imputed data sets or by calling logistf on the imputed data sets with dataout=TRUE.\n")
      formula<-as.formula(fits[[1]]$call$formula)
      nimp<-length(data)
    } 
  }
  else {
    nimp<-length(data)
    stop(paste("Please provide a list of fit objects or a mira object with fits on the ",nimp," imputed data sets.\n"))
  }     
  
  
  imputations<-nimp
  if (missing(which) & missing(variable))
    stop("You must specify which (a one-sided formula) or variable (the name of the variable).")
  if (missing(control))
    control <- logistf.control()
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
  if (!missing(which)) {
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
  
  if (missing(weightvar)) {
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
  
  if(is.null(from)){
    lower.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.lower[pos])))
    from<-min(lower.collect) ###von dem den index nehmen und davon das PL CI ausrechnen
  } 
  
  if(is.null(to)){
    upper.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.upper[pos])))
    to<-max(upper.collect)
  } 
  
  
  
  estimate<-mean(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients[pos])))
  
  iter<-numeric(0)
  
  loglik<-unlist(lapply(1:imputations, function(x) fits[[x]]$loglik['full']))
  beta<-t(matrix(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients)),k,imputations))
  
  
  lpdf<-function(zz,z) logistf.pdf(x=xyw[imputation.indicator==zz,1:k], y=xyw[imputation.indicator==zz,k+1], 
                                   weight=xyw[imputation.indicator==zz,k+2], beta=beta[zz,],loglik=loglik[zz],
                                   pos=pos, firth=firth, offset=offset, control=control, b=z, old=old)$pdf
  
  z_seq<-seq(from, to, (to-from)/steps)
  
  if(keep==FALSE){
    f=function(z)  mean(unlist(lapply(1:imputations, function(zz) lpdf(zz,z))))
    ### lasse z laufen von from nach to
    ### evaluiere f an allen z's und errechne daraus profile
    pdf_mat<-NULL
    profile.mat<-NULL
    pdf_seq<-unlist(lapply(z_seq, function(Z) f(Z)))
  } else {
    f=function(z)  unlist(lapply(1:imputations, function(zz) lpdf(zz,z)))
    #   pdf_mat<-matrix(0,steps+1,imputations)
    pdf_mat<-matrix(unlist(lapply(z_seq, function(Z) f(Z))), imputations, steps+1)   
    profile.mat<- -qnorm(pdf_mat)**2
    pdf_seq<-apply(pdf_mat,2,mean)
  }
  
  
  chisq_seq<- -qnorm(pdf_seq)**2
  res<-list(beta=z_seq, cdf=pdf_seq, profile=chisq_seq, cdf.matrix=pdf_mat, profile.matrix=profile.mat, call=match.call())
  attr(res,"class")<-c("logistf.CLIP.profile","logistf.profile")
  res
}

