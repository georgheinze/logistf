#' Add or Drop All Possible Single Terms to/from a \code{logistf} Model
#'
#'Compute all the single terms in the scope argument that can be added to or dropped 
#'from the model, fit those models and compute a table of the changes in fit.
#'
#'\code{drop1} and \code{add1} generate a table where for each variable the penalized 
#'likelihood ratio chi-squared, the degrees of freedom, and the p-value for dropping/adding this variable are given.
#'
#' @param object A fitted \code{logistf} object
#' @param scope The scope of variables considered for adding or dropping. Should be a 
#' vector of variable names. Can be left missing; the method will then use all variables 
#' in the object‘s data slot which are not identified as the response variable.
#' @param test The type of test statistic. Currently, only the PLR test (penalized likelihood 
#' ratio test) is allowed for logistf fits.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A matrix with \code{nvar} rows and 3 columns (Chisquared, degrees of freedom, p-value).
#' @export
#'
#' @examples
#' data(sex2) 
#' fit<-logistf(data=sex2, case~1, pl=FALSE) 
#' add1(fit)
#'  
#' fit2<-logistf(data=sex2, case~age+oc+dia+vic+vicl+vis) 
#' drop1(fit2)
#' 
#' @rdname add1
#' @method add1 logistf
#' @exportS3Method add1 logistf
add1.logistf<-function(object, scope, test="PLR", ...){
  if(missing(scope)) scope<-colnames(object$data)
  if(is.numeric(scope)) scope<-colnames(object$data)[scope]
  whereisresponse<-(match(colnames(model.frame(object$formula,data=object$data))[1],scope))
  if(!is.na(whereisresponse)) scope<-scope[-whereisresponse]
  scope<-scope[is.na(match(scope, attr(terms(object),"term.labels")))]
  variables<-scope
  
  nvar<-length(variables)
  mat<-matrix(0,nvar,3)
  for(i in 1:nvar){
    newform<-as.formula(paste("~.+",variables[i]))
    res<-anova(object, update(object,formula=newform))
    mat[i,1]<-res$chisq
    mat[i,2]<-res$df
    mat[i,3]<-res$pval
  }
  rownames(mat)<-variables
  colnames(mat)<-c("ChiSq","df","P-value")
  return(mat)
}

#' @method drop1 logistf
#' @rdname add1
#' @exportS3Method drop1 logistf
drop1.logistf<-function(object, scope, test="PLR", full.penalty.vec=NULL, ...){
  variables<-attr(terms(object$formula, data = object$data),"term.labels")
  if(!is.null(full.penalty.vec)){ #exclude already removed variables - see backward
    matched <- match(full.penalty.vec, variables)+1 #+1: to include intercept
    ind <- (1:7)[-matched] #for col.fit.object
    variables <- variables[-(matched-1)]
  }
  nvar<-length(variables)
  mat<-matrix(0,nvar,3) #initialise output
  for(i in 1:nvar){ #for every variable: omit from the object model fitin anova
    #full.penalty option of backward function
    if (!is.null(full.penalty.vec)){
      newform <- as.formula(paste("~", paste(variables[i], paste(full.penalty.vec,collapse="+"), sep="+")))
      res<-anova(object, formula=newform, method="nested", col.fit.object=ind)
    }
    else {
      newform<-as.formula(paste("~",variables[i]))
      res<-anova(object, formula=newform, method="nested")
    }
    mat[i,1]<-res$chisq
    mat[i,2]<-res$df
    mat[i,3]<-res$pval
  }
  rownames(mat)<-variables
  colnames(mat)<-c("ChiSq","df","P-value")
  return(mat)
}