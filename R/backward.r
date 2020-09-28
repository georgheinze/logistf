#' Backward Elimination/Forward Selection of Model Terms in logistf Models
#'
#' These functions provide simple backward elimination/forward selection procedures for logistf models.
#' 
#' The variable selection is simply performed by repeatedly calling add1 or drop1 methods for logistf, 
#' and is based on penalized likelihood ratio test. It can also properly handle variables that were 
#' defined as factors in the original data set.
#' 
#' @param object A fitted logistf model object. To start with an empty model, create a model fit 
#' with a formula=<y>~1, pl=FALSE. (Replace <y> by your response variable.)
#' @param scope The scope of variables to add/drop from the model. Can be missing for backward, backward will use 
#' the terms of the object fit. Alternatively, an arbitrary vector of variable names can be given, to allow 
#' that only some of the variables will be competitively selected or dropped. Has to be provided for forward.
#' @param steps The number of forward selection/backward elimination steps.
#' @param slstay For \code{backward}, the significance level to stay in the model.
#' @param slentry For \code{forward}, the significance level to enter the model.
#' @param trace If \code{TRUE}, protocols selection steps.
#' @param printwork If \code{TRUE}, prints each working model that is visited by the selection procedure.
#' @param full.penalty If \code{TRUE} penalty is not taken from current model but from start model.
#' @param pl For forward, computes profile likelihood confidence intervals for the final model if \code{TRUE}.
#' @param ... Further arguments to be passed to methods.
#'
#' @return An updated \code{logistf, flic} or \code{flac} fit with the finally selected model.
#'
#' @examples 
#' data(sex2) 
#' fit<-logistf(data=sex2, case~1, pl=FALSE) 
#' fitf<-forward(fit, scope = c("dia", "age")) 
#' 
#' fit2<-logistf(data=sex2, case~age+oc+vic+vicl+vis+dia) 
#' fitb<-backward(fit2)
#' 
#' @importFrom utils tail
#' 
#' @rdname backward
#' @export backward
backward <- function(object,...){
  UseMethod("backward",object)
}
#' @exportS3Method backward logistf
#' @method backward logistf
#' @rdname backward
backward.logistf <- function(object, scope, steps=1000, slstay=0.05, trace=TRUE, printwork=FALSE,full.penalty=FALSE, ...){
  istep<-0 #index of steps
  
  mf <- match.call(expand.dots =FALSE)
  m <- match("object", names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  variables <- object$terms[-1]
  
  terms.fit <- object$call$terms.fit
  if(!is.null(terms.fit)) stop("Please call backward on a logistf-object with all terms fitted.")

  working<-object
  if(trace){
    cat("Step ", istep, ": starting model\n")
    if(printwork) {
      print(working)
      cat("\n\n")
    }
  }
  if(missing(scope)) scope<-variables #scope missing - use terms of object fit
  if (full.penalty) { #if TRUE, use start model and save all removals in vector removal
    removal <- vector()
  }
  while(istep<steps & working$df>1){
    if(full.penalty && istep!=0){ #check with istep!=0 if removal is an empty vector - not possible to pass such to drop1
      mat <- drop1(object, full.penalty.vec=removal)
    }
    else {
      mat<-drop1(working)
    }
    istep<-istep+1
    if(all(mat[,3]<slstay)) {
      break
    } #check p-value of variables in scope
    inscope<-match(scope,rownames(mat))
    inscope<-inscope[!is.na(inscope)]
    if (full.penalty){
      removal <- c(removal, rownames(mat)[mat[,3]==max(mat[inscope,3])])
      curr_removal <- removal[istep] #if two pvalues are the same, both are taken 
    }
    else { #if full.penalty = FALSE: save only current removal
      removal<-rownames(mat)[mat[,3]==max(mat[inscope,3])] #remove highest pvalue - if p-values are the same fot two variables: delete both
      curr_removal <- removal
    }
    #check if object$formula contains a dot shortcut i.e. last character: 
    if (sapply(grepl("\\.", working$formula), tail, 1)){
      char <- paste(working$terms[-1], collapse=" + ") #get variables from data without intercept
      variables <- paste(char, paste(removal, collapse="-"), sep="-") #remove target variable
      newform <- as.formula(paste("", variables, sep="~")) #coerce to formula
    }
    else {
        newform=as.formula(paste("~.-",paste(curr_removal, collapse = "-")))
    }
    if(!full.penalty){ #update working only if full.penalty==FALSE
      if(working$df==2 | working$df==mat[mat[,3]==max(mat[,3]),2]){
        working<-update(working, formula=newform, pl=FALSE, evaluate = FALSE)
        working <- eval.parent(working)
      }
      else {
        working<-update(working, formula. =newform, evaluate = FALSE)
        working <- eval.parent(working)
        
      }
    }
    else {
      working$df <- working$df-1 #update only degrees of freedom in case of full.penalty 
    }
    if(trace){
      cat("Step ", istep, ": removed ", curr_removal, " (P=", max(mat[,3]),")\n")
      if(printwork) {
        print(working)
        cat("\n\n")
      }
    }
  }
  if(trace) cat("\n")
  if(full.penalty){
    tmp <- match(removal, variables)
    tofit <- object$terms[-(tmp+1)]
    working<-update(working,terms.fit=tofit, evaluate = FALSE)
    working <- eval.parent(working)
  }
  return(working)
}

#' @exportS3Method backward flic
#' @method backward flic
#' @rdname backward
backward.flic<-function(object, scope, steps=1000, slstay=0.05, trace=TRUE, printwork=FALSE,full.penalty=FALSE, ...){
   message("It is intended to call backward() on a logistf-object and afterwards to call flic() on the reduced model.")
}

#' @describeIn backward Forward Selection 
#' @export forward
forward <- function(object,...){
  UseMethod("forward",object)
}
#' @exportS3Method forward logistf
#' @method forward logistf
#' @rdname backward
forward.logistf<-function(object, scope, steps=1000, slentry=0.05, trace=TRUE, printwork=FALSE, pl=TRUE, ...){
  istep<-0
  
  mf <- match.call(expand.dots =FALSE)
  m <- match("object", names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  variables <- object$terms[-1]
  
  terms.fit <- object$call$terms.fit
  if(!is.null(terms.fit)) stop("Please call forward on a logistf-object with all terms fitted.")
  
  working<-object
  if(missing(scope)) {
    stop("Please provide scope (vector of variable names).\n")
  }
  else if(is.numeric(scope)) {
    scope<-attr(terms(object),"term.labels")[scope]
  }
  else {
    if(!sum(is.na(match(variables, scope)))==length(variables)) scope<-scope[-(match(variables,scope))]
  }

  if(trace){
    cat("Step ", istep, ": starting model\n")
    if(printwork) {
        print(working)
        cat("\n\n")
        }
  }

  inscope<-scope
  while(istep<steps & length(inscope)>=1){
    istep<-istep+1
    mat<-add1(working, scope=inscope)
    if(all(mat[,3]>slentry)) break
    index<-(1:nrow(mat))[mat[,3]==min(mat[,3])]
    if(length(index)>1) index<-index[mat[index,1]==max(mat[index,1])]
    addvar<-rownames(mat)[index]
    newform=as.formula(paste("~.+",addvar))
    working<-update(working, formula=newform, pl=FALSE)
    newindex<-is.na(match(inscope, attr(terms(object),"term.labels")))
    if(all(is.na(newindex)==TRUE)) break
    #inscope<-inscope[-index]
    inscope <- inscope[-match(rownames(mat)[index],inscope)]
    if(trace){
      cat("Step ", istep, ": added ", addvar, " (P=", mat[addvar,3],")\n")
      if(printwork) {
        print(working)
        cat("\n\n")
      }
    }
  }
   if(pl) working<-update(working, pl=TRUE)
   if(trace) cat("\n")
   return(working)
}


#' @exportS3Method backward flac
#' @method backward flac
#' @rdname backward
backward.flac<-function(object, steps=1000, slstay=0.05, trace=TRUE, printwork=FALSE,full.penalty=FALSE,...){
  stop("Currently not implemented.")
  istep<-0 #index of steps
  mf <- match.call(expand.dots =FALSE)
  m <- match("object", names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  variables <- object$terms[-1]
  
  terms.fit <- object$call$terms.fit
  if(!is.null(terms.fit)) stop("Please call backward on a logistf-object with all terms fitted.")
  
  working<-object
  if(trace){
    cat("Step ", istep, ": starting model\n")
    if(printwork){
      print(working)
      cat("\n\n")
    }
  }
  scope<-attr(terms(working$formula),"term.labels") #scope missing - use terms of object fit
  if(full.penalty) { #if TRUE, use start model and save all removals in vector removal
    removal <- vector()
  }
  while(istep<steps & working$df>=1){
    if(full.penalty && istep!=0){ #check with istep!=0 if removal is an empty vector - not possible to pass such to drop1
      mat <- drop1(object, full.penalty.vec=removal)
    }
    else {
      mat<-drop1(working)
    }
    istep<-istep+1
    if(all(mat[,3]<slstay)) {
      break
    } #check p-value of variables in scope
    inscope<-match(scope,rownames(mat))
    inscope<-inscope[!is.na(inscope)]
    if (full.penalty){
      removal <- c(removal, rownames(mat)[mat[,3]==max(mat[inscope,3])])
      curr_removal <- removal[istep] #if two pvalues are the same, both are taken 
    }
    else { #if full.penalty = FALSE: save only current removal
      removal<-rownames(mat)[mat[,3]==max(mat[inscope,3])] #remove highest pvalue
      curr_removal <- removal
    }
    #check if object$formula contains a dot shortcut i.e. last character: 
    if (sapply(grepl("\\.", working$formula), tail, 1)){
      char <- paste(working$terms[-1], collapse=" + ") #get variables from data without intercept
      variables <- paste(char, paste(removal, collapse="-"), sep="-") #remove target variable
      newform <- as.formula(paste("", variables, sep="~")) #coerce to formula
    }
    else {
      newform=as.formula(paste("~.-",paste(curr_removal, collapse = "-")))
    }
    if(!full.penalty){ #udate working only if full.penalty==FALSE
      if(working$df==1 | working$df==mat[mat[,3]==max(mat[,3]),2]){
        working<-update(working, formula=newform, pl=FALSE)
      }
      else {
        working<-update(working, formula=newform)
      }
    }
    if(trace){
      cat("Step ", istep, ": removed ", curr_removal, " (P=", max(mat[,3]),")\n")
      if(printwork) {
        print(working)
        cat("\n\n")
      }
    }
  }
  if(trace) {
    cat("\n")
  }
  if(full.penalty) {
    tmp <- match(removal, variables)
    tofit <- object$terms[-(tmp+1)]
    working<-update(working,formula=working$formula, terms.fit=tofit)
  }
  return(working)
}