#' Analysis of Penalized Deviance for \code{logistf} Models
#'
#'This method compares hierarchical and non-hierarchical logistf models using 
#'penalized likelhood ratio tests. It replaces the function logistftest of former 
#'versions of logistf.
#'
#'Comparing models fitted by penalized methods, one must consider that the penalized likelihoods 
#'are not directly comparable, since a penalty is involved. Or in other words, inserting zero for 
#'some regression coefficients will not lead to the same penalized likelihood as if the corresponding 
#'variables are simply "unknown" to a model. The anova method takes care that the same penalty is 
#'used for two hierarchically nested models, and if the models are not hierarchically nested, it will 
#'first relate each penalized likelihood to its null penalized likelihood, and only compare the resulting 
#'penalized likelihod ratio statistics. The chi-squared approximation for this latter method (PLR) is 
#'considered less accurate than that of the nested method. Nevertheless, it is the only way to go for 
#'comparison of non-nested models.
#'
#' @param object A fitted \code{logistf} model object
#' @param fit2 Another fitted \code{logistf} model object, to be compared with \code{object}
#' @param formula Alternatively to \code{fit2}, a formula which specifies terms to omit from the object model fit.
#' @param method One of c("nested","PLR"). nested is the default for hierarchically nested 
#' models, and will compare the penalized likelihood ratio statistics (minus twice
#' the difference between maximized penalized log likelihood and null penalized
#' log likelihood), where the null penalized log likelihood is computed from the 
#' same, hierarchically superior model. Note that unlike in maximum likelihood 
#' analysis, the null penalized likelihood depends on the penalty (Jeffreys prior) 
#' which itself depends on the scope of variables of the hierarchically superior 
#' model. PLR compares the difference in penalized likelihood ratio between the 
#' two models, where for each model the null penalized likelihood is computed 
#' within the scope of variables in that model. For PLR, the models need not be 
#' hierarchically nested.
#' @param ... Further arguments passed to the method.
#'
#' @return An object of class \code{anova.logistf} with items
#'   \item{chisq}{the chisquared statistic for the model comparison}
#'   \item{df}{The degrees of freedom}
#'   \item{pval}{The p-value}
#'   \item{call}{The function call}
#'   \item{method}{The method of comparison (input)}
#'   \item{model1}{The first model}
#'   \item{model2}{The second model which was compared to the first model}
#'   \item{PLR1}{The PLR statistic of the first model}
#'   \item{PLR2}{the PLR statistic of the second model; for the nested method, this will be the drop in chi-squared due to setting the coefficients to zero}
#' @export
#'
#' @examples
#' data(sex2) 
#' fit<-logistf(data=sex2, case~age+oc+dia+vic+vicl+vis)
#' 
#' #simultaneous test of variables vic, vicl, vis:
#' anova(fit, formula=~vic+vicl+vis)
#' 
#' #test versus a simpler model
#' fit2<-logistf(data=sex2, case~age+oc+dia)
#' # or: fit2<-update(fit, case~age+oc+dia)
#' anova(fit,fit2)
#' 
#' # comparison of non-nested models (with different df):
#' fit3<-logistf(data=sex2, case~age+vic+vicl+vis)
#' anova(fit2,fit3, method="PLR")
#' 
#' 
#' @rdname anova
#' @method anova logistf
#' @exportS3Method anova logistf
anova.logistf<-function(object,  fit2, formula, method="nested", ...){
 # methods: "PLR": take difference in PLR, "nested": proper method for nested models
 # needed in logistf class: $firth, $data
  
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("object","fit2","formula","method"), names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  fit2 <- eval(mf$fit2, parent.frame())
  
  fit1<-object
  if(missing(formula)){
    out2 <- fit2$formula
    if(fit1$df<fit2$df) {
    ff0<-fit1 #swap
    fit1<-fit2
    fit2<-ff0
    }
  } else {
    out2 <- formula
  }
   
  if(method=="PLR"){
    if(fit1$df==fit2$df) stop("Models not comparable (equal df).\n")
    df<-abs(fit1$df-fit2$df)
    PLR1<- -2 * (fit1$loglik['null']-fit1$loglik['full'])
    PLR2<- -2 * (fit2$loglik['null']-fit2$loglik['full'])
    chisq<-PLR1-PLR2
    if (chisq<0) chisq<-0
    pval<-1-pchisq(chisq,df)
    model2<-as.character(fit2$formula)
  } else if(method=="nested"){
    f1<-fit1$formula
    a<-attr(terms(fit1),"term.labels")
    #check if just intercept fitted: 
    if(missing(formula)){
      f2<-fit2$formula
      b<-attr(terms(fit2),"term.labels")
      # find out about which model is nested in the other
      upper<-f1
      lower<-f2
      if(!any(is.na(match(a,b)))) {
        lower<-f1
        upper<-f2
        ab<-a
        a<-b
        b<-ab
      } else if(any(is.na(match(b,a)))) {
        stop("Models are not nested. Try method=PLR on non-nested models.\n")
      }
      aug<-a[is.na(match(a,b))]
      f3<-paste("~",aug[1])
      if(length(aug)>1) {
        for(j in 2:length(aug)) {
          f3<-paste(f3, aug[j], sep="+")
        }
      }
      f3<-paste(f3, "-1")
      f3<-as.formula(f3, env = environment(object$formula))
    } else {
      f3<-as.formula(paste(paste(as.character(formula), collapse=""),"-1",collapse=""), env = environment(object$formula))
    }
      
    test<-logistftest(object=fit1, test = f3, firth=fit1$firth, weights=fit1$weights,...)
      
    chisq<- -2 * (test$loglik['null']-test$loglik['full'])
    PLR1<- -2 * (fit1$loglik['null']-fit1$loglik['full'])
    PLR2<-PLR1-chisq
    df<-test$df
    pval<-test$prob
    model2<-as.character(f3)
  }
  res<-list(chisq=chisq, df=df, pval=pval, call=match.call(), method=method, model1=as.character(fit1$formula), model2=out2, PLR1=PLR1, PLR2=PLR2)
  attr(res,"class")<-"anova.logistf"
  return(res)
}

#' @method anova flic
#' @exportS3Method anova flic
anova.flic<-function(object,  fit2, formula, method="nested", ...){
  anova.logistf(object,  fit2, formula, ...)
}

#' @exportS3Method print anova.logistf
print.anova.logistf <-function(x,...){
    cat("Comparison of logistf models:\n")
    obj<-x
    pdat<-data.frame(Formula=c(as.character(obj$model1),as.character(obj$model2)), ChiSquared=c(obj$PLR1,obj$PLR2))
    pdatrows<-capture.output(print(pdat))
    for(i in 1:3) cat(pdatrows[i],"\n")
    cat("\nMethod: ", obj$method, "\n")
    cat("Chi-Squared: ", obj$chisq, "  df=",obj$df,"  P=", obj$pval,"\n")
}

