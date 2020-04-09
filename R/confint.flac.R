#' @exportS3Method confint flac
confint.flac <- function(object,parm, ...){
  level<-1-object$alpha
  levstr<-paste(level*100,"%",sep="")
  cimat<-cbind(object$ci.lower,object$ci.upper)
  rownames(cimat)<-names(object$ci.lower)
  colnames(cimat)<-c(paste("Lower ",levstr,sep=""),paste("Upper ",levstr,sep=""))
  return(cimat)
}