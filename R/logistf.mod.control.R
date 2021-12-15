#' Controls additional parameters for \code{logistf}
#'
#' Sets parameters for \code{logistf} calls.
#'
#' @param tau  Penalization parameter (default = 0.5)
#' @param terms.fit A numeric vector of terms to fit. Intercept has to be included if needed.
#'
#' @return
#'    \item{tau}{Penalization parameter (default = 0.5)}
#'    \item{terms.fit}{A numeric vector of terms to fit. Intercept has to be included if needed.}
#' @export
#' 
#' @encoding UTF-8
#' @examples
#' data(sexagg)
#' fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, 
#' modcontrol=logistf.mod.control(terms.fit=c(1,2)))
#' summary(fit2)
#' 
logistf.mod.control <- function(tau=0.5, terms.fit = NULL){
      if(!is.null(terms.fit) & !is.numeric(terms.fit)){
          stop("Please provide a numeric vector of terms to fit.")
      }
    res<-list(tau=tau,
              terms.fit = terms.fit, 
              call=match.call())
    attr(res, "class")<-"logistf.mod.control"
    return(res)
  }