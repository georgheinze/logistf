#' Generate a Basis Matrix for Natural Cubic Splines
#'
#' Generate the B-spline basis matrix for a natural cubic spline. This function calls [splines::ns()] but sets the boundary knots, 
#' if not given, according to the number of knots and to the guidelines of [Hmisc::rcspline.eval()]. 
#' For further information please see the individual documentations.
#' 
#' @param x the predictor variable
#' @param df degrees of freedom. If supplied, df - 1 - intercept knots at suitably chosen quantiles of x are set.
#' @param knots knots of the spline.
#' @param intercept if \code{TRUE}, an intercept is included in the basis
#' @param Boundary.knots boundary points at which to impose the natural boundary conditions and anchor the B-spline basis
#'
#' @seealso [splines::ns()],[Hmisc::rcspline.eval()]
#' @return A matrix of dimension \code{length(x) * df} where either \code{df} was supplied or if knots were supplied, 
#' \code{df = length(knots) + 1 + intercept}.
#' @export
#'
#' @examples
#' 
#' data(sex2)
#' sex2$agec <- runif(nrow(sex2), 18, 25)
#' rcslf(sex2$agec, df = 3)
#' summary(fit <- logistf(case ~ rcslf(sex2$agec, df = 3), data = sex2))
rcslf <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots) {
    if(missing(Boundary.knots)){
      xx <- x[!is.na(x)]
      n <- length(xx)
      nknots <- if (is.null(df) && length(knots)>0) length(knots)+2 else df - 1 - intercept + 2
      outer <- if (nknots > 3) 0.05 else 0.1 
      if (nknots > 6) {
        outer <- 0.025
      }
      if(is.null(knots)){
        p <- seq(outer, 1 - outer, length = nknots)
        all <- quantile(xx, p)
        Boundary.knots <- c(all[1], all[nknots])
        knots <- all[-c(1,nknots)]
      }
      else {
        p <- seq(outer, 1 - outer, length = 2)
        Boundary.knots <- quantile(xx, p)
      }
    }
  splines::ns(x=x, df=df, knots=knots, intercept=intercept, Boundary.knots = Boundary.knots)
}
