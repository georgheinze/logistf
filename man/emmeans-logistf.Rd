\name{emmeans-logistf}
\alias{emmeans-logistf}

\title{Emmeans support for logistf}

\description{
Support for the \code{emmeans} package is available. See below for an example of using \code{emmeans::emmeans()} with a \code{logistf} object.
}


\examples{

data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)

emmeans::emmeans(fit, ~age+dia)

}
