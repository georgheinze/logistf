library(logistf)
form <- formula('case ~ age+oc+vic+vicl+vis+dia')
df <- sex2

tol_coef <- 1e-1
tol <- 1e-5

L0 <- -170


#Basic: Newton Raphson
##Coefficients
expect_equal(
  suppressWarnings(logistf(form, df)$coefficients), 
  c('(Intercept)'=0.12025404, 'age'=-1.10598130, 'oc'=-0.06881673, 'vic'=2.26887464, 'vicl'=-2.11140816, 'vis'=-0.78831694, 'dia'=3.09601263), 
  tolerance = tol_coef
)
##Controls: 
###maxit:default
expect_lte(
  suppressWarnings(logistf(form, df)$iter), 
  logistf.control()$maxit
)
###maxit:100
expect_lte(
  suppressWarnings(logistf(form, df, control=logistf.control(maxit=100))$iter), 
  100
)
###Convergence lconv, gconv, xconv
expect_lte(
  suppressWarnings(logistf(form, df)$conv[1]), 
  logistf.control()$lconv
)
expect_lte(
  suppressWarnings(logistf(form, df)$conv[2]), 
  logistf.control()$gconv
)
expect_lte(
  suppressWarnings(logistf(form, df)$conv[3]), 
  logistf.control()$xconv
)

##Loglik of Null model
expect_gt(
  suppressWarnings(logistf(form, df)$loglik[1]), 
  L0
)
##Loglik of full model: greater than loglik of null model
expect_gt(
  suppressWarnings(logistf(form, df)$loglik[2]), 
  L0
)





