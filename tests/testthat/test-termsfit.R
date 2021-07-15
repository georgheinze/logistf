library(logistf)
form <- formula('case ~ age+oc+vic+vicl+vis+dia')
terms.fit <- c(1,2)
df <- sex2

fitcontrol = logistf.fit.control(terms.fit = terms.fit)

tol_coef <- 1e-1
tol <- 1e-5

L0 <- -170


#Basic: Newton Raphson
##Coefficients
expect_equal(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol)$coefficients), 
  c('(Intercept)'= 0.3023695, 'age'=-0.8159534, 'oc'=0, 'vic'=0, 'vicl'=0, 'vis'=0, 'dia'=0), 
  tolerance = tol_coef
)
##Controls: 
###maxit:default
expect_lte(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol)$iter), 
  logistf.control()$maxit
)


##Loglik of Null model
expect_gt(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol)$loglik[1]), 
  L0
)
##Loglik of full model: greater than loglik of null model
expect_gt(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol)$loglik[2]), 
  L0
)

#Basic:IRLS
##Coefficients
expect_equal(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol, control = logistf.control(fit = 'IRLS'))$coefficients), 
  c('(Intercept)'= 0.3023695, 'age'=-0.8159534, 'oc'=0, 'vic'=0, 'vicl'=0, 'vis'=0, 'dia'=0), 
  tolerance = tol_coef
)
##Controls: 
###maxit:default
expect_lte(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol, control = logistf.control(fit = 'IRLS'))$iter), 
  logistf.control()$maxit
)

##Loglik of Null model
expect_gt(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol, control = logistf.control(fit = 'IRLS'))$loglik[1]), 
  L0
)
##Loglik of full model: greater than loglik of null model
expect_gt(
  suppressWarnings(logistf(form, df, fitcontrol = fitcontrol, control = logistf.control(fit = 'IRLS'))$loglik[2]), 
  L0
)