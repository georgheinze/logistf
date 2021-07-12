library(logistf)
form <- formula('case ~ age+oc+vic+vicl+vis+dia')
terms.fit <- c(1,2)
df <- sex2

fitcontrol = logistf.fit.control(terms.fit = terms.fit)

tol_coef <- 1e-1
tol <- 1e-5

L0 <- -169.6506

expect_equal(
  suppressWarnings(flic(logistf(form, df, fitcontrol = fitcontrol))$coefficients),
  suppressWarnings(flic(form, df, fitcontrol = fitcontrol)$coefficients)
)

expect_equal(
  suppressWarnings(flac(logistf(form, df, fitcontrol = fitcontrol), data = df)$coefficients),
  suppressWarnings(flac(form, df, fitcontrol = fitcontrol)$coefficients)
)
