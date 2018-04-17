## ------------------------------------------------------------------------
# Linear regression
data(Prostate)
head(Prostate$X)
head(Prostate$y)

## ------------------------------------------------------------------------
fit <- ncvreg(Prostate$X, Prostate$y)

## ----opts.label='fig'----------------------------------------------------
plot(fit)

## ----opts.label='fig'----------------------------------------------------
coef(fit, lambda=0.1)

## ----opts.label='fig'----------------------------------------------------
cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
plot(cvfit)

## ------------------------------------------------------------------------
coef(cvfit)

## ------------------------------------------------------------------------
predict(cvfit, X=head(Prostate$X))
predict(cvfit, type="nvars")

