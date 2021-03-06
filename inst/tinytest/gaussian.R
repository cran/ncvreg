library(glmnet, quietly=TRUE)

# ncvreg works for linear regression
n <- 50
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
y <- rnorm(n, X%*%b)
beta <- lm(y~X)$coef
scad <- coef(ncvreg(X, y, lambda=1:0, penalty="SCAD", eps=.0001), which=2)
mcp <- coef(ncvreg(X, y, lambda=1:0, penalty="MCP", eps=.0001), which=2)
expect_equivalent(scad, beta, tolerance=.01)
expect_equivalent(mcp, beta, tolerance=.01)

# logLik() is correct: gaussian
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
y <- rnorm(n)
fit.mle <- lm(y~X)
fit <- ncvreg(X, y, lambda.min=0)
expect_equivalent(logLik(fit)[100], logLik(fit.mle)[1], tol= .001)
expect_equivalent(AIC(fit)[100], AIC(fit.mle), tol= .001)

# ncvreg reproduces lasso: gaussian
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
par(mfrow=c(3,2))
y <- rnorm(n)
nlasso <- coef(fit <- ncvreg(X, y, penalty="lasso"))
plot(fit, log=TRUE)
glasso <- as.matrix(coef(fit <- glmnet(X, y, lambda=fit$lambda)))
plot(fit, "lambda")
expect_equivalent(nlasso,  glasso, tolerance=.01)

# cv.ncvreg() options work for gaussian
n <- 50
p <- 10
X <- matrix(rnorm(n*p), ncol=p)
b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=1)

par(mfrow=c(2,2))
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")
summary(cvfit)
predict(cvfit, type="coefficients")
predict(cvfit, type="vars")
predict(cvfit, type="nvars")

b <- c(-3, 3, rep(0, 8))
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")

b <- rep(0, 10)
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.ncvreg(X, y)
plot(cvfit, type="all")

###############################################
# ncvreg dependencies work: gaussian
###############################################

# Predict
fit <- ncvreg(X, y, lambda.min=0)
p <- predict(fit, X, 'link', lambda=0.1)
p <- predict(fit, X, 'link')
p <- predict(fit, X, 'response')
p <- predict(fit, X, 'coef')
p <- predict(fit, X, 'vars')
p <- predict(fit, X, 'nvars')

# Integers
X <- matrix(rpois(500, 1), 50, 10)
y <- rpois(50, 1)
fit <- ncvreg(X, y)

# Data frame
fit <- ncvreg(as.data.frame(X), y)

# Penalty factor
fit <- ncvreg(as.data.frame(X), y, penalty.factor=c(0:9))

# User lambdas
fit <- ncvreg(as.data.frame(X), y, lambda=c(1, 0.1, 0.01))

# ReturnX
fit <- ncvreg(as.data.frame(X), y, returnX=TRUE)

# Constant columns
fit <- ncvreg(cbind(5, X), y)

# Plot
plot(fit)
plot(fit, log.l=TRUE)

# Summary
summary(fit, which=10)
summary(fit, lam=0.05)
