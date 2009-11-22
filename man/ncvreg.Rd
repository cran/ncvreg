\name{ncvreg}
\alias{ncvreg}
\title{Fit an MCP- or SCAD-penalized regression path}
\description{Fit coefficients paths for MCP- or SCAD-penalized
  regression models over a grid of values for the regularization
  parameter lambda.  Fits linear and logistic regression models.}
\usage{
ncvreg(X, y, family=c("gaussian","binomial"), penalty=c("MCP","SCAD"),
a=3, lambda.min=ifelse(n>p,.001,.05), n.lambda=100, eps=.001,
max.iter=500, convex=TRUE)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{ncvreg}
  standardizes the data and includes an intercept by default.}
  \item{y}{The response vector.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{penalty}{The penalty to be applied to the model.  Either "MCP"
    (the default) or "SCAD".}
  \item{a}{The tuning parameter of the MCP/SCAD penalty (see details).}
  \item{lambda.min}{The smallest value for lambda, as a fraction of
    lambda.max.  Default is .001 if the number of observations is larger
    than the number of covariates and .05 otherwise.}
  \item{n.lambda}{The number of lambda values.  Default is 100.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{.001}.}
  \item{max.iter}{Maximum number of iterations.  Default is 500.  See
    details.}
  \item{convex}{Calculate index for which objective function ceases to
    be locally convex?  Default is TRUE.}}
\details{
  The sequence of models indexed by the regularization parameter
  \code{lambda} is fit using a coordinate descent algorithm.  For
  logistic regression models, some care is taken to avoid model
  saturation; the algorithm may exit early in this setting.  The
  objective function is defined to be
  \deqn{1/(2*n)RSS + penalty}
  for \code{"gaussian"} and
  \deqn{-(1/n) loglik + penalty}
  for \code{"binomial"}, where the likelihood is from a traditional
  generalized linear model for the log-odds of an event.\cr
  \cr
  This algorithm is stable, very efficient, and generally converges
  quite rapidly to the solution.  For logistic regression, adaptive
  rescaling (see reference) is used.\cr
  \cr
  The convexity diagnostics rely on a fine covering of
  (lambda.min,lambda.max); choosing a low value of \code{n.lambda} may
  produce unreliable results.}
\value{
  An object with S3 class \code{"ncvreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{n.lambda}.}
  \item{iter}{A vector of length \code{n.lambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
  \item{lambda}{The sequence of regularization parameter values in the
    path.}
  \item{penalty}{Same as above.}
  \item{family}{Same as above.}
  \item{a}{Same as above.}
  \item{convex.min}{The last index for which the objective function is
    locally convex.  The smallest value of lambda for which the
    objective function is convex is therefore \code{lambda[convex.min]},
    with corresponding coefficients \code{beta[,convex.min]}.}}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  Available at \url{http://web.as.uky.edu/statistics/techreports/tr403/tr403.pdf}.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{plot.ncvreg}}
\examples{
## Linear regression
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

par(mfrow=c(2,2))
fit <- ncvreg(X,y)
plot(fit,main="a=3")
fit <- ncvreg(X,y,a=10)
plot(fit,main="a=10")
fit <- ncvreg(X,y,a=1.5)
plot(fit,main="a=1.5")
fit <- ncvreg(X,y,penalty="SCAD")
plot(fit,main="SCAD")

## Logistic regression
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd

par(mfrow=c(2,2))
fit <- ncvreg(X,y,family="binomial")
plot(fit,main="a=3")
fit <- ncvreg(X,y,family="binomial",a=10)
plot(fit,main="a=10")
fit <- ncvreg(X,y,family="binomial",a=1.5)
plot(fit,main="a=1.5")
fit <- ncvreg(X,y,family="binomial",penalty="SCAD")
plot(fit,main="SCAD")
}
\keyword{models}
\keyword{regression}
