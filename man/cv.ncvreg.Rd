\name{cv.ncvreg}
\alias{cv.ncvreg}
\title{Cross-validation for ncvreg}
\description{Performs k-fold cross validation for MCP- or SCAD-penalized
  regression models over a grid of values for the regularization
  parameter lambda.}
\usage{
cv.ncvreg(X,y,family="gaussian",penalty="MCP",nfolds=10,seed,...)
}
\arguments{
  \item{X}{The design matrix, without an intercept.  \code{ncvreg}
  standardizes the data and includes an intercept by default.}
  \item{y}{The response vector.}
  \item{family}{Either "gaussian" or "binomial", depending on the
    response.}
  \item{penalty}{The penalty to be applied to the model.  Either "MCP"
    (the default) or "SCAD".}
  \item{nfolds}{The number of cross-validation folds.  Default is 10.}
  \item{seed}{You may set the seed of the random number generator in
    order to obtain reproducible results.}
  \item{...}{Additional arguments to \code{ncvreg}.}}
\details{
  The function calls \code{ncvreg} \code{nfolds} + 1 times; once with
  all of the data, and once leaving each fold out.  The cross-validation
  error is based on the loss function used to fit the model; see
  \code{ncvreg} for details.}
\value{
  An object with S3 class \code{"ncvreg"}.  In addition to the
  \code{ncvreg} fit, the object contains:
  \item{allerror}{A matrix with \code{nfolds} rows and \code{nlambda}
    columns, containing the cross-validation error for each fit.}
  \item{error}{The error for each value of \code{lambda}, averaged
    across the cross-validation folds.}
  \item{cv}{The index of \code{lambda} corresponding to the minimum
    cross-validation error.}}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  University of
  Kentucky Technical Report \#403, available at \url{http://web.as.uky.edu/statistics/users/pbreheny/ncvreg.pdf}.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

fit <- cv.ncvreg(X,y)
plot(fit$lambda,fit$error)
beta <- fit$beta[,fit$cv]
}
\keyword{models}
\keyword{regression}
