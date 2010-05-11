\name{ncvreg-package}
\alias{ncvreg-package}
\docType{package}
\title{Regularization paths for SCAD- and MCP-penalized regression models}
\description{Efficient algorithms for fitting regularization paths for
  linear or logistic regression models penalized by MCP or SCAD, with
  optional additional L2 penalty ("Mnet").}
\details{
  \tabular{ll}{
    Package: \tab ncvreg\cr
    Type: \tab Package\cr
    Version: \tab 2.1\cr
    Date: \tab 2010-05-11\cr
    License: \tab GPL-2\cr}
  Accepts a design matrix \code{X} and vector of responses \code{y},
  produces the regularization path over a grid of values for the tuning
  parameter \code{lambda}.  Also provides methods for plotting and for
  determining locally convex regions of the coefficients paths.}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  University of
  Kentucky Technical Report \#403, available at \url{http://web.as.uky.edu/statistics/users/pbreheny/ncvreg.pdf}.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

fit <- ncvreg(X,y)
plot(fit)
}
