\name{ncvreg-package}
\alias{ncvreg-package}
\docType{package}
\title{Regularization paths for SCAD- and MCP-penalized regression models}
\description{Efficient algorithms for fitting the regularization path for linear or logistic regression models penalized by MCP or SCAD.}
\details{
  \tabular{ll}{
    Package: \tab ncvreg\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2009-11-20\cr
    License: \tab GPL-2\cr}
  Accepts a design matrix \code{X} and vector of responses \code{y},
  produces the regularization path over a grid of values for the tuning
  parameter \code{lambda}.  Also provides methods for plotting and for
  determining locally convex regions of the coefficients paths.}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  Available at \url{http://web.as.uky.edu/statistics/techreports/tr403/tr403.pdf}.}
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
