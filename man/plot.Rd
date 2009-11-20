\name{plot.ncvreg}
\alias{plot.ncvreg}
\title{Plot coefficients from a "ncvreg" object}
\description{Produces a plot of the coefficient paths for a fitted
  \code{"ncvreg"} object.}
\usage{
\method{plot}{ncvreg}(x, col, alpha=1, xlab, ylab=expression(hat(beta)), log.l=FALSE, shade=TRUE, ylim, ...)
}
\arguments{
  \item{x}{Fitted \code{"ncvreg"} model.}
  \item{col}{Color of coefficient lines.  User can specify either a
    single color or a vector of colors, one for each coefficient.
    Defaults to an evenly spaced rainbow of colors.}
  \item{alpha}{Controls alpha-blending.  Default is alpha=1.}
  \item{xlab}{Horizontal axis label.}
  \item{ylab}{Vertical axis label.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    FALSE.} 
  \item{shade}{Should nonconvex region be shaded?  Default is TRUE.}
  \item{ylim}{Vertical axis limits.}
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  Available at \url{http://web.as.uky.edu/statistics/techreports/tr403/tr403.pdf}.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(prostate)
X <- as.matrix(prostate[,1:8])
y <- prostate$lpsa

fit <- ncvreg(X,y)
plot(fit)
plot(fit,color="black")
plot(fit,log=TRUE)
}
\keyword{models}
\keyword{regression}
