\name{predict.ncvreg}
\alias{predict.ncvreg}
\alias{coef.ncvreg}
\title{Model predictions based on a fitted "ncvreg" object.}
\description{
Similar to other predict methods, this function returns predictions from
a fitted \code{"ncvreg"} object.
}
\usage{
\method{predict}{ncvreg}(object,X,which=1:length(object$lambda),type=c("link","response","class","coefficients"),...)
\method{coef}{ncvreg}(object,which=1:length(object$lambda),...)
}
\arguments{
  \item{object}{Fitted \code{"ncvreg"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"}}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients.}
  \item{\dots}{Not used.}}
\value{The object returned depends on type.}
\references{Breheny, P. and Huang, J. (2009)  Coordinate descent
  algorithms for nonconvex penalized regression methods.  University of
  Kentucky Technical Report \#403, available at \url{http://web.as.uky.edu/statistics/users/pbreheny/ncvreg.pdf}.}
\author{Patrick Breheny <patrick.breheny@uky.edu>}
\seealso{\code{ncvreg}}
\examples{
data(heart)
X <- as.matrix(heart[,1:9])
y <- heart$chd

par(mfrow=c(2,2))
fit <- ncvreg(X,y,family="binomial")
coef(fit)
predict(fit,X,type="link")
predict(fit,X,type="response")
predict(fit,X,type="class")
}
\keyword{models}
\keyword{regression}
