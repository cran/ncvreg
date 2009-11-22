setupLambda <- function(X,y,family,penalty,a,gamma,lambda.min,n.lambda)
  {
    n <- nrow(X)
    p <- ncol(X)

    ## Determine lambda.max
    if (family=="gaussian")
      {
        r <- y - mean(y)
        lambda.max <- max(abs(crossprod(X,r)/n))    
      }
    if (family=="binomial")
      {
        fit <- glm(y~0,family="binomial")
        pi. <- fit$fitted.values
        w <- pi.*(1-pi.)
        r = (y - pi.)/w
        lambda.max <- max(abs(crossprod(X,w*r)/n))
      }
    
    if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=n.lambda-1)),0)
    else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=n.lambda))
    return(lambda)
  }
