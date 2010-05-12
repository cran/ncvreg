cv.ncvreg <- function(X,y,family="gaussian",penalty="MCP",nfolds=10,seed,...)
  {
    if (!missing(seed)) set.seed(seed)
    n <- length(y)
    fit <- ncvreg(X,y,family=family,penalty=penalty,...)
    error <- array(NA,dim=c(nfolds,length(fit$lambda)))
    
    if (family=="gaussian")
      {
        cv.ind <- ceiling((1:n)/n*nfolds)
      }
    else if (family=="binomial")
      {
        ind1 <- which(y==1)
        ind0 <- which(y==0)
        n1 <- length(ind1)
        n0 <- length(ind0)
        cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
        cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
        cv.ind <- numeric(n)
        cv.ind[y==1] <- cv.ind1
        cv.ind[y==0] <- cv.ind0
      }

    for (i in 1:nfolds)
      {
        X1 <- X[cv.ind!=i,]
        y1 <- y[cv.ind!=i]
        X2 <- X[cv.ind==i,]
        y2 <- y[cv.ind==i]

        fit.i <- ncvreg(X1,y1,family=family,penalty=penalty,lambda=fit$lambda,...)
        yhat <- predict(fit.i,X2,type="response")
        error[i,1:ncol(yhat)] <- loss.ncvreg(y2,yhat,family)
      }
    fit$allerror <- error
    fit$error <- apply(error,2,mean,na.rm=TRUE)
    fit$cv <- which.min(fit$error)
    return(fit)
  }

