ncvreg <- function(X, y, family=c("gaussian","binomial"), penalty=c("MCP","SCAD"), a=3, lambda.min=ifelse(n>p,.001,.05), n.lambda=100, eps=.001, max.iter=500, convex=TRUE)
  {
    ## Error checking
    family <- match.arg(family)
    penalty <- match.arg(penalty)
    if (a <= 1 & penalty=="MCP") stop("a must be greater than 1 for the MC penalty")
    if (a <= 2 & penalty=="SCAD") stop("a must be greater than 2 for the SCAD penalty")
    if (n.lambda < 2) stop("n.lambda must be at least 2")

    ## Set up XX, yy, lambda
    n <- length(y)
    p <- ncol(X)
    meanx <- apply(X,2,mean)
    normx <- sqrt(apply((t(X)-meanx)^2,1,sum)/n)
    if (any(normx < 0.0001)) stop("X contains columns which are numerically constant, please remove them; an intercept is included automatically")
    XX <- scale(X,meanx,normx)
    if (family=="gaussian") yy <- y - mean(y)
    else yy <- y
    lambda <- setupLambda(XX,yy,family,penalty,a,gamma,lambda.min,n.lambda)

    ## Fit
    if (family=="gaussian")
      {
        fit <- .C("cdfit_gaussian",double(p*n.lambda),integer(n.lambda),as.double(XX),as.double(yy),as.integer(n),as.integer(p),penalty,as.double(lambda),as.integer(n.lambda),as.double(eps),as.integer(max.iter),as.double(a))
        beta <- rbind(0,matrix(fit[[1]],nrow=p))
        iter <- fit[[2]]
      }
    if (family=="binomial")
      {
        fit <- .C("cdfit_binomial",double(n.lambda),double(p*n.lambda),integer(n.lambda),as.double(XX),as.double(yy),as.integer(n),as.integer(p),penalty,as.double(lambda),as.integer(n.lambda),as.double(eps),as.integer(max.iter),as.double(a))
        beta <- rbind(fit[[1]],matrix(fit[[2]],nrow=p))
        iter <- fit[[3]]
        
        ## Eliminate saturated lambda values, if any
        ind <- !is.na(beta[1,])
        beta <- beta[,ind]
        iter <- iter[ind]
        lambda <- lambda[ind]
      }
    if (any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda")

    if (convex) convex.min <- convexMin(beta,XX,penalty,a,family)

    ## Unstandardize
    beta[-1,] <- beta[-1,]/normx
    if (family=="gaussian") beta[1,] <- mean(y) - crossprod(meanx,beta[-1,,drop=FALSE])
    if (family=="binomial") beta[1,] <- beta[1,] - crossprod(meanx,beta[-1,,drop=FALSE])

    ## Names
    if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
    else varnames <- colnames(X)
    varnames <- c("(Intercept)",varnames)
    dimnames(beta) <- list(varnames,round(lambda,digits=4))

    ## Output
    val <- list(beta=beta,iter=iter,lambda=lambda,penalty=penalty,family=family,a=a,convex.min=convex.min)
    class(val) <- "ncvreg"
    return(val)
  }
