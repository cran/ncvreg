convexMin <- function(beta,X,penalty,a,family)
  {
    n <- nrow(X)
    p <- ncol(X)
    l <- ncol(beta)
    
    if (penalty=="MCP") p.. <- 1/a
    else if (penalty=="SCAD") p.. <- 1/(a-1)

    val <- NULL
    for (i in 2:(l-1))
      {
        if (is.na(beta[1,i+1])) break
        A1 <- beta[-1,i]==0
        A2 <- beta[-1,i+1]==0
        if (all(A1==A2)) next
        U <- A1&A2
        Xu <- X[,!U]
        if (family=="gaussian")
          {
            cmin.i <- min(eigen(crossprod(Xu)/n-diag(rep(p..,sum(!U))))$values)
          }
        if (family=="binomial")
          {
            eta <- beta[1,i+1] + X%*%beta[-1,i+1]
            pi. <- exp(eta)/(1+exp(eta))
            w <- as.numeric(pi.*(1-pi.))
            w[eta > log(.9999/.0001)] <- .0001
            w[eta < log(.0001/.9999)] <- .0001
            Xu <- sqrt(w) * cbind(1,Xu)
            xwxn <- crossprod(Xu)/n
            cmin.i <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*p..)))$values)
          }
        if (cmin.i < 0)
          {
            val <- i
            break
          }
      }
    return(val)
  }
