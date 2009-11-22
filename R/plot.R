plot.ncvreg <- function(x, col, alpha=1, xlab, ylab=expression(hat(beta)), log.l=FALSE, shade=TRUE, ylim, ...)
  {
    zeros <- which(apply(abs(x$beta),1,sum)==0)
    beta <- x$beta[-c(1,zeros),,drop=FALSE]
    p <- nrow(beta)
    if (missing(ylim)) ylim <- range(beta)
    if (missing(col)) col <- hsv(seq(0,1,len=(p+1)),alpha=alpha)[1:p]
    l <- x$lambda
    if (log.l)
      {
        l <- log(l)
        if (missing(xlab)) xlab <- expression(log(lambda))
      }
    else if (missing(xlab)) xlab <- expression(lambda)

    plot(l,rep(0,length(l)),ylim=ylim,xlim=rev(range(l)),xlab=xlab,ylab=ylab,type="l",...)
    
    if (shade & !is.null(x$convex.min))
      {
        l1 <- l[x$convex.min]
        l2 <- min(l)
        polygon(x=c(l1,l2,l2,l1),y=c(ylim[1],ylim[1],ylim[2],ylim[2]),col="gray85",border=FALSE) 
      }
    for (i in 1:p)
      {
        lines(l,beta[i,],col=col[i])
      }
    abline(h=0)
  }
