plot.ncvreg <- function(x, col, alpha=1, log.l=FALSE, shade=TRUE, ...)
  {
    zeros <- which(apply(abs(x$beta),1,sum)==0)
    beta <- x$beta[-c(1,zeros),,drop=FALSE]
    p <- nrow(beta)
    l <- x$lambda
    if (log.l)
      {
        l <- log(l)
        xlab <- expression(log(lambda))
      }
    else xlab <- expression(lambda)
    ylim <- range(beta)

    plot.args = list(x=l, y=rep(0,length(l)), ylim=ylim, xlab=xlab, ylab=expression(hat(beta)), type="n", xlim=rev(range(l)))
    new.args = list(...)
    if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)

    if (shade & !is.null(x$convex.min))
      {
        l1 <- l[x$convex.min]
        l2 <- min(l)
        polygon(x=c(l1,l2,l2,l1),y=c(ylim[1],ylim[1],ylim[2],ylim[2]),col="gray85",border=FALSE) 
      }
    if (missing(col)) col <- hsv(seq(0,1,len=(p+1)),alpha=alpha)[1:p]
    n.col <- length(col)
    for (i in 1:p)
      {
        lines(l,beta[i,],col=col[(i-1)%%n.col+1])
      }
    abline(h=0)
  }
