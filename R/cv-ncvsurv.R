#' @rdname cv.ncvreg
#' @export

cv.ncvsurv <- function(X, y, ..., cluster, nfolds=10, seed, fold, se=c('quick', 'bootstrap'), returnY=FALSE, trace=FALSE) {
  se <- match.arg(se)

  # Coersion
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (!inherits(y, "matrix")) {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be a matrix or able to be coerced to a matrix", call.=FALSE)
    if (ncol(y) != 2) stop("y must have two columns for survival data: time-on-study and a censoring indicator", call.=FALSE)
  }
  if (typeof(y) == "integer") storage.mode(y) <- "double"

  # Complete data fit
  fit <- ncvsurv(X=X, y=y, ...)

  # Set up folds
  n <- nrow(X)
  sde <- sqrt(.Machine$double.eps)
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  if (missing(fold)) {
    ind1 <- which(fit$fail==1)
    ind0 <- which(fit$fail==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds
    fold0[fold0==0] <- nfolds
    fold <- integer(n)
    fold[fit$fail==1] <- sample(fold1)
    fold[fit$fail==0] <- sample(fold0)
  } else {
    fold <- fold[fit$order]
    nfolds <- max(fold)
  }

  Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$penalty.factor <- fit$penalty.factor
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("fold","fit","X", "y", "cv.args"), envir=environment())
    parallel::clusterCall(cluster, function() library(ncvreg))
    fold.results <- parallel::parLapply(cl=cluster, X=1:nfolds, fun=cvf.surv, XX=X, y=y, fold=fold, cv.args=cv.args)
  }

  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #", i, sep="","\n")
      res <- cvf.surv(i, X, y, fold, cv.args)
    }
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(Y), 2, all))
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  # Return
  if (se == "quick") {
    L <- loss.ncvsurv(y, Y, total=FALSE)
    cve <- apply(L, 2, sum)/sum(fit$fail)
    cvse <- apply(L, 2, stats::sd)*sqrt(nrow(L))/sum(fit$fail)
  } else {
    cve <- as.double(loss.ncvsurv(y, Y))/sum(fit$fail)
    cvse <- se.ncvsurv(y, Y)/sum(fit$fail)
  }
  min <- which.min(cve)

  val <- list(cve=cve, cvse=cvse, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=cve[1])
  if (returnY) val$Y <- Y
  structure(val, class=c("cv.ncvsurv", "cv.ncvreg"))
}
cvf.surv <- function(i, XX, y, fold, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i,]
  fit.i <- do.call("ncvsurv", cv.args)

  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i,]
  nl <- length(fit.i$lambda)
  yhat <- predict(fit.i, X2)

  list(nl=length(fit.i$lambda), yhat=yhat)
}
