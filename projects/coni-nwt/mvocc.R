if (FALSE) {
## plain mv occupancy model
## Y: M x J matrix
## X: M x nOP matrix
## Z: MxJ x nDP matrix (visits stacked)
.mvocc <- function(Y, X, Z, offsetx, offsetz, weights,
method="Nelder-Mead", inits, control=list(), hessian=FALSE, ...) {

    good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    Y <- ifelse(Y > 0, 1, 0)

    J <- ncol(Y)
    M <- nrow(Y)
    nDP <- ncol(Z)
    nOP <- ncol(X)
    nP <- nDP + nOP

    if (missing(offsetx))
        offsetx <- rep(0, M)
    if (missing(offsetz))
        offsetz <- rep(0, M)
    if (missing(weights))
        weights <- rep(1, M)
    if (missing(inits))
        inits <- rep(0, nP)

    yvec <- as.numeric(Y)
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(Y, na.rm = TRUE) == 0, 1, 0)
    rowProds <- function (x, na.rm = FALSE)
        exp(rowSums(log(x), na.rm = na.rm))

    nll <- function(params) {
        psi <- plogis(X %*% params[1:nOP] + offsetx)
        pvec <- plogis(Z %*% params[(nOP + 1):nP] + offsetz)
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1
        cpmat <- matrix(cp, M, J)
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        loglik <- sum(loglik * weights)
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    o <- optim(inits, nll,
        method = method, hessian = hessian, control = control, ...)
    o
}

summary(t(replicate(100, {
phi <- 0.5
pd <- 0.7
M <- 1000
J <- 10

W <- rbinom(M, 1, phi)
Y <- matrix(0, M, J)
for (j in 1:J)
    Y[,j] <- rbinom(M, W, pd)

X <- matrix(1, M, 1)
Z <- matrix(1, M*J, 1)

o <- .mvocc(Y, X, Z)
plogis(o$par)
})))

## modified mv occupancy model
## Y: M x J matrix
## X: M x nOP matrix
## Z: MxJ x nDP matrix (visits stacked)
.mvocc2 <- function(Y, X, Z, p, A, weights,
method="Nelder-Mead", inits, control=list(), hessian=FALSE, ...) {

    good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    Y <- ifelse(Y > 0, 1, 0)

    J <- ncol(Y)
    M <- nrow(Y)
    nDP <- ncol(Z)
    nOP <- ncol(X)
    nP <- nDP + nOP

    if (missing(p))
        p <- rep(1, M)
    if (missing(A))
        A <- rep(1, M)
    if (missing(weights))
        weights <- rep(1, M)
    if (missing(inits))
        inits <- rep(0, nP)

    yvec <- as.numeric(Y)
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(Y, na.rm = TRUE) == 0, 1, 0)
    rowProds <- function (x, na.rm = FALSE)
        exp(rowSums(log(x), na.rm = na.rm))

    nll <- function(params) {
        D <- exp(X %*% params[1:nOP])
        psi <- 1-exp(-D * A)
        pvec <- plogis(Z %*% params[(nOP + 1):nP]) * p
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1
        cpmat <- matrix(cp, M, J)
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        loglik <- sum(loglik * weights)
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    o <- optim(inits, nll,
        method = method, hessian = hessian, control = control, ...)
    o
}

summary(t(replicate(100, {
    M <- 1000
    J <- 10

    p <- 0.8 # prob availability
    q <- 0.55 # prob use
    D <- 0.1 # density
    A <- 2^2*pi # effective area

    W <- rbinom(M, 1, 1-exp(-D*A))
    table(W)
    Y <- matrix(0, M, J)
    for (j in 1:J)
        Y[,j] <- rbinom(M, W, p*q)

    X <- matrix(1, M, 1)
    Z <- matrix(1, M*J, 1)

    o <- .mvocc2(Y, X, Z, p, A)
    c(exp(o$par[1]), plogis(o$par[2]))
})))


## CONI specific model
## p: availability from survival model
## lambda: Poisson mean from conditional likelihood step
.mvocc3 <- function(Y, X, Z, p, lambda, weights,
method="Nelder-Mead", inits, control=list(), hessian=FALSE, ...) {

    good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    Y <- ifelse(Y > 0, 1, 0)

    J <- ncol(Y)
    M <- nrow(Y)
    nDP <- ncol(Z)
    nOP <- ncol(X)
    nP <- nDP + nOP

    if (missing(p))
        p <- rep(1, M)
    if (missing(lambda))
        lambda <- rep(1, M)
    if (missing(weights))
        weights <- rep(1, M)
    if (missing(inits))
        inits <- rep(0, nP)

    yvec <- as.numeric(Y)
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(Y, na.rm = TRUE) == 0, 1, 0)
    rowProds <- function (x, na.rm = FALSE)
        exp(rowSums(log(x), na.rm = na.rm))

    nll <- function(params) {
        delta <- plogis(X %*% params[1:nOP])
        psi <- delta * (1-exp(-lambda))
        pvec <- plogis(Z %*% params[(nOP + 1):nP]) * p
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1
        cpmat <- matrix(cp, M, J)
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        loglik <- sum(loglik * weights)
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    o <- optim(inits, nll,
        method = method, hessian = hessian, control = control, ...)
    o
}

summary(t(replicate(100, {
    M <- 1000
    J <- 10

    p <- 0.8 # prob availability
    q <- 0.55 # prob use
    D <- 0.1 # density
    A <- 2^2*pi # effective area
    lambda <- A*D
    delta <- 0.9 # ZI is 1-delta

    N <- rpois(M, lambda) * rbinom(M, 1, delta)
    W <- ifelse(N > 0, 1, 0)
    table(W, N)
    Y <- matrix(0, M, J)
    for (j in 1:J)
        Y[,j] <- rbinom(M, W, p*q)

    X <- matrix(1, M, 1)
    Z <- matrix(1, M*J, 1)

    o <- .mvocc3(Y, X, Z, p, lambda)
    plogis(o$par)
})))
}

library(DEoptim)
## improved version: nice output provided
## CONI specific model
## p: availability from survival model
## lambda: Poisson mean from conditional likelihood step
mvocc <- function(Y, X, Z, p, lambda, weights,
method="Nelder-Mead", init, control=list(), hessian=TRUE, DElimit=10, ...) {

    good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
    .solvenear <- function(x) {
        xinv <- try(solve(x), silent = TRUE)
        if (inherits(xinv, "try-error"))
            xinv <- as.matrix(solve(Matrix::nearPD(x)$mat))
        xinv
    }

    Y <- ifelse(Y > 0, 1, 0)

    J <- ncol(Y)
    M <- nrow(Y)
    nDP <- ncol(Z)
    nOP <- ncol(X)
    nP <- nDP + nOP

    if (missing(p))
        p <- matrix(1, M, J)
    if (missing(lambda))
        lambda <- rep(1, M)
    if (missing(weights))
        weights <- rep(1, M)
    if (missing(init))
        init <- rep(0, nP)

    yvec <- as.numeric(Y)
    p <- as.numeric(p)
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(Y, na.rm = TRUE) == 0, 1, 0)
    rowProds <- function (x, na.rm = FALSE)
        exp(rowSums(log(x), na.rm = na.rm))

    nll <- function(params) {
        delta <- plogis(X %*% params[1:nOP])
        psi <- delta * (1-exp(-lambda))
        pvec <- plogis(Z %*% params[(nOP + 1):nP]) * p
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1
        cpmat <- matrix(cp, M, J)
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        loglik <- sum(loglik * weights)
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    if (method == "DE") {
        up <- rep(DElimit, nP)
        lo <- -up
        opt <- DEoptim(fn=nll, lower=lo, upper=up,
            control=list(trace=FALSE, itermax=nP*200), ...)
        cf <- opt$optim$bestmem
        ll <- -opt$optim$bestval
        if (hessian) {
            hess <- optimHess(opt$optim$bestmem, nll)
            S <- .solvenear(hess)
        } else {
            S <- NULL
        }
    } else {
        opt <- optim(init, nll,
            hessian=hessian, method=method, control = control, ...)
        cf <- opt$par
        ll <- -opt$value
        S <- if (hessian)
            .solvenear(opt$hessian) else NULL
    }

    if (is.null(colnames(X)))
        colnames(X) <- paste0("X", seq_len(ncol(X))-1L)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("Z", seq_len(ncol(Z))-1L)
    colnames(X) <- gsub("[[:punct:]]", "", colnames(X))
    colnames(Z) <- gsub("[[:punct:]]", "", colnames(Z))
    NAM <- c(colnames(X), colnames(Z))
    names(cf) <- NAM
    if (is.null(S))
        S <- matrix(NA, np, nP)
    dimnames(S) <- list(NAM, NAM)

    se <- sqrt(diag(S))
    tstat <- cf/se
    pval <- 2 * pnorm(-abs(tstat))
    coefs <- cbind(cf, se, tstat, pval)
    colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(coefs) <- NAM
    out <- list(coef=cf, loglik=ll, vcov=S, nobs=nrow(Y), df=length(cf),
        summary=coefs, init=init, method=method, lambda=lambda, results=opt,
        nobs=M,
        nvisits=J)
    class(out) <- "mvocc"
    out
}

logLik.mvocc <- function (object, ...)
    structure(object$loglik, df = object$df, nobs=object$nobs, class = "logLik")

coef.mvocc <- function(object, ...)
    object$coef

vcov.mvocc <- function(object, ...)
    object$vcov

summary.mvocc <- function(object, ...)
    printCoefmat(object$summary, signif.legend = TRUE, ...)


