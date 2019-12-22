## plain mv occupancy model
## Y: M x J matrix
## X: M x nOP matrix
## Z: MxJ x nDP matrix (visits stacked)
mvocc <- function(Y, X, Z, offsetx, offsetz, weights,
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

o <- mvocc(Y, X, Z)
plogis(o$par)

## modified mv occupancy model
## Y: M x J matrix
## X: M x nOP matrix
## Z: MxJ x nDP matrix (visits stacked)
mvocc2 <- function(Y, X, Z, p, A, weights,
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

o <- mvocc2(Y, X, Z, p, A)
c(exp(o$par[1]), plogis(o$par[2]))

