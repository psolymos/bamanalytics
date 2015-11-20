## simulate the data for single location

simulateData <- 
function(r, d, sigma, D0, T, A=1, p=1, seed = NULL) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    T <- T[1]
    if (T < 2)
        stop("T<2 not allowed")
    if (any(A < 0))
        stop("A<0 not allowed")
    if (any(p < 0))
        stop("p<0 not allowed")
    if (any(p > 1))
        stop("p>1 not allowed")
    if (any(d <= 0))
        stop("d<=0 not allowed")
    r <- rep(r, T)[1:T]
    d <- rep(d, T)[1:T]
    A <- rep(A, T)[1:T]
    p <- rep(p, T)[1:T]
    #N <- numeric(T)
    #N[1] <- D0 * A[1]
    #D <- numeric(T)
    #D[1] <- D0
    #x <- log(N)
    logD <- numeric(T)
    logD[1] <- log(D0)
    for(t in 2:T) {
        logD[t] <- logD[t-1] + log(A[t-1]) - log(A[t]) +
            r[t] * (1 - exp(logD[t-1]) / d[t-1]) + rnorm(1, 0, sigma)
        #x[t] <- x[t-1] + r[t] * (1 - exp(x[t-1]) / K[t-1]) + rnorm(1, 0, sigma)
    }
    D <- exp(logD)
    N <- D * A
    Y <- rpois(T, N*p)
    out <- data.frame(year=1:T, Y=Y, N=N, K=d*A, A=A, p=p, r=r, d=d, sigma=sigma)
    class(out) <- c("ricker", "data.frame")
    out
}

x <- simulateData(r=-0.08, d=7.23, sigma=0.66, D0=1, A=1, p=1, T=25)
plot(as.ts(x$Y), col=2)
lines(as.ts(x$N), col=4)

#xx <- replicate(10, simulateData(r=0.5, d=5, sigma=0.1, D0=1, A=1, p=1, T=25), simplify=FALSE)
xx <- replicate(10, simulateData(r=-0.1, d=5, sigma=0.2, D0=1, A=1, p=1, T=25), simplify=FALSE)
xxy <- rowSums(sapply(xx, "[[", "Y"))
xxn <- rowSums(sapply(xx, "[[", "N"))
xxk <- rowSums(sapply(xx, "[[", "K"))

plot(as.ts(xxy), col=2)
lines(as.ts(xxn), col=4)


## Ricker logistic
library(dclone)
model1 <- structure(
    c(" model ",
    " { ",
    "     logD[1] <- log(Y[1]/(p[1]*A[1])) ",
    "     for (t in 2:T) { ",
    "         Y[t] ~ dpois(exp(logD[t])*A[t]*p[t]) ",
    "         logD[t] ~ dnorm(mu[t], 1/sigma^2) ",
    "         mu[t] <-  logD[t-1] + log(A[t-1]) - log(A[t]) + ",
    "             r * (1 - exp(logD[t-1]) / d) ",
    "     } ",
    "     r ~ dnorm(0.00000E+00, 0.001) ",
    "     logd ~ dnorm(0.00000E+00, 0.1) ",
    "     logsigma ~ dnorm(0.00000E+00, 1) ",
    "     sigma <- exp(logsigma) ",
    "     d <- exp(logd) ",
    " } "), class = "custommodel")
dat <- list(Y=xxy, T=length(xxy), A=rep(1,length(xxy)), p=rep(1,length(xxy)))
m <- jags.fit(dat, c("r","d","sigma"), model1)
coef(m)

## covariate effects on r

T <- 25
#clim <- sin(seq(-5,5,len=T))
#r <- 0.5 + 0.5*clim
r <- 0.5
x <- simulateData(r=r, d=5, sigma=0.1, D0=4, A=1, p=1, T=T)
plot(as.ts(x$Y), col=2)
lines(as.ts(x$N), col=4)

xx <- replicate(10, simulateData(r=r, d=5, sigma=0.1, D0=1, A=1, p=1, T=T), simplify=FALSE)
xxy <- rowSums(sapply(xx, "[[", "Y"))
xxn <- rowSums(sapply(xx, "[[", "N"))
xxk <- rowSums(sapply(xx, "[[", "K"))

plot(as.ts(xxy), col=2)
lines(as.ts(xxn), col=4)


model2 <- structure(
    c(" model ",
    " { ",
    "     logD[1] <- log(Y[1]/(p[1]*A[1])) ",
#    "     r[1] <- r0 + r1*x[1] ",
    "     for (t in 2:T) { ",
    "         Y[t] ~ dpois(exp(logD[t])*A[t]*p[t]) ",
    "         logD[t] ~ dnorm(mu[t], 1/sigma^2) ",
    "         mu[t] <-  logD[t-1] + log(A[t-1]) - log(A[t]) + ",
    "             r * (1 - exp(logD[t-1]) / exp(logd[t-1])) ",
#    "         r[t] ~ dnorm(r0 + r1*x[t], 1/taur^2) ",
#    "         r[t] <- r0 + r1*x[t] ",
#    "         logd[t-1] <- min(logk + loglambda[t-1], log(Kmax)) ",
    "         logd[t-1] <- logk + loglambda[t-1] ",
    "     } ",
    "     r ~ dnorm(0.00000E+00, 0.001) ",
#    "     r1 ~ dnorm(0.00000E+00, 0.001) ",
    "     logk ~ dnorm(0.00000E+00, 0.001) ",
    "     logsigma ~ dnorm(0.00000E+00, 1) ",
    "     sigma <- exp(logsigma) ",
    "     k <- exp(logk) ",
#"              N[j] <- min(exp(x[j]), 10000) "
    " } "), class = "custommodel")
dat <- list(Y=xxy, T=length(xxy), #Kmax=10^4,
    A=rep(1,length(xxy)), p=rep(1,length(xxy)),
    #x=clim, 
    loglambda=rep(0, length(xxy)))
m <- jags.fit(dat, c("r","k","logk","logd","sigma"), model2)
#m <- jags.fit(dat, "logd", model2)
coef(m)

## problems with estimating K when r < 0 (decline)
## possibly keep only change thing without density dependence

## covariate for r is hard stuff
## - check if this works w/o density dependence?

## change in HS might also be hard to quantify
## (except maybe for forest harvest)
## this requires a carefully selected data set (CL?)

## spatially varying K

#r_fixed <- 0.5
#D0_fixed <- 1
#T_fixed <- 25

r_fixed <- -0.1
D0_fixed <- 300
d_fixed <- 500
T_fixed <- 25

x <- simulateData(r=r_fixed, d=d_fixed, sigma=0.1, D0=D0_fixed, A=1, p=1, T=T_fixed)
plot(as.ts(x$Y), col=2)
lines(as.ts(x$N), col=4)

dd <- seq(1, 5, len=10)
xx <- lapply(dd, function(d) {
  simulateData(r=r_fixed, d=d_fixed, sigma=0.1, D0=D0_fixed, A=1, p=1, T=T_fixed)
  })
xxy <- rowSums(sapply(xx, "[[", "Y"))
xxn <- rowSums(sapply(xx, "[[", "N"))
xxk <- rowSums(sapply(xx, "[[", "K"))

plot(as.ts(xxy), col=2)
lines(as.ts(xxn), col=4)

## dense
dat <- list(Y=xxy, T=length(xxy), A=rep(1,length(xxy)), p=rep(1,length(xxy)))
m1 <- jags.fit(dat, c("r","logd","sigma"), model1)
coef(m1)

## spatially varying K and resampling

zz <- t(replicate(nrow(xx[[1]]), rbinom(10, 1, 0.5)))
xxy2 <- rowSums(sapply(xx, "[[", "Y")*zz)
xxn2 <- rowSums(sapply(xx, "[[", "N")*zz)
xxk2 <- rowSums(sapply(xx, "[[", "K")*zz)

aaa <- rowSums(zz)
ddd <- matrix(dd, nrow(zz), ncol(zz), byrow=TRUE)
ddd <- rowSums(ddd*zz) / aaa


plot(as.ts(xxy2/aaa), col=2, lwd=2)
lines(as.ts(xxn2/aaa), col=4, lwd=2)
lines(as.ts(xxy/10), col=2, lty=2)
lines(as.ts(xxn/10), col=4, lty=2)
lines(as.ts(10*ddd), col=3)

## standardized sparse
dat <- list(Y=round(xxy2/aaa), T=length(xxy2), A=rep(1,length(xxy2)), p=rep(1,length(xxy2)))
#dat <- list(Y=xxy2, T=length(xxy2), A=aaa, p=rep(1,length(xxy2)))
m2 <- jags.fit(dat, c("r","logd","sigma"), model1)
coef(m2)

model2 <- structure(
    c(" model ",
    " { ",
    "     logD[1] <- log(Y[1]/(p[1]*A[1])) ",
#    "     r[1] <- r0 + r1*x[1] ",
    "     for (t in 2:T) { ",
    "         Y[t] ~ dpois(exp(logD[t])*A[t]*p[t]) ",
    "         logD[t] ~ dnorm(mu[t], 1/sigma^2) ",
    "         mu[t] <-  logD[t-1] + log(A[t-1]) - log(A[t]) + ",
    "             r * (1 - exp(logD[t-1]) / exp(logd[t-1])) ",
#    "         r[t] ~ dnorm(r0 + r1*x[t], 1/taur^2) ",
#    "         r[t] <- r0 + r1*x[t] ",
#    "         logd[t-1] <- min(logk + loglambda[t-1], log(Kmax)) ",
    "         logd[t-1] <- logk + loglambda[t-1] ",
    "     } ",
    "     r ~ dnorm(0.00000E+00, 0.001) ",
#    "     r1 ~ dnorm(0.00000E+00, 0.001) ",
    "     logk ~ dnorm(0.00000E+00, 1) ",
    "     logsigma ~ dnorm(0.00000E+00, 1) ",
    "     sigma <- exp(logsigma) ",
    "     k <- exp(logk) ",
#"              N[j] <- min(exp(x[j]), 10000) "
    " } "), class = "custommodel")

## model based sparse
dat <- list(Y=xxy2, T=length(xxy2), #Kmax=10^4,
    A=aaa, p=rep(1,length(xxy2)),
    #x=clim, 
    loglambda=log(ddd/2))
m3 <- jags.fit(dat, c("r","k","sigma"), model2)
#m <- jags.fit(dat, "logd", model2)
coef(m3)

cbind(coef(m1),coef(m2),coef(m3))

## explore DC, prior on K makes a lot of difference

## 

## standardized sparse for BAM data

load("c:/bam/Feb2014/data_ec_2014.Rdata")
library(mefa4)
d <- DAT[,c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND","YEAR",
    "BCR","PROV","CAWA","QPAD_CAWA")]
d$bcrprov <- interaction(d$BCR, d$PROV, drop=TRUE)
rm(DAT, xy)
d$BBS <- substr(as.character(d$PCODE),1,3) == "BBS"
d1 <- droplevels(d[d$BBS,])
d2 <- droplevels(d[!d$BBS,])

#d11 <- d1[d1$bcrprov=="6.AB",]
#d21 <- d2[d2$bcrprov=="6.AB",]

y1 <- Xtab(CAWA~bcrprov+YEAR,d1)
a1 <- Xtab(exp(QPAD_CAWA)~bcrprov+YEAR,d1)
y2 <- Xtab(CAWA~bcrprov+YEAR,d2, cdrop="2012")
a2 <- Xtab(exp(QPAD_CAWA)~bcrprov+YEAR,d2, cdrop="2012")

y11 <- as.matrix(y1)["6.AB",]
a11 <- as.matrix(a1)["6.AB",]
y21 <- as.matrix(y2)["6.AB",]
a21 <- as.matrix(a2)["6.AB",]

y12 <- colSums(y1)
a12 <- colSums(a1)
y22 <- colSums(y2)
a22 <- colSums(a2)

plot(as.ts(y11/a11))

model3 <- structure(
    c(" model ",
    " { ",
    "     logD[1] <- log(Y[1]/(p[1]*A[1])) ",
    "     for (t in 2:T) { ",
    "         Y[t] ~ dpois(exp(logD[t])*A[t]*p[t]) ",
    "         logD[t] ~ dnorm(mu[t], 1/sigma^2) ",
    "         mu[t] <-  logD[t-1] + log(A[t-1]) - log(A[t]) + ",
    "             r * (1 - exp(logD[t-1]) / d) ",
    "     } ",
    "     r ~ dnorm(0.00000E+00, 0.001) ",
    "     logd ~ dnorm(0.00000E+00, 0.1) ",
    "     logsigma ~ dnorm(0.00000E+00, 1) ",
    "     sigma <- exp(logsigma) ",
    "     d <- exp(logd) ",
    " } "), class = "custommodel")

fun <- function(y, a) {
    dat <- list(Y=round(y), T=length(y), A=a, p=rep(1,length(y)))
    jags.fit(dat, c("r","d","sigma"), model3)
}

m11 <- fun(y11,a11)
m21 <- fun(y21,a21)
m12 <- fun(y12,a12)
m22 <- fun(y22,a22)

cbind(coef(m11),coef(m21),coef(m12),coef(m22))

par(mfrow=c(2,2))
plot(as.ts(y11/a11), main="BBS, BCR6-AB", ylab="Y/A")
plot(as.ts(y21/a21), main="BAM, BCR6-AB", ylab="Y/A")
plot(as.ts(y12/a12), main="BBS, Canada", ylab="Y/A")
plot(as.ts(y22/a22), main="BAM, Canada", ylab="Y/A")


library(PVAClone)
m11 <- pva(y11/a11, gompertz("normal"), c(5,10))
m21 <- pva(ifelse(y21/a21==0, NA, y21/a21), gompertz("normal"), c(5,10))
m12 <- pva(y12/a12, gompertz("normal"), c(5,10))
m22 <- pva(y22/a22, gompertz("normal"), c(5,10))

cl <- makeCluster(3)
k <- c(5,10,50)
m11 <- pva(y11/a11, ricker("normal"), k, cl=cl)
m21 <- pva(ifelse(y21/a21==0, NA, y21/a21), ricker("normal"), k, cl=cl)
m12 <- pva(y12/a12, ricker("normal"), k, cl=cl)
m22 <- pva(y22/a22, ricker("normal"), k, cl=cl)

cbind(coef(m11),coef(m21),coef(m12),coef(m22))


## using the forecast package

library(forecast)

r_fixed <- -0.05
D0_fixed <- 300
d_fixed <- 500
T_fixed <- 25
sigma_fixed <- 0.1
nloc <- 1000


## dense data, r, K
xx <- replicate(nloc, simulateData(r=r_fixed, d=d_fixed, sigma=sigma_fixed, 
    D0=D0_fixed, A=1, p=1, T=T_fixed), simplify=FALSE)
xxy <- rowSums(sapply(xx, "[[", "Y"))
xxn <- rowSums(sapply(xx, "[[", "N"))
xxk <- rowSums(sapply(xx, "[[", "K"))
plot(as.ts(xxy))

(drK <- auto.arima(xxy))
(drK <- ets(log(xxy)))
plot(forecast(drK))

## sparse data, r, K
zz <- t(replicate(nrow(xx[[1]]), rbinom(nloc, 1, 0.5)))
xxy2 <- rowSums(sapply(xx, "[[", "Y")*zz)
xxn2 <- rowSums(sapply(xx, "[[", "N")*zz)
xxk2 <- rowSums(sapply(xx, "[[", "K")*zz)
plot(as.ts(xxy2))

(srK <- auto.arima(xxy2))
(srK <- ets(log(xxy2)))
plot(forecast(srK))

## dense data, r, Ki
dd <- seq(100, 300, len=nloc)
xx2 <- lapply(dd, function(d) {
  simulateData(r=r_fixed, d=d_fixed, sigma=sigma_fixed, 
        D0=D0_fixed, A=1, p=1, T=T_fixed)
  })
xxy3 <- rowSums(sapply(xx2, "[[", "Y"))
xxn3 <- rowSums(sapply(xx2, "[[", "N"))
xxk3 <- rowSums(sapply(xx2, "[[", "K"))
plot(as.ts(xxy3))

(drKi <- auto.arima(xxy3))
(drKi <- ets(log(xxy3)))
plot(forecast(drKi))

## sparse data, r, Ki
xxy4 <- rowSums(sapply(xx2, "[[", "Y")*zz)
xxn4 <- rowSums(sapply(xx2, "[[", "N")*zz)
xxk4 <- rowSums(sapply(xx2, "[[", "K")*zz)
plot(as.ts(xxy4))

(srKi <- auto.arima(xxy4))
(srKi <- ets(log(xxy4)))
plot(forecast(srKi))

