## preliminaries for MPI
.Last <- function() {
    if (getOption("CLUSTER_ACTIVE")) {
        stopCluster(cl)
        cat("active cluster stopped by .Last\n")
    } else {
        cat("no active cluster found\n")
    }
}
options("CLUSTER_ACTIVE" = FALSE)
library(snow)
library(Rmpi)
library(rlecuyer)


set.seed(1234)
yr <- 1995:2014
yr0 <- yr[1]
n <- 1000
T <- length(yr)
nh <- 10
ah <- rexp(nh)
ah <- ah/sum(ah)
h <- t(rmultinom(n, 1, ah))
h <- h[,colSums(h) > 0]
x <- apply(h, 1, which.max)
h[,1] <- 1
nh <- ncol(h)
beta <- rnorm(nh, -1, 0.5) # habitat effect
#sigma <- 0#0.5
#tau <- 0#0.5
#theta <- -0.03/10
A <- 1
p <- 1
#sp <- 0.1
B <- 100

## Delta: % yearly change
## sp: spatial sampling density
## an: anchor site proportion
## hb: habitat bias
fun <-
function(Delta=0, sp=1, an=0, hb=FALSE, sigma=0, tau=0)
{
    theta <- log(Delta/100 + 1)

    eps <- rnorm(n, 0, sigma)
    omega <- rnorm(T, 0, tau)

    spm <- matrix(0, n, T)
    D <- Y <- matrix(NA, n, T)
    for (i in 1:T) {
        D[,i] <- exp(drop(h %*% beta) + theta*(yr[i] - yr0) + eps + omega[i])
        Y[,i] <- rpois(n, D[,i] * A * p)
        w <- if (hb)
            rexp(nh) else rep(1, nh)
        spm[sample.int(n, round(n*sp), prob=w[as.integer(x)]),i] <- 1
    }
    if (an > 0) {
        spm[sample.int(n, round(n*an)),] <- 1
    }

    #### point level data
    ## dense
    dat <- data.frame(Y = as.numeric(Y),
        sp=as.logical(spm),
        x=as.factor(x), yr=rep(yr, each=n))
    ## sparse
    Dat <- dat[dat$sp,]

    m1 <- glm(Y ~ x, dat, family="poisson", offset=rep(log(A*p), nrow(dat)))
    m2 <- glm(Y ~ x, Dat, family="poisson", offset=rep(log(A*p), nrow(Dat)))
    f1 <- matrix(exp(drop(h %*% coef(m1))), n, T)
    f2 <- matrix(exp(drop(h %*% coef(m2))), n, T)

    #cbind(truth=beta, m1=coef(m1), m2=coef(m2))

    #### s2d pooled data
    pdat <- data.frame(Y=colSums(Y), yr=yr, A=n*A, sD=colSums(f1))
    pDat <- data.frame(Y=colSums(Y*spm), yr=yr, A=colSums(spm), sD=colSums(f2*spm))

    #with(pdat, plot(yr, Y/A, type="l"))
    #with(pDat, lines(yr, Y/A, col=2))

    #coef(lm(log(Y) ~ yr - 1, pdat, offset=log(pdat$sD*pdat$A*p)))
    #coef(lm(log(Y) ~ yr - 1, pDat, offset=log(pDat$sD*pDat$A*p)))

    mods <- list(
        ## pt level fwd
        d_pnt_fwd = glm(Y ~ yr, dat, family="poisson",
            offset=rep(log(A*p), nrow(dat)) + log(fitted(m1))),
        s_pnt_fwd = glm(Y ~ yr, Dat, family="poisson",
            offset=rep(log(A*p), nrow(Dat)) + log(fitted(m2))),
        ## pt level joint
        d_pnt_jnt = glm(Y ~ x + yr, dat, family="poisson",
            offset=rep(log(A*p), nrow(dat))),
        s_pnt_jnt = glm(Y ~ x + yr, Dat, family="poisson",
            offset=rep(log(A*p), nrow(Dat))),
        ## s2d pooled
        d_pld = glm(Y ~ yr, pdat, family="poisson",
            offset=log(pdat$sD*pdat$A*p)),
        s_pld = glm(Y ~ yr, pDat, family="poisson",
            offset=log(pDat$sD*pDat$A*p)))
    rval <- cbind(Delta=sapply(mods, function(z) coef(z)["yr"]),
        t(sapply(mods, function(z) confint.default(z)["yr",])))
    rval <- 100*(exp(rval) - 1)
    rval <- cbind(rval,
        H=ifelse(sign(rval[,2])*sign(rval[,3]) > 0, sign(rval[,1]), 0))
    rval
}


bfun <-
function(Delta=0, sp=1, an=0, hb=FALSE, B=100, ...)
{
    vv <- lapply(seq_len(B),
        function(i, ...) fun(Delta=Delta, sp=sp, an=an, hb=hb, ...))
    Del <- sapply(vv, function(z) z[,"Delta"])
    Pow <- rowMeans(ifelse(sapply(vv,
        function(z) z[,"H"]) != 0, 1, 0))
    Pow2 <- rowMeans(ifelse(sapply(vv,
        function(z) z[,"H"]) == sign(Delta), 1, 0))
    MSE <- rowMeans((Del - Delta)^2)
    #Var <- rowMeans((Del - rowMeans(Del))^2)
    Bias <- rowMeans(Del - Delta)
    Var <- MSE - Bias^2
    out <- cbind(Mean=rowMeans(Del), Pow=Pow, Pow2=Pow2, MSE=MSE, Bias=Bias, Var=Var)
    attr(out, "settings") <- c(Delta=Delta, sp=sp, hb=as.integer(hb))
    out
}

vals1 <- expand.grid(
    Delta=c(-10, -5, -3, -1, 0, 1, 3, 5, 10),
    sp=c(0.1, 0.5, 0.8),
    an=c(0.1, 0.5, 0.8),
    tau=0,
    sigma=0,
    hb=c(FALSE, TRUE))
vals2 <- expand.grid(
    Delta=c(-10, -5, -3, -1, 0, 1, 3, 5, 10),
    sp=1,
    an=0.3,
    tau=c(0, 0.5, 1),
    sigma=c(0, 0.5, 1),
    hb=FALSE)

vals <- vals2

#system.time(aa <- bfun(Delta=vals[1,1], sp=vals[1,2],
#    an=vals[1,3], hb=vals[1,4], B=10))

(args <- commandArgs(trailingOnly = TRUE))
nodes <- as.numeric(args[1])
ncl <- nodes * 12

## parallel stuff
cl <- makeMPIcluster(ncl)
options("CLUSTER_ACTIVE" = TRUE)

clusterSetupRNG (cl, type = "RNGstream")
clusterExport(cl, ls())

res <- parLapply(cl, 1:nrow(vals),
    function(i) bfun(Delta=vals[i,"Delta"], sp=vals[i,"sp"],
        an=vals[i,"an"], hb=vals[i,"hb"], B=B, 
        sigma=vals[i,"sigma"], tau=vals[i,"tau"]))

#res <- lapply(1:2,
#    function(i) bfun(Delta=vals[i,1], sp=vals[i,2],
#        an=vals[i,3], hb=vals[i,4], B=2))

## shutting down safely
stopCluster(cl)
options("CLUSTER_ACTIVE" = FALSE)

## save results
save.image(file="glmm-res3-sigmatau.Rdata")

mpi.quit("no")


## plots

## graphics
library(lattice)
load("~/Dropbox/bam/full_life_cycle/glmm-res3.Rdata")
rm(.Last)

pow <- t(sapply(res, function(z) z[,"Pow"]))

vals$Beta <- log(vals$Delta / 100 + 1)
vals$Decadal <- 100 * (exp(vals$Beta * 10) - 1)

op <- par(mfrow=c(3,3), las=1, mar=c(4,4,2,1)+0.1)
for (i in 1:3) {
    for (j in 1:3) {
        sp <- unique(vals$sp)[i]
        an <- unique(vals$an)[j]
        ss <- pow[!vals$hb & vals$sp == sp & vals$an == an,]
        matplot(sort(unique(vals$Delta)), ss[,5:6],
            type="l", ylim=c(0,1),
            xlab="", ylab="",
            main=paste0("sp=", sp, ", an=", an),
            lty=1, lwd=2, col=c(2,4))
        if (j==1 && i==2)
            title(ylab="Power")
        if (i==3 && j==2)
            title(xlab="Annual trend (%)")
    }
}
legend("bottomright", col=c(2,4), lty=1, lwd=2,
       legend=c("dense","sparse"), bty="n")
par(op)



load("~/Dropbox/bam/full_life_cycle/glmm-res3-sigmatau.Rdata")

pow <- t(sapply(res, function(z) z[,"Pow"]))

vals$Beta <- log(vals$Delta / 100 + 1)
vals$Decadal <- 100 * (exp(vals$Beta * 10) - 1)

op <- par(mfrow=c(3,3), las=1, mar=c(4,4,2,1)+0.1)
for (i in 1:3) {
    for (j in 1:3) {
        sigma <- unique(vals$sigma)[i]
        tau <- unique(vals$tau)[j]
        ss <- pow[!vals$hb & vals$sigma == sigma & vals$tau == tau,]
        matplot(sort(unique(vals$Delta)), ss[,5:6],
            type="l", ylim=c(0,1),
            xlab="", ylab="",
            main=paste0("sig=", sigma, ", tau=", tau),
            lty=1, lwd=2, col=c(2,4))
        if (j==1 && i==2)
            title(ylab="Power")
        if (i==3 && j==2)
            title(xlab="Annual trend (%)")
    }
}
legend("bottomright", col=c(2,4), lty=1, lwd=2,
       legend=c("dense","sparse"), bty="n")
par(op)



i <- 4
(Which <- rownames(res[[1]])[i])
dat <- data.frame(vals, t(sapply(res, function(z) z[Which,])))

with(dat, plot(Delta, Mean, cex=0.1+Pow))
abline(0,1)

contourplot(Bias ~ sp * Delta | hb, data = dat)

Dat <- lapply(1:6, function(j)
    data.frame(vals, t(sapply(res, function(z) z[j,]))))

j <- 5
k <- "Pow"
par(mfrow=c(1,3))
spx <- 0.1
xm <- cbind(Dense_Rand=Dat[[j]][vals$sp==spx & !vals$hb,k],
            Sparse_Rand=Dat[[j+1]][vals$sp==spx & !vals$hb,k],
            Dense_Pref=Dat[[j]][vals$sp==spx & vals$hb,k],
            Sparse_Pref=Dat[[j+1]][vals$sp==spx & vals$hb,k])
matplot(unique(vals$Delta), xm, col=c(4,4,2,2), lty=c(1,2,1,2), type="l", lwd=2,
        xlab="True trend", ylab="P(reject H0 | H1 true)")
spx <- 0.5
xm <- cbind(Dense_Rand=Dat[[j]][vals$sp==spx & !vals$hb,k],
            Sparse_Rand=Dat[[j+1]][vals$sp==spx & !vals$hb,k],
            Dense_Pref=Dat[[j]][vals$sp==spx & vals$hb,k],
            Sparse_Pref=Dat[[j+1]][vals$sp==spx & vals$hb,k])
matplot(unique(vals$Delta), xm, col=c(4,4,2,2), lty=c(1,2,1,2), type="l", lwd=2,
        xlab="True trend", ylab="P(reject H0 | H1 true)")
spx <- 0.8
xm <- cbind(Dense_Rand=Dat[[j]][vals$sp==spx & !vals$hb,k],
            Sparse_Rand=Dat[[j+1]][vals$sp==spx & !vals$hb,k],
            Dense_Pref=Dat[[j]][vals$sp==spx & vals$hb,k],
            Sparse_Pref=Dat[[j+1]][vals$sp==spx & vals$hb,k])
matplot(unique(vals$Delta), xm, col=c(4,4,2,2), lty=c(1,2,1,2), type="l", lwd=2,
        xlab="True trend", ylab="P(reject H0 | H1 true)")



par(mfrow=c(3,2))
for (j in 1:6) {
    with(Dat[[j]][vals$sp==0.1,], plot(Delta, Pow, cex=1+Pow))
    #abline(0,1)
    abline(h=0)
}



