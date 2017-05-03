library(knitr)
library(mefa4)
library(pbapply)
library(pROC)
ROOT <- "e:/peter/bam/Apr2016/out"
source("~/repos/bamanalytics/R/makingsense_functions.R")

PROJECT <- "bam"
Date <- "2016-12-01"
level <- 0.9

fstat <- function(x, level=0.95) {
    out <- quantile(x, c(0.5, (1-level)/2, 1 - (1-level)/2))
    names(out) <- c("Median", "LCL", "UCL")
    out
}
chfun <- function(Na, Nb, ta, tb) {
    100 * ((Nb/Na)^(1/(tb-ta)) - 1)
}

e <- new.env()
load(file.path(ROOT, "data", "pack_2016-12-01.Rdata"), envir=e)
#load(file.path(ROOT, "data", "pack_2017-04-19.Rdata"), envir=e)

mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,c(Terms, "Units")]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
xn <- xn[rownames(Xn),]
off <- e$OFF[rownames(xn),]
bbb <- unique(e$BB)
bb <- e$BB

rm(e)
INTERNAL <- 1:nrow(xn) %in% bbb
ss1 <- which(!INTERNAL) # test portion
ss2 <- which(INTERNAL) # non-test portion

modTab <- getFancyModsTab(mods)
xnh <- nonDuplicated(xn, HABTR, TRUE)[,c("HAB","HABTR","isNF","isDev",
    "isWet","isOpn","isDM","isDec","isMix")]
xnh <- xnh[c("ConifDense", "ConifSparse","ConifOpen",
    "DecidDense", "DecidSparse", "DecidOpen",
    "MixedDense", "MixedSparse", "MixedOpen",
    "WetDense", "WetSparse", "WetOpen",
    "Shrub", "Grass", "Barren", "Agr", "Devel"),]

spp <- "CAWA"

Y <- yy[,spp]
Y1 <- ifelse(yy[,spp]>0, 1, 0)
off1 <- off[,spp]

fn <- file.path(ROOT, "results", "cawa",
    paste0(PROJECT, "_", spp, "_", Date, ".Rdata"))
load(fn)
#100 * sum(getOK(res)) / length(res)
est_hab <- getEst(res, stage = 2, X=Xn)
est_habhgt <- getEst(res, stage = 3, X=Xn)
est_dtb <- getEst(res, stage = 4, X=Xn)
est_wet <- getEst(res, stage = 5, X=Xn)
est <- getEst(res, stage = length(mods)-1, X=Xn)
est_yr <- getEst(res, stage = length(mods), X=Xn)
#dcawa <- data.frame(PKEY=rownames(xn), cawa=Y)
#dcawa$validation <- 0
#dcawa$validation[ss1] <- 1
#write.csv(dcawa, row.names=FALSE, file="w:/bam-cawa/cawa-pkeys-2016-12-13.csv")

## pseudo R^2

## predictions for internal points
mn_in <- matrix(0, length(ss2), length(mods)+1)
colnames(mn_in) <- c("NULL", names(mods))
rownames(mn_in) <- rownames(xn)[ss2]
for (i in 0:length(mods)) {
    est_i <- getEst(res, stage = i, X=Xn)
    col_keep <- colSums(abs(est_i) > 0) != 0
    pr <- exp(sapply(1:nrow(est_i), function(j)
        Xn[ss2,colnames(est_i[,col_keep,drop=FALSE]),drop=FALSE] %*%
        est_i[j,col_keep]))
    mn_in[,i+1] <- rowMeans(pr)
}

## POISSON
## Null model: intercept and offset
ll0 <- sum(dpois(Y[ss2], mn_in[,"NULL"] * exp(off1[ss2]), log=TRUE))
## Saturated: one parameter per observation
lls <- sum(dpois(Y[ss2], Y[ss2], log=TRUE))
## Full: our smoothed prediction
llf <- sapply((1:length(mods))+1, function(i)
    sum(dpois(Y[ss2], mn_in[,i] * exp(off1[ss2]), log=TRUE)))
R2D <- 1-(lls-llf)/(lls-ll0)
names(R2D) <- colnames(mn_in)[-1]
plot(R2D,type="b",ylim=c(0,1))
mean(Y[ss2])
colMeans(mn_in)
colMeans(mn_in)/mean(Y[ss2])

## BINOMIAL
## Null model: intercept and offset
ll0 <- sum(dbinom(Y1[ss2], 1, 1-exp(-(mn_in[,"NULL"] * exp(off1[ss2]))), log=TRUE))
## Saturated: one parameter per observation
lls <- sum(dbinom(Y1[ss2], 1, Y1[ss2], log=TRUE)) # should be 0
## Full: our smoothed prediction
llf <- sapply((1:length(mods))+1, function(i)
    sum(dbinom(Y1[ss2], 1, 1-exp(-(mn_in[,i] * exp(off1[ss2]))), log=TRUE)))
R2Db <- 1-(lls-llf)/(lls-ll0)


## ROC curves and AUC

## predictions for external points
mn <- matrix(0, length(ss1), length(mods)+1)
colnames(mn) <- c("NULL", names(mods))
rownames(mn) <- rownames(xn)[ss1]
for (i in 0:length(mods)) {
    est_i <- getEst(res, stage = i, X=Xn)
    col_keep <- colSums(abs(est_i) > 0) != 0
    pr <- exp(sapply(1:nrow(est_i), function(j)
        Xn[ss1,colnames(est_i[,col_keep,drop=FALSE]),drop=FALSE] %*%
        est_i[j,col_keep]))
    mn[,i+1] <- rowMeans(pr)
}
## ROC/AUC
rocAll1 <- lapply(1:ncol(mn), function(i) {
      pp <- mn[,i] * exp(off1[ss1])
      roc(Y1[ss1], pp)
  })
names(rocAll1) <- c("NULL", names(mods))
auc <- sapply(rocAll1, function(z) as.numeric(z$auc))
k <- (auc-auc[1])/(1-auc[1])

op <- par(mfrow=c(1,2))
tmp <- barplot(k, ylim=c(0,1), space=0.2, ylab="AUC", xlab="Stages",
    col="lightgrey", border=NA, width=1)
text(tmp, k, round(k, 3), col=4, cex=0.75, pos=3)
#tmp <- barplot(auc, ylim=c(0,1), space=0.2, ylab="AUC", xlab="Stages",
#    col="lightgrey", border=NA, width=1)
#text(tmp, auc, round(auc, 3), col=4, cex=0.75)
plot(rocAll1[["NULL"]], col=4, lty=2)
lines(rocAll1[["Clim"]], col=4)
legend("bottomright", col=4, lty=c(2,1), legend=c("NULL", "Clim"), bty="n")
par(op)

## regional ROC
## need to use all the data: only 2 regions with enough validation points

table(xn$Units,Y1)
ftable(INTERNAL, xn$Units, Y1)
# AK: Alaska, BCR 2+4
# N: Arctic/Shield, BCR3+7
# Mt: Rockies, BCR 5+9+10
# W: West, BCR 6+11 and BCR 8 in MB & SK
# 8E: BCR 8 east (not MB or SK)
# HW: Hardwood, BCR 12+13+23
# At: Atlantic, BCR 14

i <- 6
rocReg <- lapply(c("W", "8E", "HW", "At"), function(REG) {
      ppp <- c(mn_in[,i] * exp(off1[ss2]), mn[,i] * exp(off1[ss1]))
      yyy <- c(Y1[ss2], Y1[ss1])
      sss <- xn$Units[c(ss2, ss1)] == REG
      roc(yyy[sss], ppp[sss])
  })
names(rocReg) <- c("W", "8E", "HW", "At")

sapply(rocReg, function(z) as.numeric(z$auc))

plot(rocAll1[["Clim"]])
lines(rocReg[["W"]], col=2)
lines(rocReg[["8E"]], col=3)
lines(rocReg[["HW"]], col=4)
lines(rocReg[["At"]], col=5)


## QQ plots
vals <- seq(min(Y[ss1]), max(Y[ss1]), 1)
pobs <- numeric(length(vals))
names(pobs) <- vals
tab <- table(Y[ss1])
for (k in vals) {
    pobs[as.character(k)] <- if (as.character(k) %in% names(tab))
        tab[as.character(k)]/length(ss1) else 0
}
pobs <- cumsum(pobs)

pexp <- matrix(0, length(vals), ncol(mn))
rownames(pexp) <- vals
colnames(pexp) <- colnames(mn)
for (i in 1:ncol(mn)) {
    m_pexp <- matrix(0, length(ss1), length(vals))
    colnames(m_pexp) <- vals
    for (k in vals) {
        m_pexp[,as.character(k)] <- dpois(x=rep(k, length(ss1)),
            lambda=mn[,i] * exp(off1[ss1]))
    }
    m_pexp[,ncol(m_pexp)] <- 1 - rowSums(m_pexp[,-ncol(m_pexp)])
    pexp[,i] <- cumsum(colMeans(m_pexp))
}

p_min <- min(cbind(pobs,pexp))
op <- par(mfrow=c(1,2))
plot(pobs, pexp[,"NULL"], type="b", pch=rownames(pexp), main="NULL",
    xlab="Observed", ylab="Expected", xlim=c(p_min, 1), ylim=c(p_min, 1), col=4)
abline(0, 1, col="grey")
plot(pobs, pexp[,"Clim"], type="b", pch=rownames(pexp), main="Clim",
    xlab="Observed", ylab="Expected", xlim=c(p_min, 1), ylim=c(p_min, 1), col=4)
abline(0, 1, col="grey")
par(op)

## Ranking

pro <- table(Y[ss1])/length(ss1)
op <- par(mfrow=c(1,2))
tmp <- boxplot(mn[,"NULL"] * exp(off1[ss1]) ~ Y[ss1], range=0,
    at=cumsum(2*pro^0.2), width=2*pro^0.2, main="NULL",
    xlab="Observed count", ylab="Corrected density", col="grey",
    ylim=c(0, max(mn)))
boxplot(mn[,"Clim"] * exp(off1[ss1]) ~ Y[ss1], range=0,
    at=cumsum(2*pro^0.2), width=2*pro^0.2, main="Clim",
    xlab="Observed count", ylab="Corrected density", col="grey",
    ylim=c(0, max(mn)))
par(op)


op <- par(mfrow=c(1,2))
ResourceSelection:::.mep(Y[ss1], mn[,"NULL"] * exp(off1[ss1]),
    link="log", type="unique", level=0.9, main="NULL",
    xlab="Observed count", ylab="Corrected density", ylim=c(0,max(mn)), xlim=c(-0.2,4.2))
ResourceSelection:::.mep(Y[ss1], mn[,"Clim"] * exp(off1[ss1]),
    link="log", type="unique", level=0.9, main="Clim",
    xlab="Observed count", ylab="Corrected density", ylim=c(0,max(mn)), xlim=c(-0.2,4.2))
par(op)

#op <- par(mfrow=c(1,2))
#ResourceSelection:::.mep(Y[ss1], mn[,"NULL"] * exp(off1[ss1]),
#    level = 0.9, link = "log", type = "unique",
#    main="NULL", xlim=c(0,4), ylim=c(0, max(mn)))
#ResourceSelection:::.mep(Y[ss1], mn[,"Clim"] * exp(off1[ss1]),
#    level = 0.9, link = "log", type = "unique",
#    main="Clim", xlim=c(0,4), ylim=c(0, max(mn)))
#par(op)

## include OCCC

