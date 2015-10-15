##---
##title: "QPAD estimation"
##author: "Peter Solymos"
##date: "May 7, 2015"
##output:
##  pdf_document:
##    toc: true
##    toc_depth: 2
##---

## this is modified version for estimating leave-one-out XV for projects

### Preliminaries

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"

## Load required packages
library(mefa4)
library(pbapply)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2015-10-08.Rdata"))

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","PCODE","X","Y")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
sort(table(pkDur$PCODE))
## only 1 obs in project (all others are >20)
pkDur <- droplevels(pkDur[pkDur$PCODE != "GLDRPCCL01",])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])

if (FALSE) {
library(rworldmap)
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(pkDur[,c("X","Y")], pch=19,
    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255), cex=0.2)

}

## models to consider
NAMES <- list("0"="INTERCEPT")
ff <- list(~ 1)

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL


## define pcode, and excl (excl=T to exclude pcode, excl=F to use only pcode
fitDurFun2 <- function(spp, fit=TRUE, type=c("rem","mix"), pcode=NULL, excl=TRUE) {
    rn <- intersect(rownames(pkDur), rownames(xtDur[[spp]]))
    X0 <- pkDur[rn,]
    if (!is.null(pcode)) {
        X0 <- if (excl)
            X0[X0$PCODE != pcode,,drop=FALSE] else X0[X0$PCODE == pcode,,drop=FALSE]
    }
    rn <- rownames(X0)
    Y0 <- as.matrix(xtDur[[spp]][rn,])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdur$x)))
    stopifnot(length(setdiff(levels(X0$DURMETH), rownames(ltdur$end))) == 0)
    ## interval end matrix
    D <- ltdur$end[match(X0$DURMETH, rownames(ltdur$end)),]
    ## exclude 0 sum and <1 interval rows
    iob <- rowSums(Y0) > 0 & rowSums(!is.na(D)) > 1
    if (sum(iob)==0)
        return(structure("0 observation with multiple duration (1)",
            class="try-error"))
    if (sum(iob)==1)
        return(structure("1 observation with multiple duration (2)",
            class="try-error"))
    X <- droplevels(X0[iob,])
    Y0 <- Y0[iob,]
    D <- D[iob,]
    n <- nrow(D)
    ## arranging counts into position
    Yid <- ltdur$id[match(X$DURMETH, rownames(ltdur$id)),]
    Y <- matrix(NA, n, ncol(ltdur$id))
    for (i in seq_len(n)) {
        w <- Yid[i,]
        w <- w[!is.na(w)]
        Y[i,seq_len(length(w))] <- Y0[i,w]
    }
    if (fit) {
        res <- list()
        for (i in seq_len(length(ff))) {
            f <- as.formula(paste0("Y | D ", paste(as.character(ff[[i]]), collapse=" ")))
            mod <- try(cmulti(f, X, type=type))
            if (!inherits(mod, "try-error")) {
                rval <- mod[c("coefficients","vcov","nobs","loglik")]
                rval$p <- length(coef(mod))
                rval$names <- NAMES[[i]]
            } else {
                rval <- mod
            }
            res[[names(NAMES)[i]]] <- rval
        }
        ## number of observations
        #res$n <- n
    } else {
        res <- list(Y=Y, D=D, n=n)
    }
    res
}

SPP <- names(xtDur)
Pcodes <- levels(pkDur$PCODE)

resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(Pcodes)) {
        cat("Singing rate (rem) estimation for", SPP[i], "in", Pcodes[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem",
            pcode=Pcodes[j], excl=TRUE))
        resDurPcode1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem",
            pcode=Pcodes[j], excl=FALSE))
    }
}
str(resDurBAMless1)
str(resDurPcode1)
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem.Rdata")

resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(Pcodes)) {
        cat("Singing rate (mix) estimation for", SPP[i], "in", Pcodes[j], "\n")
        flush.console()
        resDurBAMless1_mix[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="mix",
            pcode=Pcodes[j], excl=TRUE))
        resDurPcode1_mix[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="mix",
            pcode=Pcodes[j], excl=FALSE))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix.Rdata")

if (FALSE) {
xtfunproj <- function(spp) {
    rn <- intersect(rownames(pkDur), rownames(xtDur[[spp]]))
    Y <- rowSums(xtDur[[spp]][rn,])
    tb <- table(PCODE=pkDur[rn,"PCODE"], n=ifelse(Y > 0, 1, 0))
    cbind(ntot=rowSums(tb), n1=tb[,"1"])
}
nn <- list()
for (spp in SPP)
    nn[[spp]] <- try(xtfunproj(spp))
nn <- nn[!sapply(nn, inherits, "try-error")]
}

## bias/variance/MSE for m0 and mbs

ROOT2 <- "~/Dropbox/bam/duration_ms/revisionOct2015"

library(MASS)
library(detect)
load_BAM_QPAD(1)
.BAMCOEFS$version
load(file.path(ROOT2, "BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_QPAD_v3_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS
.BAMCOEFSmix$version

aic0 <- .BAMCOEFS$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))

cfall0 <- sapply(SPP, function(spp) exp(coefBAMspecies(spp, 0, 0)$sra))
dim(cfall0) <- c(length(SPP), 1)
dimnames(cfall0) <- list(SPP, "phi_0")
cfallb <- t(sapply(SPP, function(spp) {
    tmp <- unname(.BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients)
    c(phi_b=exp(tmp[1]), c=plogis(tmp[2]))
    }))

pkDur <- droplevels(pkDur[pkDur$DURMETH %in% c("G","R","Y","Ya","Z"),])
rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
pkDur <- droplevels(pkDur[rn,])

dput(colnames(xtDur[[1]]))

ff <- list(
    ~ 1,
    ~ JDAY,
    ~ TSSR,
    ~ JDAY + I(JDAY^2),
    ~ TSSR + I(TSSR^2),
    ~ JDAY + TSSR,
    ~ JDAY + I(JDAY^2) + TSSR,
    ~ JDAY + TSSR + I(TSSR^2),
    ~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2),
    ~ TSLS,
    ~ TSLS + I(TSLS^2),
    ~ TSLS + TSSR,
    ~ TSLS + I(TSLS^2) + TSSR,
    ~ TSLS + TSSR + I(TSSR^2),
    ~ TSLS + I(TSLS^2) + TSSR + I(TSSR^2))
names(ff) <- 0:14
NAMES <- list(
    "0"="INTERCEPT",
    "1"=c("INTERCEPT", "JDAY"),
    "2"=c("INTERCEPT", "TSSR"),
    "3"=c("INTERCEPT", "JDAY", "JDAY2"),
    "4"=c("INTERCEPT", "TSSR", "TSSR2"),
    "5"=c("INTERCEPT", "JDAY", "TSSR"),
    "6"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR"),
    "7"=c("INTERCEPT", "JDAY", "TSSR", "TSSR2"),
    "8"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR", "TSSR2"),
    "9"=c("INTERCEPT", "TSLS"),
    "10"=c("INTERCEPT", "TSLS", "TSLS2"),
    "11"=c("INTERCEPT", "TSLS", "TSSR"),
    "12"=c("INTERCEPT", "TSLS", "TSLS2", "TSSR"),
    "13"=c("INTERCEPT", "TSLS", "TSSR", "TSSR2"),
    "14"=c("INTERCEPT", "TSLS", "TSLS2", "TSSR", "TSSR2"))

ymin <- 1
B <- 1000

spp <- "OVEN"


xtfun <- function(spp, ymin=1) {
    Y <- groupSums(xtDur[[spp]][rn,], 2,
        c("0-3", "xxx", "0-3", "0-3", "xxx", "xxx", "0-3", "0-3",
        "xxx", "3-5", "3-5", "xxx", "xxx", "3-5", "5-10", "5-10",
        "5-10", "5-10", "xxx", "5-10", "5-10", "5-10", "5-10", "5-10"))
    Y <- as.matrix(Y)[,c("0-3","3-5","5-10")]
    YY <- cbind("0-3"=Y[,1], "0-5"=Y[,1]+Y[,2], "0-10"=rowSums(Y))
    YY <- YY[YY[,3] >= ymin,,drop=FALSE]
    nrow(YY)
}

xvfun <- function(spp, B=1000, ymin=1) {
    Y <- groupSums(xtDur[[spp]][rn,], 2,
        c("0-3", "xxx", "0-3", "0-3", "xxx", "xxx", "0-3", "0-3",
        "xxx", "3-5", "3-5", "xxx", "xxx", "3-5", "5-10", "5-10",
        "5-10", "5-10", "xxx", "5-10", "5-10", "5-10", "5-10", "5-10"))
    Y <- as.matrix(Y)[,c("0-3","3-5","5-10")]
    YY <- cbind("0-3"=Y[,1], "0-5"=Y[,1]+Y[,2], "0-10"=rowSums(Y))
    YY <- YY[YY[,3] >= ymin,,drop=FALSE]
    pkDurS <- pkDur[rownames(YY),]


    cfi00 <- coefBAMspecies(spp, "0", "0")$sra
    vci00 <- drop(vcovBAMspecies(spp, "0", "0")$sra)

    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov

    best0 <- colnames(aic0)[which.min(aic0[spp,])]
    bestb <- colnames(aicb)[which.min(aicb[spp,])]
    X0 <- model.matrix(ff[[best0]], pkDurS)
    Xb <- model.matrix(ff[[bestb]], pkDurS)

    cf0 <- .BAMCOEFS$sra_estimates[[spp]][[best0]]$coefficients
    vc0 <- .BAMCOEFS$sra_estimates[[spp]][[best0]]$vcov
    cfb <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$coefficients
    vcb <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$vcov

    b_out <- list()
    for (j in 1:B) {
        phi00 <- exp(rnorm(1, cfi00, sqrt(vci00)))
        pcf0b <- mvrnorm(1, cfi0b, vci0b)

        pp_0 <- 1-exp(-c(3,5,10)*phi00)
        pp_b <- 1-plogis(pcf0b[2])*exp(-c(3,5,10)*exp(pcf0b[1]))

        YC3_0 <- (YY[,1] * pp_0[3]) / (pp_0[1] * YY[,3])
        YC5_0 <- (YY[,2] * pp_0[3]) / (pp_0[2] * YY[,3])
        YC3_b <- (YY[,1] * pp_b[3]) / (pp_b[1] * YY[,3])
        YC5_b <- (YY[,2] * pp_b[3]) / (pp_b[2] * YY[,3])

        ## time varying models

        cft0 <- mvrnorm(1, cf0, vc0)
        cftb <- mvrnorm(1, cfb, vcb)
        phi0i <- exp(X0 %*% cft0)
        phib <- exp(cftb[1])
        cbi <- plogis(Xb %*% cftb[-1])
        #phi0i <- exp(X0 %*% cf0)
        #phib <- exp(cfb[1])
        #scbi <- plogis(Xb %*% cfb[-1])

        YC3i_0 <- (YY[,1] * (1-exp(-10*phi0i))) / ((1-exp(-3*phi0i)) * YY[,3])
        YC5i_0 <- (YY[,2] * (1-exp(-10*phi0i))) / ((1-exp(-5*phi0i)) * YY[,3])
        YC3i_b <- (YY[,1] * (1-cbi*exp(-10*phib))) / ((1-cbi*exp(-3*phib)) * YY[,3])
        YC5i_b <- (YY[,2] * (1-cbi*exp(-10*phib))) / ((1-cbi*exp(-5*phib)) * YY[,3])

        out <- list()
        out$Y10 <- mean(YY[,3])
        out$m0 <- c(YC3=mean(YC3_0), YC5=mean(YC5_0))
        out$mb <- c(YC3=mean(YC3_b), YC5=mean(YC5_b))
        out$mt <- c(YC3=mean(YC3i_0), YC5=mean(YC5i_0))
        out$mbt <- c(YC3=mean(YC3i_b), YC5=mean(YC5i_b))

        b_out[[j]] <- unlist(out)
    }

    b_out <- do.call(rbind, b_out)

    theta <- b_out[,1]
    theta <- 1
    theta_hat <- b_out[,-1]

    MSE <- colSums((theta_hat - theta)^2) / B
    Var <- colSums(t(t(theta_hat) - colMeans(theta_hat))^2) / B
    Bias <- sqrt(MSE - Var)
    data.frame(MSE=MSE, Var=Var, Bias=Bias)
}

res <- list()
for (spp in SPP) {
cat(spp, "\n");flush.console()
res[[spp]] <- try(xvfun(spp))
}
nob <- sapply(SPP, xtfun, ymin=1)

save(res, nob, file="~/Dropbox/bam/duration_ms/revisionOct2015/var-bias-res.Rdata")


res2 <- res[!sapply(res, inherits, "try-error")]

aaa <- data.frame(Var=as.numeric(sapply(res2, "[[", "Var")),
    MSE=as.numeric(sapply(res2, "[[", "MSE")),
    Bias=as.numeric(sapply(res2, "[[", "Bias")),
    Model=rep(c("0","b","t","bt"), each=2),
    Duration=c("3","5"),
    Species=rep(names(res2), each=8))
aaa$n <- nob[match(aaa$Species, names(nob))]

m1 <- lm(Var ~ Model + Duration + n, aaa[aaa$n >= 1,])
m2 <- lm(Bias ~ Model + Duration + n, aaa[aaa$n >= 1,])
a1 <- anova(update(m1, . ~ . + Species))
a1$Perc <- 100 * a1[,"Sum Sq"] / sum(a1[,"Sum Sq"])
a2 <- anova(update(m2, . ~ . + Species))
a2$Perc <- 100 * a2[,"Sum Sq"] / sum(a2[,"Sum Sq"])
summary(m1)
a1
summary(m2)
a2

ng <- 1:200
sf <- function(x) quantile(x, 0.9, na.rm=TRUE)
maxVar <- cbind(max3_0 = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "t" & aaa$n >= z])),
    max5_t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxBias <- cbind(max3_0 = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "t" & aaa$n >= z])),
    max5_t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))

par(mfrow=c(4,2))
for (i in c(1,3,5,7)+1) {
plot(ng, maxVar[,i], main=paste("Var", colnames(maxVar)[i]),
    type="l", lwd=2, col=2, ylim=c(0, max(maxVar)))
lines(ng, maxVar[,i-1], lty=2, col=2, lwd=2)
plot(ng, maxBias[,i], main=paste("Bias", colnames(maxVar)[i]),
    type="l", lwd=2, col=2, ylim=c(min(maxBias), max(maxBias)))
lines(ng, maxBias[,i-1], lty=2, col=2, lwd=2)
}



## sra 0 vs b models

aic0 <- .BAMCOEFS$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
colnames(aic0) <- paste0("m0_", colnames(aic0))
colnames(aicb) <- paste0("mb_", colnames(aicb))
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))

aic0 <- aic0[SPP,]
aicb <- aicb[SPP,]
aic <- cbind(aic0[SPP,], aicb[SPP,])

np <- sapply(SPP, function(z) selectmodelBAMspecies(z)$sra$nobs[1])

waic0 <- t(apply(aic0, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}))
waicb <- t(apply(aicb, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}))
waic <- t(apply(aic, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}))

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
best <- colnames(aic)[apply(aic, 1, which.min)]

par(mfrow=c(1,3))
plot(np, waic0[,1], log="x", ylim=c(0,1), pch=ifelse(best0=="0", "o", "+"))
plot(np, waicb[,1], log="x", ylim=c(0,1), pch=ifelse(bestb=="0", "o", "+"))
plot(np, waic[,"m0_0"] + waic[,"mb_0"], log="x", ylim=c(0,1),
    pch=ifelse(best %in% c("m0_0", "mb_0"), "o", "+"))

best0
 0  1 10 11 12 13 14  2  3  4  5  6  7  8  9 
23 17 16  6 19 24 22 14 15 23 10 12 33  4 15 
> table(bestb)
bestb
 0  1 10 11 12 13 14  2  3  4  5  7  8  9 
50 18 12 14 11 34 20 15 12 20 10 23  4 10 
> table(best)
best
 m0_0  m0_1 m0_10 m0_11 m0_12 m0_14  m0_2  m0_3  m0_4  m0_5  m0_6  m0_7  m0_8  m0_9 
   11     8     7     3     2     2     8     8    10     4     5     8     1     8 
 mb_0  mb_1 mb_10 mb_11 mb_12 mb_13 mb_14  mb_2  mb_3  mb_4  mb_5  mb_7  mb_8  mb_9 
    9     8    10    12    10    32    19     8     7    17     7    18     4     7 
