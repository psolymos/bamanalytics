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
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionOct2015"

B <- 1000

## Load required packages
library(mefa4)
library(pbapply)
library(detect)
library(MASS)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2015-10-08.Rdata"))

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","PCODE")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
sort(table(pkDur$PCODE))
## only 1 obs in project (all others are >20)
pkDur <- droplevels(pkDur[pkDur$PCODE != "GLDRPCCL01",])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])


## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL

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

aic0 <- aic0[SPP,]
aicb <- aicb[SPP,]
aic <- cbind(aic0[SPP,], aicb[SPP,])

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
best <- colnames(aic)[apply(aic, 1, which.min)]
names(best) <- names(best0) <- names(bestb) <- SPP

rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
pkDur <- droplevels(pkDur[rn,])

## random PCODE
set.seed(1000)
pkDur$PCODErnd <- pkDur$PCODE[sample.int(nrow(pkDur))]

D <- ltdur$end[match(pkDur$DURMETH, rownames(ltdur$end)),]
## exclude 0 sum and <1 interval rows
nOK <- table(PC=pkDur$PCODE, n_int=rowSums(!is.na(D)))

nPC <- rowSums(nOK) - nOK[,"1"]
nPC <- sort(nPC[nPC > 0])
PC <- names(nPC)

## define pcode, and excl (excl=T to exclude pcode, excl=F to use only pcode)
## ff_id is the subset of ff list to be used
fitDurFun3 <- function(spp, fit=TRUE, type=c("rem","mix"), 
pcode=NULL, excl=TRUE, rnd=FALSE) 
{
    ff <- list(
        "0"=~ 1,
        "1"=~ JDAY,
        "2"=~ TSSR,
        "3"=~ JDAY + I(JDAY^2),
        "4"=~ TSSR + I(TSSR^2),
        "5"=~ JDAY + TSSR,
        "6"=~ JDAY + I(JDAY^2) + TSSR,
        "7"=~ JDAY + TSSR + I(TSSR^2),
        "8"=~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2),
        "9"=~ TSLS,
        "10"=~ TSLS + I(TSLS^2),
        "11"=~ TSLS + TSSR,
        "12"=~ TSLS + I(TSLS^2) + TSSR,
        "13"=~ TSLS + TSSR + I(TSSR^2),
        "14"=~ TSLS + I(TSLS^2) + TSSR + I(TSSR^2))
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
    if (excl) {
        if (type == "rem")
            ff_id <- unique(c("0", best0[spp]))
        if (type == "mix")
            ff_id <- unique(c("0", bestb[spp]))
        ff <- ff[ff_id]
        NAMES <- NAMES[ff_id]
    }

    rn <- intersect(rownames(pkDur), rownames(xtDur[[spp]]))
    X0 <- pkDur[rn,]
    if (!is.null(pcode)) {
        xPCODE <- if (rnd)
            X0$PCODErnd else X0$PCODE
        keep <- if (excl)
            xPCODE != pcode else xPCODE == pcode
        X0 <- X0[keep,,drop=FALSE]

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
            #mod <- try(cmulti.fit(Y, D, model.matrix(ff[[i]], X), type=type))
            if (!inherits(mod, "try-error")) {
                rval <- mod[c("coefficients","vcov","nobs","loglik")]
                rval$p <- length(coef(mod))
                rval$names <- NAMES[[i]]
            } else {
                rval <- mod
                attributes(rval) <- NULL
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

RND <- TRUE

SPP0 <- SPP
SPP1 <- SPP[1:60]
SPP2 <- SPP[61:120]
SPP3 <- SPP[121:180]
SPP4 <- SPP[181:length(SPP)]

## ------------------- rem

SPP <- SPP1
resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (rem) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-1rnd.Rdata")

SPP <- SPP2
resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (rem) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-2rnd.Rdata")

SPP <- SPP3
resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (rem) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-3rnd.Rdata")

SPP <- SPP4
resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (rem) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-4rnd.Rdata")

## ------------------- mix

SPP <- SPP1
resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (mix) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-1rnd.Rdata")

SPP <- SPP2
resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (mix) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-2rnd.Rdata")

SPP <- SPP3
resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (mix) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-3rnd.Rdata")

SPP <- SPP4
resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(PC)) {
        cat("Singing rate (mix) estimation for", SPP[i], "excluding", PC[j], "\n")
        flush.console()
        resDurBAMless1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=TRUE, rnd=RND))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE, rnd=RND))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-4rnd.Rdata")


## read all the suff, read in lines 1-86 (RND defs there!)

isRND <- TRUE

## rem
e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-1.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-1rnd.Rdata"), envir=e)
resDurBAMless1 <- e$resDurBAMless1
resDurPcode1 <- e$resDurPcode1

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-2.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-2rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-3.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-3rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-4.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-4rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

## mix
e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-1.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-1rnd.Rdata"), envir=e)
resDurBAMless1_mix <- e$resDurBAMless1_mix
resDurPcode1_mix <- e$resDurPcode1_mix

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-2.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-2rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-3.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-3rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-4.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-4rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

ff <- list(
        "0"=~ 1,
        "1"=~ JDAY,
        "2"=~ TSSR,
        "3"=~ JDAY + I(JDAY^2),
        "4"=~ TSSR + I(TSSR^2),
        "5"=~ JDAY + TSSR,
        "6"=~ JDAY + I(JDAY^2) + TSSR,
        "7"=~ JDAY + TSSR + I(TSSR^2),
        "8"=~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2),
        "9"=~ TSLS,
        "10"=~ TSLS + I(TSLS^2),
        "11"=~ TSLS + TSSR,
        "12"=~ TSLS + I(TSLS^2) + TSSR,
        "13"=~ TSLS + TSSR + I(TSSR^2),
        "14"=~ TSLS + I(TSLS^2) + TSSR + I(TSSR^2))
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

aic_fun <- function(x) {
#    if (inherits(x, "try-error"))
    if (is.character(x))
        Inf else -2*x$loglik + 2*x$p
}
coef_fun <- function(x, id) {
#    if (inherits(x, "try-error"))
    if (is.character(x[[id]]))
        NA else x[[id]]$coefficients
}
vcov_fun <- function(x, id) {
#    if (inherits(x, "try-error"))
    if (is.character(x[[id]]))
        NA else x[[id]]$vcov
}
pred_fun <- function(X, c, v, type, ptonly=FALSE) {
    tmp <- data.frame(q3.mu=NA,q3.var=NA,
        q5.mu=NA,q5.var=NA,
        q10.mu=NA,q10.var=NA)

    if (any(is.na(c)))
        return(tmp)
    if (any(is.na(v)))
        return(tmp)

    if (ptonly) {
        if (type=="rem") {
            p3 <- qlogis(1-exp(-3*exp(X %*% c)))
            p5 <- qlogis(1-exp(-5*exp(X %*% c)))
            p10 <- qlogis(1-exp(-10*exp(X %*% c)))
        } else {
            p3 <- qlogis(1-plogis(X %*% c[-1])*exp(-3*exp(c[1])))
            p5 <- qlogis(1-plogis(X %*% c[-1])*exp(-5*exp(c[1])))
            p10 <- qlogis(1-plogis(X %*% c[-1])*exp(-10*exp(c[1])))
        }
        out <- data.frame(q3.mu=p3, q3.var=NA, 
            q5.mu=p5, q5.var=NA, q10.mu=p10, q10.var=NA)
    } else {
        v <- try(as.matrix(Matrix::nearPD(v)$mat))
        if (inherits(v, "try-error"))
            return(tmp)
        cfs <- try(mvrnorm(B, c, v))
        if (inherits(cfs, "try-error"))
            return(tmp)
        if (type=="rem") {
            p3 <- apply(cfs, 1, function(z) qlogis(1-exp(-3*exp(X %*% z))))
            p5 <- apply(cfs, 1, function(z) qlogis(1-exp(-5*exp(X %*% z))))
            p10 <- apply(cfs, 1, function(z) qlogis(1-exp(-10*exp(X %*% z))))
        } else {
            p3 <- apply(cfs, 1, function(z) qlogis(1-plogis(X %*% z[-1])*exp(-3*exp(z[1]))))
            p5 <- apply(cfs, 1, function(z) qlogis(1-plogis(X %*% z[-1])*exp(-5*exp(z[1]))))
            p10 <- apply(cfs, 1, function(z) qlogis(1-plogis(X %*% z[-1])*exp(-10*exp(z[1]))))
        }
        q3 <- t(apply(p3, 1, function(z) c(mean(z), sd(z)^2)))
        q5 <- t(apply(p5, 1, function(z) c(mean(z), sd(z)^2)))
        q10 <- t(apply(p10, 1, function(z) c(mean(z), sd(z)^2)))
        colnames(q3) <- colnames(q5) <- colnames(q10) <- c("mu","var")
        out <- data.frame(q3=q3, q5=q5, q10=q10)
    }
    out
}

logphi_fun <- function(X, c, v, level=0.95) {

    tmp <- matrix(NA, nrow(X), 3)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pt <- drop(X %*% c)
    tmp[,1] <- pt
    v <- try(as.matrix(Matrix::nearPD(v)$mat))
    if (inherits(v, "try-error"))
        return(tmp)
    cfs <- try(mvrnorm(B, c, v))
    if (inherits(cfs, "try-error"))
        return(tmp)
    pr <- apply(cfs, 1, function(z) drop(X %*% z))
    ci <- t(apply(pr, 1, quantile, a))
    tmp[,2:3] <- ci
    tmp
}
pred_fun_old <- function(X, c, v, type, ptonly=FALSE) {
    tmp <- data.frame(q3.med=NA,q3.lo=NA,q3.up=NA,
        q5.med=NA,q5.lo=NA,q5.up=NA, q10.med=NA,q10.lo=NA,q10.up=NA)
    if (ptonly)
        tmp <- data.frame(q3=NA,q5=NA,q10=NA)
    if (any(is.na(c)))
        return(tmp)
    if (any(is.na(v)))
        return(tmp)
    v <- try(as.matrix(Matrix::nearPD(v)$mat))
    if (inherits(v, "try-error"))
        return(tmp)
    cfs <- try(mvrnorm(B, c, v))
    if (inherits(cfs, "try-error"))
        return(tmp)
    if (ptonly) {
        if (type=="rem") {
            p3 <- 1-exp(-3*exp(X %*% c))
            p5 <- 1-exp(-5*exp(X %*% c))
            p10 <- 1-exp(-10*exp(X %*% c))
        } else {
            p3 <- 1-plogis(X %*% c[-1])*exp(-3*exp(c[1]))
            p5 <- 1-plogis(X %*% c[-1])*exp(-5*exp(c[1]))
            p10 <- 1-plogis(X %*% c[-1])*exp(-10*exp(c[1]))
        }
    } else {
        if (type=="rem") {
            p3 <- apply(cfs, 1, function(z) 1-exp(-3*exp(X %*% z)))
            p5 <- apply(cfs, 1, function(z) 1-exp(-5*exp(X %*% z)))
            p10 <- apply(cfs, 1, function(z) 1-exp(-10*exp(X %*% z)))
        } else {
            p3 <- apply(cfs, 1, function(z) 1-plogis(X %*% z[-1])*exp(-3*exp(z[1])))
            p5 <- apply(cfs, 1, function(z) 1-plogis(X %*% z[-1])*exp(-5*exp(z[1])))
            p10 <- apply(cfs, 1, function(z) 1-plogis(X %*% z[-1])*exp(-10*exp(z[1])))
        }
        q3 <- t(apply(p3, 1, quantile, c(0.5, 0.05, 0.95)))
        q5 <- t(apply(p5, 1, quantile, c(0.5, 0.05, 0.95)))
        q10 <- t(apply(p10, 1, quantile, c(0.5, 0.05, 0.95)))
        colnames(q3) <- colnames(q5) <- colnames(q10) <- c("med","lo","up")
    }
    out <- if (ptonly)
        data.frame(q3=p3, q5=p5, q10=p10) else data.frame(q3=q3, q5=q5, q10=q10)
    out
}
compare_ci <- function(lo1, up1, lo2, up2) {
    if (any(is.na(c(lo1, up1, lo2, up2))))
        return(NA)
    x1lessthan2lo <- lo1 < lo2 & up1 < lo2
    x1greaterthan2up <- lo1 > up2 & up1 > up2
    !(x1lessthan2lo | x1greaterthan2up)
}
waic_fun <- function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}

compare_distr <- function(mu1, var1, mu2, var2, n) {
    tval <- (mu1-mu2) / sqrt(var1/n + var2/n)
    df <- ((var1/n) + (var2/n))^2 / (((var1/n)^2 / (n-1)) + ((var2/n)^2 / (n-1)))
    pval <- 2*pt(-abs(tval), df)
    c(t=tval, df=df, p=pval)
}
if (FALSE) { 
#https://gist.github.com/psolymos/2890e9ab87946163242f#file-welch-t-test-r
## Welch t-test

### Simple Normal-Normal example
mu1 <- 1
mu2 <- 2
var1 <- 1
var2 <- 0.5
n <- 10^3
y1 <- rnorm(n, mu1, sqrt(var1))
y2 <- rnorm(n, mu2, sqrt(var2))

## base implementation
(tt <- t.test(y1, y2, var.equal=FALSE))

## by hand
mu1 <- mean(y1)
mu2 <- mean(y2)
var1 <- sd(y1)^2
var2 <- sd(y2)^2
(tval <- (mu1-mu2) / sqrt(var1/n + var2/n))
(df <- ((var1/n) + (var2/n))^2 / (((var1/n)^2 / (n-1)) + ((var2/n)^2 / (n-1))))
(pval <- 2*pt(-abs(tval), df))

## a function
compare_distr <- function(mu1, var1, mu2, var2) {
    tval <- (mu1-mu2) / sqrt(var1/n + var2/n)
    df <- ((var1/n) + (var2/n))^2 / (((var1/n)^2 / (n-1)) + ((var2/n)^2 / (n-1)))
    pval <- 2*pt(-abs(tval), df)
    c(t=tval, df=df, p=pval)
}
compare_distr(mu1, var1, mu2, var2)
}


ptonly <- FALSE # point pred (TRUE) or CI as well (FALSE)
external <- TRUE # TRUE: (project)(the rest of BAM), FALSE: ((proj) all BAM)

SPP0 <- SPP
SPP1 <- SPP[1:60]
SPP2 <- SPP[61:120]
SPP3 <- SPP[121:180]
SPP4 <- SPP[181:length(SPP)]

#spp <- "OVEN"
#pc <- "LMWELL"
xvres <- list()

for (spp in SPP) {
#for (spp in SPP1) {
#for (spp in SPP2) {
#for (spp in SPP3) {
#for (spp in SPP4) {

tab <- list()
for (pc in PC) {

gc()
OK <- c(!inherits(resDurPcode1[[spp]][[pc]], "try-error"),
    !inherits(resDurBAMless1[[spp]][[pc]], "try-error"),
    !inherits(resDurPcode1_mix[[spp]][[pc]], "try-error"),
    !inherits(resDurBAMless1_mix[[spp]][[pc]], "try-error"))

if (all(OK)) {

    cat(spp, pc, "\n");flush.console()
    ## project pc exluded
    best_99_m0 <- unname(best0[spp])
    best_99_mb <- unname(bestb[spp])
    type_99_best <- if (which.min(aic[spp,]) > 15) "mix" else "rem"
    best_99 <- unname(best[spp])

    if (external) {
        c_99_m0 <- coef_fun(resDurBAMless1[[spp]][[pc]], "0")
        v_99_m0 <- vcov_fun(resDurBAMless1[[spp]][[pc]], "0")
        c_99_mb <- coef_fun(resDurBAMless1_mix[[spp]][[pc]], "0")
        v_99_mb <- vcov_fun(resDurBAMless1_mix[[spp]][[pc]], "0")

        c_99_best0 <- coef_fun(resDurBAMless1[[spp]][[pc]], best_99_m0)
        v_99_best0 <- vcov_fun(resDurBAMless1[[spp]][[pc]], best_99_m0)
        c_99_bestb <- coef_fun(resDurBAMless1_mix[[spp]][[pc]], best_99_mb)
        v_99_bestb <- vcov_fun(resDurBAMless1_mix[[spp]][[pc]], best_99_mb)
    } else {
        c_99_m0 <- coef_fun(.BAMCOEFS$sra_estimates[[spp]], "0")
        v_99_m0 <- vcov_fun(.BAMCOEFS$sra_estimates[[spp]], "0")
        c_99_mb <- coef_fun(.BAMCOEFSmix$sra_estimates[[spp]], "0")
        v_99_mb <- vcov_fun(.BAMCOEFSmix$sra_estimates[[spp]], "0")

        c_99_best0 <- coef_fun(.BAMCOEFS$sra_estimates[[spp]], best_99_m0)
        v_99_best0 <- vcov_fun(.BAMCOEFS$sra_estimates[[spp]], best_99_m0)
        c_99_bestb <- coef_fun(.BAMCOEFSmix$sra_estimates[[spp]], best_99_mb)
        v_99_bestb <- vcov_fun(.BAMCOEFSmix$sra_estimates[[spp]], best_99_mb)
    }
    if (type_99_best == "rem") {
        c_99_best <- c_99_best0
        v_99_best <- v_99_best0
    } else {
        c_99_best <- c_99_bestb
        v_99_best <- v_99_bestb
    }

    ## only project pc
    aic_1_m0 <- sapply(resDurPcode1[[spp]][[pc]], aic_fun)
    aic_1_mb <- sapply(resDurPcode1_mix[[spp]][[pc]], aic_fun)

    waic_1_m0 <- waic_fun(aic_1_m0)
    waic_1_mb <- waic_fun(aic_1_mb)
    waic_1    <- waic_fun(c(aic_1_m0, aic_1_mb))
    w0_1_m0 <- waic_1_m0["0"]
    w0_1_mb <- waic_1_mb["0"]
    w0_1 <- sum(waic_1[which(names(waic_1)=="0")])

    waic_99_m0 <- waic_fun(aic0[spp,])
    waic_99_mb <- waic_fun(aicb[spp,])
    waic_99    <- waic_fun(aic[spp,])
    w0_99_m0 <- waic_99_m0["0"]
    w0_99_mb <- waic_99_mb["0"]
    w0_99 <- sum(waic_99[which(names(waic_99)=="0")])

    best_1_m0 <- names(aic_1_m0)[which.min(aic_1_m0)]
    best_1_mb <- names(aic_1_mb)[which.min(aic_1_mb)]
    type_1_best <- if (which.min(c(aic_1_m0, aic_1_mb)) > 15) "mix" else "rem"
    best_1 <- names(c(aic_1_m0, aic_1_mb))[which.min(c(aic_1_m0, aic_1_mb))]

    c_1_m0 <- coef_fun(resDurPcode1[[spp]][[pc]], "0")
    v_1_m0 <- vcov_fun(resDurPcode1[[spp]][[pc]], "0")
    c_1_mb <- coef_fun(resDurPcode1_mix[[spp]][[pc]], "0")
    v_1_mb <- vcov_fun(resDurPcode1_mix[[spp]][[pc]], "0")

    c_1_best0 <- coef_fun(resDurPcode1[[spp]][[pc]], best_1_m0)
    v_1_best0 <- vcov_fun(resDurPcode1[[spp]][[pc]], best_1_m0)
    c_1_bestb <- coef_fun(resDurPcode1_mix[[spp]][[pc]], best_1_mb)
    v_1_bestb <- vcov_fun(resDurPcode1_mix[[spp]][[pc]], best_1_mb)
    if (type_1_best == "rem") {
        c_1_best <- c_1_best0
        v_1_best <- v_1_best0
    } else {
        c_1_best <- c_1_bestb
        v_1_best <- v_1_bestb
    }

    df <- pkDur
    ## random or not
    xPCODE <- if (isRND)
        df$PCODErnd else df$PCODE
    df$y <- rowSums(xtDur[[spp]][rownames(df),])
    df$pc <- ifelse(xPCODE == pc, 1L, 0L)
    df99 <- df[xPCODE != pc,,drop=FALSE]
    df1 <- df[xPCODE == pc,,drop=FALSE]

    n99 <- nrow(df99)
    n1 <- nrow(df1)
    det99 <- sum(df99$y > 0)
    det1 <- sum(df1$y > 0)

#    mod0 <- glm(pc ~ 1, df, family=binomial("logit"))
#    mod1 <- glm(pc ~ JDAY + TSSR + TSLS, df, family=binomial("logit"))
#    logLR <- as.numeric(logLik(mod1) - logLik(mod0))


    ## model matrix for project based prediction
#    X1_1 <- if (ptonly && best_1 == "0")
#        matrix(1, 1, 1) else model.matrix(ff[[best_1]], df1)
#    X1_99 <- if (ptonly && best_99 == "0")
#        matrix(1, 1, 1) else model.matrix(ff[[best_99]], df1)
#    X1 <- if (ptonly)
#        matrix(1, 1, 1) else matrix(1, n1, 1)

    X1_1 <- if (ptonly && best_1_m0 == "0")
        matrix(1, 1, 1) else model.matrix(ff[[best_1_m0]], df1)
    X1_99 <- if (ptonly && best_99_m0 == "0")
        matrix(1, 1, 1) else model.matrix(ff[[best_99_m0]], df1)
    X1 <- if (ptonly)
        matrix(1, 1, 1) else matrix(1, n1, 1)

    LEVEL <- 0.9

    pr_1_m0 <- logphi_fun(X1, c_1_m0, v_1_m0, level=LEVEL)
    pr_99_m0 <- logphi_fun(X1, c_99_m0, v_99_m0, level=LEVEL)
    pr_1_best0 <- logphi_fun(X1_1, c_1_best0, v_1_best0, level=LEVEL)
    pr_99_best0 <- logphi_fun(X1_99, c_99_best0, v_99_best0, level=LEVEL)
    OK0 <- as.integer(pr_1_m0[,2] <= pr_99_m0[,1] & pr_1_m0[,3] >= pr_99_m0[,1])
    OKbest0 <- as.integer(pr_1_best0[,2] <= pr_99_best0[,1] & pr_1_best0[,3] >= pr_99_best0[,1])

    tab[[pc]] <- data.frame(spp=spp,
        pc=pc,
        n99=n99, n1=n1, det99=det99, det1=det1, 
        best99=best_99, best1=best_1, type99=type_99_best, type1=type_1_best,
        w1_m0=w0_1_m0, w1_mb=w0_1_mb, w1=w0_1, w99_m0=w0_99_m0, w99_mb=w0_99_mb, w99=w0_99,
        OK0=mean(OK0), OKbest0=mean(OKbest0))
    }
}
tab <- do.call(rbind, tab)

xvres[[spp]] <- tab
}
#save(xvres, file=file.path(ROOT2, "xval-summary.Rdata"))
#save(xvres, file=file.path(ROOT2, "xval-summary-rnd.Rdata"))

if (!isRND) {
    #save(xvres, file=file.path(ROOT2, "xval-summary-1.Rdata"))
    #save(xvres, file=file.path(ROOT2, "xval-summary-2.Rdata"))
    #save(xvres, file=file.path(ROOT2, "xval-summary-3.Rdata"))
    save(xvres, file=file.path(ROOT2, "xval-summary-4.Rdata"))
}
if (isRND) {
    #save(xvres, file=file.path(ROOT2, "xval-summary-1rnd.Rdata"))
    #save(xvres, file=file.path(ROOT2, "xval-summary-2rnd.Rdata"))
    #save(xvres, file=file.path(ROOT2, "xval-summary-3rnd.Rdata"))
    #save(xvres, file=file.path(ROOT2, "xval-summary-4rnd.Rdata"))
}








load(file.path(ROOT2, "xval-dis.Rdata"))
#load(file.path(ROOT2, "xval-summary-new.Rdata"))


## load stuff
ROOT <- "c:/bam/May2015"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionOct2015"
library(mefa4)

isRND <- TRUE

e <- new.env()
if (isRND)
    load(file.path(ROOT2, "xval-summary-1rnd.Rdata"), envir=e)
if (!isRND)
    load(file.path(ROOT2, "xval-summary-1.Rdata"), envir=e)
xvres <- e$xvres
e <- new.env()
if (isRND)
    load(file.path(ROOT2, "xval-summary-2rnd.Rdata"), envir=e)
if (!isRND)
    load(file.path(ROOT2, "xval-summary-2.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)
e <- new.env()
if (isRND)
    load(file.path(ROOT2, "xval-summary-3rnd.Rdata"), envir=e)
if (!isRND)
    load(file.path(ROOT2, "xval-summary-3.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)
e <- new.env()
if (isRND)
    load(file.path(ROOT2, "xval-summary-4rnd.Rdata"), envir=e)
if (!isRND)
    load(file.path(ROOT2, "xval-summary-4.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)

c1 <- c("pr_1_m0.q3.mu", "pr_1_m0.q3.var", 
    "pr_1_m0.q5.mu", "pr_1_m0.q5.var", 
    "pr_1_m0.q10.mu", "pr_1_m0.q10.var", 
    "pr_1_mb.q3.mu", "pr_1_mb.q3.var", 
    "pr_1_mb.q5.mu", "pr_1_mb.q5.var", 
    "pr_1_mb.q10.mu", "pr_1_mb.q10.var", 
    "pr_1_best.q3.mu", "pr_1_best.q3.var", 
    "pr_1_best.q5.mu", "pr_1_best.q5.var", 
    "pr_1_best.q10.mu", "pr_1_best.q10.var") 
c99 <- c("pr_99_m0.q3.mu", "pr_99_m0.q3.var", 
    "pr_99_m0.q5.mu", "pr_99_m0.q5.var", 
    "pr_99_m0.q10.mu", "pr_99_m0.q10.var", 
    "pr_99_mb.q3.mu", "pr_99_mb.q3.var", 
    "pr_99_mb.q5.mu", "pr_99_mb.q5.var", 
    "pr_99_mb.q10.mu", "pr_99_mb.q10.var", 
    "pr_99_best.q3.mu", "pr_99_best.q3.var", 
    "pr_99_best.q5.mu", "pr_99_best.q5.var", 
    "pr_99_best.q10.mu", "pr_99_best.q10.var")
mc1 <- matrix(c1, ncol=2, byrow=TRUE)
mc99 <- matrix(c99, ncol=2, byrow=TRUE)
mc <- cbind(mc1, mc99)

compare_distr2 <- function(x, j, level=0.95) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    Cols <- mc[j,]
    ci <- x[, Cols[1]] + x[, Cols[2]] %o% fac
    as.integer(ci[,1] <= x[, Cols[3]] & ci[,2] >= x[, Cols[3]])
}
compare_distr2(x, 1, 0.9)
aa <- t(sapply(1:nrow(x), compare_distr2, j=1, x=x))
j <- 8
xv <- list()
for (spp in names(xvres)) {
    
    compare_distr2(xvres[[spp]], 1, 0.9)
    
    
}

spp <- "OVEN"
x <- xvres[[spp]]
x <- x[rowSums(is.na(x)) == 0,]
#lm(pr_1_best.q5.mu ~ pr_99_best.q5.mu + det1, x)
d <- x[,"pr_1_best.q5.mu"] - x[,"pr_99_best.q5.mu"]
plot(jitter(x$det1), abs(d), main=spp, pch=21, log="x")

plot(x[,"pr_1_best.q5.mu"], x[,"pr_99_best.q5.mu"],
    ylim=c(-10,10), xlim=c(-10,10))
abline(0,1,col=2)



xv <- do.call(rbind, xvres)
xv$dis <- dis[match(xv$pc, rownames(dis)), "MEAN"]
xv$occ <- (xv$det99 + xv$det1) / (xv$n99 + xv$n1)
xv$occ1 <- xv$det1 / xv$n1
xv$best1is0 <- ifelse(xv$best1==0, 1, 0)
xv$best99is0 <- ifelse(xv$best99==0, 1, 0)

keep <- rep(TRUE, nrow(xv))
keep[xv$pr_1_m0.q10 < 0.05 | xv$pr_1_mb.q10 < 0.05 | xv$pr_1_best.q10 < 0.05] <- FALSE
keep[xv$pr_1_m0.q3 > 0.95 | xv$pr_1_mb.q3 > 0.95 | xv$pr_1_best.q3 > 0.95] <- FALSE
keep[is.na(xv$pr_1_m0.q3) | is.na(xv$pr_1_mb.q3) | is.na(xv$pr_1_best.q3)] <- FALSE
xvx <- xv[keep,]

xvx$w1 <- xvx$best_10 / xvx$n1

plot(w1 ~ log(n1), xvx, pch=ifelse(xvx$best1is0>0, "x", "."))
abline(lm(w1 ~ log(n1), xvx), col=2)

plot(w1 ~ det1, xvx, pch=21, col=ifelse(xvx$best1is0>0, "blue", "grey"), xlim=c(0,1000))
plot(w1 ~ log(det1), xvx, pch=21, col=ifelse(xvx$best1is0>0, "blue", "grey"))
abline(lm(w1 ~ log(det1), xvx), col=2)
#boxplot(log(det1) ~ best1is0, xvx)

nn <- 0:5000#max(xvx$det1)
#xx <- sapply(nn, function(z) mean(xvx$best1is0[xvx$det1 >= z]))
xx <- sapply(nn, function(z) mean(c(xv$best1is0[xv$det1 >= z], xv$best99is0[xv$det99 >= z]) ))
plot(nn, xx, type="l", ylim=c(0, max(xx)), xlab="Number of detections",
    ylab="P(constant m0 or mb)")
rug(xvx$det1)

## -- old

plot(pr_1_best$q3.med, pr_99_best$q3.med, ylim=c(0,1), xlim=c(0,1))
abline(0,1,col=2)

plot_fun <- function() {
op <- par(mfrow=c(3,3))

plot(1:n1, pr_1_m0$q3.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_m0$q3.lo, y1=pr_1_m0$q3.up, col=ifelse(sdiff3_m0, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_m0$q3.lo, y1=pr_99_m0$q3.up, col=ifelse(sdiff3_m0, 3, 4), lwd=3)

plot(1:n1, pr_1_mb$q3.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_mb$q3.lo, y1=pr_1_mb$q3.up, col=ifelse(sdiff3_mb, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_mb$q3.lo, y1=pr_99_mb$q3.up, col=ifelse(sdiff3_mb, 3, 4), lwd=3)

plot(1:n1, pr_1_best$q3.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_best$q3.lo, y1=pr_1_best$q3.up, col=ifelse(sdiff3_best, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_best$q3.lo, y1=pr_99_best$q3.up, col=ifelse(sdiff3_best, 3, 4), lwd=3)

plot(1:n1, pr_1_m0$q5.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_m0$q5.lo, y1=pr_1_m0$q5.up, col=ifelse(sdiff5_m0, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_m0$q5.lo, y1=pr_99_m0$q5.up, col=ifelse(sdiff5_m0, 3, 4), lwd=3)

plot(1:n1, pr_1_mb$q5.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_mb$q5.lo, y1=pr_1_mb$q5.up, col=ifelse(sdiff5_mb, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_mb$q5.lo, y1=pr_99_mb$q5.up, col=ifelse(sdiff5_mb, 3, 4), lwd=3)

plot(1:n1, pr_1_best$q5.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_best$q5.lo, y1=pr_1_best$q5.up, col=ifelse(sdiff5_best, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_best$q5.lo, y1=pr_99_best$q5.up, col=ifelse(sdiff5_best, 3, 4), lwd=3)

plot(1:n1, pr_1_m0$q10.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_m0$q10.lo, y1=pr_1_m0$q10.up, col=ifelse(sdiff10_m0, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_m0$q10.lo, y1=pr_99_m0$q10.up, col=ifelse(sdiff10_m0, 3, 4), lwd=3)

plot(1:n1, pr_1_mb$q10.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_mb$q10.lo, y1=pr_1_mb$q10.up, col=ifelse(sdiff10_mb, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_mb$q10.lo, y1=pr_99_mb$q10.up, col=ifelse(sdiff10_mb, 3, 4), lwd=3)

plot(1:n1, pr_1_best$q10.med, ylim=c(0,1), type="n")
segments(x0=(1:n1)-0.2, y0=pr_1_best$q10.lo, y1=pr_1_best$q10.up, col=ifelse(sdiff10_best, 3, 2), lwd=3)
segments(x0=(1:n1)+0.2, y0=pr_99_best$q10.lo, y1=pr_99_best$q10.up, col=ifelse(sdiff10_best, 3, 4), lwd=3)
par(op)
}

## write a function that calculates point pred & CI depending on rem/mix
## calculates CIs and compare -- write compare fun.

#logLR <- list()
dis <- list()
for (pc in PC) {
    cat(pc, "\n");flush.console()

    df <- pkDur
    #df$y <- rowSums(xtDur[[spp]][rownames(df),])
    df$pc <- ifelse(df$PCODE == pc, 1L, 0L)
    df99 <- df[df$PCODE != pc,,drop=FALSE]
    df1 <- df[df$PCODE == pc,,drop=FALSE]

    jd0 <- range(df$JDAY)
    jd99 <- range(df99$JDAY)
    jd1 <- range(df1$JDAY)
    ts0 <- range(df$TSSR)
    ts99 <- range(df99$TSSR)
    ts1 <- range(df1$TSSR)
    ls0 <- range(df$TSLS)
    ls99 <- range(df99$TSLS)
    ls1 <- range(df1$TSLS)

    djd <- if (compare_ci(jd1[1], jd1[2], jd99[1], jd99[2]))
        (min(jd1[2], jd99[2]) - max(jd1[1], jd99[1])) / diff(jd0) else 0
    dts <- if (compare_ci(ts1[1], ts1[2], ts99[1], ts99[2]))
        (min(ts1[2], ts99[2]) - max(ts1[1], ts99[1])) / diff(ts0) else 0
    dls <- if (compare_ci(ls1[1], ls1[2], ls99[1], ls99[2]))
        (min(ls1[2], ls99[2]) - max(ls1[1], ls99[1])) / diff(ls0) else 0
    dis[[pc]] <- c(MEAN=mean(c(djd, dts, dls)), JDAY=djd,
        TSSR=dts, TSLS=dls)

    #mod0 <- glm(pc ~ 1, df, family=binomial("logit"))
    #mod1 <- glm(pc ~ JDAY + TSSR + TSLS, df, family=binomial("logit"))
    #logLR[[pc]] <- as.numeric(logLik(mod1) - logLik(mod0))
}
#save(logLR, file=file.path(ROOT2, "xval-logLR.Rdata"))
dis <- do.call(rbind, dis)
save(dis, file=file.path(ROOT2, "xval-dis.Rdata"))

e <- new.env()
load(file.path(ROOT2, "xval-summary-int-1.Rdata"), envir=e)
xvres <- e$xvres
e <- new.env()
load(file.path(ROOT2, "xval-summary-int-2.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)
e <- new.env()
load(file.path(ROOT2, "xval-summary-int-3.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)
e <- new.env()
load(file.path(ROOT2, "xval-summary-int-4.Rdata"), envir=e)
xvres <- c(xvres, e$xvres)

for (i in 1:length(xvres)) {
    xvres[[i]]$pc <- rownames(xvres[[i]])
}

xv <- do.call(rbind, xvres)
xv$dis <- dis[match(xv$pc, rownames(dis)), "MEAN"]
xv$occ <- (xv$det99 + xv$det1) / (xv$n99 + xv$n1)
xv$occ1 <- xv$det1 / xv$n1
xv$best1is0 <- ifelse(xv$best1==0, 1, 0)

#xv <- xv[xv$det1 >= 200,]

xv1 <- xv[,c("m0_3","m0_5","m0_10","mb_3","mb_5","mb_10","best_3","best_5","best_10")] / xv$n1
summary(xv1)
xv$ok <- xv1$best_3 >= 0.9

mod <- glm(cbind(m0_3, n1-m0_3) ~ occ + dis, xv, family=binomial)

summary(mm <- glm(best1is0 ~ n1, xv, family=binomial))

boxplot(n1 ~ ok, xv)

plot(xv$det1, xv1$best_3)

plot(xv1$best_3 ~ sqrt(xv$n1))
abline(lm(xv1$best_3 ~ sqrt(xv$n1)), col=2)
