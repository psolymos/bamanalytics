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

D <- ltdur$end[match(pkDur$DURMETH, rownames(ltdur$end)),]
## exclude 0 sum and <1 interval rows
nOK <- table(PC=pkDur$PCODE, n_int=rowSums(!is.na(D)))

nPC <- rowSums(nOK) - nOK[,"1"]
nPC <- sort(nPC[nPC > 0])
PC <- names(nPC)

## define pcode, and excl (excl=T to exclude pcode, excl=F to use only pcode)
## ff_id is the subset of ff list to be used
fitDurFun3 <- function(spp, fit=TRUE, type=c("rem","mix"), pcode=NULL, excl=TRUE) 
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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-1.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-2.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-3.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (rem) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="rem",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1, resDurPcode1,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-rem-4.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-1.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-2.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-3.Rdata")

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
            pcode=PC[j], excl=TRUE))
        cat("Singing rate (mix) estimation for", SPP[i], "in", PC[j], "\n")
        flush.console()
        resDurPcode1_mix[[SPP[i]]][[PC[j]]] <- try(fitDurFun3(SPP[i], TRUE, type="mix",
            pcode=PC[j], excl=FALSE))
    }
}
save(resDurBAMless1_mix, resDurPcode1_mix,
    file="~/Dropbox/bam/duration_ms/revisionOct2015/xval-mix-4.Rdata")


## read all the suff

## rem
e <- new.env()
load(file.path(ROOT2, "xval-rem-1.Rdata"), envir=e)
resDurBAMless1 <- e$resDurBAMless1
resDurPcode1 <- e$resDurPcode1

e <- new.env()
load(file.path(ROOT2, "xval-rem-2.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
load(file.path(ROOT2, "xval-rem-3.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
load(file.path(ROOT2, "xval-rem-4.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

## mix
e <- new.env()
load(file.path(ROOT2, "xval-mix-1.Rdata"), envir=e)
resDurBAMless1_mix <- e$resDurBAMless1_mix
resDurPcode1_mix <- e$resDurPcode1_mix

e <- new.env()
load(file.path(ROOT2, "xval-mix-2.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
load(file.path(ROOT2, "xval-mix-3.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
load(file.path(ROOT2, "xval-mix-4.Rdata"), envir=e)
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
    if (inherits(x, "try-error"))
        Inf else -2*x$loglik + 2*x$p
}
coef_fun <- function(x, id) {
    if (inherits(x, "try-error"))
        rep(NA, NAMES[id]) else x[[id]]$coefficients
}
vcov_fun <- function(x, id) {
    if (inherits(x, "try-error"))
        matrix(NA, NAMES[id], NAMES[id]) else x[[id]]$vcov
}

spp <- "OVEN"
pc <- "LMWELL"

## project pc exluded
best_99_m0 <- unname(best0[spp])
best_99_mb <- unname(bestb[spp])
type_99_best <- if (which.min(aic[spp,]) > 15) "mix" else "rem"
best_99 <- unname(best[spp])

c_99_m0 <- coef_fun(resDurBAMless1[[spp]][[pc]], "0")
v_99_m0 <- vcov_fun(resDurBAMless1[[spp]][[pc]], "0")
c_99_mb <- coef_fun(resDurBAMless1_mix[[spp]][[pc]], "0")
v_99_mb <- vcov_fun(resDurBAMless1_mix[[spp]][[pc]], "0")

c_99_best0 <- coef_fun(resDurBAMless1[[spp]][[pc]], best_99_m0)
v_99_best0 <- vcov_fun(resDurBAMless1[[spp]][[pc]], best_99_m0)
c_99_bestb <- coef_fun(resDurBAMless1_mix[[spp]][[pc]], best_99_mb)
v_99_bestb <- vcov_fun(resDurBAMless1_mix[[spp]][[pc]], best_99_mb)
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
df$y <- rowSums(xtDur[[spp]][rownames(df),])
df$pc <- ifelse(df$PCODE == pc, 1L, 0L)
df99 <- df[df$PCODE != pc,,drop=FALSE]
df1 <- df[df$PCODE == pc,,drop=FALSE]

n99 <- nrow(df99)
n1 <- nrow(df1)
det99 <- sum(df99$y > 0)
det1 <- sum(df1$y > 0)

mod0 <- glm(pc ~ 1, df, family=binomial("logit"))
mod1 <- glm(pc ~ JDAY + TSSR + TSLS, df, family=binomial("logit"))
logLR <- as.numeric(logLik(mod1) - logLik(mod0))

## model matrix for project based prediction
X1_1 <- model.matrix(ff[[best_1]], df1)
X1_99 <- model.matrix(ff[[best_99]], df1)


## write a function that calculates point pred & CI depending on rem/mix
pred_fun(X, c, v, type) {
    if (type=="rem") {
        cfs <- mvrnorm(B, c, v)
    } else {
        
    }
}
## calculates CIs and compare -- write compare fun.