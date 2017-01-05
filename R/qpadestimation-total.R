##---
##title: "QPAD estimation"
##author: "Peter Solymos"
##date: "May 7, 2015"
##output:
##  pdf_document:
##    toc: true
##    toc_depth: 2
##---

### Preliminaries

## Define root folder where data are stored
ROOT <- "e:/peter/bam/May2015"

## Load required packages
library(mefa4)
library(pbapply)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

sppt <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")
rownames(sppt) <- sppt$Species_ID
spp_singing <- sort(rownames(sppt)[sppt$Singing_birds])

## Load preprocesses data
#load(file.path(ROOT, "out", "new_offset_data_package_2016-03-21.Rdata"))
load(file.path(ROOT, "out", "new_offset_data_package_2016-12-01.Rdata"))

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","YEAR","PCODE")]
pkDur$TSSR2 <- pkDur$TSSR^2
pkDur$JDAY2 <- pkDur$JDAY^2
pkDur$DSLS <- pkDur$TSLS
pkDur$DSLS2 <- pkDur$DSLS^2
## JDAY outliers
pkDur$JDAY[pkDur$JDAY > 0.55] <- NA
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
pkDur$TSLS <- NULL

## save this and spatial/yearly variation for next version
## make sure to leave all TSSR (i.e. after noon) in the data
if (FALSE) {
pkDur$TSSRsin <- sin(pkDur$TSSR_orig * 2 * pi)
pkDur$TSSRcos <- cos(pkDur$TSSR_orig * 2 * pi)
pkDur$TSSRsin2 <- sin(pkDur$TSSR_orig * 2 * pi)^2
pkDur$TSSRcos2 <- cos(pkDur$TSSR_orig * 2 * pi)^2

NAMES <- list(
    "1"=c("INTERCEPT", "JDAY", "JDAY2", "TSSRsin", "TSSRcos"),
    "2"=c("INTERCEPT", "JDAY", "JDAY2", "TSSRsin", "TSSRcos", "TSSRsin2", "TSSRcos2"),
    "3"=c("INTERCEPT", "JDAY", "JDAY2", "JDAY3", "TSSRsin", "TSSRcos"),
    "4"=c("INTERCEPT", "JDAY", "JDAY2", "JDAY3", "TSSRsin", "TSSRcos", "TSSRsin2", "TSSRcos2"),
    "5"=c("INTERCEPT", "JDAY", "JDAY2", "JDAY3"),
    "6"=c("INTERCEPT", "TSSRsin", "TSSRcos", "TSSRsin2", "TSSRcos2"))
ff <- list(
    "1"= ~ JDAY + JDAY2 + TSSRsin + TSSRcos,
    "2"= ~ JDAY + JDAY2 + TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2,
    "3"= ~ JDAY + JDAY2 + JDAY3 + TSSRsin + TSSRcos,
    "4"= ~ JDAY + JDAY2 + JDAY3 + TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2,
    "5"= ~ JDAY + JDAY2 + JDAY3,
    "6"= ~ TSSRsin + TSSRcos + TSSRsin2 + TSSRcos2)
}

## models to consider
NAMES <- list(
    "0"="(Intercept)",
    "1"=c("(Intercept)", "JDAY"),
    "2"=c("(Intercept)", "TSSR"),
    "3"=c("(Intercept)", "JDAY", "JDAY2"),
    "4"=c("(Intercept)", "TSSR", "TSSR2"),
    "5"=c("(Intercept)", "JDAY", "TSSR"),
    "6"=c("(Intercept)", "JDAY", "JDAY2", "TSSR"),
    "7"=c("(Intercept)", "JDAY", "TSSR", "TSSR2"),
    "8"=c("(Intercept)", "JDAY", "JDAY2", "TSSR", "TSSR2"),
    "9"=c("(Intercept)", "DSLS"),
    "10"=c("(Intercept)", "DSLS", "DSLS2"),
    "11"=c("(Intercept)", "DSLS", "TSSR"),
    "12"=c("(Intercept)", "DSLS", "DSLS2", "TSSR"),
    "13"=c("(Intercept)", "DSLS", "TSSR", "TSSR2"),
    "14"=c("(Intercept)", "DSLS", "DSLS2", "TSSR", "TSSR2"))
ff <- list(
    ~ 1,
    ~ JDAY,
    ~ TSSR,
    ~ JDAY + JDAY2,
    ~ TSSR + TSSR2,
    ~ JDAY + TSSR,
    ~ JDAY + JDAY2 + TSSR,
    ~ JDAY + TSSR + TSSR2,
    ~ JDAY + JDAY2 + TSSR + TSSR2,
    ~ DSLS,
    ~ DSLS + DSLS2,
    ~ DSLS + TSSR,
    ~ DSLS + DSLS2 + TSSR,
    ~ DSLS + TSSR + TSSR2,
    ~ DSLS + DSLS2 + TSSR + TSSR2)
names(ff) <- 0:14

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
#xtDur[["NONE"]] <- NULL
xtDur <- xtDur[intersect(names(xtDur), spp_singing)]

tot <- xtDur[[1]]
for (i in 2:length(xtDur)) {
    tot <- tot + xtDur[[i]]
}
xtDur <- list(TOTA=tot)

fitDurFun <- function(spp, fit=TRUE, type=c("rem","mix")) {
    rn <- intersect(rownames(pkDur), rownames(xtDur[[spp]]))
    X0 <- pkDur[rn,]
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
        res <- list(Y=Y, D=D, n=n, pkey=rownames(X))
    }
    res
}


SPP <- names(xtDur)
resDur <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("Singing rate data check for", SPP[i], "\n")
    flush.console()
    resDur[[i]] <- try(fitDurFun(SPP[i], FALSE))
}
names(resDur) <- SPP
resDurOK <- resDur[!sapply(resDur, inherits, "try-error")]
c(OK=length(resDurOK), failed=length(resDur)-length(resDurOK), all=length(resDur))
t(sapply(resDurOK, "[[", "n"))
resDurData <- resDurOK

## estimate species with data
SPP <- names(resDurOK)
resDur <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("Singing rate estimation for", SPP[i], date(), "\n")
    flush.console()
    resDur[[i]] <- try(fitDurFun(SPP[i], TRUE, type="rem"))
}
names(resDur) <- SPP
resDurOK <- resDur[!sapply(resDur, inherits, "try-error")]
c(OK=length(resDurOK), failed=length(resDur)-length(resDurOK), all=length(resDur))
resDur <- resDurOK

save(resDur, resDurData,
    file=file.path(ROOT, "out", "TOTA_estimates_SRA_QPAD_v2016.Rdata"))

### Distance sampling

## non NA subset for distance related estimates
pkDis <- dat[,c("PKEY","TREE","TREE3","HAB_NALC1","HAB_NALC2","DISMETH")]
pkDis <- droplevels(pkDis[rowSums(is.na(pkDis)) == 0,])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDis <- droplevels(pkDis[pkDis$DISMETH != "W",])
pkDis$CTREE <- pkDis$TREE3

pkDis$WNALC <- pkDis$HAB_NALC2
levels(pkDis$WNALC)[levels(pkDis$WNALC) %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
pkDis$NALC <- pkDis$WNALC
levels(pkDis$NALC)[levels(pkDis$NALC) %in% c("Wet")] <- "Open"

pkDis$WNALCTREE <- pkDis$HAB_NALC1
levels(pkDis$WNALCTREE)[levels(pkDis$WNALCTREE) %in% c("Agr","Barren","Devel",
    "Grass", "Shrub","ConifOpen","DecidOpen","MixedOpen")] <- "Open"
pkDis$NALCTREE <- pkDis$WNALCTREE
levels(pkDis$NALCTREE)[levels(pkDis$NALCTREE) %in% c("Wet")] <- "Open"

pkDis$LCC2 <- as.factor(ifelse(pkDis$WNALC %in% c("Open", "Wet"), "OpenWet", "Forest"))
pkDis$LCC4 <- pkDis$WNALC
levels(pkDis$LCC4) <- c(levels(pkDis$LCC4), "DecidMixed")
pkDis$LCC4[pkDis$WNALC %in% c("Decid", "Mixed")] <- "DecidMixed"
pkDis$LCC4 <- droplevels(pkDis$LCC4)
pkDis$LCC4 <- relevel(pkDis$LCC4, "DecidMixed")

## models to consider

if (FALSE) {
NAMES <- list(
    "0"="INTERCEPT",
    "1"=c("INTERCEPT", "TREE"),
    "2"=c("INTERCEPT", "NALCOpen", "NALCDecid", "NALCMixed"),
    "3"=c("INTERCEPT", "WNALCOpen", "WNALCDecid", "WNALCMixed", "WNALCWet"),
    "4"=c("INTERCEPT", "NALCOpen", "NALCDecid", "NALCMixed", "TREE"),
    "5"=c("INTERCEPT", "WNALCOpen", "WNALCDecid", "WNALCMixed", "WNALCWet",
        "TREE"))
#    "2"=c("INTERCEPT", "CTREESparse", "CTREEDense"),
#    "7"=c("INTERCEPT", "NALCTREEOpen", "NALCTREEConifSparse", "NALCTREEDecidDense",
#        "NALCTREEDecidSparse", "NALCTREEMixedDense", "NALCTREEMixedSparse"),
#    "8"=c("INTERCEPT", "WNALCTREEOpen", "WNALCTREEConifSparse", "WNALCTREEDecidDense",
#        "WNALCTREEDecidSparse", "WNALCTREEMixedDense", "WNALCTREEMixedSparse",
#        "WNALCTREEWet"))
## keep: 0, 1, 3, 4, 5, 6
ff <- list(
    ~ 1, # 0
    ~ TREE, # 1
    ~ NALC, # 2
    ~ WNALC, # 3
    ~ NALC + TREE, # 4
    ~ WNALC + TREE) # 5
#    ~ CTREE, # 2
#    ~ NALCTREE, # 7
#    ~ WNALCTREE) # 8
}

NAMES <- list(
    "0"="(Intercept)",
    "1"=c("(Intercept)", "TREE"),
    "2"=c("(Intercept)", "LCC2OpenWet"),
    "3"=c("(Intercept)", "LCC4Conif", "LCC4Open", "LCC4Wet"),
    "4"=c("(Intercept)", "LCC2OpenWet", "TREE"),
    "5"=c("(Intercept)", "LCC4Conif", "LCC4Open", "LCC4Wet", "TREE"))
ff <- list(
    ~ 1, # 0
    ~ TREE, # 1
    ~ LCC2, # 2
    ~ LCC4, # 3
    ~ LCC2 + TREE, # 4
    ~ LCC4 + TREE) # 5
names(ff) <- 0:5

## crosstab for species
xtDis <- Xtab(ABUND ~ PKEY + dis + SPECIES, pc)
#xtDis[["NONE"]] <- NULL
xtDis <- xtDis[intersect(names(xtDis), spp_singing)]

tot <- xtDis[[1]]
for (i in 2:length(xtDis)) {
    tot <- tot + xtDis[[i]]
}
xtDis <- list(TOTA=tot)

fitDisFun <- function(spp, fit=TRUE) {
    rn <- intersect(rownames(pkDis), rownames(xtDis[[spp]]))
    X0 <- pkDis[rn,]
    Y0 <- as.matrix(xtDis[[spp]][rn,])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdis$x)))
    stopifnot(length(setdiff(levels(X0$DISMETH), rownames(ltdis$end))) == 0)
    ## interval end matrix
    D <- ltdis$end[match(X0$DISMETH, rownames(ltdis$end)),]
#    D <- D / 100 # 100 m units
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
    Yid <- ltdis$id[match(X$DISMETH, rownames(ltdis$id)),]
    Y <- matrix(NA, n, ncol(ltdis$id))
    for (i in seq_len(n)) {
        w <- Yid[i,]
        w <- w[!is.na(w)]
        Y[i,seq_len(length(w))] <- Y0[i,w]
    }
    if (fit) {
        res <- list()
        for (i in seq_len(length(ff))) {
            f <- as.formula(paste0("Y | D ", paste(as.character(ff[[i]]), collapse=" ")))
            mod <- try(cmulti(f, X, type="dis"))
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
        res <- list(Y=Y, D=D, n=n, pkey=rownames(X))
    }
    res
}


SPP <- names(xtDis)
resDis <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("EDR data check for", SPP[i], "\n")
    flush.console()
    resDis[[i]] <- try(fitDisFun(SPP[i], FALSE))
}
names(resDis) <- SPP
resDisOK <- resDis[!sapply(resDis, inherits, "try-error")]
c(OK=length(resDisOK), failed=length(resDis)-length(resDisOK), all=length(resDis))
t(sapply(resDisOK, "[[", "n"))
resDisData <- resDisOK

## estimate species with data
SPP <- names(resDisOK)
resDis <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("EDR estimation for", SPP[i], date(), "\n")
    flush.console()
    resDis[[i]] <- try(fitDisFun(SPP[i], TRUE))
}
names(resDis) <- SPP
resDisOK <- resDis[!sapply(resDis, inherits, "try-error")]
c(OK=length(resDisOK), failed=length(resDis)-length(resDisOK), all=length(resDis))
resDis <- resDisOK

save(resDis, resDisData,
    file=file.path(ROOT, "out", "TOTA_estimates_EDR_QPAD_v2016.Rdata"))


### Putting things together
## need to load/recreate pkDis!

ROOT <- "e:/peter/bam/May2015"

## n.min is threshold above which all models are considered
## n.con is threshold above which the 0 constant model is considered
n.con <- 25#50
n.min <- 75#50
n.min.class <- 5 # min number of detections within each class

type <- "rem"
load(file.path(ROOT, "out", "TOTA_estimates_SRA_QPAD_v2016.Rdata"))
load(file.path(ROOT, "out", "TOTA_estimates_EDR_QPAD_v2016.Rdata"))
e <- new.env()
load(file.path(ROOT, "out", "new_offset_data_package_2016-12-01.Rdata"), envir=e)

## 0/1 table for successful model fit
edr_mod <- t(sapply(resDis, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))
sra_mod <- t(sapply(resDur, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))
tmp <- union(rownames(edr_mod), rownames(sra_mod))
edr_models <- matrix(0L, length(tmp), ncol(edr_mod))
dimnames(edr_models) <- list(tmp, colnames(edr_mod))
edr_models[rownames(edr_mod),] <- edr_mod
sra_models <- matrix(0L, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod
## data deficient cases: dropped factor levels
## need to check length of NAMES and length of COEF --> correct factor level
## designation is possible, if != bump up AIC/BIC
for (spp in rownames(sra_mod)) {
    for (mid in colnames(sra_models)) {
        if (!inherits(resDur[[spp]][[mid]], "try-error")) {
            if (type == "rem")
                lcf <- length(resDur[[spp]][[mid]]$coefficients)
            if (type == "mix")
                lcf <- length(resDur[[spp]][[mid]]$coefficients) - 1
            lnm <- length(resDur[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("SRA conflict for", spp, "model", mid,
                "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                sra_models[spp,mid] <- 0
            }
        } else {
            resDur[[spp]][[mid]] <- structure("Error", class = "try-error")
            #attributes(resDur[[spp]][[mid]]) <- NULL
            #class(resDur[[spp]][[mid]]) <- "try-error"
        }
    flush.console()
    }
}
for (spp in rownames(edr_mod)) {
    ## data for checking detections in classes
    Dat <- pkDis[resDisData[[spp]]$pkey,c("LCC2","LCC4")]
    for (mid in colnames(edr_models)) {
        if (!inherits(resDis[[spp]][[mid]], "try-error")) {
            lcf <- length(resDis[[spp]][[mid]]$coefficients)
            lnm <- length(resDis[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("EDR conflict for", spp, "model", mid,
                "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                edr_models[spp,mid] <- 0
            } else {
                if (mid %in% c("2", "4") && min(table(Dat$LCC2)) < n.min.class) {
                    cat("EDR LCC2 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("3", "5") && min(table(Dat$LCC4)) < n.min.class) {
                    cat("EDR LCC4 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("1", "4", "5") &&
                    resDis[[spp]][[mid]]$coefficients["log.tau_TREE"] > 0) {
                    cat("EDR TREE > 0 issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
            }
        } else {
            resDis[[spp]][[mid]] <- structure("Error", class = "try-error")
            attributes(resDis[[spp]][[mid]]) <- NULL
            class(resDis[[spp]][[mid]]) <- "try-error"
        }
    flush.console()
    }
}
## no constant model means exclude that species
edr_models[edr_models[,1]==0,] <- 0L
sra_models[sra_models[,1]==0,] <- 0L

phi0 <- sapply(resDur, function(z) exp(z[["0"]]$coef))
names(phi0) <- names(resDur)
sra_models[names(phi0)[phi0 < 0.01],] <- 0L

sum(edr_models)
sum(sra_models)

colSums(edr_models)
colSums(sra_models)

edr_nmod <- ncol(edr_mod)
sra_nmod <- ncol(sra_mod)

## sample sizes
edr_n <- sra_n <- numeric(length(tmp))
names(edr_n) <- names(sra_n) <- tmp
edr_nn <- sapply(resDis, function(z) ifelse(inherits(z[["0"]], "try-error"),
    NA, z[["0"]]$nobs))
edr_n[names(edr_nn)] <- edr_nn
sra_nn <- sapply(resDur, function(z) ifelse(inherits(z[["0"]], "try-error"),
    NA, z[["0"]]$nobs))
sra_n[names(sra_nn)] <- sra_nn

## exclude all models for species with < n.con observations
sra_models[sra_n < n.con, ] <- 0L
edr_models[edr_n < n.con, ] <- 0L
## exclude non-constant model for species with n.con < observations < n.min
sra_models[sra_n < n.min & sra_n >= n.con, 2:ncol(sra_models)] <- 0L
edr_models[edr_n < n.min & edr_n >= n.con, 2:ncol(edr_models)] <- 0L
table(rowSums(sra_models))
table(rowSums(edr_models))
table(rowSums(sra_models), rowSums(edr_models))

## spp to keep
spp <- tmp[rowSums(edr_models) > 0 & rowSums(sra_models) > 0]
length(spp)

edr_models <- edr_models[spp,,drop=FALSE]
sra_models <- sra_models[spp,,drop=FALSE]
edr_n <- edr_n[spp]
sra_n <- sra_n[spp]

## no. of parameters

edr_df <- sapply(resDis[["TOTA"]][1:edr_nmod], "[[", "p")
sra_df <- sapply(resDur[["TOTA"]][1:sra_nmod], "[[", "p")

## estimates
edr_estimates <- resDis[spp]
sra_estimates <- resDur[spp]

## species table
spp_table <- data.frame(spp="TOTA",
    scientific_name="Multicantare aviarium",
    common_name="Total Singing Bird Abundance")
rownames(spp_table) <- "TOTA"

## get variable names for different models
sra_list <- sapply(sra_estimates[["TOTA"]], function(z) paste(z$names, collapse=" + "))
edr_list <- sapply(edr_estimates[["TOTA"]], function(z) paste(z$names, collapse=" + "))

## loglik values
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}
edr_loglik <- edr_models
edr_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:edr_nmod) { # models
        if (edr_models[i,j] > 0)
            edr_loglik[i,j] <- resDis[[i]][[j]]$loglik
    }
}

## AIC/BIC
sra_aic <- sra_aicc <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_aicc[] <- Inf
sra_bic[] <- Inf
edr_aic <- edr_aicc <- edr_bic <- edr_loglik
edr_aic[] <- Inf
edr_aicc[] <- Inf
edr_bic[] <- Inf
for (i in spp) {
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_aicc[i,] <- sra_aic[i,] + (2*sra_df*(sra_df+1)) / (sra_n[i]-sra_df-1)
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
    edr_aic[i,] <- -2*edr_loglik[i,] + 2*edr_df
    edr_aicc[i,] <- edr_aic[i,] + (2*edr_df*(edr_df+1)) / (edr_n[i]-edr_df-1)
    edr_bic[i,] <- -2*edr_loglik[i,] + log(edr_n[i])*edr_df
}

## model ranking
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA
edr_aicrank <- t(apply(edr_aic, 1, rank))*edr_models
edr_aicrank[edr_aicrank==0] <- NA

sra_aiccrank <- t(apply(sra_aicc, 1, rank))*sra_models
sra_aiccrank[sra_aiccrank==0] <- NA
edr_aiccrank <- t(apply(edr_aicc, 1, rank))*edr_models
edr_aiccrank[edr_aiccrank==0] <- NA

sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA
edr_bicrank <- t(apply(edr_bic, 1, rank))*edr_models
edr_bicrank[edr_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_aiccbest <- apply(sra_aiccrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])
edr_aicbest <- apply(edr_aicrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_aiccbest <- apply(edr_aiccrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_bicbest <- apply(edr_bicrank, 1, function(z) colnames(edr_models)[which.min(z)])

table(edr_aicbest, sra_aicbest)
rowSums(table(edr_aicbest, sra_aicbest))
colSums(table(edr_aicbest, sra_aicbest))

table(edr_aiccbest, sra_aiccbest)
rowSums(table(edr_aiccbest, sra_aiccbest))
colSums(table(edr_aiccbest, sra_aiccbest))

table(edr_bicbest, sra_bicbest)
rowSums(table(edr_bicbest, sra_bicbest))
colSums(table(edr_bicbest, sra_bicbest))

version <- "3_tota"

bamcoefs <- list(spp=spp,
    spp_table=spp_table,
    edr_list=edr_list,
    sra_list=sra_list,
    edr_models=edr_models,
    sra_models=sra_models,
    edr_n=edr_n,
    sra_n=sra_n,
    edr_df=edr_df,
    sra_df=sra_df,
    edr_loglik=edr_loglik,
    sra_loglik=sra_loglik,
    edr_aic=edr_aic,
    sra_aic=sra_aic,
    edr_aicc=edr_aicc,
    sra_aicc=sra_aicc,
    edr_bic=edr_bic,
    sra_bic=sra_bic,
    edr_aicrank=edr_aicrank,
    sra_aicrank=sra_aicrank,
    edr_aiccrank=edr_aiccrank,
    sra_aiccrank=sra_aiccrank,
    edr_bicrank=edr_bicrank,
    sra_bicrank=sra_bicrank,
    edr_aicbest=edr_aicbest,
    sra_aicbest=sra_aicbest,
    edr_aiccbest=edr_aiccbest,
    sra_aiccbest=sra_aiccbest,
    edr_bicbest=edr_bicbest,
    sra_bicbest=sra_bicbest,
    edr_estimates=edr_estimates,
    sra_estimates=sra_estimates,
    version=version)
.BAMCOEFS <- list2env(bamcoefs)

save(.BAMCOEFS, file=file.path(ROOT, "out", "TOTA_BAMCOEFS_QPAD_v3.rda"))
#toDump <- as.list(.BAMCOEFS)
#dump("toDump", file=file.path(ROOT, "out", "BAMCOEFS_QPAD_v3.Rdump"))


### Plot species specific results

library(MASS)
library(QPAD)
load(file.path("e:/peter/bam/May2015", "out", "TOTA_BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version # should be 3_tota
plot_BAM_QPAD("TOTA")


## calculating offsets

#### Calculate the offsets (optional) -----------------------------------------
## QPADv3

## Define root folder where data are stored
ROOT <- "e:/peter/bam/Apr2016"
## Load required packages
library(mefa4)
library(QPAD)
## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

load(file.path(ROOT, "out", paste0("data_package_2016-12-01.Rdata")))

library(QPAD)
#load_BAM_QPAD(3)
load(file.path("e:/peter/bam/May2015", "out", "TOTA_BAMCOEFS_QPAD_v3.rda"))
getBAMversion()
sppp <- getBAMspecieslist()

offdat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS","JULIAN")],
    SS[match(PKEY$SS, rownames(SS)),c("TREE","TREE3","HAB_NALC1","HAB_NALC2", "SPRNG")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
offdat$DSLS <- (offdat$JULIAN - offdat$SPRNG) / 365
summary(offdat)

offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
offdat$DSLS2 <- offdat$DSLS^2
offdat$LCC4 <- as.character(offdat$HAB_NALC2)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4,
        c("DecidMixed", "Conif", "Open", "Wet"))
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
offdat$LCC2 <- factor(offdat$LCC2, c("Forest", "OpenWet"))
table(offdat$LCC4, offdat$LCC2)
offdat$MAXDIS <- offdat$MAXDIS / 100

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","DSLS","TSSR2","JDAY2","DSLS2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
#offdat$OKp <- rowSums(is.na(offdat[,c("TSSR","JDAY","DSLS")])) == 0
#offdat$OKq <- rowSums(is.na(offdat[,c("TREE","LCC4")])) == 0
#Xp <- model.matrix(~TSSR+TSSR2+JDAY+JDAY2+DSLS+DSLS2, offdat[offdat$OKp,])
#Xq <- model.matrix(~LCC2+LCC4+TREE, offdat[offdat$OKq,])

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
cat(spp, "\n");flush.console()
p <- rep(NA, nrow(offdat))
A <- q <- p

## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
#vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0

p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))

phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))

ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

OFF[,spp] <- log(p) + log(A) + log(q)

}

(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,])) # BARS GCSP
which(!is.finite(Ra[2,]))

save(OFF,
    file=file.path(ROOT, "out", "TOTA_offsets-v3_2016-12-01.Rdata"))

