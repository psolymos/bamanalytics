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

ROOT <- "c:/bam/May2015"

R <- 200
#spp <- "OVEN"
level <- 0.9
version <- 3
prob <- c(0, 1) + c(1, -1) * ((1-level)/2)

library(MASS)
library(QPAD)
load(file.path(ROOT, "out", "TOTA_BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version # should be 3_tota

jd <- seq(0.35, 0.55, 0.01) # TSSR
ts <- seq(-0.25, 0.5, 0.01) # JDAY
ls <- seq(0, 0.25, len=length(jd)) # DSLS

xp1 <- expand.grid(JDAY=jd, # ---------- CHECK !!!
    TSSR=ts)
xp1$JDAY2 <- xp1$JDAY^2
xp1$TSSR2 <- xp1$TSSR^2
xp1$Jday <- xp1$JDAY * 365
xp1$Tssr <- xp1$TSSR * 24

xp2 <- expand.grid(DSLS=ls, # ---------- CHECK !!!
    TSSR=ts)
xp2$DSLS2 <- xp2$DSLS^2
xp2$TSSR2 <- xp2$TSSR^2
xp2$Dsls <- xp2$DSLS * 365
xp2$Tssr <- xp2$TSSR * 24

Xp1 <- model.matrix(~., xp1)
#colnames(Xp1)[1] <- "INTERCEPT"
Xp2 <- model.matrix(~., xp2)
#colnames(Xp2)[1] <- "INTERCEPT"

if (getBAMversion() < 3) {
    lc <- seq(1, 5, 1)
    tr <- seq(0, 1, 0.01)
    xq <- expand.grid(LCC=as.factor(lc),
        TREE=tr)
} else {
#    lc2 <- factor(c("Forest", "OpenWet"), c("Forest", "OpenWet"))
    lc <- factor(c("DecidMixed", "Conif", "Open", "Wet"),
        c("DecidMixed", "Conif", "Open", "Wet"))
    tr <- seq(0, 1, 0.1)
    xq <- expand.grid(LCC4=lc, TREE=tr)
    xq$LCC2 <- as.character(xq$LCC4)
    xq$LCC2[xq$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
    xq$LCC2[xq$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
    xq$LCC2 <- factor(xq$LCC2, c("Forest", "OpenWet"))
}
Xq0 <- model.matrix(~., xq)
#colnames(Xq)[1] <- "INTERCEPT"

SPP <- getBAMspecieslist()
cfall <- exp(t(sapply(SPP, function(spp)
    unlist(coefBAMspecies(spp, 0, 0)))))
t <- seq(0, 10, 0.1)
r <- seq(0, 4, 0.05)
pp <- sapply(SPP, function(spp) sra_fun(t, cfall[spp,1]))
qq <- sapply(SPP, function(spp) edr_fun(r, cfall[spp,2]))

pdf(paste0("~/Dropbox/bam/qpad_v3/QPAD_res_v",
    getBAMversion(),".pdf"), onefile=TRUE, width=9, height=12)
for (spp in SPP) {

cat(spp, "\n");flush.console()

## model weights
wp <- selectmodelBAMspecies(spp)$sra$wBIC
wq <- selectmodelBAMspecies(spp)$edr$wBIC
names(wp) <- rownames(selectmodelBAMspecies(spp)$sra)
names(wq) <- rownames(selectmodelBAMspecies(spp)$edr)
nsra <- selectmodelBAMspecies(spp)$sra$nobs[1]
nedr <- selectmodelBAMspecies(spp)$edr$nobs[1]

## constant model
#lphi0 <- coefBAMspecies(spp, 0, 0)$sra
#lphise0 <- sqrt(vcovBAMspecies(spp, 0, 0)$sra[1])
#lphipi0 <- quantile(rnorm(R, lphi0, lphise0), prob)
#p <- cbind(est=sra_fun(t, exp(lphi0)))#,
#           pi1=sra_fun(t, exp(lphipi0[1])),
#           pi2=sra_fun(t, exp(lphipi0[2])))
#ltau0 <- coefBAMspecies(spp, 0, 0)$edr
#ltause0 <- sqrt(vcovBAMspecies(spp, 0, 0)$edr[1])
#ltaupi0 <- quantile(rnorm(R, ltau0, ltause0), prob)
#q <- cbind(est=edr_fun(r, exp(ltau0)))#,
#           pi1=edr_fun(r, exp(ltaupi0[1])),
#           pi2=edr_fun(r, exp(ltaupi0[2])))

## covariate effects
mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

Xp <- if (mi$sra %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
Xp <- Xp[,names(cfi$sra),drop=FALSE]
lphi1 <- drop(Xp %*% cfi$sra)
#pcf1 <- mvrnorm(R, cfi$sra, vci$sra)
#lphipi1 <- t(apply(apply(pcf1, 1, function(z) drop(Xp %*% z)),
#    1, quantile, prob))
#px <- cbind(est=exp(lphi1), exp(lphipi1))
pmat <- matrix(exp(lphi1), length(jd), length(ts))
pmax <- sra_fun(10, max(exp(lphi1)))
pmat <- sra_fun(3, pmat)
pmax <- 1

Xq <- Xq0[,names(cfi$edr),drop=FALSE]
ltau1 <- drop(Xq %*% cfi$edr)
#pcf1 <- mvrnorm(R, cfi$edr, vci$edr)
#ltaupi1 <- t(apply(apply(pcf1, 1, function(z) drop(Xq %*% z)),
#    1, quantile, prob))
#qx <- cbind(est=exp(ltau1), exp(ltaupi1))
qmat <- matrix(exp(ltau1), length(lc), length(tr))
qmax <- edr_fun(0.5, max(exp(ltau1)))
qmat <- edr_fun(1, qmat)
qmax <- 1

#library(lattice)
#levelplot(px[,1] ~ Jday * Tssr, xp)

op <- par(las=1, mfrow=c(3,2))

barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", nsra, ") v", getBAMversion()),
    ylab="Model weight", xlab="Model ID")
barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", nedr, ") v", getBAMversion()),
    ylab="Model weight", xlab="Model ID")

plot(t, pp[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
matlines(t, pp, col="grey", lwd=1, lty=1)
lines(t, pp[,spp], col=1, lwd=2)

plot(r*100, qq[,spp], type="n", ylim=c(0,1),
     xlab="Point count radius (m)",
     ylab="Probability of detection")
#polygon(100*c(r, rev(r)), c(q[,2], rev(q[,3])),
#        col="grey", border="grey")
matlines(r*100, qq, col="grey", lwd=1, lty=1)
lines(r*100, qq[,spp], col=1, lwd=2)
abline(v=cfall[spp,2]*100, lty=2)
rug(cfall[,2]*100, side=1, col="grey")
box()

xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mi$sra))
box()
image(1:nlevels(xq$LCC4), tr*100, qmat,
      col = rev(grey(seq(0, qmax, len=12))), axes=FALSE,
      xlab="Land cover types", ylab="Percent tree cover",
      main=paste("Best model:", mi$edr))
if (version < 3)
    axis(1, 1:5, c("DConif","DDecid","SConif","SDecid","Open"))
if (version > 2)
    axis(1, 1:nlevels(xq$LCC4), levels(xq$LCC4))
axis(2)
box()

par(op)
}
dev.off()

## comparing v2 and v3 results

library(QPAD)

load_BAM_QPAD(2)
getBAMversion()
spp2 <- getBAMspecieslist()
n2 <- cbind(sra.n=.BAMCOEFS$sra_n, edr.n=.BAMCOEFS$edr_n)
est2 <- exp(t(sapply(spp2, function(i) unlist(coefBAMspecies(i)))))

load_BAM_QPAD(3)
getBAMversion()
spp3 <- getBAMspecieslist()
n3 <- cbind(sra.n=.BAMCOEFS$sra_n, edr.n=.BAMCOEFS$edr_n)
est3 <- exp(t(sapply(spp3, function(i) unlist(coefBAMspecies(i)))))

## "YWAR" is "YEWA"
setdiff(spp2,spp3)
setdiff(spp3,spp2)

spp2[spp2=="YWAR"] <- "YEWA"
rownames(n2)[rownames(n2)=="YWAR"] <- "YEWA"
rownames(est2)[rownames(est2)=="YWAR"] <- "YEWA"

setdiff(spp2,spp3)

n2 <- n2[spp2,]
est2 <- est2[spp2,]
n3v <- n3[spp2,]
est3v <- est3[spp2,]

par(mfrow=c(2,2))
plot(n2[,1], n3v[,1], main="n sra", xlab="v2", ylab="v3");abline(0,1)
abline(0,2,lty=2)
abline(0,4,lty=2)
abline(0,8,lty=2)
plot(n2[,2], n3v[,2], main="n edr", xlab="v2", ylab="v3");abline(0,1)
abline(0,2,lty=2)
abline(0,4,lty=2)
abline(0,8,lty=2)
plot(est2[,1], est3v[,1], main="srate", xlab="v2", ylab="v3",
    cex=sqrt((n3v[,1]-n2[,1])/n2[,1]))
abline(0,1)
plot(est2[,2], est3v[,2], main="EDR", xlab="v2", ylab="v3",
    cex=sqrt((n3v[,2]-n2[,2])/n2[,2]))
abline(0,1)

e <- new.env()
load(file.path(ROOT, "out", "new_offset_data_package_2015-05-14.Rdata"), envir=e)
TAX <- e$TAX
rownames(TAX) <- TAX$Species_ID
TAX <- droplevels(TAX[spp3,])
rm(e)

tab <- data.frame(getBAMspeciestable(),
    v2=n2[match(spp3,spp2),], v2=est2[match(spp3,spp2),],
    v3=n3, v3=est3,
    Order=TAX$Order, Family=TAX$Family_Sci)
write.csv(tab, file=file.path(ROOT, "out", "offsets-spp-list.csv"), row.names=FALSE)

data.frame(Order=sort(table(tab$Order)))
data.frame(Family=sort(table(tab$Family)))

cbind(table(tab$Order, tab$v3.sra.n > 74), table(tab$Order))


## exploration

wp <- t(sapply(SPP, function(spp) selectmodelBAMspecies(spp)$sra$weights))
colnames(wp) <- 0:14
wq <- t(sapply(SPP, function(spp) selectmodelBAMspecies(spp)$edr$weights))
colnames(wq) <- 0:5
np <- sapply(SPP, function(spp) selectmodelBAMspecies(spp)$sra$nobs[1])
nq <- sapply(SPP, function(spp) selectmodelBAMspecies(spp)$edr$nobs[1])
Hp <- apply(wp, 1, function(z) sum(z^2))
Hq <- apply(wq, 1, function(z) sum(z^2))

par(mfrow=c(1,2))
plot(np, wp[,"0"], pch=21, cex=0.5+Hp, log="x")
summary(np[wp[,"0"] > 0.95])
summary(np[wp[,"0"] <= 0.95])

plot(nq, wq[,"0"], pch=21, cex=0.5+Hq, log="x")
summary(nq[wq[,"0"] > 0.95])
summary(nq[wq[,"0"] <= 0.95])


## sra 0 vs b models

e <- new.env()
load(file.path(ROOT, "out", "BAMCOEFS_QPAD_v3_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS
.BAMCOEFSmix$version

aic0 <- .BAMCOEFS$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
colnames(aic0) <- paste0("m0_", colnames(aic0))
colnames(aicb) <- paste0("mb_", colnames(aicb))
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))
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

best <- colnames(aic)[apply(aic, 1, which.min)]
data.frame(table(best))
by(np, best %in% c("m0_0", "mb_0"), summary)


cfall0 <- sapply(SPP, function(spp) exp(coefBAMspecies(spp, 0, 0)$sra))
dim(cfall0) <- c(length(SPP), 1)
dimnames(cfall0) <- list(SPP, "phi_0")
cfallb <- t(sapply(SPP, function(spp) {
    tmp <- unname(.BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients)
    c(phi_b=exp(tmp[1]), c=plogis(tmp[2]))
    }))

t <- seq(0, 10, 0.1)
pp0 <- sapply(SPP, function(spp) 1-exp(-t*cfall0[spp,"phi_0"]))
ppb <- sapply(SPP, function(spp) 1-cfallb[spp,"c"]*exp(-t*cfallb[spp,"phi_b"]))

f <- function(spp) {
    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    !inherits(try(mvrnorm(R, cfi0b, vci0b)), "try-error")
}
OK <- sapply(SPP, f)
table(OK)
SPP <- SPP[OK]

w00 <- waic0[SPP,"m0_0"]
w0b <- waicb[SPP,"mb_0"]
w0 <- pmax(w00, w0b)
H0 <- apply(waic0[SPP,], 1, function(z) sum(z^2))
Hb <- apply(waicb[SPP,], 1, function(z) sum(z^2))
H <- apply(waic[SPP,], 1, function(z) sum(z^2))
n <- np[SPP]

par(mfrow=c(1,3))
plot(n, w00, pch=21, cex=1+2*H0, log="x")
lines(n[order(n)], fitted(glm(w00 ~ n, family=binomial, weights=H0))[order(n)], col=2)
plot(n, w0b, pch=21, cex=1+2*Hb, log="x")
lines(n[order(n)], fitted(glm(w0b ~ n, family=binomial, weights=Hb))[order(n)], col=2)
plot(n, w0, pch=21, cex=1+2*H, log="x")
lines(n[order(n)], fitted(glm(w0 ~ n, family=binomial, weights=H))[order(n)], col=2)

## best model is 0
table((1:ncol(waic))[apply(waic[SPP,], 1, which.max)] <= 15)


pdf(paste0("~/Dropbox/bam/qpad_v3/QPAD_res_v",
    getBAMversion(),"_mix.pdf"), onefile=TRUE, width=10, height=12)
for (spp in SPP) {

cat(spp, "\n");flush.console()

## model weights
wp <- waic0[spp,]
wq <- waicb[spp,]
names(wp) <- colnames(waic0)
names(wq) <- colnames(waicb)
wph <- waic[spp,names(wp)]
wqh <- waic[spp,names(wq)]
names(wp) <- 0:14
names(wq) <- 0:14


## covariate effects
mi <- bestmodelBAMspecies(spp, type="AIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

mib <- as.character(0:14)[which.min(aicb[spp,])]
cfib <- .BAMCOEFSmix$sra_estimates[[spp]][[mib]]$coefficients
vcib <- .BAMCOEFSmix$sra_estimates[[spp]][[mib]]$vcov

#     TSSR             JDAY            TSLS
# Min.   :-0.315   Min.   :0.351   Min.   :-0.101
# 1st Qu.: 0.063   1st Qu.:0.433   1st Qu.: 0.103
# Median : 0.149   Median :0.455   Median : 0.131
# Mean   : 0.141   Mean   :0.455   Mean   : 0.133
# 3rd Qu.: 0.234   3rd Qu.:0.479   3rd Qu.: 0.164
# Max.   : 0.520   Max.   :0.641   Max.   : 0.442
# NA's   :8455     NA's   :5804    NA's   :17255
jd <- seq(0.35, 0.55, 0.01)
ts <- seq(-0.3, 0.5, 0.01)
ls <- seq(-0.1, 0.4, len=length(jd))

xp1 <- expand.grid(JDAY=jd, # ---------- CHECK !!!
    TSSR=ts)
xp1$JDAY2 <- xp1$JDAY^2
xp1$TSSR2 <- xp1$TSSR^2
xp1$Jday <- xp1$JDAY * 365
xp1$Tssr <- xp1$TSSR * 24

xp2 <- expand.grid(TSLS=ls, # ---------- CHECK !!!
    TSSR=ts)
xp2$TSLS2 <- xp2$TSLS^2
xp2$TSSR2 <- xp2$TSSR^2
xp2$Tsls <- xp2$TSLS * 365
xp2$Tssr <- xp2$TSSR * 24

Xp1 <- model.matrix(~., xp1)
colnames(Xp1)[1] <- "INTERCEPT"
Xp2 <- model.matrix(~., xp2)
colnames(Xp2)[1] <- "INTERCEPT"

Xp <- if (mi$sra %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
Xp <- Xp[,names(cfi$sra),drop=FALSE]
lphi1 <- drop(Xp %*% cfi$sra)
pmat <- matrix(exp(lphi1), length(jd), length(ts))
pmax <- sra_fun(10, max(exp(lphi1)))
pmat <- sra_fun(3, pmat)
pmax <- 1

Xpb <- if (mib %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
cfib1 <- exp(cfib[1])
cfib2 <- cfib[-1]
names(cfib2) <- sapply(strsplit(names(cfib2), "_"), "[[", 2)
names(cfib2)[1] <- "INTERCEPT"
names(cfib2)[names(cfib2) == "I(TSSR^2)"] <- "TSSR2"
names(cfib2)[names(cfib2) == "I(TSLS^2)"] <- "TSLS2"
names(cfib2)[names(cfib2) == "I(JDAY^2)"] <- "JDAY2"

Xpb <- Xpb[,names(cfib2),drop=FALSE]
lphi1b <- 1-plogis(drop(Xpb %*% cfib2))*exp(-3*cfib1)
pmatb <- matrix(lphi1b, length(jd), length(ts))

## CI for m0, m0
cfi00 <- coefBAMspecies(spp, "0", mi$edr)$sra
vci00 <- drop(vcovBAMspecies(spp, "0", mi$edr)$sra)
phi00 <- exp(rnorm(R, cfi00, sqrt(vci00)))
ci00 <- sapply(phi00, function(z) 1-exp(-t*z))

cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
pcf1b <- mvrnorm(R, cfi0b, vci0b)
ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))


op <- par(las=1, mfrow=c(3,2))

barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", np[spp], ") v", getBAMversion(),
    " wb=", round(sum(wph),2)),
    ylab="Model weight", xlab="Model ID")
barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", np[spp], ") v", getBAMversion(),
    " wb=", round(sum(wqh),2)),
    ylab="Model weight", xlab="Model ID")

plot(t, pp0[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, pp0, col="grey", lwd=1, lty=1)
matlines(t, ci00, col="grey", lwd=1, lty=1)
lines(t, pp0[,spp], col=1, lwd=2)

plot(t, ppb[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, ppb, col="grey", lwd=1, lty=1)
matlines(t, ci0b, col="grey", lwd=1, lty=1)
lines(t, ppb[,spp], col=1, lwd=2)

xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mi$sra))
box()

xval <- if (mib %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmatb,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mib %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mib))
box()

par(op)
}
dev.off()
