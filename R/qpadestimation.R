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
ROOT <- "c:/bam/May2015"

## Load required packages
library(mefa4)
library(pbapply)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2015-05-07.Rdata"))

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","DURMETH")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])

## models to consider
NAMES <- list(
    "0"="INTERCEPT",
    "1"=c("INTERCEPT", "JDAY"),
    "2"=c("INTERCEPT", "TSSR"),
    "3"=c("INTERCEPT", "JDAY", "JDAY2"),
    "4"=c("INTERCEPT", "TSSR", "TSSR2"),
    "5"=c("INTERCEPT", "JDAY", "TSSR"),
    "6"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR"),
    "7"=c("INTERCEPT", "JDAY", "TSSR", "TSSR2"),
    "8"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR", "TSSR2"))
ff <- list(
    ~ 1,
    ~ JDAY,
    ~ TSSR,
    ~ JDAY + I(JDAY^2),
    ~ TSSR + I(TSSR^2),
    ~ JDAY + TSSR,
    ~ JDAY + I(JDAY^2) + TSSR,
    ~ JDAY + TSSR + I(TSSR^2),
    ~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2))

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)

fitDurFun <- function(spp, fit=TRUE) {
    ## get nonzero PCs
    Y0 <- as.matrix(xtDur[[spp]][rownames(pkDur),])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdur$x)))
    ## interval end matrix
    D <- ltdur$end[match(pkDur$DURMETH, rownames(ltdur$end)),]
    ## exclude 0 sum and <1 interval rows
    iob <- rowSums(Y0) > 0 & rowSums(!is.na(D)) > 1
    if (sum(iob)==0)
        return(structure("0 observation with multiple duration (1)", 
            class="try-error"))
    if (sum(iob)==1)
        return(structure("1 observation with multiple duration (2)", 
            class="try-error"))
    X <- droplevels(pkDur[iob,])
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
    
    ## integer mode -- faster, but DO NOT use for intervals (<1)
    if (fit) {
        res <- list()
        for (i in seq_len(length(ff))) {
            f <- as.formula(paste0("Y | D ", paste(as.character(ff[[i]]), collapse=" ")))
            mod <- try(cmulti(f, X, type="rem"))
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
    resDur[[i]] <- try(fitDurFun(SPP[i], TRUE))
}
names(resDur) <- SPP
resDurOK <- resDur[!sapply(resDur, inherits, "try-error")]
c(OK=length(resDurOK), failed=length(resDur)-length(resDurOK), all=length(resDur))
resDur <- resDurOK

save(resDur, resDurData, 
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2015.Rdata"))

### Distance sampling

## non NA subset for distance related estimates
pkDis <- dat[,c("PKEY","TREE","TREE3","HAB_NALC1","HAB_NALC2","DISMETH")]
pkDis <- droplevels(pkDis[rowSums(is.na(pkDis)) == 0,])
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

## models to consider

NAMES <- list(
    "0"="INTERCEPT",
    "1"=c("INTERCEPT", "TREE"),
    "2"=c("INTERCEPT", "CTREESparse", "CTREEDense"),
    "3"=c("INTERCEPT", "NALCOpen", "NALCDecid", "NALCMixed"),
    "4"=c("INTERCEPT", "WNALCOpen", "WNALCDecid", "WNALCMixed", "WNALCWet"),
    "5"=c("INTERCEPT", "NALCOpen", "NALCDecid", "NALCMixed", "TREE"),
    "6"=c("INTERCEPT", "WNALCOpen", "WNALCDecid", "WNALCMixed", "WNALCWet", 
        "TREE"),
    "7"=c("INTERCEPT", "NALCTREEOpen", "NALCTREEConifSparse", "NALCTREEDecidDense", 
        "NALCTREEDecidSparse", "NALCTREEMixedDense", "NALCTREEMixedSparse"),
    "8"=c("INTERCEPT", "WNALCTREEOpen", "WNALCTREEConifSparse", "WNALCTREEDecidDense", 
        "WNALCTREEDecidSparse", "WNALCTREEMixedDense", "WNALCTREEMixedSparse", 
        "WNALCTREEWet"))
ff <- list(
    ~ 1, # 0
    ~ TREE, # 1
    ~ CTREE, # 2
    ~ NALC, # 3
    ~ WNALC, # 4
    ~ NALC + TREE, # 5
    ~ WNALC + TREE, # 6
    ~ NALCTREE, # 7
    ~ WNALCTREE) # 8

## crosstab for species
xtDis <- Xtab(ABUND ~ PKEY + dis + SPECIES, pc)

fitDisFun <- function(spp, fit=TRUE) {
    ## get nonzero PCs
    Y0 <- as.matrix(xtDis[[spp]][rownames(pkDis),])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdis$x)))
    ## interval end matrix
    D <- ltdis$end[match(pkDis$DISMETH, rownames(ltdis$end)),]
    D <- D / 100 # 100 m units
    ## exclude 0 sum and <1 interval rows
    iob <- rowSums(Y0) > 0 & rowSums(!is.na(D)) > 1
    if (sum(iob)==0)
        return(structure("0 observation with multiple duration (1)", 
            class="try-error"))
    if (sum(iob)==1)
        return(structure("1 observation with multiple duration (2)", 
            class="try-error"))
    X <- droplevels(pkDis[iob,])
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
    
    ## integer mode -- faster, but DO NOT use for intervals (<1)
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
        res <- list(Y=Y, D=D, n=n)
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
    file=file.path(ROOT, "out", "estimates_EDR_QPAD_v2015.Rdata"))


### Putting things together

ROOT <- "c:/bam/May2015"

## n.min is threshold above which all models are considered
## n.con is threshold above which the 0 constant model is considered
n.con <- 25
n.min <- 75

load(file.path(ROOT, "out", "estimates_SRA_QPAD_v2015.Rdata"))
load(file.path(ROOT, "out", "estimates_EDR_QPAD_v2015.Rdata"))

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
for (spp in tmp) {
    for (mid in colnames(sra_models)) {
        if (!inherits(resDur[[spp]][[mid]], "try-error")) {
            lcf <- length(resDur[[spp]][[mid]]$coefficients)
            lnm <- length(resDur[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("SRA conflict for", spp, "model", mid, 
                "( lcf =", lcf, ", lnm =", lnm, "\n")
                sra_models[spp,mid] <- 0
            }
        }
    }
    for (mid in colnames(edr_models)) {
        if (!inherits(resDis[[spp]][[mid]], "try-error")) {
            lcf <- length(resDis[[spp]][[mid]]$coefficients)
            lnm <- length(resDis[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("EDR conflict for", spp, "model", mid, 
                "( lcf =", lcf, ", lnm =", lnm, "\n")
                edr_models[spp,mid] <- 0
            }
        }
    }
}
## no constant model means exclude that species
edr_models[edr_models[,1]==0,] <- 0L
sra_models[sra_models[,1]==0,] <- 0L
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

edr_models <- edr_models[spp,]
sra_models <- sra_models[spp,]
edr_n <- edr_n[spp]
sra_n <- sra_n[spp]

## no. of parameters

edr_df <- sapply(resDis[["OVEN"]][1:edr_nmod], "[[", "p")
sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

## estimates
edr_estimates <- resDis[spp]
sra_estimates <- resDur[spp]

## this needs to be updated !!! ---------------------------------------- FIXME
spp_table <- data.frame(spp=spp, 
    scientific_name=NA,
    common_name=NA)
rownames(spp_table) <- spp
if (FALSE) {
load("lifehist_final.Rdata")
rownames(taxo2) <- taxo2$SPECIES
spp_table <- data.frame(spp=spp, 
    scientific_name=taxo2[spp, "SCIENTIFIC NAME"],
    common_name=taxo2[spp, "ENGLISH NAME"])
rownames(spp_table) <- spp
}

## get variable names for different models
sra_list <- sapply(sra_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))
edr_list <- sapply(edr_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

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
sra_aic <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_bic[] <- Inf
edr_aic <- edr_bic <- edr_loglik
edr_aic[] <- Inf
edr_bic[] <- Inf
for (i in spp) {
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
    edr_aic[i,] <- -2*edr_loglik[i,] + 2*edr_df
    edr_bic[i,] <- -2*edr_loglik[i,] + log(edr_n[i])*edr_df
}

if (FALSE) { ## --------------------------------------------------- CHECK !!!
## constrain TREE (-) estimates
for (i in spp) {
    ## TREE
    if (edr_models[i, "1"] == 1) {
        if (edr_estimates[[i]][["1"]]$coef[2] > 0) {
            edr_aic[i, "1"] <- Inf
            edr_bic[i, "1"] <- Inf
            cat(i, "tree", round(edr_estimates[[i]][["1"]]$coef[2], 4), "\n")
        }
    }
}
}

## model ranking
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA
edr_aicrank <- t(apply(edr_aic, 1, rank))*edr_models
edr_aicrank[edr_aicrank==0] <- NA
sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA
edr_bicrank <- t(apply(edr_bic, 1, rank))*edr_models
edr_bicrank[edr_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])
edr_aicbest <- apply(edr_aicrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_bicbest <- apply(edr_bicrank, 1, function(z) colnames(edr_models)[which.min(z)])

table(edr_aicbest, sra_aicbest)
rowSums(table(edr_aicbest, sra_aicbest))
colSums(table(edr_aicbest, sra_aicbest))

table(edr_bicbest, sra_bicbest)
rowSums(table(edr_bicbest, sra_bicbest))
colSums(table(edr_bicbest, sra_bicbest))

version <- "3"

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
    edr_bic=edr_bic, 
    sra_bic=sra_bic,
    edr_aicrank=edr_aicrank, 
    sra_aicrank=sra_aicrank, 
    edr_bicrank=edr_bicrank, 
    sra_bicrank=sra_bicrank,
    edr_aicbest=edr_aicbest, 
    sra_aicbest=sra_aicbest, 
    edr_bicbest=edr_bicbest, 
    sra_bicbest=sra_bicbest,
    edr_estimates=edr_estimates, 
    sra_estimates=sra_estimates,
    version=version)
.BAMCOEFS <- list2env(bamcoefs)

save(.BAMCOEFS, file=file.path(ROOT, "out", "BAMCOEFS_QPAD_v3.rda"))
