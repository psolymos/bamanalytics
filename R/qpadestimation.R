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
load(file.path(ROOT, "out", "new_offset_data_package_2015-05-11.Rdata"))

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
xtDur[["NONE"]] <- NULL

fitDurFun <- function(spp, fit=TRUE) {
    rn <- intersect(rownames(pkDur), rownames(xtDur[[spp]]))
    X0 <- pkDur[rn,]
    Y0 <- as.matrix(xtDur[[spp]][rn,])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdur$x)))
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

## crosstab for species
xtDis <- Xtab(ABUND ~ PKEY + dis + SPECIES, pc)
xtDis[["NONE"]] <- NULL

fitDisFun <- function(spp, fit=TRUE) {
    rn <- intersect(rownames(pkDis), rownames(xtDis[[spp]]))
    X0 <- pkDis[rn,]
    Y0 <- as.matrix(xtDis[[spp]][rn,])
    ## make sure that columns (intervals) match up
    stopifnot(all(colnames(Y0) == colnames(ltdis$x)))
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
        } else {
            resDur[[spp]][[mid]] <- structure("Error", class = "try-error")
            #attributes(resDur[[spp]][[mid]]) <- NULL
            #class(resDur[[spp]][[mid]]) <- "try-error"
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
        } else {
            resDis[[spp]][[mid]] <- structure("Error", class = "try-error")
            attributes(resDis[[spp]][[mid]]) <- NULL
            class(resDis[[spp]][[mid]]) <- "try-error"
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
e <- new.env()
load(file.path(ROOT, "out", "new_offset_data_package_2015-05-11.Rdata"), envir=e)
tax <- e$TAX
rownames(tax) <- tax$Species_ID
spp_table <- data.frame(spp=spp,
    scientific_name=tax[spp, "Scientific_Name"],
    common_name=tax[spp, "English_Name"])
rownames(spp_table) <- spp
spp_table <- droplevels(spp_table)

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
toDump <- as.list(.BAMCOEFS)
dump("toDump", "v3.R")


### Plot species spacific results

R <- 1000
#spp <- "OVEN"
level <- 0.9
version <- 3

prob <- c(0, 1) + c(1, -1) * ((1-level)/2)
library(MASS)
library(detect)
load_BAM_QPAD(1)
.BAMCOEFS$version
if (version > 2)
    load("~/Dropbox/bam/qpad_v3/BAMCOEFS_QPAD_v3.rda")
.BAMCOEFS$version

SPP <- getBAMspecieslist()
cfall <- exp(t(sapply(SPP, function(spp)
    unlist(coefBAMspecies(spp, 0, 0)))))
t <- seq(0, 10, 0.1)
r <- seq(0, 4, 0.05)
pp <- sapply(SPP, function(spp) sra_fun(t, cfall[spp,1]))
qq <- sapply(SPP, function(spp) edr_fun(r, cfall[spp,2]))

pdf(paste0("~/Dropbox/bam/qpad_v3/QPAD_res_v",
    getBAMversion(),".pdf"), onefile=TRUE, width=8, height=12)
for (spp in SPP) {

cat(spp, "\n");flush.console()

## model weights
wp <- selectmodelBAMspecies(spp)$sra$weights
wq <- selectmodelBAMspecies(spp)$edr$weights
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
mi <- bestmodelBAMspecies(spp)
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

##      TSSR             JDAY
## Min.   :-0.315   Min.   :0.351
## 1st Qu.: 0.071   1st Qu.:0.433
## Median : 0.159   Median :0.458
## Mean   : 0.148   Mean   :0.456
## 3rd Qu.: 0.241   3rd Qu.:0.479
## Max.   : 0.520   Max.   :0.548
## NA's   :13971    NA's   :10679
jd <- seq(0.35, 0.55, 0.01)
ts <- seq(-0.3, 0.5, 0.01)
xp <- expand.grid(JDAY=jd, # ---------- CHECK !!!
    TSSR=ts)
xp$JDAY2 <- xp$JDAY^2
xp$TSSR2 <- xp$TSSR^2
xp$Jday <- xp$JDAY * 365
xp$Tssr <- xp$TSSR * 24

Xp <- model.matrix(~., xp)
colnames(Xp)[1] <- "INTERCEPT"
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

if (getBAMversion() < 3) {
    lc <- seq(1, 5, 1)
    tr <- seq(0, 1, 0.01)
    xq <- expand.grid(LCC=as.factor(lc),
        TREE=tr)
} else {
    lc <- seq(1, 5, 1) # WNALC
    tr <- seq(0, 1, 0.1)
    xq <- expand.grid(WNALC=as.factor(lc),
        TREE=tr)
    xq$CTREE <- factor("Open", c("Open","Sparse","Dense"))
    xq$CTREE[xq$TREE > 0.25] <- "Sparse"
    xq$CTREE[xq$TREE > 0.60] <- "Dense"
    levels(xq$WNALC) <- c("Conif",
                          "Open",
                          "Decid",
                          "Mixed",
                          "Wet")
    xq$NALC <- xq$WNALC
    levels(xq$NALC)[levels(xq$NALC) == "Wet"] <- "Open"
    xq$WNALCTREE <- factor(NA, c("ConifDense",
                                 "Open",
                                 "ConifSparse",
                                 "DecidDense",
                                 "DecidSparse",
                                 "MixedDense",
                                 "MixedSparse",
                                 "Wet"))
    xq$WNALCTREE[xq$WNALC == "Open"] <- "Open"
    xq$WNALCTREE[xq$WNALC == "Wet"] <- "Wet"
    xq$WNALCTREE[xq$WNALC != "Wet" & xq$TREE <= 0.25] <- "Open"
    xq$WNALCTREE[xq$WNALC == "Conif" & xq$TREE > 0.25] <- "ConifSparse"
    xq$WNALCTREE[xq$WNALC == "Decid" & xq$TREE > 0.25] <- "DecidSparse"
    xq$WNALCTREE[xq$WNALC == "Mixed" & xq$TREE > 0.25] <- "MixedSparse"
    xq$WNALCTREE[xq$WNALC == "Conif" & xq$TREE > 0.60] <- "ConifDense"
    xq$WNALCTREE[xq$WNALC == "Decid" & xq$TREE > 0.60] <- "DecidDense"
    xq$WNALCTREE[xq$WNALC == "Mixed" & xq$TREE > 0.60] <- "MixedDense"
    xq$NALCTREE <- xq$WNALCTREE
    levels(xq$NALCTREE)[levels(xq$NALCTREE) == "Wet"] <- "Open"

}
Xq <- model.matrix(~., xq)
colnames(Xq)[1] <- "INTERCEPT"
Xq <- Xq[,names(cfi$edr),drop=FALSE]
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
     ylab="Probability of singing")
#polygon(100*c(r, rev(r)), c(q[,2], rev(q[,3])),
#        col="grey", border="grey")
matlines(r*100, qq, col="grey", lwd=1, lty=1)
lines(r*100, qq[,spp], col=1, lwd=2)
abline(v=cfall[spp,2]*100, lty=2)
rug(cfall[,2]*100, side=1, col="grey")

image(jd*365, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab="Julian days", ylab="Hours since sunrise",
    main=paste("Best model:", mi$sra))
image(lc, tr*100, qmat,
      col = rev(grey(seq(0, qmax, len=12))), axes=FALSE,
      xlab="Land cover types", ylab="Percent tree cover",
      main=paste("Best model:", mi$edr))
axis(1, 1:5, levels(xq$WNALC))
axis(2)
box()

par(op)
}
dev.off()

