## =============================================================================
## load stuff from processing -----------------------------------------

### Preliminaries

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2016"

## Load required packages
library(MASS)
library(mefa4)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-02.Rdata"))

## =============================================================================
## BAM-wise estimation -------------------------------------------

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","YEAR","PCODE")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])

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

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL

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
    resDur[[i]] <- try(fitDurFun(SPP[i], TRUE, type="rem"))
    #resDur[[i]] <- try(fitDurFun(SPP[i], TRUE, type="mix"))
}
names(resDur) <- SPP
resDurOK <- resDur[!sapply(resDur, inherits, "try-error")]
c(OK=length(resDurOK), failed=length(resDur)-length(resDurOK), all=length(resDur))
resDur <- resDurOK

save(resDur, resDurData,
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2016.Rdata"))
    #file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2016_mix.Rdata"))

## =============================================================================
## summarize results a la QPAD (no distance sampling this time) --------------

### Putting things together

type <- "rem"

## n.min is threshold above which all models are considered
## n.con is threshold above which the 0 constant model is considered
n.con <- 2
n.min <- 25 # (max df = 5, 5*5=25)

cat(type, "\n")
if (type == "rem")
    load(file.path(ROOT, "out", "estimates_SRA_QPAD_v2016.Rdata"))
if (type == "mix")
    load(file.path(ROOT, "out", "estimates_SRA_QPAD_v2016_mix.Rdata"))

## 0/1 table for successful model fit
sra_mod <- t(sapply(resDur, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))
tmp <- rownames(sra_mod)
sra_models <- matrix(0L, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod
## data deficient cases: dropped factor levels
## need to check length of NAMES and length of COEF --> correct factor level
## designation is possible, if != bump up AIC/BIC
for (spp in tmp) {
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
    }
}
## no constant model means exclude that species
sra_models[sra_models[,1]==0,] <- 0L
colSums(sra_models)

sra_nmod <- ncol(sra_mod)

## sample sizes
sra_n <- numeric(length(tmp))
names(sra_n) <- tmp
sra_nn <- sapply(resDur, function(z) ifelse(inherits(z[["0"]], "try-error"),
    NA, z[["0"]]$nobs))
sra_n[names(sra_nn)] <- sra_nn

## exclude all models for species with < n.con observations
sra_models[sra_n < n.con, ] <- 0L
## exclude non-constant model for species with n.con < observations < n.min
sra_models[sra_n < n.min & sra_n >= n.con, 2:ncol(sra_models)] <- 0L
table(rowSums(sra_models))

## spp to keep
spp <- tmp[rowSums(sra_models) > 0]

if (FALSE) {
spplist <- read.csv("c:/Users/Peter/Dropbox/bam/duration_ms/revisionMarch2016/duration-ms-species_2Mar2016_SM.csv")
rownames(spplist) <- spplist$Species_ID
spplist[spplist$Same_taxa != "",]
spplist <- spplist[spplist$Same_taxa == "",]
compare_sets(spp, rownames(spplist))

e <- new.env()
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-02.Rdata"), envir=e)
TAX <- e$TAX
#TAX <- TAX[TAX$Species_ID %in% spp,]
#rownames(TAX) <- TAX$Species_ID
TAX <- nonDuplicated(TAX, Species_ID, TRUE)
TAX <- TAX[spp,c("Species_ID", "English_Name", "Family_Sci")]
rm(e)

aa <- with(spplist, table(Family_Sci, Singing_birds == "TRUE"))
aa[aa[,1]>0,]
aa[aa[,2]>0,]
TAX2 <- TAX[!(rownames(TAX) %in% rownames(spplist)), ]
TAX2$Singing_birds <- TRUE
TAX2$Singing_birds[!(TAX2$Family_Sci %in% rownames(aa[aa[,2]>0,]))] <- NA
TAX2$Singing_birds[TAX2$Family_Sci %in% rownames(aa[aa[,1]>0,])] <- FALSE

tax <- rbind(spplist[,colnames(TAX2)], TAX2)
write.csv(tax, file="c:/Users/Peter/Dropbox/bam/duration_ms/revisionMarch2016/duration-ms-species_3Mar2016.csv")
}
tax <- read.csv("c:/Users/Peter/Dropbox/bam/duration_ms/revisionMarch2016/duration-ms-species_3Mar2016.csv")
rownames(tax) <- tax$Species_ID
compare_sets(spp, rownames(tax))
spp <- sort(intersect(spp, rownames(tax[tax$Singing_birds,])))

sra_models <- sra_models[spp,]
sra_n <- sra_n[spp]

## no. of parameters

sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

## estimates
sra_estimates <- resDur[spp]

## species table
e <- new.env()
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-02.Rdata"), envir=e)
tax2 <- e$TAX
#rownames(tax2) <- tax2$Species_ID
tax2 <- nonDuplicated(tax2, Species_ID, TRUE)
spp_table <- data.frame(spp=spp,
    scientific_name=tax2[spp, "Scientific_Name"],
    common_name=tax2[spp, "English_Name"])
rownames(spp_table) <- spp
spp_table <- droplevels(spp_table)

## get variable names for different models
sra_list <- sapply(sra_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

## loglik values
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}

## AIC/BIC
sra_aic <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_bic[] <- Inf
for (i in spp) {
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
}

## model ranking
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA
sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])

version <- "3"

bamcoefs <- list(spp=spp,
    spp_table=spp_table,
    sra_list=sra_list,
    sra_models=sra_models,
    sra_n=sra_n,
    sra_df=sra_df,
    sra_loglik=sra_loglik,
    sra_aic=sra_aic,
    sra_bic=sra_bic,
    sra_aicrank=sra_aicrank,
    sra_bicrank=sra_bicrank,
    sra_aicbest=sra_aicbest,
    sra_bicbest=sra_bicbest,
    sra_estimates=sra_estimates,
    version=version,
    type=type)
.BAMCOEFS <- list2env(bamcoefs)

if (type == "rem") {
    save(.BAMCOEFS, file=file.path(ROOT, "out", "BAMCOEFS_duration_rem.rda"))
}
if (type == "mix") {
    save(.BAMCOEFS, file=file.path(ROOT, "out", "BAMCOEFS_duration_mix.rda"))
}


## =============================================================================
## project-wise estimation -------------------------------------------

## models to consider
NAMES <- list("0"="INTERCEPT")
ff <- list(~ 1)

## define pcode, and excl (excl=T to exclude pcode, excl=F to use only pcode)
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

Pcodes <- levels(pkDur$PCODE)

#resDurBAMless1 <- list()
resDurPcode1 <- list()
for (i in 1:length(SPP)) {
#    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(Pcodes)) {
        cat("Singing rate (rem) estimation for", SPP[i], "in", Pcodes[j], "\n")
        flush.console()
#        resDurBAMless1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem",
#            pcode=Pcodes[j], excl=TRUE))
        resDurPcode1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem",
            pcode=Pcodes[j], excl=FALSE))
    }
}
save(resDurPcode1, #resDurBAMless1, 
    file="~/Dropbox/bam/duration_ms/revisionMarch2016/xval-Pcode1-rem.Rdata")

#resDurBAMless1_mix <- list()
resDurPcode1_mix <- list()
for (i in 1:length(SPP)) {
#    resDurBAMless1_mix[[SPP[i]]] <- list()
    resDurPcode1_mix[[SPP[i]]] <- list()
    for (j in 1:length(Pcodes)) {
        cat("Singing rate (mix) estimation for", SPP[i], "in", Pcodes[j], "\n")
        flush.console()
#        resDurBAMless1_mix[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="mix",
#            pcode=Pcodes[j], excl=TRUE))
        resDurPcode1_mix[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="mix",
            pcode=Pcodes[j], excl=FALSE))
    }
}
save(resDurPcode1_mix, #resDurBAMless1_mix, 
    file="~/Dropbox/bam/duration_ms/revisionMarch2016/xval-Pcode1-mix.Rdata")


## =============================================================================
## 3-5-10 validation variance / bias -----------------------------------------
## bias/variance/MSE for m0/m0t and mb/mbt

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_rem.rda"), envir=e)
.BAMCOEFSrem <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS

aic0 <- .BAMCOEFSrem$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))

#cfall0 <- sapply(SPP, function(spp) exp(coefBAMspecies(spp, 0, 0)$sra))
cfall0 <- sapply(SPP, function(spp) {
    exp(unname(.BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients))
    })
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

xvfun <- function(spp, B=1000, ymin=1) {
    Y <- groupSums(xtDur[[spp]][rn,], 2,
        c("0-3", "xxx", "0-3", "0-3", "xxx", "xxx", "0-3", "0-3",
        "xxx", "3-5", "3-5", "xxx", "xxx", "3-5", "5-10", "5-10",
        "5-10", "5-10", "xxx", "5-10", "5-10", "5-10", "5-10", "5-10"))
    Y <- as.matrix(Y)[,c("0-3","3-5","5-10")]
    YY <- cbind("0-3"=Y[,1], "0-5"=Y[,1]+Y[,2], "0-10"=rowSums(Y))
    YY <- YY[YY[,3] >= ymin,,drop=FALSE]
    pkDurS <- pkDur[rownames(YY),]

    cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov


    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov

    best0 <- colnames(aic0)[which.min(aic0[spp,])]
    bestb <- colnames(aicb)[which.min(aicb[spp,])]
    X0 <- model.matrix(ff[[best0]], pkDurS)
    Xb <- model.matrix(ff[[bestb]], pkDurS)

    cf0 <- .BAMCOEFSrem$sra_estimates[[spp]][[best0]]$coefficients
    vc0 <- .BAMCOEFSrem$sra_estimates[[spp]][[best0]]$vcov
    cfb <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$coefficients
    vcb <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$vcov

    vci0b <- Matrix::nearPD(vci0b)$mat
    vc0 <- Matrix::nearPD(vc0)$mat
    vcb <- Matrix::nearPD(vcb)$mat

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
        out$m0t <- c(YC3=mean(YC3i_0), YC5=mean(YC5i_0))
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

nob <- pbsapply(SPP, xtfun, ymin=1)
res <- list()
for (spp in SPP) {
cat(spp, date(), "\n");flush.console()
res[[spp]] <- try(xvfun(spp))
}

save(res, nob, file=file.path(ROOT2, "var-bias-res.Rdata"))

load(file.path(ROOT2, "var-bias-res.Rdata"))
res2 <- res[!sapply(res, inherits, "try-error")]

aaa <- data.frame(Var=as.numeric(sapply(res2, "[[", "Var")),
    MSE=as.numeric(sapply(res2, "[[", "MSE")),
    Bias=as.numeric(sapply(res2, "[[", "Bias")),
    Model=rep(c("0","b","0t","bt"), each=2),
    Duration=c("3","5"),
    Species=as.factor(rep(names(res2), each=8)))
aaa$n <- nob[match(aaa$Species, names(nob))]
aaa$logn <- log(aaa$n)
dim(aaa)
aaa <- aaa[!is.na(aaa$Var) & aaa$n >= 2,]
summary(aaa)
dim(aaa)
aaa$Mixture <- ifelse(aaa$Model %in% c("b","bt"), 1, 0)
aaa$Timevar <- ifelse(aaa$Model %in% c("0t","bt"), 1, 0)

library(lme4)
m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1m <- lmer(Var ~ (Mixture + Timevar + Duration + logn)^2 + (1 | Species), aaa)
m2m <- lmer(Bias ~ (Mixture + Timevar + Duration + logn)^2 + (1 | Species), aaa)
cbind(fixef(m1m), coef(m1))
cbind(fixef(m2m), coef(m2))
## minor diffs: use lm

#m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
#m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1 <- step(m1)
m2 <- step(m2)

#m1 <- lm(Var ~ Model + Duration + logn, aaa)
#m2 <- lm(Bias ~ Model + Duration + logn, aaa)
summary(m1)
summary(m2)

a1 <- anova(update(m1, . ~ . + Species))
a1$Perc <- 100 * a1[,"Sum Sq"] / sum(a1[,"Sum Sq"])
a2 <- anova(update(m2, . ~ . + Species))
a2$Perc <- 100 * a2[,"Sum Sq"] / sum(a2[,"Sum Sq"])
a1
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
    max3_0t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
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
    max3_0t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
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

## =============================================================================
## species specific predictions plots



## =============================================================================
## when does model fail?

load("~/Dropbox/bam/duration_ms/revisionMarch2016/xval-Pcode1-rem.Rdata")
load("~/Dropbox/bam/duration_ms/revisionMarch2016/xval-Pcode1-mix.Rdata")

SPP <- names(resDurPcode1)
problem <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        ii <- intersect(rownames(xtDur[[spp]]), ii)
        yy <- rowSums(xtDur[[spp]][ii, ,drop=FALSE])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                ymean=mean(yy),
                logphi=NA, se_logphi=NA,
                    nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    ymean=mean(yy),
                    logphi=NA, se_logphi=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    ymean=mean(yy),
                    logphi=unname(tmp3$coefficients), 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problem0 <- do.call(rbind, problem)

problem <- list()
for (spp in names(resDurPcode1_mix)) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1_mix[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        ii <- intersect(rownames(xtDur[[spp]]), ii)
        yy <- rowSums(xtDur[[spp]][ii, ,drop=FALSE])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                ymean=mean(yy),
                logphi=NA, 
                logitc=NA, 
                se_logphi=NA,
                se_logitc=NA,
                cor=NA,
                nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    ymean=mean(yy),
                    logphi=NA, 
                    logitc=NA, 
                    se_logphi=NA,
                    se_logitc=NA,
                    cor=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    ymean=mean(yy),
                    logphi=unname(tmp3$coefficients)[1], 
                    logitc=unname(tmp3$coefficients)[2], 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    se_logitc=sqrt(tmp3$vcov[2,2]),
                    cor=cov2cor(tmp3$vcov)[1,2],
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problemb <- do.call(rbind, problem)

save(problem0, problemb, file=file.path(ROOT2, "problem.Rdata"))

ss <- c("0 observation with multiple duration (1)", "1 observation with multiple duration (2)")

dat0 <- problem0[!(problem0$msg %in% ss),]
dat0$okfit <- ifelse(is.na(dat0$logphi), 0, 1)
dat0$okse <- ifelse(is.na(dat0$se_logphi), 0, 1)
dat0$msg <- NULL
mod00 <- glm(okfit ~ log(ndet+1), dat0, family=binomial)
mod0se0 <- glm(okse ~ log(ndet+1), dat0, family=binomial)

datb <- problemb[!(problemb$msg %in% ss),]
datb$okfit <- ifelse(is.na(datb$logphi), 0, 1)
datb$okse <- ifelse(is.na(datb$se_logphi), 0, 1)
datb$msg <- NULL
modb0 <- glm(okfit ~ log(ndet+1), datb, family=binomial)
modbse0 <- glm(okse ~ log(ndet+1), datb, family=binomial)

mod01 <- glm(okfit ~ log(n1+1), dat0, family=binomial)
mod0se1 <- glm(okse ~ log(n1+1), dat0, family=binomial)
mod02 <- glm(okfit ~ log(n2+1), dat0, family=binomial)
mod0se2 <- glm(okse ~ log(n2+1), dat0, family=binomial)

modb1 <- glm(okfit ~ log(n1+1), datb, family=binomial)
modbse1 <- glm(okse ~ log(n1+1), datb, family=binomial)
modb2 <- glm(okfit ~ log(n2+1), datb, family=binomial)
modbse2 <- glm(okse ~ log(n2+1), datb, family=binomial)


op <- par(mfrow=c(1,3))
for (i in 1:3) {
if (i == 1) {
    mod0 <- mod00
    modb <- modb0
    mod0se <- mod0se0
    modbse <- modbse0
    xlab <- "Number of >0 survey counts"
}
if (i == 2) {
    mod0 <- mod01
    modb <- modb1
    mod0se <- mod0se1
    modbse <- modbse1
    xlab <- "Number of >1 survey counts"
}
if (i == 3) {
    mod0 <- mod02
    modb <- modb2
    mod0se <- mod0se2
    modbse <- modbse2
    xlab <- "Number of >2 survey counts"
}
nmax <- 100
ndat <- data.frame(n=2:nmax, 
    p0=plogis(coef(mod0)[1] + coef(mod0)[2]*log(1 + 2:nmax)),
    pb=plogis(coef(modb)[1] + coef(modb)[2]*log(1 + 2:nmax)),
    p0se=plogis(coef(mod0se)[1] + coef(mod0se)[2]*log(1 + 2:nmax)),
    pbse=plogis(coef(modbse)[1] + coef(modbse)[2]*log(1 + 2:nmax)))

plot(ndat[,1], ndat[,2], col=2, type="l", ylim=c(0,1), lwd=2,
    xlab=xlab, ylab="Probability",
    xlim=c(0,nmax))
lines(ndat[,1], ndat[,3], col=4, lwd=2)
lines(ndat[,1], ndat[,4], col=2, lwd=2, lty=2)
lines(ndat[,1], ndat[,5], col=4, lwd=2, lty=2)
abline(h=0.9, lty=1)
abline(v=ndat[which.min(abs(ndat[,2]-0.9)),1], col=2)
abline(v=ndat[which.min(abs(ndat[,3]-0.9)),1], col=4)
legend("bottomright", col=c(2,2,4,4), lty=c(1,2,1,2), lwd=2, 
    legend=c("m0 fit", "m0 SE", "mb fit", "mb SE"), bty="n")
}
par(op)

dat0se <- dat0[!is.na(dat0$se_logphi),]
datbse <- datb[!is.na(datb$se_logphi),]

points(se_logphi ~ nobs, dat0se, ylim=c(0,50), xlim=c(0, 100),
    pch=19, col=rgb(50,50,50,50,maxColorValue=255))

plot(se_logphi ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(se_logitc ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(abs(cor) ~ nobs, datbse, xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
hist(datbse$cor)





## =============================================================================
## figures etc
library(rworldmap)
rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
X0 <- pkDur[rn,]
D <- ltdur$end[match(X0$DURMETH, rownames(ltdur$end)),]
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(pkDur[, c("X","Y")], pch=".",
    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255))
points(X0[rowSums(!is.na(D)) > 1, c("X","Y")], pch=19,
    col="red", cex=0.3)


## =============================================================================
## Model selection results

## sra 0 vs b models

aic0 <- .BAMCOEFSrem$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
colnames(aic0) <- paste0("m0_", colnames(aic0))
colnames(aicb) <- paste0("mb_", colnames(aicb))
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))

aic0 <- aic0[SPP,]
aicb <- aicb[SPP,]
aic <- cbind(aic0[SPP,], aicb[SPP,])

np <- .BAMCOEFSrem$sra_n[SPP]

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
names(best0) <- names(bestb) <- names(best) <- SPP

par(mfrow=c(1,3))
plot(np, waic0[,1], log="x", ylim=c(0,1), pch=ifelse(best0=="0", "o", "+"))
plot(np, waicb[,1], log="x", ylim=c(0,1), pch=ifelse(bestb=="0", "o", "+"))
plot(np, waic[,"m0_0"] + waic[,"mb_0"], log="x", ylim=c(0,1),
    pch=ifelse(best %in% c("m0_0", "mb_0"), "o", "+"))

MAX <- 5000
nn <- 25:MAX
ww <- sapply(nn, function(z) mean((rowSums(waic[,grepl("_0", colnames(waic))]))[np >= z]))
ww2 <- sapply(nn, function(z) mean((rowSums(waic[,grepl("mb_", colnames(waic))]))[np >= z]))

par(mfrow=c(1,2))
plot(nn, 100*(1-ww), type="l", ylim=100*c(0.75, 1), xlab="Number of detections",
    ylab="% time varying", xlim=c(0,MAX))
rug(np)
plot(nn, 100*ww2, type="l", ylim=100*c(0.75, 1), xlab="Number of detections",
    ylab="% mixture", xlim=c(0,MAX))
rug(np)

table(Timevar=!grepl("_0", best), Mixture=grepl("mb_", best))

hasJD <- paste0(rep(c("m0_", "mb_"), each=6), rep(c(1, 3, 5:8), 2))
hasSR <- paste0(rep(c("m0_", "mb_"), each=10), rep(c(2, 4:8, 11:14), 2))
hasLS <- paste0(rep(c("m0_", "mb_"), each=6), rep(c(9:14), 2))

SPPx <- names(np)[np >= 25]
waicx <- waic[SPPx,]
bestx <- best[SPPx]
best0x <- best0[SPPx]
bestbx <- bestb[SPPx]
best0bx <- sapply(strsplit(bestx, "_"), "[[", 2)
Support <- t(sapply(as.character(0:14), function(z) {
    c(m0=sum(best0x == z), mb=sum(bestbx == z), Combined=sum(best0bx == z))
}))
Support

## m0 and mb combined here, only spp >=25 det
## constant
sum(grepl("_0", bestx)) * 100 / length(SPPx)
## JDAY
sum(bestx %in% hasJD) * 100 / length(SPPx)
## TSLS
sum(bestx %in% hasLS) * 100 / length(SPPx)
## TSSR
sum(bestx %in% hasSR) * 100 / length(SPPx)
## JDAY & TSSR
sum(bestx %in% intersect(hasJD, hasSR)) * 100 / length(SPPx)
## TSLS & TSSR
sum(bestx %in% intersect(hasLS, hasSR)) * 100 / length(SPPx)
## JDAY or TSLS
sum(bestx %in% c(hasJD, hasLS)) * 100 / length(SPPx)
## JDAY or TSLS & TSSR
sum(bestx %in% intersect(c(hasJD, hasLS), hasSR)) * 100 / length(SPPx)












## -- old

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

## read all the suff, read in lines 1-86 (RND defs there!)

## random PCODE IDs or not (projects as they are)
isRND <- FALSE

## rem -- subsets are for sets of species to be combined
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

## mix -- subsets are for sets of species to be combined
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
waic_fun <- function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}



## when does model fail?

problem <- list()
for (spp in names(resDurPcode1)) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        yy <- rowSums(xtDur[[spp]][ii, ])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                logphi=NA, se_logphi=NA,
                    nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=NA, se_logphi=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=unname(tmp3$coefficients), 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problem0 <- do.call(rbind, problem)

problem <- list()
for (spp in names(resDurPcode1_mix)) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1_mix[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        yy <- rowSums(xtDur[[spp]][ii, ])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                logphi=NA, 
                logitc=NA, 
                se_logphi=NA,
                se_logitc=NA,
                cor=NA,
                nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=NA, 
                    logitc=NA, 
                    se_logphi=NA,
                    se_logitc=NA,
                    cor=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=unname(tmp3$coefficients)[1], 
                    logitc=unname(tmp3$coefficients)[2], 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    se_logitc=sqrt(tmp3$vcov[2,2]),
                    cor=cov2cor(tmp3$vcov)[1,2],
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problemb <- do.call(rbind, problem)

save(problem0, problemb, file=file.path(ROOT2, "problem.Rdata"))

ss <- c("0 observation with multiple duration (1)", "1 observation with multiple duration (2)")

dat0 <- problem0[!(problem0$msg %in% ss),]
dat0$okfit <- ifelse(is.na(dat0$logphi), 0, 1)
dat0$okse <- ifelse(is.na(dat0$se_logphi), 0, 1)
dat0$msg <- NULL
mod00 <- glm(okfit ~ log(ndet+1), dat0, family=binomial)
mod0se0 <- glm(okse ~ log(ndet+1), dat0, family=binomial)

datb <- problemb[!(problemb$msg %in% ss),]
datb$okfit <- ifelse(is.na(datb$logphi), 0, 1)
datb$okse <- ifelse(is.na(datb$se_logphi), 0, 1)
datb$msg <- NULL
modb0 <- glm(okfit ~ log(ndet+1), datb, family=binomial)
modbse0 <- glm(okse ~ log(ndet+1), datb, family=binomial)

mod01 <- glm(okfit ~ log(n1+1), dat0, family=binomial)
mod0se1 <- glm(okse ~ log(n1+1), dat0, family=binomial)
mod02 <- glm(okfit ~ log(n2+1), dat0, family=binomial)
mod0se2 <- glm(okse ~ log(n2+1), dat0, family=binomial)

modb1 <- glm(okfit ~ log(n1+1), datb, family=binomial)
modbse1 <- glm(okse ~ log(n1+1), datb, family=binomial)
modb2 <- glm(okfit ~ log(n2+1), datb, family=binomial)
modbse2 <- glm(okse ~ log(n2+1), datb, family=binomial)


op <- par(mfrow=c(1,3))
for (i in 1:3) {
if (i == 1) {
    mod0 <- mod00
    modb <- modb0
    mod0se <- mod0se0
    modbse <- modbse0
    xlab <- "Number of >0 survey counts"
}
if (i == 2) {
    mod0 <- mod01
    modb <- modb1
    mod0se <- mod0se1
    modbse <- modbse1
    xlab <- "Number of >1 survey counts"
}
if (i == 3) {
    mod0 <- mod02
    modb <- modb2
    mod0se <- mod0se2
    modbse <- modbse2
    xlab <- "Number of >2 survey counts"
}
nmax <- 100
ndat <- data.frame(n=2:nmax, 
    p0=plogis(coef(mod0)[1] + coef(mod0)[2]*log(1 + 2:nmax)),
    pb=plogis(coef(modb)[1] + coef(modb)[2]*log(1 + 2:nmax)),
    p0se=plogis(coef(mod0se)[1] + coef(mod0se)[2]*log(1 + 2:nmax)),
    pbse=plogis(coef(modbse)[1] + coef(modbse)[2]*log(1 + 2:nmax)))

plot(ndat[,1], ndat[,2], col=2, type="l", ylim=c(0,1), lwd=2,
    xlab=xlab, ylab="Probability",
    xlim=c(0,nmax))
lines(ndat[,1], ndat[,3], col=4, lwd=2)
lines(ndat[,1], ndat[,4], col=2, lwd=2, lty=2)
lines(ndat[,1], ndat[,5], col=4, lwd=2, lty=2)
abline(h=0.9, lty=1)
abline(v=ndat[which.min(abs(ndat[,2]-0.9)),1], col=2)
abline(v=ndat[which.min(abs(ndat[,3]-0.9)),1], col=4)
legend("bottomright", col=c(2,2,4,4), lty=c(1,2,1,2), lwd=2, 
    legend=c("m0 fit", "m0 SE", "mb fit", "mb SE"), bty="n")
}
par(op)

dat0se <- dat0[!is.na(dat0$se_logphi),]
datbse <- datb[!is.na(datb$se_logphi),]

points(se_logphi ~ nobs, dat0se, ylim=c(0,50), xlim=c(0, 100),
    pch=19, col=rgb(50,50,50,50,maxColorValue=255))

plot(se_logphi ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(se_logitc ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(abs(cor) ~ nobs, datbse, xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
hist(datbse$cor)





## sra m0 vs mb models


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


R <- 1000

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

rownames(TAX) <- TAX$Species_ID
tab <- droplevels(TAX[SPP,c("Species_ID", "English_Name", "Family_Sci")])
tab$n <- np
tab$OK <- OK
tab$phi_0 <- round(cfall0[,1], 4)
tab$phi_b <- round(cfallb[,1], 4)
tab$c <- round(cfallb[,2], 4)

write.csv(tab, row.names=FALSE, file=file.path(ROOT2, "duration-ms-species.csv"))

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


pdf(paste0(ROOT2, "/m0-vs-mb.pdf"), onefile=TRUE, width=10, height=12)
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
    " w0=", round(sum(wph),2)),
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





JDAY,
TSSR,
JDAY + I(JDAY^2)
TSSR + I(TSSR^2)
JDAY + TSSR
JDAY + I(JDAY^2) + TSSR
JDAY + TSSR + I(TSSR^2)
JDAY + I(JDAY^2) + TSSR + I(TSSR^2)
TSLS
TSLS + I(TSLS^2)
TSLS + TSSR
TSLS + I(TSLS^2) + TSSR
TSLS + TSSR + I(TSSR^2)
TSLS + I(TSLS^2) + TSSR + I(TSSR^2)
