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
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-21.Rdata"))

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
SPP <- names(xtDur)

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
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2016_rem.Rdata"))
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
    load(file.path(ROOT2, "estimates_SRA_QPAD_v2016_rem.Rdata"))
if (type == "mix")
    load(file.path(ROOT2, "estimates_SRA_QPAD_v2016_mix.Rdata"))

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
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-21.Rdata"), envir=e)
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
df0 <- matrix(.BAMCOEFSrem$sra_df, nrow(aic0), 15, byrow=TRUE)
dfb <- matrix(.BAMCOEFSmix$sra_df, nrow(aicb), 15, byrow=TRUE)
n0 <- .BAMCOEFSrem$sra_n
nb <- .BAMCOEFSmix$sra_n
## this is AICc
aic0 <- aic0 + (2*df0*(df0+1)) / (n0-df0-1)
aicb <- aicb + (2*dfb*(dfb+1)) / (nb-dfb-1)

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
    #Bias <- sqrt(MSE - Var)
    Bias <- colSums(theta_hat - theta) / B
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
                ymean=mean(yy[yy>0]),
                ymedian=unname(median(yy[yy>0])),
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
                    ymean=mean(yy[yy>0]),
                    ymedian=unname(median(yy[yy>0])),
                    logphi=NA, se_logphi=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    ymean=mean(yy[yy>0]),
                    ymedian=unname(median(yy[yy>0])),
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
                ymean=mean(yy[yy>0]),
                ymedian=unname(median(yy[yy>0])),
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
                    ymean=mean(yy[yy>0]),
                    ymedian=unname(median(yy[yy>0])),
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
                    ymean=mean(yy[yy>0]),
                    ymedian=unname(median(yy[yy>0])),
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

