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
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","PCODE")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
sort(table(pkDur$PCODE))
## only 1 obs in project (all others are >20)
pkDur <- droplevels(pkDur[pkDur$PCODE != "GLDRPCCL01",])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])

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
#resDurBAMless1_mix <- vector("list", length(SPP))
#resDurPcode1_mix <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    resDurBAMless1[[SPP[i]]] <- list()
    resDurPcode1[[SPP[i]]] <- list()
    for (j in 1:length(Pcodes)) {
        cat("Singing rate estimation for", SPP[i], "in", Pcodes[j], "\n")
        flush.console()
        resDurBAMless1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem", 
            pcode=Pcodes[j], excl=TRUE))
        resDurPcode1[[SPP[i]]][[Pcodes[j]]] <- try(fitDurFun2(SPP[i], TRUE, type="rem", 
            pcode=Pcodes[j], excl=FALSE))
    }
}
str(resDurBAMless1)
str(resDurPcode1)

