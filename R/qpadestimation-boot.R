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
#load(file.path(ROOT, "out", "new_offset_data_package_2016-12-01.Rdata"))
load(file.path(ROOT, "out", "new_offset_data_package_2017-03-01.Rdata")) # this has ABMI interval fix
#load(file.path(ROOT, "out", "new_offset_data_package_2017-02-13-all-species.Rdata")) # all spp

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

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
#xtDur[["NONE"]] <- NULL
xtDur <- xtDur[intersect(names(xtDur), spp_singing)]

fitDurFun0 <- function(spp, fit=TRUE, type=c("rem","mix"), B=0) {
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
    iseq <- seq_len(n)
    is1 <- rowSums(Y, na.rm=TRUE) > 0
    if (fit) {
        res <- list()
        for (i in seq_len(B+1)) {
            ii <- if (i == 1)
                iseq[is1] else sample(iseq[is1], sum(is1), replace=TRUE)
            Yb <- Y[ii,,drop=FALSE]
            Db <- D[ii,,drop=FALSE]
            Xb <- X[ii,,drop=FALSE]
            mod <- try(cmulti(Yb | Db ~ 1, Xb, type=type))
            if (!inherits(mod, "try-error")) {
                rval <- mod[c("coefficients","vcov","nobs","loglik")]
                rval$p <- length(coef(mod))
                rval$names <- "(Intercept)"
            } else {
                rval <- mod
            }
            res[[i]] <- rval
        }
    } else {
        res <- list(Y=Y, D=D, n=n, pkey=rownames(X))
    }
    res
}

B <- 199

library(QPAD)
load_BAM_QPAD(3)
SPP <- getBAMspecieslist()

resDur <- vector("list", length(SPP))
names(resDur) <- SPP
for (i in SPP) {
    cat("Singing rate estimation for", i, date(), "\n")
    flush.console()
    resDur[[i]] <- try(fitDurFun0(i, TRUE, type="rem", B=B))
}
resDur_rem <- resDur

save(resDur_rem,
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2016-ABMIfix-boot0rem.Rdata"))

resDur <- vector("list", length(SPP))
names(resDur) <- SPP
for (i in SPP) {
    cat("Singing rate estimation for", i, date(), "\n")
    flush.console()
    resDur[[i]] <- try(fitDurFun0(i, TRUE, type="mix", B=B))
}
resDur_mix <- resDur

save(resDur_mix,
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2016-ABMIfix-boot0mix.Rdata"))
