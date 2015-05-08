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

load(file.path(ROOT, "out", "new_offset_data_package_2015-05-07.Rdata"))

save(resDur, resDurData, 
    file=file.path(ROOT, "out", "estimates_SRA_QPAD_v2015.Rdata"))





## -- old

## might worth to fit intercept model for full data (and not only PKEYs where covariates are available) -- maybe for next ms or just for select spp
## compare sample sizes

### !!!!!!!!!!! Analysis for DApq ms and offset/density db

USE_ROAD <- FALSE

library(mefa4)
source("c:/Dropbox/bam/D.fit.R")

## load pre-processed stuff
setwd("c:/Dropbox/bam/DApq3")
load("lifehist_final.Rdata")
load("BAM_V4.Rdata")

## keep songbird records only
spp <- as.character(taxo2$SPECIES[taxo2$Order %in% "PASSERIFORMES"])
PCTBL <- PCTBL[PCTBL$BEH %in% c(1, 6),] # 1=heard, 6=seen and heard
levels(taxo2$Order) <- toupper(levels(taxo2$Order))
PCTBL <- PCTBL[PCTBL$SPECIES %in% 
  taxo2$SPECIES[taxo2$Order=="PASSERIFORMES"],]
PCTBL$SPECIES <- PCTBL$SPECIES[drop=TRUE]

## Restrict to Boreal
PKEY <- PKEY[PKEY$BOR_LOC != "OUT",] ## incl B-Alpine and Hemiboreal

## no roadside data
if (!USE_ROAD)
    PKEY <- PKEY[PKEY$OnRoad != "Y",]

with(PKEY, table(DURMETH, OnRoad))
with(PKEY, table(DISTMETH, OnRoad))

#### **********************************************************************
#### **** DISTANCE EFFECT ESTIMATION
#### **********************************************************************

## exclude unknown and single distance bands, and 10+ DURATION
pctblDist <- PCTBL[!(PCTBL$DISTMETH %in% c("D","F","G","K","O")) & # single dist band
    !is.na(PCTBL$DISTMETH) & # Dist code is NA
    !(PCTBL$DURATION %in% c(10,11,3,8,9)),] # outside of 10 min or unknown max duration

xtDist <- Xtab(ABUND ~ PKEY + DISTANCE + SPECIES, pctblDist,
    cdrop=c("4","5","6","9")) # drop 

distclasses <- nonDuplicated(PCTBL[,c("DISTANCE",
    "DistanceDescrip", "DIST_START", "DIST_END")], DISTANCE, TRUE)
distclasses[colnames(xtDist[[1]]),]

tmp <- nonDuplicated(pctblDist, PKEY, TRUE)
ii <- intersect(rownames(tmp), rownames(PKEY))
covarDist <- data.frame(TREE=PKEY[ii,"TREE"], 
    ROAD=ifelse(PKEY[ii,"OnRoad"]=="Y",1,0), 
    LCC=PKEY[ii,"lccs"], 
    PCODE=PKEY[ii,"PCODE"], 
    DISTMETH=PKEY[ii,"DISTMETH"])
rownames(covarDist) <- PKEY[ii,"PKEY"]
covarDist <- covarDist[rowSums(is.na(covarDist))==0,]
#covarDist <- covarDist[covarDist$DISTMETH != "U",]

aa <- distclasses[colnames(xtDist[[1]]),]
covarDist$ROAD <- as.factor(covarDist$ROAD) ## this will change to road type
table(covarDist$DISTMETH)

distintervals3 <- matrix(c(
     75,Inf, NA, NA, NA, NA, NA,    # A
     50,100,Inf, NA, NA, NA, NA,    # B
     50,Inf, NA, NA, NA, NA, NA,    # C
     50, 75,100,Inf, NA, NA, NA,    # H
     25, 50, 75,Inf, NA, NA, NA,    # I
     50,100,150, NA, NA, NA, NA,    # L
     25, 50, 75,100,125,150,Inf,    # M
    100,Inf, NA, NA, NA, NA, NA,    # N
     30, 50, 75, NA, NA, NA, NA,    # P
     50, 75,100, NA, NA, NA, NA,    # R
     25, 50, 75,100,Inf, NA, NA,    # S
     50,100, NA, NA, NA, NA, NA,    # T
     ## note: U is faked: intervals are pooled according to B
     50,100,Inf, NA, NA, NA, NA),    # U
     13, 7, byrow=TRUE)
rownames(distintervals3) <- c("A","B","C","H","I","L","M","N", "P", "R", "S", "T","U")
distintervals3 <- distintervals3 / 100

if (!USE_ROAD) {
    NAMES <- list("1"=c("INTERCEPT", "TREE"),
        "2"=c("INTERCEPT", "LCC2", "LCC3", "LCC4", "LCC5"))
    ff <- list(~ TREE,
        ~ LCC)
} else {
    NAMES <- list("1"=c("INTERCEPT", "TREE"),
        "2"=c("INTERCEPT", "LCC2", "LCC3", "LCC4", "LCC5"),
        "3"=c("INTERCEPT", "ROAD"),
        "4"=c("INTERCEPT", "TREE", "ROAD"),
        "5"=c("INTERCEPT", "LCC2", "LCC3", "LCC4", "LCC5", "ROAD"))
    ff <- list(~ TREE,
        ~ LCC,
        ~ ROAD,
        ~ TREE + ROAD,
        ~ LCC + ROAD)
}

#SPP <- "OVEN"
fitDistFun <- function(SPP, fit=TRUE) {
    ## get nonzero PCs
    iob <- rowSums(xtDist[[SPP]])>0
    if (sum(iob)==0)
        return(structure("0 observation with multiple distance (1)", 
            class="try-error"))
    if (sum(iob)==1)
        return(structure("1 observation with multiple distance (2)", 
            class="try-error"))
    xt <- xtDist[[SPP]][iob,]
    ## inner join, so only PKEYs with all available covariates are used
    xDist <- Mefa(xt, covarDist, join="inner")
    if (dim(xDist)[1]<=1)
        return(structure("<2 observation with multiple distance (3)", 
            class="try-error"))
    ## sampling protocol IDs
    ii <- samp(xDist)$DISTMETH
    ii <- as.character(ii)
    ## pool columns where end times are same
    Y <- xtab(xDist)
    ## some how dense matrix is quicker
    Y <- as.matrix(Y)
    ## end distances
    D <- distintervals3[as.character(samp(xDist)$DISTMETH),]
    ## pieces for each protocol type
    nans <- rowSums(is.na(distintervals3))
    if (any(ii == "A")) {
            WHICH <- ii == "A"
            YA <- cbind(matrix(Y[WHICH,c("25","14")], sum(WHICH)))
            YA <- cbind(YA, matrix(NA, nrow(YA), nans["A"]))
            iA <- which(WHICH)
        } else {
            YA <- iA <- NULL
        }
    if (any(ii == "B")) {
            WHICH <- ii == "B"
            YB <- cbind(matrix(Y[WHICH,c("1","2","3")], sum(WHICH)))
            YB <- cbind(YB, matrix(NA, nrow(YB), nans["B"]))
            iB <- which(WHICH)
        } else {
            YB <- iB <- NULL
        }
    if (any(ii == "C")) {
            WHICH <- ii == "C"
            YC <- cbind(matrix(Y[WHICH,c("1","15")], sum(WHICH)))
            YC <- cbind(YC, matrix(NA, nrow(YC), nans["C"]))
            iC <- which(WHICH)
        } else {
            YC <- iC <- NULL
        }
    if (any(ii == "H")) {
            WHICH <- ii == "H"
            YH <- cbind(matrix(Y[WHICH,c("1","12","13","3")], sum(WHICH)))
            iH <- which(WHICH)
            YH <- cbind(YH, matrix(NA, nrow(YH), nans["H"]))
        } else {
            YH <- iH <- NULL
        }
    if (any(ii == "I")) {
            WHICH <- ii == "I"
            YI <- cbind(matrix(Y[WHICH,c("10","11","12","14")], sum(WHICH)))
            iI <- which(WHICH)
            YI <- cbind(YI, matrix(NA, nrow(YI), nans["I"]))
        } else {
            YI <- iI <- NULL
        }
    if (any(ii == "L")) {
            WHICH <- ii == "L"
            YL <- cbind(matrix(Y[WHICH,c("1","2","21")], sum(WHICH)))
            iL <- which(WHICH)
            YL <- cbind(YL, matrix(NA, nrow(YL), nans["L"]))
        } else {
            YL <- iL <- NULL
        }
    if (any(ii == "M")) {
            WHICH <- ii == "M"
            YM <- cbind(matrix(Y[WHICH,c("10","11","12","13","16","17","20")], sum(WHICH)))
            iM <- which(WHICH)
            YM <- cbind(YM, matrix(NA, nrow(YM), nans["M"]))
        } else {
            YM <- iM <- NULL
        }
    if (any(ii == "N")) {
            WHICH <- ii == "N"
            YN <- cbind(matrix(Y[WHICH,c("8","3")], sum(WHICH)))
            iN <- which(WHICH)
            YN <- cbind(YN, matrix(NA, nrow(YN), nans["N"]))
        } else {
            YN <- iN <- NULL
        }
    if (any(ii == "P")) {
            WHICH <- ii == "P"
            YP <- cbind(matrix(Y[WHICH,c("23","24","12")], sum(WHICH)))
            iP <- which(WHICH)
            YP <- cbind(YP, matrix(NA, nrow(YP), nans["P"]))
        } else {
            YP <- iP <- NULL
        }
    if (any(ii == "R")) {
            WHICH <- ii == "R"
            YR <- cbind(matrix(Y[WHICH,c("1","12","13")], sum(WHICH)))
            iR <- which(WHICH)
            YR <- cbind(YR, matrix(NA, nrow(YR), nans["R"]))
        } else {
            YR <- iR <- NULL
        }
    if (any(ii == "S")) {
            WHICH <- ii == "S"
            YS <- cbind(matrix(Y[WHICH,c("10","11","12","13","3")], sum(WHICH)))
            iS <- which(WHICH)
            YS <- cbind(YS, matrix(NA, nrow(YS), nans["S"]))
        } else {
            YS <- iS <- NULL
        }
    if (any(ii == "T")) {
            WHICH <- ii == "T"
            YT <- cbind(matrix(Y[WHICH,c("1","2")], sum(WHICH)))
            iT <- which(WHICH)
            YT <- cbind(YT, matrix(NA, nrow(YT), nans["T"]))
        } else {
            YT <- iT <- NULL
        }
    if (any(ii == "U")) {
            WHICH <- ii == "U"
            tmp1 <- rowSums(matrix(Y[WHICH,c("27","28","29","30","31")], sum(WHICH)))
            tmp2 <- rowSums(matrix(Y[WHICH,c("32","33","34","35","36")], sum(WHICH)))
            tmp3 <- rowSums(matrix(Y[WHICH,c("16","17","20")], sum(WHICH)))
            YU <- cbind(tmp1,tmp2,tmp3)
            YU <- cbind(YU, matrix(NA, nrow(YU), nans["U"]))
            iU <- which(WHICH)
        } else {
            YU <- iU <- NULL
        }
    ## putting pieces together
    Y <- rbind(YA,YB,YC,YH,YI,YL,YM,YN,YP,YR,YS,YT,YU)
    if (is.null(Y))
        return(structure("0 observation with multiple distance (4)", 
            class="try-error"))
    iall <- c(iA,iB,iC,iH,iI,iL,iM,iN,iP,iR,iS,iT,iU)
    D <- D[iall,]
    ## integer mode -- faster, but DO NOT use for intervals (<1)
    storage.mode(Y) <- "integer"
    YY <- unname(Y)
    DD <- unname(D)
    X <- samp(xDist)[iall,]
    nn <- table(ii, samp(xDist)$ROAD)
    if (fit) {
        ## allocate return value
        res <- if (!is.null(ff))
            vector("list", 1+length(ff)) else vector("list", 1)
        ## intercept model
        res[[1]] <- try(D.fit(YY, DD, NULL, type="dist"))
        if (!inherits(res[[1]], "try-error")) {
            res[[1]]$p <- 1
            res[[1]]$names <- "INTERCEPT"
        }
        ## covariates models
        if (!is.null(ff)) {
            for (i in 1:length(ff)) {
                XX <- model.matrix(ff[[i]], X)
                XX <- unname(XX)
                res[[(i+1)]] <- try(D.fit(YY, DD, XX, type="dist"))
                if (!inherits(res[[(i+1)]], "try-error")) {
                    res[[(i+1)]]$p <- length(res[[(i+1)]]$coef)
                    res[[(i+1)]]$names <- NAMES[[i]]
                }
            }
        }
        names(res) <- if (!is.null(ff))
            0:length(ff) else 0
        ## number of observations
        res$n <- nn # attr(res, "nobs") <- nn
    } else {
        res <- list(Y=YY, D=DD, i=iall, n=nn)
    }
    res
}

#str(fitDistFun("OVEN", FALSE))
#str(fitDistFun("OVEN", TRUE))

## see which species has enough data

SPP <- names(xtDist)
resDist <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("EDR estimation for", SPP[i], "\n")
    flush.console()
    resDist[[i]] <- try(fitDistFun(SPP[i], FALSE))
}
names(resDist) <- SPP
resDistOK <- resDist[!sapply(resDist, inherits, "try-error")]
c(OK=length(resDistOK), failed=length(resDist)-length(resDistOK), all=length(resDist))
t(sapply(resDistOK, function(z) colSums(z$n)))
resDistData <- resDistOK

## estimate species with data
SPP <- names(resDistOK)
resDist <- vector("list", length(SPP))
date()
for (i in 1:length(SPP)) {
    cat("EDR estimation for", SPP[i], "\n")
    flush.console()
    resDist[[i]] <- try(fitDistFun(SPP[i], TRUE))
}
date()
names(resDist) <- SPP
resDistOK <- resDist[!sapply(resDist, inherits, "try-error")]
c(OK=length(resDistOK), failed=length(resDist)-length(resDistOK), all=length(resDist))

if (!USE_ROAD)
    save(resDist, resDistData, file="estimates_EDR_QPADpaper.Rdata")

if (USE_ROAD)
    save(resDist, resDistData, file="estimates_EDR_QPADpaper_withRoad.Rdata")

#### **********************************************************************
#### **** DURATION
#### **********************************************************************

load("ABMI_V2011.Rdata")

## ABMI: keep only aerial detections and Passerines
keep <- PCTBL_abmi$BEHAVIOUR %in% c("Counter-singing","Singing","Calling")#levels(res$BEHAVIOUR)
PCTBL_abmi <- PCTBL_abmi[keep,]
spp1 <- as.character(taxo$COMMON_NAME[taxo$ORDER_NAME=="Passeriformes"])
spp1 <- spp1[!is.na(spp1)]
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$COMMON_NAME %in% spp1,]
PCTBL_abmi$COMMON_NAME <- PCTBL_abmi$COMMON_NAME[drop=TRUE]

PCTBL_abmi$TBB_TIME_1ST_DETECTED[PCTBL_abmi$TBB_TIME_1ST_DETECTED %in% 
    c("DNC", "NONE", "VNA")] <- NA
PCTBL_abmi$TBB_TIME_1ST_DETECTED <- as.numeric(as.character(PCTBL_abmi$TBB_TIME_1ST_DETECTED))
PCTBL_abmi$period1st <- as.numeric(cut(PCTBL_abmi$TBB_TIME_1ST_DETECTED, c(-1, 200, 400, 600)))

PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_1 %in% c("0","1"),]
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_2 %in% c("0","1"),]
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_3 %in% c("0","1"),]
PCTBL_abmi$TBB_INTERVAL_1 <- as.integer(PCTBL_abmi$TBB_INTERVAL_1) - 1L
PCTBL_abmi$TBB_INTERVAL_2 <- as.integer(PCTBL_abmi$TBB_INTERVAL_2) - 1L
PCTBL_abmi$TBB_INTERVAL_3 <- as.integer(PCTBL_abmi$TBB_INTERVAL_3) - 1L

tmp <- col(matrix(0,nrow(PCTBL_abmi),3)) * 
    PCTBL_abmi[,c("TBB_INTERVAL_1","TBB_INTERVAL_2","TBB_INTERVAL_3")]
tmp[tmp==0] <- NA
tmp <- cbind(999,tmp)
PCTBL_abmi$period123 <- apply(tmp, 1, min, na.rm=TRUE)
with(PCTBL_abmi, table(period1st, period123))
PCTBL_abmi$period1 <- pmin(PCTBL_abmi$period1st, PCTBL_abmi$period123)
with(PCTBL_abmi, table(period1st, period1))
with(PCTBL_abmi, table(period123, period1))

PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$boreal,]

## Species names in ABMI
levels(taxo2[["ENGLISH NAME"]])[levels(taxo2[["ENGLISH NAME"]])=="Blue-headed Vireo"] <- "Blue-headed (solitary) Vireo"
levels(taxo2[["ENGLISH NAME"]])[levels(taxo2[["ENGLISH NAME"]])=="Black-and-white Warbler"] <- "Black and White Warbler"
levels(taxo2[["ENGLISH NAME"]])[levels(taxo2[["ENGLISH NAME"]])=="Baltimore Oriole"] <- "Baltimore (northern) Oriole"
levels(taxo2[["ENGLISH NAME"]])[levels(taxo2[["ENGLISH NAME"]])=="Western Wood-Pewee"] <- "Western Wood Pewee"
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$COMMON_NAME %in% 
    intersect(PCTBL_abmi$COMMON_NAME, taxo2[["ENGLISH NAME"]]),]
PCTBL_abmi$SPECIES <- taxo2$SPECIES[match(PCTBL_abmi$COMMON_NAME, taxo2[["ENGLISH NAME"]])]
PCTBL_abmi$SPECIES <- PCTBL_abmi$SPECIES[drop=TRUE]

PCTBL2 <- with(PCTBL_abmi, data.frame(
    PKEY=as.factor(Label),
    SS=as.factor(Label2),
    PERIOD=as.integer(period1+15),
    DISTANCE=7L,
    SPECIES=SPECIES,
    ABUND=1,
    DISTANCECODE="D",
    DURATIONCODE="X",
    source="ABMI"))
PCTBL3 <- with(PCTBL, data.frame(
    PKEY=PKEY,
    SS=SS,
    PERIOD=DURATION,
    DISTANCE=DISTANCE,
    SPECIES=SPECIES,
    ABUND=ABUND,
    DISTANCECODE=DISTMETH,
    DURATIONCODE=DURMETH,
    source="BAM"))
PCTBL4 <- rbind(PCTBL3, PCTBL2)

if (FALSE) {
PCTBLx <- PCTBL2
PCTBLx$SS <- with(PCTBL_abmi, paste(SiteLabel, "PC", TBB_POINT_COUNT, sep="_"))
PCTBLx$PKEY <- paste(PCTBLx$SS, PCTBL_abmi$YEAR, 1, sep=":")
PCTBLx$PCODE <- "ABMI"
PCTBLx$BEH <- 1
write.csv(PCTBLx, file="ABMI_PCTBL_2003-11.csv")

PCTBL_abmi$SS2 <- with(PCTBL_abmi, paste(SiteLabel, "PC", TBB_POINT_COUNT, sep="_"))
PCTBL_abmi$PKEY2 <- paste(PCTBL_abmi$SS2, PCTBL_abmi$YEAR, 1, sep=":")
PKEYx2 <- nonDuplicated(PCTBL_abmi, PKEY2)
write.csv(PKEYx2, file="ABMI_PKEY_2003-11.csv")
}

## exclude unknown and single duration interval
pctblDur <- PCTBL4[!(PCTBL4$DURATIONCODE %in% c("A","B","D","E","J")) &
    !is.na(PCTBL4$DURATIONCODE),]
rm(PCTBL2, PCTBL3)
gc()

tmp <- PCTBL[,c("DURATION", "DUR_Descrip", "Dur_Start", "DUR_end")]
tmp <- tmp[!is.na(tmp$DURATION),]
durclasses <- nonDuplicated(tmp, DURATION, TRUE)

## drop unknown and >10 min
xtDur <- Xtab(ABUND~PKEY + PERIOD + SPECIES, pctblDur, 
    cdrop=c("3","4","8","9","10")) # 9-10 is for >10 min of K
durclasses[colnames(xtDur[[1]]),]

PKEY2 <- with(PKEY_abmi, data.frame(
    PKEY=Label,
    SS=Label2,
    PCODE="ABMI",
    JDAY=JDAY,
    MaxDuration=10,
    MaxDistance=Inf,
    TSSR=TSSR,
    TSLS=TSLS,
    DURMETH="X"))

## ABMI TSSR is OK since all is relative to Edmonton DST
PKEY3 <- PKEY[,c("PKEY","SS",
    "PCODE","JDAY","MaxDuration",
    "MaxDistance","TSSR","TSLS","DURMETH")]
PKEY4 <- rbind(PKEY3, PKEY2)
rownames(PKEY4) <- PKEY4$PKEY

tmp <- nonDuplicated(pctblDur, PKEY, TRUE)
ii <- intersect(rownames(tmp), rownames(PKEY4))
covarDur <- PKEY4[ii,c("TSSR","DURMETH","JDAY")]
#covarDur <- PKEY4[ii,c("TSLS","TSSR","DURMETH","JDAY")]
#PKEY4$ROAD <- ifelse(PKEY4$NEAR_DIST>=0 &!is.na(PKEY4$NEAR_DIST),1,0)
#covarDur <- PKEY4[ii,c("TSLS","TSSR","DURMETH","ROAD")]
covarDur <- covarDur[rowSums(is.na(covarDur))==0,]
aa <- durclasses[colnames(xtDur[[1]]),]
rm(PKEY2, PKEY3)
gc()

table(covarDur$DURMETH)

durintervals3 <- matrix(c(
    5, 10,NA,NA,            # C
    5, 10,NA,NA,            # F
    3, 5, 10,NA,            # G
    3,10, NA,NA,            # H
    3, 5, NA,NA,            # I
    5, 10,NA,NA,            # K >10 min excluded (PERIOD 9 and 10)
    10,15,NA,NA,            # Q
    3, 5, 8, 10,            # R
    3, 6, 10,NA,            # T (???)
    3.333, 6.666, 10, NA),  # X
    10, 4, byrow=TRUE)
rownames(durintervals3) <- c("C","F","G","H","I","K","Q","R","T","X")
#storage.mode(durintervals3) <- "integer"

NAMES <- list("1"=c("INTERCEPT", "JDAY"),
        "2"=c("INTERCEPT", "TSSR"),
        "3"=c("INTERCEPT", "JDAY", "JDAY2"),
        "4"=c("INTERCEPT", "TSSR", "TSSR2"),
        "5"=c("INTERCEPT", "JDAY", "TSSR"),
        "6"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR"),
        "7"=c("INTERCEPT", "JDAY", "TSSR", "TSSR2"),
        "8"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR", "TSSR2"))
ff <- list(~ JDAY,
    ~ TSSR,
    ~ JDAY + I(JDAY^2),
    ~ TSSR + I(TSSR^2),
    ~ JDAY + TSSR,
    ~ JDAY + I(JDAY^2) + TSSR,
    ~ JDAY + TSSR + I(TSSR^2),
    ~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2))

#SPP <- "OVEN"
fitDurFun <- function(SPP, fit=TRUE) {
    ## get nonzero PCs
    iob <- rowSums(xtDur[[SPP]])>0
    if (sum(iob)==0)
        return(structure("0 observation with multiple duration (1)", 
            class="try-error"))
    if (sum(iob)==1)
        return(structure("1 observation with multiple duration (2)", 
            class="try-error"))
    xt <- xtDur[[SPP]][iob,]
    ## inner join, so only PKEYs with all available covariates are used
    xDur <- Mefa(xt, covarDur, join="inner")
    if (dim(xDur)[1]<=1)
        return(structure("<2 observation with multiple duration (3)", 
            class="try-error"))
    ## sampling protocol IDs
    ii <- samp(xDur)$DURMETH
    ii <- as.character(ii)
#    ii[ii=="F"] <- "C"
    ## pool columns where end times are same
#    Y <- groupSums(xtab(xDur), 2, c(3,5,5,10,10,10,15,20,3,5,10))
    Y <- xtab(xDur)
    ## some how dense matrix is quicker
    Y <- as.matrix(Y)
    ## end times
    D <- durintervals3[as.character(samp(xDur)$DURMETH),]
    ## pieces for each protocol type
    if (any(ii %in% c("C","F"))) {
            WHICH <- ii %in% c("C","F")
            YC <- cbind(matrix(Y[WHICH,c("1","2")], sum(WHICH)), NA, NA)
            iC <- which(WHICH)
        } else {
            YC <- iC <- NULL
        }
    if (any(ii == "G")) {
            WHICH <- ii == "G"
            YG <- cbind(matrix(Y[WHICH,c("5","7","2")], sum(WHICH)), NA)
            iG <- which(WHICH)
        } else {
            YG <- iG <- NULL
        }
    if (any(ii == "H")) {
            WHICH <- ii == "H"
            YH <- cbind(matrix(Y[WHICH,c("5","6")], sum(WHICH)), NA, NA)
            iH <- which(WHICH)
        } else {
            YH <- iH <- NULL
        }
    if (any(ii == "I")) {
            WHICH <- ii == "I"
            YI <- cbind(matrix(Y[WHICH,c("5","7")], sum(WHICH)), NA, NA)
            iI <- which(WHICH)
        } else {
            YI <- iI <- NULL
        }
    if (any(ii == "K")) {
            WHICH <- ii == "K"
            YK <- cbind(matrix(Y[WHICH,c("1","2")], sum(WHICH)), NA, NA)
            iK <- which(WHICH)
        } else {
            YK <- iK <- NULL
        }
#    if (any(ii == "Q")) {
#            WHICH <- ii == "Q"
#            YQ <- cbind(matrix(Y[WHICH,c("10","15")], sum(WHICH)), NA, NA)
#            iQ <- which(WHICH)
#        } else {
#            YQ <- iQ <- NULL
#        }
#    if (any(ii == "R")) {
#            WHICH <- ii == "R"
#            YR <- matrix(Y[WHICH,c("3","5","8","10")], sum(WHICH))
#            iR <- which(WHICH)
#        } else {
#            YR <- iR <- NULL
#        }
    if (any(ii == "X")) {
            WHICH <- ii == "X"
            Yabmi <- cbind(matrix(Y[WHICH,c("16","17","18")], sum(WHICH)), NA)
            iabmi <- which(WHICH)
        } else {
            Yabmi <- iabmi <- NULL
        }
    ## putting pieces together
#    Y <- rbind(YC, YG, YH, YI, YK, YQ, YR, Yabmi)
    Y <- rbind(YC, YG, YH, YI, YK, Yabmi)
    if (is.null(Y))
        return(structure("0 observation with multiple duration (4)", 
            class="try-error"))
#    iall <- c(iC, iG, iH, iI, iK, iQ, iR, iabmi)
    iall <- c(iC, iG, iH, iI, iK, iabmi)
    D <- D[iall,]
    ## integer mode -- faster, but DO NOT use for intervals (<1)
    storage.mode(Y) <- "integer"
    YY <- unname(Y)
    DD <- unname(D)
    X <- samp(xDur)[iall,]
    nn <- table(X$DURMETH)
#    nn <- table(X$DURMETH, X$ROAD)
    if (fit) {
        ## allocate return value
        res <- if (!is.null(ff))
            vector("list", 1+length(ff)) else vector("list", 1)
        ## intercept model
        res[[1]] <- try(D.fit(YY, DD, NULL, type="dur"))
        if (!inherits(res[[1]], "try-error")) {
            res[[1]]$p <- 1
            res[[1]]$names <- "INTERCEPT"
        }
        ## covariates models
        if (!is.null(ff)) {
            for (i in 1:length(ff)) {
                XX <- model.matrix(ff[[i]], X)
                XX <- unname(XX)
                res[[(i+1)]] <- try(D.fit(YY, DD, XX, type="dur"))
                if (!inherits(res[[(i+1)]], "try-error")) {
                    res[[(i+1)]]$p <- length(res[[(i+1)]]$coef)
                    res[[(i+1)]]$names <- NAMES[[i]]
                }
            }
        }
        names(res) <- if (!is.null(ff))
            0:length(ff) else 0
        ## number of observations
        res$n <- nn # attr(res, "nobs") <- nn
    } else {
        res <- list(Y=YY, D=DD, i=iall, n=nn)
    }
    res
}

#str(fitDurFun("OVEN", FALSE))
#str(fitDurFun("OVEN", TRUE))

SPP <- names(xtDur)
resDur <- vector("list", length(SPP))
for (i in 1:length(SPP)) {
    cat("Singing rate estimation for", SPP[i], "\n")
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
date()
for (i in 1:length(SPP)) {
    cat("Singing rate estimation for", SPP[i], "\n")
    flush.console()
    resDur[[i]] <- try(fitDurFun(SPP[i], TRUE))
}
date()
names(resDur) <- SPP
resDurOK <- resDur[!sapply(resDur, inherits, "try-error")]
c(OK=length(resDurOK), failed=length(resDur)-length(resDurOK), all=length(resDur))

#save(resDur, resDurData, file="estimates_SRA_QPADpaper.Rdata")

if (!USE_ROAD)
    save(resDur, resDurData, file="estimates_SRA_QPADpaper.Rdata")

if (USE_ROAD)
    save(resDur, resDurData, file="estimates_SRA_QPADpaper_withRoad.Rdata")

## putting things together

USE_ROAD <- FALSE
n.min <- 75

setwd("c:/Dropbox/bam/DApq3")
if (USE_ROAD) {
    load("estimates_SRA_QPADpaper_withRoad.Rdata")
    load("estimates_EDR_QPADpaper_withRoad.Rdata")
} else {
    load("estimates_SRA_QPADpaper.Rdata")
    load("estimates_EDR_QPADpaper.Rdata")
}

## 0/1 table for successful model fit
edr_mod <- t(sapply(resDist, function(z) 
    ifelse(sapply(z, inherits, what="try-error"), 0, 1)))
edr_mod <- edr_mod[,-ncol(edr_mod)]
sra_mod <- t(sapply(resDur, function(z) 
    ifelse(sapply(z, inherits, what="try-error"), 0, 1)))
sra_mod <- sra_mod[,-ncol(sra_mod)]
tmp <- union(rownames(edr_mod), rownames(sra_mod))
edr_models <- matrix(0, length(tmp), ncol(edr_mod))
dimnames(edr_models) <- list(tmp, colnames(edr_mod))
edr_models[rownames(edr_mod),] <- edr_mod
sra_models <- matrix(0, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod
edr_models[edr_models[,1]==0,] <- 0
sra_models[sra_models[,1]==0,] <- 0
colSums(edr_models)
colSums(sra_models)

edr_nmod <- ncol(edr_mod)
sra_nmod <- ncol(sra_mod)

## sample sizes
edr_n <- sra_n <- numeric(length(tmp))
names(edr_n) <- names(sra_n) <- tmp
edr_nn <- sapply(resDist, function(z) sum(z$n))
edr_n[names(edr_nn)] <- edr_nn
sra_nn <- sapply(resDur, function(z) sum(z$n))
sra_n[names(sra_nn)] <- sra_nn

## spp to keep
#spp <- edr_n >= n.min & sra_n >= n.min &
#    rowSums(edr_models) > 0 & rowSums(sra_models) > 0
spp <- edr_n >= n.min & sra_n >= n.min &
    rowSums(edr_models) > 0 & rowSums(sra_models) > 0
spp <- tmp[spp]

## add in EDR+ species
if (FALSE) {
spp1 <- edr_n >= n.min & rowSums(edr_models)
spp1 <- tmp[spp1]

spp2 <- sra_n >= n.min & rowSums(sra_models) > 0
spp2 <- tmp[spp2]

spp3 <- edr_n >= n.min
spp3 <- tmp[spp3]

spp4 <- sra_n >= n.min
spp4 <- tmp[spp4]

spp5 <- edr_n >= 25
spp5 <- tmp[spp5]

spp6 <- sra_n >= 25
spp6 <- tmp[spp6]

spp7 <- edr_n >= n.min & sra_n >= n.min
spp7 <- tmp[spp7]

spp8 <- edr_n >= 25 & sra_n >= 25
spp8 <- tmp[spp8]

setdiff(spp1,spp)
setdiff(spp2,spp)
setdiff(spp3,spp)
setdiff(spp4,spp)
setdiff(spp5,spp)
setdiff(spp6,spp)
setdiff(spp7,spp)
setdiff(spp8,spp)
}

edr_models <- edr_models[spp,]
sra_models <- sra_models[spp,]
edr_n <- edr_n[spp]
sra_n <- sra_n[spp]

## no. of parameters

edr_df <- sapply(resDist[["OVEN"]][1:edr_nmod], "[[", "p")
sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

## estimates
edr_estimates <- resDist[spp]
for (i in spp)
    edr_estimates[[i]]$n <- NULL
sra_estimates <- resDur[spp]
for (i in spp)
    sra_estimates[[i]]$n <- NULL

load("lifehist_final.Rdata")
rownames(taxo2) <- taxo2$SPECIES
spp_table <- data.frame(spp=spp, 
    scientific_name=taxo2[spp, "SCIENTIFIC NAME"],
    common_name=taxo2[spp, "ENGLISH NAME"])
rownames(spp_table) <- spp

## get variable names for different models
sra_list <- sapply(sra_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))
edr_list <- sapply(edr_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

## loglik values
edr_loglik <- edr_models
edr_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:edr_nmod) { # models
        if (edr_models[i,j] > 0)
            edr_loglik[i,j] <- resDist[[i]][[j]]$loglik
    }
}
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}


## AIC/BIC
edr_aic <- edr_bic <- edr_loglik
edr_aic[] <- Inf
edr_bic[] <- Inf
sra_aic <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_bic[] <- Inf
for (i in spp) {
    edr_aic[i,] <- -2*edr_loglik[i,] + 2*edr_df
    edr_bic[i,] <- -2*edr_loglik[i,] + log(edr_n[i])*edr_df
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
}

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

if (USE_ROAD) {
    ## constrain TREE (-) and ROAD (+) estimates
    for (i in spp) {
        ## ROAD
        if (edr_models[i, "3"] == 1) {
            if (edr_estimates[[i]][["3"]]$coef[2] < 0) {
                edr_aic[i, "3"] <- Inf
                edr_bic[i, "3"] <- Inf
                cat(i, "road\n")
            }
        }
        ## TREE + ROAD
        if (edr_models[i, "4"] == 1) {
            if (edr_estimates[[i]][["4"]]$coef[2] > 0 ||
                edr_estimates[[i]][["4"]]$coef[3] < 0) {
                edr_aic[i, "4"] <- Inf
                edr_bic[i, "4"] <- Inf
                cat(i, "tree+road\n")
            }
        }
        ## LCC + ROAD
        if (edr_models[i, "5"] == 1) {
            if (edr_estimates[[i]][["5"]]$coef[6] < 0) {
                edr_aic[i, "5"] <- Inf
                edr_bic[i, "5"] <- Inf
                cat(i, "lcc+road\n")
            }
        }
    }
}

## model ranking
edr_aicrank <- t(apply(edr_aic, 1, rank))*edr_models
edr_aicrank[edr_aicrank==0] <- NA
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA
edr_bicrank <- t(apply(edr_bic, 1, rank))*edr_models
edr_bicrank[edr_bicrank==0] <- NA
sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA

edr_aicbest <- apply(edr_aicrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_bicbest <- apply(edr_bicrank, 1, function(z) colnames(edr_models)[which.min(z)])
sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])

table(edr_aicbest, sra_aicbest)
rowSums(table(edr_aicbest, sra_aicbest))
colSums(table(edr_aicbest, sra_aicbest))

table(edr_bicbest, sra_bicbest)
rowSums(table(edr_bicbest, sra_bicbest))
colSums(table(edr_bicbest, sra_bicbest))

if (!USE_ROAD)
    version <- "2"
if (USE_ROAD)
    version <- "2"

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

if (USE_ROAD)
    save(.BAMCOEFS, file="BAMCOEFS_QPADpaper_withRoad.rda")

if (!USE_ROAD)
    save(.BAMCOEFS, file="BAMCOEFS_QPADpaper.rda")
