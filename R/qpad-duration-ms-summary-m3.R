## this version incorporates 3 model types
## M0, Mf, Mb
## Define root folder where data are stored
ROOT <- "e:/peter/bam/May2015"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2017"

## Load required packages
library(MASS)
library(mefa4)
library(detect)

## Load functions kept in separate file
#source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2017-03-01.Rdata"))

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","YEAR","PCODE","X","Y","SS")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])

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
names(ff) <- names(NAMES)

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL

## map
library(rworldmap)
rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
X0 <- pkDur[rn,]
D <- ltdur$end[match(X0$DURMETH, rownames(ltdur$end)),]
png(file.path(ROOT2, "tabfig", "Fig1_map.png"), height=600, width=800)
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(pkDur[, c("X","Y")], pch=".",
    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255))
points(X0[rowSums(!is.na(D)) > 1, c("X","Y")], pch=19,
    col="red", cex=0.3)
dev.off()

ld <- data.frame(table(pkDur$DURMETH))
rownames(ld) <- ld[,1]
ld <- data.frame(ld, ltdur$end[rownames(ld),])
ld <- ld[rowSums(!is.na(ltdur$end[rownames(ld),])) > 1,]

datx <- droplevels(dat[rownames(X0)[rowSums(!is.na(D)) > 1],])

## data summaries
pkDurOK <- droplevels(pkDur[pkDur$DURMETH %in% rownames(ld),])
str(pkDurOK)

if (FALSE) {
## export data for Fig 1
ss <- droplevels(nonDuplicated(pkDurOK, SS))
#ss <- droplevels(pkDurOK)
ss <- ss[,c("SS","PCODE","DURMETH","X","Y")]

ll <- apply(ltdur$end[rownames(ld),], 1, function(z) {
    paste(c(0, z[!is.na(z)]), collapse="-")
    })
ss$INTERVALS <- as.factor(ll[match(ss$DURMETH, names(ll))])
l <- data.frame(table(ss$INTERVALS))
l$Legend <- paste0(as.character(l$Var1), " (n=", l$Freq, ")")

write.csv(ss, row.names=FALSE, file="~/Downloads/removal-ms-points-for-fig-1.csv")
}

## species
e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_rem.rda"), envir=e)
.BAMCOEFSrem <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_fmix.rda"), envir=e)
.BAMCOEFSfmix <- e$.BAMCOEFS

## species where rem model sample size is at least NMIN

## check here !!!
sb <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")

rownames(sb) <- sb$Species_ID

compare_sets(names(.BAMCOEFSrem$sra_n), names(.BAMCOEFSmix$sra_n))

NMIN <- 75

SPPfull <- sort(names(.BAMCOEFSrem$sra_n)[.BAMCOEFSrem$sra_n >= NMIN])
table(sb[SPPfull, "Singing_birds"])
## this comes from checking M0/Mb estimates (all BAM scale)
EXC1 <- c("CBCH", "CORE", "PAWR", "PSFL", "RBSA")
SPPfull <- SPPfull[!(SPPfull %in% EXC1)]


sptab <- .BAMCOEFSrem$spp_table[SPPfull,]
sptab$nfull <- .BAMCOEFSrem$sra_n[SPPfull]

SPPmix <- sort(names(.BAMCOEFSmix$sra_n)[.BAMCOEFSmix$sra_n >= NMIN])
EXC2 <- c("BEKI", "BOBO",
    "BRBL", "CLSW", "DUFL", "PAWR", "PIGR", "RBWO",
    "RNPH", "ROPI", "SAPH", "STGR", "VGSW")
SPPmix <- SPPmix[!(SPPmix %in% EXC2)]

sptab$model <- factor("rem", c("rem","mix","both"))
sptab[SPPmix, "model"] <- "both"

SPP <- rownames(sptab)[sptab$model=="both"]

cfall0 <- sapply(SPPfull, function(spp) {
    exp(unname(.BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients))
    })
dim(cfall0) <- c(length(SPPfull), 1)
dimnames(cfall0) <- list(SPPfull, "phi_0")

cfallb <- t(sapply(SPP, function(spp) {
    tmp <- unname(.BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients)
    c(phi_b=exp(tmp[1]), c=plogis(tmp[2]))
    }))

sptab$M0_phi <- cfall0[,1]
sptab$Mb_phi <- cfallb[match(SPPfull, SPP),"phi_b"]
sptab$Mb_c <- cfallb[match(SPPfull, SPP),"c"]

Projects <- table(droplevels(X0$PCODE))
SppOcc <- sapply(SPPfull, function(z) sum(rowSums(xtDur[[z]][rn,], na.rm=TRUE)>0)/length(rn))
SppYmean <- sapply(SPPfull, function(z) mean(rowSums(xtDur[[z]][rn,], na.rm=TRUE)))


## sra 0 vs f or b models

## aic is really AICc to account for small sample sizes
aic0 <- .BAMCOEFSrem$sra_aic
aicf <- .BAMCOEFSfmix$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
df0 <- matrix(.BAMCOEFSrem$sra_df, nrow(aic0), 15, byrow=TRUE)
dff <- matrix(.BAMCOEFSfmix$sra_df, nrow(aicb), 15, byrow=TRUE)
dfb <- matrix(.BAMCOEFSmix$sra_df, nrow(aicb), 15, byrow=TRUE)
n0 <- .BAMCOEFSrem$sra_n
nf <- .BAMCOEFSfmix$sra_n
nb <- .BAMCOEFSmix$sra_n

aicc0 <- aic0 + (2*df0*(df0+1)) / (n0-df0-1)
aiccf <- aicf + (2*dff*(dff+1)) / (nf-dff-1)
aiccb <- aicb + (2*dfb*(dfb+1)) / (nb-dfb-1)
colnames(aicc0) <- paste0("M0_", colnames(aic0))
colnames(aiccf) <- paste0("Mf_", colnames(aicf))
colnames(aiccb) <- paste0("Mb_", colnames(aicb))

aic0 <- aicc0[SPPfull,]
aicf <- aiccf[SPP,]
aicb <- aiccb[SPP,]
aicx <- cbind(aicc0[SPP,], aiccb[SPP,], aiccf[SPP,]) # all 3 combined

waic0 <- t(apply(aic0, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicf <- t(apply(aicf, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicb <- t(apply(aicb, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicx <- t(apply(aicx, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))

H0 <- apply(waic0, 1, function(z) sum(z^2))
Hf <- apply(waicf, 1, function(z) sum(z^2))
Hb <- apply(waicb, 1, function(z) sum(z^2))
Hx <- apply(waicx, 1, function(z) sum(z^2))

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestf <- as.character(0:14)[apply(aicf, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
bestx <- as.factor(colnames(aicx))[apply(aicx, 1, which.min)]
## these 2 are identical, shouldnt have both
bestx[bestx %in% c("Mf_0", "Mb_0")]
names(best0) <- SPPfull
names(bestb) <- names(bestf) <- names(bestx) <- SPP

tmp <- strsplit(as.character(bestx), "_")
dfb <- data.frame(spp=SPP, best=bestx, type=sapply(tmp, "[[", 1),
    model=as.integer(sapply(tmp, "[[", 2)))
btab <- table(dfb$model, dfb$type)
rownames(btab) <- paste(rownames(btab), as.character(ff))
btab <- addmargins(btab)

tb <- read.csv("~/Dropbox/bam/duration_ms/revisionMarch2017/tabfig/spptab.csv")
rownames(tb) <- tb$spp
tb$Mf_best <- as.integer(bestf)[match(rownames(tb), names(bestf))]
tb$Best3 <- bestx[match(rownames(tb), names(bestx))]
write.csv(tb, row.names=FALSE, file="~/Dropbox/bam/duration_ms/revisionMarch2017/tabfig/spptab3x.csv")
