##---
##title: "Maritimes predictions"
##author: "Peter Solymos"
##date: "June 26, 2015"
##output: pdf_document
##---

## Define root folder where data are stored
ROOT <- "d:/bam/May2015"
library(mefa4)
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## CAS data is in old files
x50v <- read.csv(file.path(ROOT, "maritimes/June19", "pred50m_casfri.csv"))
x100v <- read.csv(file.path(ROOT, "maritimes/June19", "pred100m_casfri.csv"))
x250v <- read.csv(file.path(ROOT, "maritimes/June19", "pred250m_casfri.csv"))
dd <- read.csv(file.path(ROOT, "maritimes/June19", "pred250m_dist.csv"))

rownames(x50v) <- x50v$POINTID
rownames(x100v) <- x100v$POINTID
rownames(x250v) <- x250v$POINTID
rownames(dd) <- dd$POINTID
pid <- rownames(x50v)
compare.sets(pid, rownames(x100v))
compare.sets(pid, rownames(x250v))
compare.sets(pid, rownames(dd))

x100v <- x100v[pid,]
x250v <- data.frame(x250v[pid,], dd[pid,])

## Load intersection files
x50 <- read.csv(file.path(ROOT, "maritimes/June19", "pred50m.csv"))
x100 <- read.csv(file.path(ROOT, "maritimes/June19", "pred100m.csv"))
x250 <- read.csv(file.path(ROOT, "maritimes/June19", "pred250m.csv"))
rownames(x50) <- x50v$OBJECTID
rownames(x100) <- x100v$OBJECTID
rownames(x250) <- x250v$OBJECTID
x50 <- x50[pid,]
x100 <- x100[pid,]
x250 <- x250[pid,]

compare.sets(pid, rownames(x50))
compare.sets(pid, rownames(x100))
compare.sets(pid, rownames(x250))

x50 <- data.frame(x50, x50v)
x100 <- data.frame(x100, x100v)
x250 <- data.frame(x250, x250v)
rm(x50v,x100v,x250v)

## CAS based tree species
tree_250 <- c("Abie_bals", "Abie_pice", "Abie_spp", "Acer_rubr",
    "Acer_sacc", "Acer_spp", "Alnu_spp", "Betu_alle", "Betu_papy",
    "Betu_popu", "Betu_spp", "Fagu_gran", "Frax_spp", "Hard_into",
    "Hard_nonc", "Hard_tole", "Hard_unkn", "Lari_deci", "Lari_lari",
    "NOSC_HARD", "NOSC_SOFT", "Pice_abie", "Pice_glau", "Pice_mari",
    "Pice_rube", "Pice_spp", "Pinu_bank", "Pinu_resi", "Pinu_spp",
    "Pinu_stro", "Pinu_sylv", "Popu_balb", "Popu_spp", "Popu_trem",
    "Prun_sero", "Quer_rubr", "Soft_unkn", "Thuj_occi", "Tsug_cana",
    "Ulmu_amer", "Uncl_spp",
    ## new to preds
    "Lari_kaem")
colnames(x50) <- sub("\\.", "_", colnames(x50))
colnames(x100) <- sub("\\.", "_", colnames(x100))
colnames(x250) <- sub("\\.", "_", colnames(x250))
## Checking if 50 and 100 m buffer has no tree species
## that that is unique for a <250 m scale (this should not happen, and seems OK)
setdiff(colnames(x50), tree_250)
setdiff(colnames(x100), tree_250)
setdiff(colnames(x250), tree_250)
setdiff(tree_250, colnames(x250))
## CAS related fields from the 250 m buffer
#cas_250 <- c(tree_250, "CRCL_AV", "CRCL_STD", "HT_AV", "HT_STD")
## Checking 250 m based percentages for leading tree species
xp <- as.matrix(x250[,tree_250[tree_250 %in% colnames(x250)]])
xp[is.na(xp)] <- 0
range(xp, na.rm=TRUE)
## Divide areas (m^2) by buffer size (250^2*pi = 196349.5)
xp <- xp / (250^2*pi)
stopifnot(all(xp <= 1))
sum(is.na(xp))
## Percentages sometimes add up to >1, but this should not affect
## finding the dominant species (maximum do not change by standardization,
## but we might want to know percentage for comparisons)
hist(rowSums(xp))
summary(rowSums(xp))
sum(rowSums(xp) > 1)
## Standardize with row sums where it is >1
#xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
#summary(rowSums(xp))
## Leading species
tlead <- as.factor(tree_250[apply(xp, 1, which.max)])
data.frame(sort(table(tlead)))
## Percent cover of leading species
plead <- apply(xp, 1, function(z) z[which.max(z)])

## Leading human footprint type in 250 m buffer
hflab <- c("BU","CO","OT","PC","SI","WF")
h <- as.matrix(x250[,colnames(x250) %in% hflab])
h[is.na(h)] <- 0
hist(rowSums(h))
h <- cbind(Undist=1-rowSums(h),
    groupSums(h, 2, c("Burn","Cut","Other","Cut","Other","Other")))
x250$ldist <- factor(colnames(h)[apply(h, 1, which.max)], colnames(h))
table(x250$ldist)

## Standardization
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) %in% c("", "No Data")] <- "Average" # lots of cases
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Below Average"] <- "-1"
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Average"] <- "0"
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Above Average"] <- "1"
x250$COMPLEXITY <- as.integer(as.character(x250$COMPLEXITY))

x250$CONNECTEDNESS <- x250$CONNECTIVITY / 100

x50$CROWNCL_AV <- x50$CROWNCL_AV / 100
x50$CROWNCL_STD <- x50$CROWNCL_STD / 100
x100$CROWNCL_AV <- x100$CROWNCL_AV / 100
x100$CROWNCL_STD <- x100$CROWNCL_STD / 100
x250$CRCL_AV <- x250$CROWNCL_AV / 100
x250$CRCL_STD <- x250$CROWNCL_STD / 100

x250$DTW_STD <- x250$DTW_STD / 10

#x250$HT_AV <- x250$HT_AV / 10 ----------------- missing
## used the one from x100 (where there was no STD counterpart)
x250$HT_AV <- x100$HEIGHT_AV / 10
x250$HT_STD <- x250$HT_STD / 10

x250$HUMAN_FOOTPRINT <- x250$HUMAN_FOOTPRINT / 100

x50$ROAD01 <- ifelse(x50$ROAD <= 400, 1L, 0L)
x100$ROAD01 <- ifelse(x100$ROAD <= 400, 1L, 0L)

x250$WET_LENGTH <- x250$WET_LENGTH / 1000

levels(x50$WET_VEG)[levels(x50$WET_VEG) == "WATER/EXPOSED"] <- "AQUATIC"
levels(x50$WET_VEG)[levels(x50$WET_VEG) %in%
    c("GRAMINOI","GRAMINOID","GRAMNINOID")] <- "GRAMINOID"
levels(x100$WET_VEG)[levels(x100$WET_VEG) == "WATER/EXPOSED"] <- "AQUATIC"
levels(x100$WET_VEG)[levels(x100$WET_VEG) %in%
    c("GRAMINOI","GRAMINOID","GRAMNINOID")] <- "GRAMINOID"

## Reclass leading tree species
lt <- read.csv(file.path(ROOT, "maritimes", "cas-tree-lookup.csv"))
rownames(lt) <- lt$spp
lt <- lt[colnames(xp),]
xp2 <- groupSums(xp, 2, lt$spp3)
xp2 <- xp2[,colnames(xp2) != "XXX"]
data.frame(tree=sort(table(as.factor(colnames(xp2)[apply(xp2, 1, which.max)]))))
x250$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], c(colnames(xp2), "NOCAS"))
x250$ltree[!x250$keep] <- "NOCAS" # no tree data
## Reclass at 50 m scale
xp <- as.matrix(x50[,colnames(x50) %in% tree_250])
xp[is.na(xp)] <- 0
xp <- xp / (50^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x50$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
x50$ltree[!x250$keep] <- "NOCAS" # no tree data
## Reclass at 100 m scale
xp <- as.matrix(x100[,colnames(x100) %in% tree_250])
xp[is.na(xp)] <- 0
xp <- xp / (100^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x100$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
x100$ltree[!x250$keep] <- "NOCAS" # no tree data
## Compare scales
table(x50$ltree, x250$ltree)
table(x100$ltree, x250$ltree)
## Get rid of tree cover info
x50 <- x50[,!(colnames(x50) %in% tree_250)]
x100 <- x100[,!(colnames(x100) %in% tree_250)]
x250 <- x250[,!(colnames(x250) %in% tree_250)]
## Define fields that do not change with buffer
asis <- c("OBJECTID", "ORIG_FID", "PROTECT")
## Renaming columns in 50 and 100 m tables to reflect
## 'local' scale (as opposed to 'territory' scale et 250 m)
colnames(x50)[!(colnames(x50) %in% asis)] <- paste("LOC",
    colnames(x50)[!(colnames(x50) %in% asis)], sep="_")
colnames(x100)[!(colnames(x100) %in% asis)] <- paste("LOC",
    colnames(x100)[!(colnames(x100) %in% asis)], sep="_")
xx50 <- data.frame(x50, x250[,!(colnames(x250) %in% asis)])
xx100 <- data.frame(x100, x250[,!(colnames(x250) %in% asis)])

#xx50$PROTECT <- xx50$LOC_PROTECT
#xx100$PROTECT <- xx100$LOC_PROTECT
#xx50$LOC_PROTECT <- NULL
#xx100$LOC_PROTECT <- NULL

xx50$CONNECTEDNESS[is.na(xx50$CONNECTEDNESS)] <- 0
xx100$CONNECTEDNESS[is.na(xx100$CONNECTEDNESS)] <- 0
xx50$CONNECTIVITY <- NULL
xx100$CONNECTIVITY <- NULL

xx50$HT_AV[is.na(xx50$HT_AV)] <- 0
xx100$HT_AV[is.na(xx100$HT_AV)] <- 0

xx50$HT_STD[is.na(xx50$HT_STD)] <- 0
xx100$HT_STD[is.na(xx100$HT_STD)] <- 0

colSums(is.na(xx50))
colSums(is.na(xx100))

save(xx50, xx100,
    file=file.path(ROOT, "out", "maritimes_preds.Rdata"))

### Predictions

ROOT <- "d:/bam/May2015"
library(mefa4)
#load(file.path(ROOT, "out", "maritimes_3spp.Rdata"))
#load("~/Dropbox/bam/maritimes2015/maritimes_preds.Rdata")
load("d:/bam/bam-projects/maritimes2015/maritimes_preds.Rdata")
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/bamanalytics/R/maritimes_mods.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

## Check if all variables defined in the
## 3 sets of model lists can be found in data

TermsA <- getTerms(modsA, "list")
setdiff(TermsA, colnames(xx50))
setdiff(TermsA, colnames(xx100))

TermsB <- getTerms(modsB, "list")
setdiff(TermsB, colnames(xx50))
setdiff(TermsB, colnames(xx100))

TermsC <- getTerms(modsC, "list")
setdiff(TermsC, colnames(xx50))
setdiff(TermsC, colnames(xx100))

Terms <- unique(c(TermsA, TermsB, TermsC))

## Summary of variables used in modeling
summary(xx50[,sort(Terms)])
colSums(is.na(xx50[,sort(Terms)]))
summary(xx100[,sort(Terms)])
colSums(is.na(xx100[,sort(Terms)]))


load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesA_CAWA_Allrun.Rdata"))
cawa <- list(A=res)
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesB_CAWA_Allrun.Rdata"))
cawa$B <- res
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesC_CAWA_Allrun.Rdata"))
cawa$C <- res

load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesA_RUBL_Allrun.Rdata"))
rubl <- list(A=res)
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesB_RUBL_Allrun.Rdata"))
rubl$B <- res
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesC_RUBL_Allrun.Rdata"))
rubl$C <- res

load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesA_OSFL_Allrun.Rdata"))
osfl <- list(A=res)
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesB_OSFL_Allrun.Rdata"))
osfl$B <- res
load(file.path("d:/bam/bam-projects/maritimes2015",
    "BirdCoefs_MaritimesC_OSFL_Allrun.Rdata"))
osfl$C <- res

allres <- list(CAWA=cawa, RUBL=rubl, OSFL=osfl)


xn_50B <- xx50[,TermsB]
xn_100C <- xx100[,TermsC]
Xn_50B <- model.matrix(getTerms(modsB, "formula"), xn_50B)
Xn_100C <- model.matrix(getTerms(modsC, "formula"), xn_100C)

## ----- species specific

spp <- "OSFL"
for (spp in c("CAWA","RUBL","OSFL")) { # loop for spp start

if (spp == "OSFL") {
    xn <- xn_100C
    Xn <- Xn_100C
    res <- allres[[spp]]$C
    mods <- modsC
    xx <- xx100
} else {
    xn <- xn_50B
    Xn <- Xn_50B
    res <- allres[[spp]]$B
    mods <- modsB
    xx <- xx50
}
colnames(Xn) <- fixNames(colnames(Xn))


est <- getEst(res)

capture.output(printCoefmat(getSummary(res), 3),
    file=paste0("~/Dropbox/bam/maritimes2015/preds_", spp, ".txt"))

mu <- getDataPred(res)
bmu <- rowMeans(exp(mu))
SD <- apply(exp(mu), 1, sd)
ci <- t(apply(exp(mu), 1, quantile, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm=TRUE))

CL <- rgb(210, 180, 140, alpha=1*255, max=255)
CLa <- rgb(210, 180, 140, alpha=0.25*255, max=255)

pdf(paste0("~/Dropbox/bam/maritimes2015/preds_", spp, ".pdf"), onefile=TRUE)
op <- par(las=2, mar=c(5,4,2,1)+0.1)
for (i in 1:ncol(xn)) {
    if (is.factor(xn[,i])) {
        boxplot(bmu ~ xn[,i], range=0, col=CL, main=spp, xlab=colnames(xn)[i],
            ylab="Density (males / ha)")
        abline(h=median(bmu), col="red4", lty=2)
    } else {
        if (length(unique(xn[,i])) < 5) {
            boxplot(bmu ~ xn[,i], range=0, col=CL, main=spp, xlab=colnames(xn)[i],
                ylab="Density (males / ha)")
            abline(h=median(bmu), col="red4", lty=2)
        } else {
            ii <- sample.int(nrow(xn), 5000)
            plot(xn[ii,i], bmu[ii], col=CLa, pch=19,
                main=spp, xlab=colnames(xn)[i],
                ylab="Density (males / ha)")
            lines(lowess(xn[,i], bmu), col="red4")
        }
    }
}
par(op)
dev.off()

out <- data.frame(Mean=bmu, SD=SD, ci)
colnames(out) <- paste(spp, c("Mean", "SD",
        2.5, 5, 25, 50, 75, 95, 97.5), sep="_")
out <- data.frame(xx, out)
write.csv(out, paste0("~/Dropbox/bam/maritimes2015/preds_", spp, ".csv"),
    row.names=FALSE)

## OSFL interaction between ltree & DTW_PROP
if (spp == "OSFL") {

lt <- droplevels(xn$ltree)
pdf(paste0("~/Dropbox/bam/maritimes2015/preds_", spp, "_interaction.pdf"),
    width=6, height=12)
op <- par(mar=c(5,4,2,1)+0.1, mfrow=c(5,2))
for (i in levels(lt)) {
    plot(xn$DTW_PROP[lt == i], bmu[lt == i], col=CLa, pch=19,
        main=i, xlab="DTW_PROP",
        ylab="Density (males / ha)", ylim=quantile(bmu, c(0,0.999)), xlim=range(xn$DTW_PROP))
    lines(lowess(xn$DTW_PROP[lt == i], bmu[lt == i]), col="red4")
}
par(op)
dev.off()

}

} ## loop for spp end


cawaA <- read.csv("d:/bam/bam-projects/maritimes2015/res-cawa-coef-A.csv")
cawaB <- read.csv("d:/bam/bam-projects/maritimes2015/res-cawa-coef-B.csv")
cawaC <- read.csv("d:/bam/bam-projects/maritimes2015/res-cawa-coef-C.csv")
osflA <- read.csv("d:/bam/bam-projects/maritimes2015/res-osfl-coef-A.csv")
osflB <- read.csv("d:/bam/bam-projects/maritimes2015/res-osfl-coef-B.csv")
osflC <- read.csv("d:/bam/bam-projects/maritimes2015/res-osfl-coef-C.csv")

f <- function(x) colMeans(x[,colSums(abs(x)) > 0])


## Table 4 edits
getEst <-
function(res, stage=NULL, na.out=TRUE, X) {
    if (!missing(X))
        Xn <- X
    OK <- !sapply(res, inherits, "try-error")
    if (any(!OK))
        warning(paste("try-error found:", sum(!OK)))
    ii <- sapply(res[OK], "[[", "iteration")
    est <- Xn[1:length(ii),,drop=FALSE]
    rownames(est) <- ii
    est[] <- 0
    if (is.null(stage))
        stage <- length(res[[ii[1]]]$coef)
    if (stage > 0) {
        for (i in 1:length(ii)) {
            tmp <- res[[ii[i]]]$coef[[stage]]
            names(tmp) <- fixNames(names(tmp))
            sdiff <- setdiff(names(tmp), colnames(est))
            if (length(sdiff) > 0)
                stop(paste(sdiff, collapse=" "))
            est[i,match(names(tmp), colnames(est))] <- tmp
        }
    } else {
        for (i in 1:length(ii)) {
            est[i,1] <- res[[ii[i]]]$null
        }
    }
    if (any(!OK) && na.out) {
        nas <- matrix(NA, sum(!OK), ncol(est))
        rownames(nas) <- which(!OK)
        est <- rbind(est, nas)
    }
    est
}

s <- list()
for (spp in c("CAWA","OSFL")) {
    s[[spp]] <- list()
    for (M in c("A", "B", "C")) {
        res <- allres[[spp]][[M]]
        mods <- list(A=modsA, B=modsB, C=modsC)[[M]]
        Terms <- list(A=TermsA, B=TermsB, C=TermsC)[[M]]
        if (spp == "OSFL") {
            xx <- xx100
        } else {
            xx <- xx50
        }
        xn <- xx[,Terms]
        Xn <- model.matrix(getTerms(mods, "formula"), xn)
        colnames(Xn) <- fixNames(colnames(Xn))
        est <- getEst(res, X=Xn)
        s[[spp]][[M]] <- getSummary(res)
    }
}


xn_50B <- xx50[,TermsB]
xn_100C <- xx100[,TermsC]
Xn_50B <- model.matrix(getTerms(modsB, "formula"), xn_50B)
Xn_100C <- model.matrix(getTerms(modsC, "formula"), xn_100C)

