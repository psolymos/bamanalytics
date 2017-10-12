library(knitr)
library(mefa4)
library(pbapply)
library(pROC)
ROOT <- "e:/peter/bam/Apr2016/out"
source("~/repos/bamanalytics/R/makingsense_functions.R")

PROJECT <- "bam"
Date <- "2016-12-01"
#Date <- "2017-04-19"
level <- 0.9

fstat <- function(x, level=0.95) {
    out <- quantile(x, c(0.5, (1-level)/2, 1 - (1-level)/2))
    names(out) <- c("Median", "LCL", "UCL")
    out
}
chfun <- function(Na, Nb, ta, tb) {
    100 * ((Nb/Na)^(1/(tb-ta)) - 1)
}

e <- new.env()
load(file.path(ROOT, "data", paste0("pack_", Date, ".Rdata")), envir=e)

mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,c(Terms, "Units")]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
xn <- xn[rownames(Xn),]
off <- e$OFF[rownames(xn),]
bbb <- unique(e$BB)
bb <- e$BB

rm(e)
INTERNAL <- 1:nrow(xn) %in% bbb
ss1 <- which(!INTERNAL) # test portion
ss2 <- which(INTERNAL) # non-test portion

modTab <- getFancyModsTab(mods)
xnh <- nonDuplicated(xn, HABTR, TRUE)[,c("HAB","HABTR","isNF","isDev",
    "isWet","isOpn","isDM","isDec","isMix")]
xnh <- xnh[c("ConifDense", "ConifSparse","ConifOpen",
    "DecidDense", "DecidSparse", "DecidOpen",
    "MixedDense", "MixedSparse", "MixedOpen",
    "WetDense", "WetSparse", "WetOpen",
    "Shrub", "Grass", "Barren", "Agr", "Devel"),]

spp <- "CAWA"
Y <- yy[,spp]
Y1 <- ifelse(yy[,spp]>0, 1, 0)
off1 <- off[,spp]

fn <- file.path(ROOT, "results", "cawa",
    paste0(PROJECT, "_", spp, "_", Date, ".Rdata"))
load(fn)
#100 * sum(getOK(res)) / length(res)
est_hab <- getEst(res, stage = 2, X=Xn)
est_habhgt <- getEst(res, stage = 3, X=Xn)
est_dtb <- getEst(res, stage = 4, X=Xn)
est_wet <- getEst(res, stage = 5, X=Xn)
est <- getEst(res, stage = length(mods)-1, X=Xn)
est_yr <- getEst(res, stage = length(mods), X=Xn)
#dcawa <- data.frame(PKEY=rownames(xn), cawa=Y)
#dcawa$validation <- 0
#dcawa$validation[ss1] <- 1
#write.csv(dcawa, row.names=FALSE, file="w:/bam-cawa/cawa-pkeys-2016-12-13.csv")

## Population size estimates
## All numbers are given in million males.

x0 <- read.csv("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-bf-totals.csv")
x1 <- read.csv("e:/peter/bam/pred-2016/maps/CAWA-6-2002-2016-12-01-2002-totals.csv")
x2 <- read.csv("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-totals.csv")
rownames(x0) <- x0$subreg
rownames(x1) <- x1$subreg
rownames(x2) <- x2$subreg
x1 <- x1[rownames(x2),]
x0 <- x0[rownames(x2),]

ci <- grepl("CAWA", colnames(x0))
CAN <- c("ALBERTA", "BRITISH COLUMBIA", "MANITOBA",
    "NEW BRUNSWICK", "NEWFOUNDLAND",
    "NORTHWEST TERRITORIES", "NOVA SCOTIA", "NUNAVUT",
    "ONTARIO", "PRINCE EDWARD ISLAND", "QUEBEC", "SASKATCHEWAN", "YUKON")

## Full study area
## pop size in full study area
p0 <- fstat(colSums(x0[,ci])/10^6, level)
p1 <- fstat(colSums(x1[,ci])/10^6, level)
p2 <- fstat(colSums(x2[,ci])/10^6, level)
write.csv(cbind("nss"=sum(x2[,"nSSinSubreg"]), "ndet"=sum(x2[,"nDETinSubreg"]),
    rbind("Backfilled"=p0, "2002"=p1, "2012"=p2)),
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_full.csv")
fstat(100*((colSums(x2[,ci])/colSums(x1[,ci]))^(1/(2012-2002)) -1), level)

## Canada
ss <- x2$JURS %in% CAN
p0 <- fstat(colSums(x0[ss,ci])/10^6, level)
p1 <- fstat(colSums(x1[ss,ci])/10^6, level)
p2 <- fstat(colSums(x2[ss,ci])/10^6, level)
write.csv(cbind("nss"=sum(x2[ss,"nSSinSubreg"]), "ndet"=sum(x2[ss,"nDETinSubreg"]),
    rbind("Backfilled"=p0, "2002"=p1, "2012"=p2)),
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_Canada.csv")
fstat(100*((colSums(x2[ss,ci])/colSums(x1[ss,ci]))^(1/(2012-2002)) -1), level)

## Boreal
ss <- x2$Brandt %in% c("B_ALPINE","BOREAL")
p0 <- fstat(colSums(x0[ss,ci])/10^6, level)
p1 <- fstat(colSums(x1[ss,ci])/10^6, level)
p2 <- fstat(colSums(x2[ss,ci])/10^6, level)
write.csv(cbind("nss"=sum(x2[ss,"nSSinSubreg"]), "ndet"=sum(x2[ss,"nDETinSubreg"]),
    rbind("Backfilled"=p0, "2002"=p1, "2012"=p2)),
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_Boreal.csv")
fstat(100*((colSums(x2[ss,ci])/colSums(x1[ss,ci]))^(1/(2012-2002)) -1), level)

## Hemiboreal
ss <- x2$Brandt %in% c("H_ALPINE","HEMIBOREAL")
p0 <- fstat(colSums(x0[ss,ci])/10^6, level)
p1 <- fstat(colSums(x1[ss,ci])/10^6, level)
p2 <- fstat(colSums(x2[ss,ci])/10^6, level)
write.csv(cbind("nss"=sum(x2[ss,"nSSinSubreg"]), "ndet"=sum(x2[ss,"nDETinSubreg"]),
    rbind("Backfilled"=p0, "2002"=p1, "2012"=p2)),
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_Hemiboreal.csv")
fstat(100*((colSums(x2[ss,ci])/colSums(x1[ss,ci]))^(1/(2012-2002)) -1), level)

## By jurisdiction (2012 conditions, within study area)
tmp <- as.matrix(x2[,c("nSSinSubreg","nDETinSubreg")])
colnames(tmp) <- c("nss","ndet")
by_jurs <- data.frame(
    groupSums(tmp, 1, x2$JURS),
    t(apply(groupSums(as.matrix(x2[,ci])/10^6, 1, x2$JURS), 1, fstat)))
by_jurs$Proportion <- by_jurs$Median / sum(by_jurs$Median)
write.csv(by_jurs[order(rownames(by_jurs)),],
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_Jurs.csv")

## By BCR (2012 conditions, within study area)
x2$bcr <- paste0("BCR:", as.character(x2$BCR))
by_bcr <- data.frame(
    groupSums(tmp, 1, x2$bcr),
    t(apply(groupSums(as.matrix(x2[,ci])/10^6, 1, x2$bcr), 1, fstat)))
by_bcr$Proportion <- by_bcr$Median / sum(by_bcr$Median)
write.csv(by_bcr,
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_BCR.csv")

## By BCR/jurisdiction (2012 conditions, within study area)
x2$bcrjurs <- paste0(as.character(x2$JURS), "-BCR:", as.character(x2$BCR))
by_bcrjurs <- data.frame(
    groupSums(tmp, 1, x2$bcrjurs),
    t(apply(groupSums(as.matrix(x2[,ci])/10^6, 1, x2$bcrjurs), 1, fstat)))
by_bcrjurs$Proportion <- by_bcrjurs$Median / sum(by_bcrjurs$Median)
write.csv(by_bcrjurs,
    file="e:/peter/bam/Apr2016/out/cawa/cawa-pop-est_BCR-Jurs.csv")

## Residual trend estimates

### Canada

```{r echo=FALSE}
trall <- read.csv("e:/peter/bam/Apr2016/out/cawa/cawa-trend-2016-12-01.csv")
tdall <- read.csv("e:/peter/bam/Apr2016/out/cawa/cawa-det-2016-12-01.csv")
tr1 <- trall[trall$Country != "",]
tr1 <- tr1[,c("Data", "Mean", "Median", "X2.5.", "X97.5.")]
colnames(tr1) <- c("Data", "Mean", "Median", "LCL", "UCL")
tr1$ndet <- colSums(tdall[tdall$Country=="CAN", c("ndet_all", "ndet_bam", "ndet_bbs")])
kable(tr1)
```

### BCR

```{r echo=FALSE}
tr2 <- trall[!is.na(trall$BCR),]
tr2 <- data.frame(BCR=tr2$BCR, tr2[,c("Data", "Mean", "Median", "X2.5.", "X97.5.")])
colnames(tr2) <- c("BCR", "Data", "Mean", "Median", "LCL", "UCL")
nn <- groupSums(as.matrix(tdall[, c("ndet_all", "ndet_bam", "ndet_bbs")]),
    1, tdall$BCR)
colnames(nn) <- c("all", "bam", "bbs")
nn <- Melt(nn)
tr2$ndet <- nn$value[match(paste(tr2$Data, tr2$BCR), paste(nn$cols, nn$rows))]
tr2$ndet[is.na(tr2$ndet)] <- 0
tr2 <- tr2[tr2$ndet > 0,]
kable(tr2[order(tr2$BCR, tr2$Data),], row.names=FALSE)
```

### Jurisdiction

```{r echo=FALSE}
tr3 <- trall[trall$PROV != "",]
tr3 <- data.frame(Jurisdiction=tr3$PROV, tr3[,c("Data", "Mean", "Median", "X2.5.", "X97.5.")])
colnames(tr3) <- c("Jurisdiction", "Data", "Mean", "Median", "LCL", "UCL")
nn <- groupSums(as.matrix(tdall[, c("ndet_all", "ndet_bam", "ndet_bbs")]),
    1, tdall$JURS)
colnames(nn) <- c("all", "bam", "bbs")
nn <- Melt(nn)
tr3$ndet <- nn$value[match(paste(tr3$Data, tr3$Jurisdiction), paste(nn$cols, nn$rows))]
tr3$ndet[is.na(tr3$ndet)] <- 0
tr3 <- tr3[tr3$ndet > 0,]
kable(tr3[order(tr3$Jurisdiction, tr3$Data),], row.names=FALSE)
```


### BCR & Jurisdiction

```{r echo=FALSE}
tr4 <- trall[trall$BCRPROV != "",]
tr4 <- data.frame(BCRJurs=tr4$BCRPROV, tr4[,c("Data", "Mean", "Median", "X2.5.", "X97.5.")])
colnames(tr4) <- c("BCRJurs", "Data", "Mean", "Median", "LCL", "UCL")
nn <- groupSums(as.matrix(tdall[, c("ndet_all", "ndet_bam", "ndet_bbs")]),
    1, tdall$BCRJURS)
colnames(nn) <- c("all", "bam", "bbs")
nn <- Melt(nn)
tr4$ndet <- nn$value[match(paste(tr4$Data, tr4$BCRJurs), paste(nn$cols, nn$rows))]
tr4$ndet[is.na(tr4$ndet)] <- 0
tr4 <- tr4[tr4$ndet > 0,]
kable(tr4[order(tr4$BCRJurs, tr4$Data),], row.names=FALSE)
```

