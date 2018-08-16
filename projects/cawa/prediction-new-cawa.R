library(RColorBrewer)
library(mefa4)
ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/bam/pred-2015"
ROOT3 <- "e:/peter/bam/pred-2016"
source("~/repos/bamanalytics/R/makingsense_functions.R")
#source("~/repos/bamanalytics/R/analysis_mods.R")

PROJECT <- "bam"
Date <- "2016-12-01"
#Date <- "2017-01-27"
#level <- 0.9
#spp <- "CAWA"

Stage <- 6 # which(names(mods) == "Clim")
# 2001, 2005, 2009, 2013
BASE_YEAR <- 2012#2002#2012
B_use <- 240#100#240
bfill <- FALSE

e <- new.env()
load(file.path("e:/peter/bam/Apr2016/out", "data",
    paste0("pack_", Date, ".Rdata")), envir=e)
mods <- e$mods

if (FALSE) {
load(file.path("e:/peter/bam/Apr2016/out", "data", "pack_2016-12-01.Rdata"))

DAT$xX <- (DAT$Xcl - 399300) / 1419050
DAT$xY <- (DAT$Ycl - 1367000) / 693666
DAT$xX2 <- DAT$xX^2
DAT$xY2 <- DAT$xY^2

mods$Clim <- list(
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMIJJA + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMI2 +
        CMI:DD0 + CMI:DD5 + EMT:MSP,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + CMIJJA + DD0 + MSP + TD + DD02 + CMI2 + CMIJJA2 +
        CMI:DD0 + CMIJJA:DD0 + MSP:TD,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + CMIJJA + DD5 + MSP + TD + DD52 + CMI2 + CMIJJA2 +
        CMI:DD5 + CMIJJA:DD5 + MSP:TD,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMIJJA + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + MSP:TD + MSP:EMT,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMI2 +
        CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + CMIJJA + DD0 + MSP + TD + EMT + DD02 + CMI2 + CMIJJA2 +
        CMI:DD0 + CMIJJA:DD0 + MSP:TD + MSP:EMT,
    . ~ . + xX + xY + xX2 + xY2 + xX:xY +
        CMI + CMIJJA + DD5 + MSP + TD + EMT + DD52 + CMI2 + CMIJJA2 +
        CMI:DD5 + CMIJJA:DD5 + MSP:TD + MSP:EMT)

DAT$EMT2 <- DAT$EMT^2
DAT$MSP2 <- DAT$MSP^2
DAT$TD2 <- DAT$TD^2

mods$Clim <- list(
    . ~ . + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD +
        CMI2 + CMIJJA2 + DD02 + DD52 + EMT2 + MSP2 + TD2 +
        CMI:CMIJJA + CMI:DD0 + CMI:DD5 + CMI:EMT + CMI:MSP + CMI:TD +
        CMIJJA:DD0 + CMIJJA:DD5 + CMIJJA:EMT + CMIJJA:MSP + CMIJJA:TD +
        DD0:DD5 + DD0:EMT + DD0:MSP + DD0:TD +
        DD5:EMT + DD5:MSP + DD5:TD +
        EMT:MSP + EMT:TD +
        MSP:TD)

save(list=c("DAT","YY","mods","TAX","OFF","BB"),
    file=file.path("e:/peter/bam/Apr2016/out", "data", "pack_2017-01-27.Rdata"))
}

#mods$Hgt <- NULL

Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
#xn$CMI2 <- xn$CMI^2
#xn$CMIJJA2 <- xn$CMIJJA^2
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs0 <- levels(st$LEVEL3)
inout <- read.csv("~/repos/bamanalytics/lookup/EcoregionsMaskYesNo.csv")
regs <- unique(as.character(inout[inout$Mask=="No","LEVEL3"])) # 8.1.10 ???
regs <- regs[regs != "Water"]
regs <- regs[regs != "8.1.10"]
#regs <- c("2.1.9", "2.2.1", "2.2.2", "2.3.1", "2.4.1", "2.4.2", "2.4.4",
#    "3.1.1", "3.1.3", "3.2.1", "3.2.3", "3.3.1", "3.4.5", "6.1.1")

## includes land cover, height, and fire
hab_col <- c("HABAgr", "HABBarren", "HABDecid", "HABDevel",
    "HABGrass", "HABMixed", "HABShrub", "HABWet", "HABTRAgr", "HABTRBarren",
    "HABTRConifOpen", "HABTRConifSparse", "HABTRDecidDense", "HABTRDecidOpen",
    "HABTRDecidSparse", "HABTRDevel", "HABTRGrass", "HABTRMixedDense",
    "HABTRMixedOpen", "HABTRMixedSparse", "HABTRShrub", "HABTRWetDense",
    "HABTRWetOpen", "HABTRWetSparse",
    "HGT", "HGT2", "HGT05",
    "DTB", "BRN", "LSS", "YSD", "YSF", "YSL",
    "HGT:isDM", "HGT:isWet", "HGT:isDec", "HGT:isMix",
    "HGT2:isDM", "HGT2:isWet", "HGT2:isDec", "HGT2:isMix", "HGT05:isDM",
    "HGT05:isWet", "HGT05:isDec", "HGT05:isMix")

## spp ---------------------------------

spp <- "CAWA"

fn <- file.path("e:/peter/bam/Apr2016/out", "results", tolower(spp),
    paste0(PROJECT, "_", spp, "_", Date, ".Rdata"))
load(fn)
100 * sum(getOK(res)) / length(res)
est <- getEst(res, stage = Stage, X=Xn)

B_use <- min(B_use, nrow(est))

#regi <- "10.1.1"
#regi <- "8.1.10"
for (regi in regs) {

cat(spp, regi, BASE_YEAR, ifelse(bfill, "bfill", "no-bfill"), "\n");flush.console()

load(file.path(ROOT3, "chunks", paste0("pgdat-", regi, ".Rdata")))
gc()

dat$xX <- (dat$POINT_X - 399300) / 1419050
dat$xY <- (dat$POINT_Y - 1367000) / 693666
dat$xX2 <- dat$xX^2
dat$xY2 <- dat$xY^2

dat$EMT2 <- dat$EMT^2
dat$MSP2 <- dat$MSP^2
dat$TD2 <- dat$TD^2

dat$HAB <- dat$HAB_NALC2
dat$HABTR <- dat$HAB_NALC1
dat$HGT[dat$HAB %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- 0
dat$HGT2 <- dat$HGT^2
dat$HGT05 <- dat$HGT^0.5

dat$ROAD <- 0L

## YR
dat$YR <- BASE_YEAR - 2001

## disturbance
dat$YearFire[is.na(dat$YearFire)] <- BASE_YEAR - 200
dat$YearLoss[is.na(dat$YearLoss)] <- BASE_YEAR - 200
## backfill means no forest loss, only due to fire
if (bfill)
    dat$YearLoss <- BASE_YEAR - 200

## years since fire
dat$YSF <- BASE_YEAR - dat$YearFire
dat$YSF[dat$YSF < 0] <- 200
## years since loss
dat$YSL <- BASE_YEAR - dat$YearLoss
dat$YSL[dat$YSL < 0] <- 200
## years since most recent burn or loss, with backfill option
dat$YSD <- if (bfill)
    dat$YSF else pmin(dat$YSF, dat$YSL)

## cut at 10 yrs
dat$BRN <- ifelse(dat$YSF <= 10, 1L, 0L)
dat$LSS <- ifelse(dat$YSL <= 10, 1L, 0L)
#dat$LSS[dat$YEAR < 2000] <- NA
dat$DTB <- ifelse(dat$YSD <= 10, 1L, 0L)
#dat$DTB[dat$YEAR < 2000] <- NA

## refining years since variables
AGEMAX <- 50
dat$YSD <- pmax(0, 1 - (dat$YSD / AGEMAX))
dat$YSF <- pmax(0, 1 - (dat$YSF / AGEMAX))
dat$YSL <- pmax(0, 1 - (dat$YSL / AGEMAX))
levels(dat$Brandt) <- c(levels(dat$Brandt), "OUT")
dat$Brandt[is.na(dat$Brandt)] <- "OUT"

dat$subreg <- paste(regi, dat$BCR, dat$JURS, dat$Brandt, sep=" + ")

## backfill for all but land cover
if (bfill) {
    dat$POL <- 0
    dat$LIN <- 0
}

## these are not considered in the models, thus NAs not tracked down !!!
dat$LIN <- dat$POL <- NULL

(aa <- data.frame(na=sort(colSums(is.na(dat)))))
dat0 <- dat[rowSums(is.na(dat)) == 0,]
#stopifnot((nrow(dat)-nrow(dat0))/nrow(dat) < 0.05)
Xn0 <- model.matrix(getTerms(mods[1:Stage], "formula"), dat0)
colnames(Xn0) <- fixNames(colnames(Xn0))
NR <- nrow(Xn0)
is_hf <- dat0$HAB %in% c("Devel", "Agr")

mu0 <- matrix(0, NR, B_use)
if (NR > 0) {
    for (j in 1:B_use) {
        if (bfill) {
            ii <- sample(which(!is_hf), sum(is_hf), replace=TRUE)
            Xn0[is_hf,hab_col] <- Xn0[ii,hab_col]
        }
        mu0[,j] <- drop(Xn0 %*% est[j,colnames(Xn0)])
    }
    lam <- lamfun(mu0, tr=0.99)
    rownames(lam) <- rownames(dat0)
    ## km^2 vs ha diff is 100
    lam_total <- groupSums(exp(mu0) * 100, 1, dat0$subreg)
    #rm(mu0, dat0, Xn0)
    attr(lam, "spp") <- spp
    attr(lam, "stage") <- Stage
    attr(lam, "base-year") <- BASE_YEAR
    attr(lam, "level3") <- regi
    attr(lam, "pdrop") <- (nrow(dat)-nrow(dat0))/nrow(dat) # NA values
} else {
    lam <- NULL
    lam_total <- NULL
}
gc()

if (!dir.exists(file.path(ROOT3, "species", spp)))
    dir.create(file.path(ROOT3, "species", spp))
fout <- file.path(ROOT3, "species", spp,
    paste0(spp, "-", Stage, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"),
    regi, "-", Date, ".Rdata"))
save(lam, lam_total, file=fout)
#rm(lam, lam_total)

}


## mapping starts here -----------------------------


library(RColorBrewer)
library(mefa4)
ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/bam/pred-2015"
ROOT3 <- "e:/peter/bam/pred-2016"

if (FALSE) {
## need to loop over chunks and save a proper XY where rownames match etc
## current load has BCR as NA etc
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs <- levels(st$LEVEL3)

load(file.path(ROOT3, "chunks", paste0("pgdat-", regs[1], ".Rdata")))
XY <- dat[,c("POINT_X","POINT_Y","BCR","JURS","LEVEL3","Brandt")]
for (regi in regs[-1]) {
    load(file.path(ROOT3, "chunks", paste0("pgdat-", regi, ".Rdata")))
    tmp <- dat[,c("POINT_X","POINT_Y","BCR","JURS","LEVEL3","Brandt")]
    XY <- rbind(XY, tmp)
    cat(dim(XY), "\n");flush.console()
}
XY$Brandt <- as.factor(XY$Brandt)
levels(XY$Brandt) <- c(levels(XY$Brandt), "OUT")
XY$Brandt[is.na(XY$Brandt)] <- "OUT"
save(XY, file=file.path(ROOT3, "allXY.Rdata"))
}

## pred grid
## pointid is rownames
load(file.path(ROOT3, "allXY.Rdata"))
XY$subreg <- as.factor(paste(XY$LEVEL3, XY$BCR, XY$JURS, XY$Brandt, sep=" + "))
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
#regs <- levels(st$LEVEL3)
inout <- read.csv("~/repos/bamanalytics/lookup/EcoregionsMaskYesNo.csv")
regs <- unique(as.character(inout[inout$Mask=="No","LEVEL3"])) # 8.1.10 ???
regs <- regs[regs != "Water"]
regs <- regs[regs != "8.1.10"]
XY$studyarea <- XY$LEVEL3 %in% regs
#summary(XY)
gc()

## Subregions
XY3 <- nonDuplicated(XY, subreg, TRUE)
XY3$POINT_X <- NULL
XY3$POINT_Y <- NULL
rownames(XY3) <- XY3$subreg

source("~/repos/bamanalytics/R/makingsense_functions.R")

PROJECT <- "bam"
Date <- "2016-12-01"
#Date <- "2017-01-27"

## observations
e <- new.env()
#load(file.path(ROOT, "out", "data", "pack_2016-04-18.Rdata"), envir=e)
load(file.path("e:/peter/bam/Apr2016/out", "data",
    paste0("pack_", Date, ".Rdata")), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
yy <- e$YY
xy_p <- e$DAT[,c("Xcl","Ycl")]
load(file.path("e:/peter/bam/Apr2016", "out", "SS-regions-and-xy.Rdata")) # SS01
SS <- SS01[match(e$DAT$SS, SS01$SS),]
jlt <- read.csv("~/repos/bamanalytics/lookup/jurisdictions.csv")
compare_sets(SS$JURSALPHA, jlt$JURS)
SS$JURS <- reclass(SS$JURSALPHA, jlt[,2:1])
SS$LEVEL3 <- SS$CECLEVEL3
SS$Brandt <- SS$BOREALLOC

SS$subreg <- as.factor(paste(SS$LEVEL3, SS$BCR, SS$JURS, SS$Brandt, sep=" + "))
SS$studyarea <- SS$LEVEL3 %in% regs
SS2 <- nonDuplicated(SS, SS, TRUE)[,c("PCODE","SS","X_CLCC","Y_CLCC","X_GEONAD83",
    "Y_GEONAD83","JURS","LEVEL3","Brandt","subreg","studyarea")]
yyss <- groupSums(yy, 1, SS$SS)
yyss[yyss>0] <- 1
mmss <- Mefa(yyss, SS2)
#mmss <- mmss[samp(mmss)$studyarea,]
summary(samp(mmss[samp(mmss)$studyarea,]))
yyl3 <- groupSums(xtab(mmss), 1, samp(mmss)$subreg)
#yyl3 <- yyl3[match(samp(mmss)$LEVEL3, rownames(yyl3)),]
#allSSbyL3 <- table(rep(1, nrow(mmss)), samp(mmss)$LEVEL3)
allSSbySubreg <- table(rep(1, nrow(mmss)), samp(mmss)$subreg)
rm(e)

## detection map
if (FALSE) {
png(file.path(ROOT3, "maps", paste0(fo, "-det.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
points(XY[XY$studyarea,1:2], col = "tan", pch=".")
points(xy_p[yy[,spp] == 0,c("Xcl","Ycl")], pch=19, cex=0.2, col=1)
points(xy_p[yy[,spp] > 0,c("Xcl","Ycl")], pch=19, cex=0.5, col=2)
par(op)
dev.off()
}


Stage <- 6 # which(names(mods) == "Clim")
BASE_YEAR <- 2012
bfill <- FALSE

ttt <- list()
brr <- list()

spp <- "CAWA"

gc()
fo <- paste0(spp, "-", Stage, "-", BASE_YEAR, "-", Date)
cat(fo, "\n");flush.console()

fl <- paste0(spp, "-", Stage, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), regs, "-", Date, ".Rdata")

is_null <- integer(length(fl))
names(is_null) <- fl
load(file.path(ROOT3, "species", spp, fl[1]))
if (is.null(lam)) {
    is_null[1] <- 1L
} else {
    plam <- lam
    tlam <- lam_total
}
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT3, "species", spp, fn))
    if (is.null(lam)) {
        is_null[fn] <- 1L
    } else {
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, lam_total)
    }
}
dim(plam)
sum(duplicated(rownames(plam)))

if (FALSE) {
q1 <- quantile(plam[,"Mean"], 0.99)
q2 <- quantile(plam[,"Median"], 0.99)
plam2 <- plam[plam[,"Mean"] > q1 | plam[,"Median"] > q2,]
sort(unique(round(plam2[,"Mean"],0)))

plam2 <- plam[plam[,"Mean"] > 0.2 | plam[,"Median"] > 0.2,]
plam3 <- plam[!(plam[,"Mean"] > 0.2 | plam[,"Median"] > 0.2),]

ii <- sample.int(dim(plam3)[1],10^6)
with(data.frame(plam3[ii,]), plot(Mean,SD, ylim=c(0,10),xlim=c(0,10)))

}

## already done in lamfun()
if (TRUE) {
    q <- quantile(plam[,"Mean"], 0.99)
    plam[plam[,"Mean"] > q,"Mean"] <- q

    q <- quantile(plam[,"Median"], 0.99)
    plam[plam[,"Median"] > q,"Median"] <- q
}

stopifnot(all(rownames(plam) == rownames(XY)))
gc()
#compare_sets(rownames(plam), rownames(XY2))

XY2all <- as.matrix(XY[,c("POINT_X","POINT_Y")])
## ideally, this is empty
#XY2miss <- XY[rownames(XY) %notin% rownames(plam),]
#XY2miss <- XY2miss[XY2miss$studyarea,]
#dim(XY2miss)

x <- plam[,"Median"]
save(x, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "median-pred.Rdata")))
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
br <- Lc_quantile(x, probs=probs, type="L")
if (!is.finite(br[length(br)]))
    br[length(br)] <- 1.01* max(x, na.rm=TRUE)
brr[[fo]] <- br
ttt[[fo]] <- tlam

#l <- opticut::lorenz(x)
#br <- opticut:::quantile.lorenz(l, probs)
#br <- quantile(x, seq(0, 1, by=0.2))

## total million males
tlam0 <- tlam
tlam[!is.finite(tlam)] <- 0
for (i in 1:nrow(tlam)) {
    q <- quantile(tlam[i,], 0.99)
    tlam[i, tlam[i,] > q] <- q
}
tlam[!is.finite(tlam)] <- max(tlam[is.finite(tlam)])
summary(colSums(tlam))
fstat <- function(x, level=0.95) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
}
fstat(colSums(tlam/10^6), 0.9)
fstat(colSums(tlam0/10^6), 0.9)
## quick numbers
100*sum(plam[,"Mean"])/10^6
100*sum(plam[,"Median"])/10^6


XY3s <- droplevels(XY3[rownames(tlam),])
colnames(tlam) <- paste0(spp, "_run", 1:ncol(tlam))
XY3s$nSSinSubreg <- 0
for (i in colnames(allSSbySubreg))
    if (i %in% levels(XY3s$subreg))
        XY3s$nSSinSubreg[XY3s$subreg == i] <- allSSbySubreg[1,i]
XY3s$nDETinSubreg <- 0
for (i in colnames(allSSbySubreg))
    if (i %in% rownames(yyl3))
        XY3s$nDETinSubreg[XY3s$subreg == i] <- yyl3[i,spp]
ddd <- data.frame(XY3s, tlam)
rownames(ddd) <- rownames(tlam)

write.csv(ddd, row.names=FALSE, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "totals.csv")))
save(XY3s, tlam, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "totals.Rdata")))


png(file.path(ROOT3, "maps", paste0(fo, "-mean-", BASE_YEAR, ifelse(bfill, "bf", ""), ".png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
#Col <- rev(brewer.pal(6, "RdYlBu"))
zval <- if (length(unique(round(br,10))) < 5)
    as.factor(rep(1, length(x))) else cut(x, breaks=unique(br), include.lowest=TRUE)
Col <- rev(brewer.pal(nlevels(zval), "RdYlBu"))
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
points(xy_p[yy[,spp] > 0,c("Xcl","Ycl")], pch=19, cex=0.1, col=1)
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "mean abundance"))
    #title=paste(spp, "median abundance"))
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-cov.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
TEXT <- paste0(br[-length(br)], "-", br[-1])
TEXT[length(TEXT)] <- paste0(">", br[length(br)-1])
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "darkgrey", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "SD / mean"))
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-sd.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
CoV <- plam[,"SD"] / mean(plam[,"Mean"])
zval <- cut(CoV, breaks=br)
br <- round(br*mean(plam[,"Mean"]), 4)
TEXT <- paste0(br[-length(br)], "-", br[-1], "%")
TEXT[length(TEXT)] <- paste0(">", br[length(br)-1], "%")
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "darkgrey", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "SD"))
par(op)
dev.off()




CAN <- c("ALBERTA", "BRITISH COLUMBIA", "MANITOBA",
    "NEW BRUNSWICK", "NEWFOUNDLAND",
    "NORTHWEST TERRITORIES", "NOVA SCOTIA", "NUNAVUT",
    "ONTARIO", "PRINCE EDWARD ISLAND", "QUEBEC", "SASKATCHEWAN", "YUKON")

## pop size in full study area
fstat(colSums(tlam)/10^6, 0.9)
## pop size in Canada
ss <- XY3s$JURS %in% CAN
fstat(colSums(tlam[ss,])/10^6, 0.9)
## pop size in Brandt boreal
ss <- XY3s$Brandt != "OUT"
fstat(colSums(tlam[ss,])/10^6, 0.9)
## pop size in Canada/Boreal
ss <- XY3s$JURS %in% CAN & XY3s$Brandt != "OUT"
fstat(colSums(tlam[ss,])/10^6, 0.9)

## by state/prov/terr
by_jurs <- data.frame(t(apply(groupSums(tlam/10^6, 1, XY3s$JURS), 1, fstat)))
by_jurs$perc <- by_jurs[,2] * 100 / sum(by_jurs[,2])
by_jurs <- by_jurs[order(by_jurs$perc),]
round(by_jurs, 4)
## by bcr
by_bcr <- data.frame(t(apply(groupSums(tlam/10^6, 1, XY3s$BCR), 1, fstat)))
by_bcr$perc <- by_bcr[,2] * 100 / sum(by_bcr[,2])
by_bcr <- by_bcr[order(by_bcr$perc),]
round(by_bcr, 4)

chfun <- function(Na, Nb, ta, tb) {
    100 * ((Nb/Na)^(1/(tb-ta)) - 1)
}
chfun(5.729, 5.704, 2002, 2012)

## difference maps

e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-bf-median-pred.Rdata", envir=e)
x0 <- e$x
e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2002-2016-12-01-2002-median-pred.Rdata", envir=e)
x1 <- e$x
e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-median-pred.Rdata", envir=e)
x2 <- e$x
x1 <- x1[names(x2)]
x0 <- x0[names(x2)]

## thresholds for TF maps
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")

br0 <- Lc_quantile(x0, probs=probs, type="L")
if (!is.finite(br0[length(br0)]))
    br0[length(br0)] <- 1.01* max(x0, na.rm=TRUE)
br1 <- Lc_quantile(x1, probs=probs, type="L")
if (!is.finite(br1[length(br1)]))
    br1[length(br1)] <- 1.01* max(x1, na.rm=TRUE)
br2 <- Lc_quantile(x2, probs=probs, type="L")
if (!is.finite(br2[length(br2)]))
    br2[length(br2)] <- 1.01* max(x2, na.rm=TRUE)

dd <- data.frame(UpperCutpoint=probs, bf2012=br0, nobf2002=br1, nobf2012=br2)
dd[1,] <- 0
write.csv(dd, row.names=FALSE, file="CAWA-density-breaks.csv")

#pcawa <- data.frame(pointid=names(x2), median_2012=x2, median_2002=x1, median_backfilled=x0)
#write.csv(pcawa, row.names=FALSE, file="w:/bam-cawa/cawa-pred-med-2016-12-13.csv")

#br <- c(0.05, 0.25, 0.5, 0.95, 1/rev(c(0.05, 0.25, 0.5, 0.95)))
br <- c(0.5, 0.95, 1/rev(c(0.5, 0.95)))
d21 <- x2/x1
d20 <- x2/x0
c21 <- cut(d21, c(-1, br, Inf))
c20 <- cut(d20, c(-1, br, Inf))

if (FALSE) {

all(names(d21)==names(d20))
tow <- data.frame(pointid=names(d21),
    diff_2012div2002=d21,
    diff_2012div2012bf=d20,
    index_2012div2002=c21,
    index_2012div2012bf=c20)
tow$index_2012div2002 <- as.integer(tow$index_2012div2002)
tow$index_2012div2012bf <- as.integer(tow$index_2012div2012bf)
write.csv(tow, row.names=FALSE, file="w:/bam-cawa/difference-maps-20170106.csv")
}

hist(d21,freq=F,col="gold",xlab="100 * D_2012 / D_2002",main="")
hist(d20[d20<10],freq=F,col="gold",xlab="100 * D_2012 / D_ref",main="")

Col <- brewer.pal(nlevels(c21), "RdYlGn")

png(file.path(ROOT3, "maps", paste0(fo, "-diff-2002-2012_5bin.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
zval <- as.integer(c21)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "tan", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
TEXT <- paste0(round(br[-length(br)],2), "-", round(br[-1],2), "x")
TEXT <- c(paste0("0-", br[1], "x"), TEXT, paste0(">", br[length(br)], "x"))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff 2002-2012"))
par(op)
dev.off()


png(file.path(ROOT3, "maps", paste0(fo, "-diff-bfill-2012_5bin.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
zval <- as.integer(c20)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "tan", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
TEXT <- paste0(round(br[-length(br)],2), "-", round(br[-1],2), "x")
TEXT <- c(paste0("0-", br[1], "x"), TEXT, paste0(">", br[length(br)], "x"))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff bfill-2012"))
par(op)
dev.off()

## try diffs

df <- x2-x1
br1 <- seq(min(df), 0, length.out=6)
br2 <- seq(0, max(df), length.out=6)
br <- c(br1[c(1,3,5)], br2[c(2,4,6)])
br
zval <- cut(df, br, include.lowest=TRUE, labels=FALSE)
table(zval)
Col <- brewer.pal(max(zval), "RdYlGn")
png(file.path(ROOT3, "maps", paste0("NEW_diff_12-02_", fo, ".png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",type="n",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
points(XY2all[zval==3,c("POINT_X","POINT_Y")], col = Col[3], pch=".")
points(XY2all[zval!=3,c("POINT_X","POINT_Y")], col = Col[zval[zval!=3]], pch=".")
TEXT <- paste0(round(br[-length(br)],3), ", ", round(br[-1],3))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff 2002-2012"))
par(op)
dev.off()

df <- x2-x0
br1 <- seq(min(df), 0, length.out=6)
br2 <- seq(0, max(df), length.out=6)
br <- c(br1[c(1,3,5)], br2[c(2,4,6)])
br
zval <- cut(df, br, include.lowest=TRUE, labels=FALSE)
table(zval)
Col <- brewer.pal(max(zval), "RdYlGn")
png(file.path(ROOT3, "maps", paste0("NEW_diff_12-bf_", fo, ".png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",type="n",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
points(XY2all[zval==3,c("POINT_X","POINT_Y")], col = Col[3], pch=".")
points(XY2all[zval!=3,c("POINT_X","POINT_Y")], col = Col[zval[zval!=3]], pch=".")
TEXT <- paste0(round(br[-length(br)],3), ", ", round(br[-1],3))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff bf-2012"))
par(op)
dev.off()





## --

if (FALSE) {
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
TEXT <- paste0(100*br[-length(br)], "-", 100*br[-1], "%")
TEXT[length(TEXT)] <- paste0(">", 100*br[length(br)-1], "%")
#CoV <- plam[,"SD"] / plam[,"Mean"]
CoV <- plam[,"IQR"] / plam[,"Median"]
zval <- cut(CoV, breaks=br)
if (ids$SEXT[fid] == "nam") {
    plot(XYfull[rownames(plam),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
} else {
    plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
}
points(xy1, pch=19, cex=2)
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    #title=paste(spp, "SD / mean"))
    title=paste(spp, "IQR / median"))
}

tlam <- data.frame(t(tlam))
write.csv(tlam, row.names=FALSE,
    file=file.path(ROOT2, "species", "cawa-nmbca-tabs", paste0("byregion-", fo, ".csv")))
plam <- data.frame(id=rownames(plam), median=plam[,"Median"], cov=CoV)
write.csv(plam, row.names=FALSE,
    file=file.path(ROOT2, "species", "cawa-nmbca-tabs", paste0("bypoint-", fo, ".csv")))


rm(plam)

## brr: list of Lc based breaks
## ttt: list of BCR/prov x B matrices
#save(brr, ttt, file=file.path(ROOT2, "species", "tlam-CAWA.Rdata"))


## summarize density by BCR/JURS/LCC class

dat2 <- NULL
cn <- c("BCR", "JURS", "HAB_NALC1")
for (regi in regs) {
    cat(regi, "\n");flush.console()
    load(file.path(ROOT3, "chunks", paste0("pgdat-", regi, ".Rdata")))
    dat2 <- rbind(dat2, dat[,cn])
}

fl <- paste0(spp, "-", Stage, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), regs, "-", Date, ".Rdata")
DD <- NULL
for (i in 1:length(regs)) {
    cat(regi, "\n");flush.console()
    e <- new.env()
    load(file.path(ROOT3, "species", spp, fl[i]), envir=e)
    DD <- rbind(DD, e$lam[,"Median",drop=FALSE])
}
dat2$D <- DD[rownames(dat2),]
save(dat2, file=file.path(ROOT3, "CAWA-density-with-BCRprov.Rdata"))

## define mean D and p or region (a_k)
dat2$bcrjurs <- interaction(dat2$BCR, dat2$JURS, sep="_", drop=TRUE)

AA <- Xtab(~ bcrjurs + HAB_NALC1, dat2)
NN <- Xtab(D ~ bcrjurs + HAB_NALC1, dat2) * 100 # abundance
DD <- NN / AA
DD[is.na(DD)] <- 0
DD <- as(DD, "dgCMatrix")

#save(AA, NN, DD, file=file.path(ROOT3, "CAWA-AANNDD-by-BCRprov.Rdata"))

## calculate # of BAM/BBS on and off road surven and SS numbers by LCC/region

