library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/Apr2016/out"
source("~/repos/bamanalytics/R/makingsense_functions.R")

PROJECT <- "bam"
Date <- "2016-12-01"
level <- 0.9

e <- new.env()
load(file.path(ROOT2, "data", "pack_2016-12-01.Rdata"), envir=e)

mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
#xn <- e$DAT[,Terms]
xn <- e$DAT
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
#xn <- xn[rownames(Xn),]
off <- e$OFF#[rownames(xn),]
yy <- e$YY#[rownames(xn),]
#bbb <- unique(e$BB)
bb <- e$BB
rm(e)
#ss1 <- which(!(1:nrow(xn) %in% bbb)) # test portion
#ss2 <- which((1:nrow(xn) %in% bbb)) # non-test portion

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

fn <- file.path(ROOT2, "results", "cawa",
    paste0(PROJECT, "_", spp, "_", Date, ".Rdata"))
load(fn)
est <- getEst(res, stage = length(mods)-1, X=Xn)
est_yr <- getEst(res, stage = length(mods), X=Xn)

load(file.path(ROOT2, "SS-regions-and-xy.Rdata"))
SS01$COUNTRY <- as.factor(SS01$COUNTRY)
SS01 <- droplevels(SS01[rownames(SS01) %in% as.character(xn$SS),])
SS01$BCR[SS01$BCR==0] <- NA

## fill-in BCR based on closest known value
ii <- which(SS01$BCR == 0 | is.na(SS01$BCR))
for (i in ii) {
    d <- sqrt((SS01$X_GEONAD83 - SS01$X_GEONAD83[i])^2 +
        (SS01$Y_GEONAD83 - SS01$Y_GEONAD83[i])^2)
    d[ii] <- Inf
    SS01$BCR[i] <- SS01$BCR[which.min(d)]
}
ii <- which(is.na(SS01$CECLEVEL3) | SS01$CECLEVEL3 == "Water")
for (i in ii) {
    d <- sqrt((SS01$X_GEONAD83 - SS01$X_GEONAD83[i])^2 +
        (SS01$Y_GEONAD83 - SS01$Y_GEONAD83[i])^2)
    d[ii] <- Inf
    SS01$CECLEVEL3[i] <- SS01$CECLEVEL3[which.min(d)]
}

col_keep <- colSums(abs(est) > 0) != 0
pr <- exp(sapply(1:nrow(est), function(j)
    Xn[,colnames(est[,col_keep,drop=FALSE]),drop=FALSE] %*%
    est[j,col_keep]))

DAT <- droplevels(SS01[match(xn$SS, rownames(SS01)),
    c("PCODE", "SS", "X_GEONAD83", "X_CLCC", "Y_CLCC", "X_GEONAD83", "Y_GEONAD83",
    "JURSALPHA", "BOREALLOC", "BCR", "CECLEVEL3", "COUNTRY")])
rownames(DAT) <- rownames(xn)
DAT$Y <- Y
DAT$off <- off1
DAT$YR <- xn$YR
DAT$D <- rowMeans(pr)
DAT$EW <- ifelse(DAT$X_GEONAD83 < (-97), "W", "E")
DAT$BCRPROV <- interaction(DAT$BCR, DAT$JURSALPHA, drop=TRUE, sep="_")
BBS_PCODE <- levels(DAT$PCODE)[substr(levels(DAT$PCODE), 1, 3) == "BBS"]
DAT$isBBS <- DAT$PCODE %in% BBS_PCODE
DAT$ROAD <- xn$ROAD[match(rownames(DAT), rownames(xn))]

## Boreal year effect estimates

fstat <- function(x, level=0.95, digits=3) {
    round(c(Mean=mean(x, na.rm=TRUE),
    Median=median(x, na.rm=TRUE),
    quantile(x, c((1-level)/2, 1 - (1-level)/2), na.rm=TRUE)), digits)
}
hist(100 * (exp(est_yr[,"YR"]) - 1), col="grey",
     main="", xlab="% annual population change")
round(c(fstat(100 * (exp(est_yr[,"YR"]) - 1)),
    summary(100 * (exp(est_yr[,"YR"]) - 1))), 3)

## Residual trend estimates

yr_fun <- function(i, subset=NULL, part=c("all", "bbs", "bam", "off")) {
    part <- match.arg(part)
    if (is.null(subset))
        subset <- rep(TRUE, nrow(DAT))
    dat <- DAT
    dat$SUBSET <- subset
    dat <- dat[bb[,i],]
    dat <- dat[dat$SUBSET,,drop=FALSE]
    if (part=="bbs") # BBS only
        dat <- dat[dat$isBBS,,drop=FALSE]
    if (part=="bam") # non-BBS excluding roadside surveys
        dat <- dat[!dat$isBBS,,drop=FALSE]
    if (part=="off") # non-BBS excluding roadside surveys
        dat <- dat[!dat$isBBS & dat$ROAD==0,,drop=FALSE]
#    if (part=="all") # non-BBS excluding roadside surveys
#        dat <- dat[dat$isBBS | (!dat$isBBS && dat$ROAD==0),,drop=FALSE]
    if (nrow(dat) < 1)
        return(NA)
    dat$logDoff <- log(dat$D) + dat$off
    mod <- glm(Y ~ YR, data=dat, offset=dat$logDoff, family=poisson)
    out <- 100 * (exp(coef(mod)[2]) - 1)
    #attr(out, "nobs") <- nrow(dat)
    #attr(out, "part") <- part
    out
}
dat_fun <- function(i, subset=NULL, part=c("all", "bbs", "bam", "off")) {
    part <- match.arg(part)
    if (is.null(subset))
        subset <- rep(TRUE, nrow(DAT))
    dat <- DAT
    dat$SUBSET <- subset
    dat <- dat[bb[,i],]
    dat <- dat[dat$SUBSET,,drop=FALSE]
    if (part=="bbs") # BBS only
        dat <- dat[dat$isBBS,,drop=FALSE]
    if (part=="bam") # non-BBS excluding roadside surveys
        dat <- dat[!dat$isBBS,,drop=FALSE]
    if (part=="off") # non-BBS excluding roadside surveys
        dat <- dat[!dat$isBBS & dat$ROAD==0,,drop=FALSE]
    c(n=nrow(dat), det=sum(dat$Y>0))
}
#yr_res <- pbsapply(1:240, yr_fun)

## need to find BCR/JURS/Ecoreg/Brandt values
## add subset uption into yr_fun
## run for smaller regions
## and Eeast, West, Boreal, Boreal+Hemi, Canada

## subsets

## broad subsets for ms

tres_all <- list()
for (PART in c("all","bbs","bam","off")) {
cat(PART, "\tFull\n");flush.console();gc()
tres_full <- pbsapply(1:240, yr_fun, subset=NULL, part=PART)
cat(PART, "\tCanada\n");flush.console();gc()
tres_can <- pbsapply(1:240, yr_fun, subset=DAT$COUNTRY == "CAN", part=PART)
cat(PART, "\tBoreal\n");flush.console();gc()
tres_bor <- pbsapply(1:240, yr_fun,
    subset=DAT$BOREALLOC %in% c("B_ALPINE", "BOREAL"), part=PART)
cat(PART, "\tHemiboreal\n");flush.console();gc()
tres_hem <- pbsapply(1:240, yr_fun,
    subset=DAT$BOREALLOC %in% c("H_ALPINE", "HEMIBOREAL"), part=PART)
tres_all[[PART]] <- list(Full=tres_full, Canada=tres_can, Boreal=tres_bor, Hemiboreal=tres_hem)
}
dat_all <- list()
for (PART in c("all","bbs","bam","off")) {
dat_full <- dat_fun(1, subset=NULL, part=PART)
dat_can <- dat_fun(1, subset=DAT$COUNTRY == "CAN", part=PART)
dat_bor <- dat_fun(1,
    subset=DAT$BOREALLOC %in% c("B_ALPINE", "BOREAL"), part=PART)
dat_hem <- dat_fun(1,
    subset=DAT$BOREALLOC %in% c("H_ALPINE", "HEMIBOREAL"), part=PART)
dat_all[[PART]] <- list(Full=dat_full, Canada=dat_can, Boreal=dat_bor, Hemiboreal=dat_hem)
}

for (PART in c("all","bbs","bam","off")) {
cat(PART, "\tEast\n");flush.console();gc()
tres_all[[PART]][["East"]] <- pbsapply(1:240, yr_fun, subset=DAT$EW == "E", part=PART)
cat(PART, "\tWest\n");flush.console();gc()
tres_all[[PART]][["West"]] <-  pbsapply(1:240, yr_fun, subset=DAT$EW == "W", part=PART)
}
for (PART in c("all","bbs","bam","off")) {
dat_all[[PART]][["East"]] <- dat_fun(1, subset=DAT$EW == "E", part=PART)
dat_all[[PART]][["West"]] <- dat_fun(1, subset=DAT$EW == "W", part=PART)
}

save(tres_all, dat_all, file="e:/peter/bam/Apr2016/out/cawa/trend-est-broad.Rdata")

load("e:/peter/bam/Apr2016/out/cawa/trend-est-broad.Rdata")

All <- do.call(rbind, lapply(names(tres_all), function(z) {
    zz <- t(sapply(tres_all[[z]], fstat))
    data.frame(set=z, reg=rownames(zz), zz)
}))
rownames(All) <- NULL

All2 <- do.call(rbind, lapply(names(dat_all), function(z) {
    zz <- t(sapply(dat_all[[z]], function(zzz) zzz))
    data.frame(set=z, reg=rownames(zz), zz)
}))
rownames(All2) <- NULL
All <- cbind(All,All2[,c("n","det")])
All[order(All$reg),]

ii <- DAT$COUNTRY == "CAN"
ii <- DAT$BOREALLOC %in% c("B_ALPINE", "BOREAL")
ii <- DAT$BOREALLOC %in% c("B_ALPINE", "BOREAL", "H_ALPINE", "HEMIBOREAL")
ii <- DAT$BCRPROV == "6_AB" #
## these are the PIF levels for CAWA
ii <- DAT$BCR == 6 # 6 7 8 11 12 13 14 23
ii <- DAT$JURSALPHA == "AB" # AB SK MB QC ON NB NS NT PEI
## these are the BBS levels for CAWA
ii <- DAT$BCR == 6 # 6 8 12 13 14
ii <- DAT$JURSALPHA == "AB" # AB MB QC ON NB NS+PEI

DAT$JURS2 <- DAT$JURSALPHA
levels(DAT$JURS2)[levels(DAT$JURS2) %in% c("NS","PEI")] <- "NS+PEI"
DAT$BCRPROV <- interaction(DAT$BCR, DAT$JURSALPHA, drop=TRUE, sep="_")
DAT$BCRPROV2 <- interaction(DAT$BCR, DAT$JURS2, drop=TRUE, sep="_")

## NS and PEI as separate
tres_nspei <- list()
tres_nspei[["off-NS"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "NS", part="off")
tres_nspei[["all-NS"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "NS", part="all")
tres_nspei[["bbs-NS"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "NS", part="bbs")
tres_nspei[["bam-NS"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "NS", part="bam")
tres_nspei[["off-PEI"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "PEI", part="off")
tres_nspei[["all-PEI"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "PEI", part="all")
tres_nspei[["bbs-PEI"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "PEI", part="bbs")
tres_nspei[["bam-PEI"]] <- pbsapply(1:240, yr_fun, subset=DAT$JURSALPHA == "PEI", part="bam")
sapply(tres_nspei, fstat)
save(tres_nspei,
    file=paste0("e:/peter/bam/Apr2016/out/cawa/trend-est-full-NSandPEIseparate.Rdata"))

## BBS runs -- same subset as official method (A Smith)

rt <- read.csv("e:/peter/bam/cawa-ms/bbs-trend/Canada Warbler prebugsdata.csv")

## map rt$strat.name to PCODE
rt$PCODE2 <- rt$strat.name
#levels(rt$PCODE) <- c("Alberta-BCR6", "MAINE-BCR14", "Manitoba-BCR12", "Manitoba-BCR6",
#    "Manitoba-BCR8", "MICHIGAN-BCR12", "MINNESOTA-BCR12", "New Brunswick-BCR14",
#    "NEW YORK-BCR13", "Nova Scotia Prince Edward Island-BCR14", "Ontario-BCR12",
#    "Ontario-BCR13", "Ontario-BCR8", "Quebec-BCR12", "Quebec-BCR13",
#    "Quebec-BCR14", "Quebec-BCR8", "Saskatchewan-BCR6", "VERMONT-BCR14",
#    "WISCONSIN-BCR12")
levels(rt$PCODE2) <- c("BBSAB", "BBSME", "BBSMB", "BBSMB",
    "BBSMB", "BBSMI", "BBSMN", "BBSNB",
    "BBSNY", "BBSNSPEI", "BBSON",
    "BBSON", "BBSON", "BBSQC", "BBSQC",
    "BBSQC", "BBSQC", "BBSSK", "BBSVT",
    "BBSWI")
tmp <- sapply(strsplit(as.character(rt$route), "-"), function(z) z[2])
rt$BBSroute <- paste0(rt$PCODE2, ":", tmp)
## BBSPEI and BBSNS needs to be treated the same
DAT$PCODE2 <- DAT$PCODE
levels(DAT$PCODE2)[levels(DAT$PCODE2) %in% c("BBSPEI","BBSNS")] <- "BBSNSPEI"

DAT$BBSroute <- as.character("none")
tmp <- sapply(strsplit(as.character(DAT$SS), ":"), function(z) z[2])
DAT$BBSroute[DAT$isBBS] <- paste0(DAT$PCODE2[DAT$isBBS], ":", tmp[DAT$isBBS])

DAT$useBBS <- DAT$BBSroute %in% rt$BBSroute
table(BBS=DAT$isBBS, Subset=DAT$useBBS)
table(BBS=DAT$Y, Subset=DAT$useBBS)

compare_sets(DAT$BBSroute, rt$BBSroute)

tres_usebbs <- list()
tres_usebbs[["All"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS, part="bbs")
tres_usebbs[["Can"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$COUNTRY == "CAN", part="bbs")
tres_usebbs[["US"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$COUNTRY != "CAN", part="bbs")
tres_usebbs[["BCR6"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$BCR == 6, part="bbs")
tres_usebbs[["BCR14"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$BCR == 14, part="bbs")
tres_usebbs[["BCR12"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$BCR == 12, part="bbs")
tres_usebbs[["BCR8"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$BCR == 8, part="bbs")
tres_usebbs[["BCR13"]] <- pbsapply(1:240, yr_fun, subset=DAT$useBBS & DAT$BCR == 13, part="bbs")

fstat <- function(x, level=0.95, digits=3) {
    round(c(Mean=mean(x, na.rm=TRUE),
    Median=median(x, na.rm=TRUE),
    quantile(x, c((1-level)/2, 1 - (1-level)/2), na.rm=TRUE)), digits)
}
t(sapply(tres_usebbs, fstat))

tres_usebbs2 <- list()
tres_usebbs2[["All"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS, part="bbs")
tres_usebbs2[["Can"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$COUNTRY == "CAN", part="bbs")
tres_usebbs2[["US"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$COUNTRY != "CAN", part="bbs")
tres_usebbs2[["BCR6"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$BCR == 6, part="bbs")
tres_usebbs2[["BCR14"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$BCR == 14, part="bbs")
tres_usebbs2[["BCR12"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$BCR == 12, part="bbs")
tres_usebbs2[["BCR8"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$BCR == 8, part="bbs")
tres_usebbs2[["BCR13"]] <- pbsapply(1:240, yr_fun, subset=DAT$isBBS & DAT$BCR == 13, part="bbs")

t(sapply(tres_usebbs, fstat))
t(sapply(tres_usebbs2, fstat))

save(tres_usebbs, tres_usebbs2,
    file=paste0("e:/peter/bam/Apr2016/out/cawa/trend-est-BBSsubset.Rdata"))

         Mean  Median    2.5%  97.5%
All     0.147   0.208  -1.987  2.242
US      4.046   4.054  -3.201 11.597
Can    -2.728  -2.723  -4.917 -0.640
BCR6  -14.201 -14.262 -19.191 -9.164
BCR14   6.935   7.073   2.314 11.273
BCR12   2.599   2.593   0.121  5.495
BCR8   -7.468  -7.737 -12.365 -1.213
BCR13  -2.814  -2.981 -11.845  7.339

 0.150030325	-2.383032197	2.60523833
 1.285928305	-2.002833323	5.311612216
-0.255788058	-3.103125137	2.394002052
-0.682768956	-5.342919169	4.321161989
-1.540443373	-4.720538597	1.238955143
 1.467685744	-1.32241603	    4.753110659
-0.443063439	-4.652723716	3.116792189
 0.189393728	-4.562392786	4.523406124


PART <- "off" # "all","bbs","bam","off"

tres_can <- pbsapply(1:240, yr_fun, subset=DAT$COUNTRY == "CAN", part=PART)

tres_prov <- list()
#for (i in c("AB","MB","QC","ON","NB","NS+PEI"))
for (i in levels(DAT$JURS2))
{
    cat(i, "\n");flush.console()
    ii <- DAT$JURS2 == i
    tres_prov[[i]] <- pbsapply(1:240, yr_fun, subset=ii, part=PART)
    print(fstat(tres_prov[[i]]))
}

tres_bcr <- list()
#for (i in c(6, 7, 8, 11, 12, 13, 14, 23))
for (i in unique(DAT$BCR))
{
    cat(i, "\n");flush.console()
    ii <- DAT$BCR == i
    tres_bcr[[as.character(i)]] <- pbsapply(1:240, yr_fun, subset=ii, part=PART)
    print(fstat(tres_bcr[[as.character(i)]]))
}

tres_bcrprov <- list()
#for (i in c("6_AB","6_MB","8_MB","14_NB","14_NS+PEI","8_ON","12_ON","13_ON",
#"8_QC","12_QC","13_QC","14_QC"))
for (i in levels(DAT$BCRPROV2))
{
    cat(i, "\n");flush.console()
    ii <- DAT$BCRPROV2 == i
    tres_bcrprov[[as.character(i)]] <- pbsapply(1:240, yr_fun, subset=ii, part=PART)
    print(fstat(tres_bcrprov[[i]]))
}

PART
fstat(tres_can)
t(sapply(tres_prov, fstat))
t(sapply(tres_bcr, fstat))
t(sapply(tres_bcrprov, fstat))

save(DAT, tres_prov, tres_bcr, tres_can, tres_bcrprov,
    file=paste0("e:/peter/bam/Apr2016/out/cawa/trend-est-full-", PART, ".Rdata"))

## summaries

ebam <- new.env()
ebbs <- new.env()
eall <- new.env()
eoff <- new.env()
load("e:/peter/bam/Apr2016/out/cawa/trend-est-full-bam.Rdata", envir=ebam)
load("e:/peter/bam/Apr2016/out/cawa/trend-est-full-bbs.Rdata", envir=ebbs)
load("e:/peter/bam/Apr2016/out/cawa/trend-est-full-all.Rdata", envir=eall)
load("e:/peter/bam/Apr2016/out/cawa/trend-est-full-off.Rdata", envir=eoff)
DAT <- eoff$DAT
eall$DAT <- NULL
ebam$DAT <- NULL
ebbs$DAT <- NULL
eoff$DAT <- NULL
BBS_PCODE <- levels(DAT$PCODE)[substr(levels(DAT$PCODE), 1, 3) == "BBS"]
DAT$isBBS <- DAT$PCODE %in% BBS_PCODE
DAT$Y01 <- ifelse(DAT$Y>0,1,0)
DAT$ROAD <- DAT$ROAD == 1
fstat <- function(x, level=0.95, digits=3) {
    round(c(Mean=mean(x, na.rm=TRUE),
    Median=median(x, na.rm=TRUE),
    quantile(x, c((1-level)/2, 1 - (1-level)/2), na.rm=TRUE)), digits)
}

tr <- list(all=as.list(eall), bam=as.list(ebam), bbs=as.list(ebbs), off=as.list(eoff))

## Canada
dall <- with(DAT, aggregate(Y01, list(Country=COUNTRY, BCR=BCR, JURS=JURS2,
    BCRJURS=BCRPROV), sum))
dbbs <- with(DAT[DAT$isBBS,], aggregate(Y01, list(Country=COUNTRY, BCR=BCR, JURS=JURS2,
    BCRJURS=BCRPROV), sum))
dbam <- with(DAT[!DAT$isBBS,], aggregate(Y01, list(Country=COUNTRY, BCR=BCR, JURS=JURS2,
    BCRJURS=BCRPROV), sum))
doff <- with(DAT[!DAT$isBBS & !DAT$ROAD,], aggregate(Y01, list(Country=COUNTRY, BCR=BCR, JURS=JURS2,
    BCRJURS=BCRPROV), sum))
d <- dall
d$x <- NULL
d$ndet_all <- dall$x[match(d$BCRJURS, dall$BCRJURS)]
d$ndet_bam <- dbam$x[match(d$BCRJURS, dbam$BCRJURS)]
d$ndet_bam[is.na(d$ndet_bam)] <- 0
d$ndet_bbs <- dbbs$x[match(d$BCRJURS, dbbs$BCRJURS)]
d$ndet_bbs[is.na(d$ndet_bbs)] <- 0
d$ndet_off <- doff$x[match(d$BCRJURS, doff$BCRJURS)]
d$ndet_off[is.na(d$ndet_off)] <- 0

#d <- aggregate(DAT$Y01, list(Country=DAT$COUNTRY, BCR=DAT$BCR, PROV=DAT$JURS2,
#    BCRPROV=DAT$BCRPROV), sum)
#d$ndet <- d$x
#d$prop <- aggregate(DAT$Y01, list(Country=DAT$COUNTRY, BCR=DAT$BCR,
#    PROV=DAT$JURS2, BCRPROV=DAT$BCRPROV), mean)$x

tmpd <- with(DAT[DAT$Y01>0,], table(COUNTRY, isBBS))
tmpn <- with(DAT, table(COUNTRY, isBBS))
tmpp <- tmpd/tmpn
table(DAT$COUNTRY[DAT$Y01>0])
table(DAT$COUNTRY)
table(DAT$COUNTRY[DAT$Y01>0])/table(DAT$COUNTRY)

tr1 <- data.frame(ndet=c(table(DAT$COUNTRY[DAT$Y01>0])[1], tmpd[1,1:2]),
    ntot=c(table(DAT$COUNTRY)[1], tmpn[1,1:2]))
tr1$prop <- tr1$ndet/tr1$ntot
tr1 <- data.frame(t(sapply(tr, function(z) fstat(z$tres_can))), tr1)
round(tr1,3)


tab <- data.frame(Country="Canada", BCR="", PROV="", BCRPROV="",
    Data=names(tr), t(sapply(tr, function(z) fstat(z$tres_can))))

tmp <- lapply(tr, function(z) t(sapply(z$tres_bcr, fstat)))
tmp$all <- data.frame(Country="", BCR=rownames(tmp$all),
    PROV="", BCRPROV="", Data="all", tmp$all)
tmp$bam <- data.frame(Country="", BCR=rownames(tmp$bam),
    PROV="", BCRPROV="", Data="bam", tmp$bam)
tmp$bbs <- data.frame(Country="", BCR=rownames(tmp$bbs),
    PROV="", BCRPROV="", Data="bbs", tmp$bbs)
tmp$off <- data.frame(Country="", BCR=rownames(tmp$off),
    PROV="", BCRPROV="", Data="off", tmp$off)
tab <- rbind(tab, tmp$all, tmp$bam, tmp$bbs, tmp$off)

tmp <- lapply(tr, function(z) t(sapply(z$tres_prov, fstat)))
tmp$all <- data.frame(Country="", BCR="",
    PROV=rownames(tmp$all), BCRPROV="", Data="all", tmp$all)
tmp$bam <- data.frame(Country="", BCR="",
    PROV=rownames(tmp$bam), BCRPROV="", Data="bam", tmp$bam)
tmp$bbs <- data.frame(Country="", BCR="",
    PROV=rownames(tmp$bbs), BCRPROV="", Data="bbs", tmp$bbs)
tmp$off <- data.frame(Country="", BCR="",
    PROV=rownames(tmp$off), BCRPROV="", Data="off", tmp$off)
tab <- rbind(tab, tmp$all, tmp$bam, tmp$bbs, tmp$off)

tmp <- lapply(tr, function(z) t(sapply(z$tres_bcrprov, fstat)))
tmp$all <- data.frame(Country="", BCR="",
    PROV="", BCRPROV=rownames(tmp$all), Data="all", tmp$all)
tmp$bam <- data.frame(Country="", BCR="",
    PROV="", BCRPROV=rownames(tmp$bam), Data="bam", tmp$bam)
tmp$bbs <- data.frame(Country="", BCR="",
    PROV="", BCRPROV=rownames(tmp$bbs), Data="bbs", tmp$bbs)
tmp$off <- data.frame(Country="", BCR="",
    PROV="", BCRPROV=rownames(tmp$off), Data="off", tmp$off)
tab <- rbind(tab, tmp$all, tmp$bam, tmp$bbs, tmp$off)

write.csv(tab, row.names=FALSE, file="e:/peter/bam/Apr2016/out/cawa/cawa-trend-2017-07-18.csv")
write.csv(d, row.names=FALSE, file="e:/peter/bam/Apr2016/out/cawa/cawa-det-2017-07-18.csv")

