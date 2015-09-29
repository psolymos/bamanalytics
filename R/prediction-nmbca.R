library(RColorBrewer)
library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")

PROJECT <- "bam"
spp <- "CAWA"
Date <- "2015-09-02"

## SEXT: "can", "nam" # spatial extent, (canb=canadian boreal ~ eosd)
## TEXT: "gfw", "fre" # temporal extent, gfw=2001-2013, fire=1997-2014
## LCTU: "nlc", "lcc", "eos" # land cover to use
ids <- expand.grid(
    TEXT=c("gfw", "fre"),
    SEXT=c("can", "nam"),
    LCTU=c("nlc", "lcc", "eos"))
ids <- ids[ids$SEXT == "can" | (ids$SEXT == "nam" & ids$LCTU == "nlc"),]
ids <- ids[order(ids$SEXT),]
ids$fn <- with(ids, paste0("bam_", spp, "_", 
    TEXT, "_", SEXT, "_", LCTU, "_", Date, ".Rdata"))
ids$data <- with(ids, paste0("pack_", 
    TEXT, "_", SEXT, "_", LCTU, "_", Date, ".Rdata"))
rownames(ids) <- 1:8

#Stage <- which(names(mods) == "HS")
Stage <- 6 # which(names(mods) == "Dist")
# 2001, 2005, 2009, 2013
BASE_YEAR <- 2012
fid <- 1

for (fid in 1:8) {

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)
est <- getEst(res, stage = Stage, X=Xn)

if (ids$SEXT[fid] == "can")
    regs <- gsub(".Rdata", "",
        gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))
if (ids$SEXT[fid] == "nam")
    regs <- gsub(".Rdata", "",
        gsub("pgdat-full-", "", list.files(file.path(ROOT2, "chunks2"))))


for (BASE_YEAR in c(2012, 2002)) { # consistent with BBS

#regi <- "6_AB"
for (regi in regs) {

cat(spp, 
    paste0(ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date), 
    "Stg:", Stage, 
    "Bs yr:", BASE_YEAR, 
    "Reg:", regi, "\n")
flush.console()

if (ids$SEXT[fid] == "can")
    load(file.path(ROOT2, "chunks", paste0("pgdat-", regi, ".Rdata")))
if (ids$SEXT[fid] == "nam")
    load(file.path(ROOT2, "chunks2", paste0("pgdat-full-", regi, ".Rdata")))
gc()

## placeholders: HSH, HSH2, isDM, isNF
#dat$HSH <- 0
#dat$HSH2 <- 0
dat$HAB <- 0
dat$isDM <- 0
dat$isNF <- 0
dat$ROAD <- 0L
dat$ARU <- 0L

## YR
if (ids$TEXT[fid] == "gfw")
    dat$YR <- (BASE_YEAR - 2001) / 10
if (ids$TEXT[fid] == "fre")
    dat$YR <- (BASE_YEAR - 1997) / 10

## disturbance
dat$YearFire[is.na(dat$YearFire)] <- BASE_YEAR - 200
dat$YearLoss[is.na(dat$YearLoss)] <- BASE_YEAR - 200

## years since fire 
dat$YSF <- BASE_YEAR - dat$YearFire
dat$YSF[dat$YSF < 0] <- 200
## years since loss
dat$YSL <- BASE_YEAR - dat$YearLoss
dat$YSL[dat$YSL < 0] <- 200
## years since most recent burn or loss
dat$YSD <- pmin(dat$YSF, dat$YSL)

## cut at 10 yrs
dat$BRN <- ifelse(dat$YSF <= 10, 1L, 0L)
dat$LSS <- ifelse(dat$YSL <= 10, 1L, 0L)
dat$LSS[dat$YEAR < 2000] <- NA
dat$DTB <- ifelse(dat$YSD <= 10, 1L, 0L)
dat$DTB[dat$YEAR < 2000] <- NA

## refining years since variables
AGEMAX <- 50
dat$YSD <- pmax(0, 1 - (dat$YSD / AGEMAX))
dat$YSF <- pmax(0, 1 - (dat$YSF / AGEMAX))
dat$YSL <- pmax(0, 1 - (dat$YSL / AGEMAX))

AA <- nrow(dat)
if (ids$SEXT[fid] == "nam") {
    dat$HAB_EOSD2 <- NULL
    dat$HAB_LCC2 <- NULL
}

dat0 <- dat[rowSums(is.na(dat)) == 0,]
#aa <- nrow(dat0) / AA
#aa <- 1
if (ids$LCTU[fid] == "nlc") {
    dat0$HAB <- dat0$HAB_NALC
    dat0$isDM <- dat0$isDM_NALC
    dat0$isNF <- dat0$isNF_NALC
}
if (ids$LCTU[fid] == "eos") {
    dat0$HAB <- dat0$HAB_EOSD
    dat0$isDM <- dat0$isDM_EOSD
    dat0$isNF <- dat0$isNF_EOSD
}
if (ids$LCTU[fid] == "lcc") {
    dat0$HAB <- dat0$HAB_LCC
    dat0$isDM <- dat0$isDM_LCC
    dat0$isNF <- dat0$isNF_LCC
}
Xn0 <- model.matrix(getTerms(mods[1:Stage], "formula"), dat0)
colnames(Xn0) <- fixNames(colnames(Xn0))
#rm(dat0)
NR <- nrow(Xn0)

mu0 <- matrix(0, NR, 240)
if (NR > 0) {
    for (j in 1:nrow(est)) {
        mu0[,j] <- drop(Xn0 %*% est[j,colnames(Xn0)])
    }
    lam <- lamfun(mu0)
    rownames(lam) <- rownames(dat0)
    rm(mu0, dat0, Xn0)
    attr(lam, "spp") <- spp
    attr(lam, "stage") <- Stage
    attr(lam, "base-year") <- BASE_YEAR
    attr(lam, "bcr-jurs") <- regi
} else {
    lam <- NULL
}
gc()

fout <- file.path(ROOT2, "species", paste0(tolower(spp), "-nmbca"), 
    paste0(spp, "-", Stage, "-", BASE_YEAR, "-", regi, "-",
        ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date, ".Rdata"))
save(lam, file=fout)
rm(lam)

}


}
}


## mapping starts here -----------------------------

library(RColorBrewer)
library(mefa4)
library(pbapply)

ROOT <- "c:/Users/Peter/bam"
ROOT2 <- "c:/Users/Peter/bam/pred-2015"
#ROOT <- "c:/bam/May2015"
#ROOT2 <- "e:/peter/bam/pred-2015"

source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")
load(file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))
load(file.path(ROOT2, "XYeosd.Rdata"))


PROJECT <- "bam"
spp <- "CAWA"
Date <- "2015-09-02"

## SEXT: "can", "nam" # spatial extent, (canb=canadian boreal ~ eosd)
## TEXT: "gfw", "fre" # temporal extent, gfw=2001-2013, fire=1997-2014
## LCTU: "nlc", "lcc", "eos" # land cover to use
ids <- expand.grid(
    TEXT=c("gfw", "fre"),
    SEXT=c("can", "nam"),
    LCTU=c("nlc", "lcc", "eos"))
ids <- ids[ids$SEXT == "can" | (ids$SEXT == "nam" & ids$LCTU == "nlc"),]
ids <- ids[order(ids$SEXT),]
ids$fn <- with(ids, paste0("bam_", spp, "_", 
    TEXT, "_", SEXT, "_", LCTU, "_", Date, ".Rdata"))
ids$data <- with(ids, paste0("pack_", 
    TEXT, "_", SEXT, "_", LCTU, "_", Date, ".Rdata"))
rownames(ids) <- 1:8

fid <- 1
#Stage <- which(names(mods) == "HS")
Stage <- 6 # which(names(mods) == "Dist")
# 2001, 2005, 2009, 2013
BASE_YEAR <- 2013

xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]

ttt <- list()

for (fid in 1:6) {
#for (BASE_YEAR in c(2003, 2013)) {

cat(fid, BASE_YEAR, "\n");flush.console()

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)
est <- getEst(res, stage = Stage, X=Xn)

if (ids$SEXT[fid] == "can")
    regs <- gsub(".Rdata", "",
        gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))
if (ids$SEXT[fid] == "nam")
    regs <- gsub(".Rdata", "",
        gsub("pgdat-full-", "", list.files(file.path(ROOT2, "chunks2"))))


fl <- paste0(spp, "-", Stage, "-", BASE_YEAR, "-", regs, "-",
        ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date, ".Rdata")

load(file.path(ROOT2, "species", paste0(tolower(spp), "-nmbca2"), fl[1]))
plam <- lam
tlam <- attr(lam, "total")
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT2, "species", paste0(tolower(spp), "-nmbca2"), fn))
    plam <- rbind(plam, lam)
    tlam <- rbind(tlam, attr(lam, "total"))
}
rownames(tlam) <- regs
dim(plam)
sum(duplicated(rownames(plam)))

ttt[[fid]] <- tlam

fo <- paste0(spp, "-", Stage, "-", BASE_YEAR, "-", 
        ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date)

#write.csv(tlam, file=file.path(ROOT, "out", "figs", "nmbca2", paste0(fo, ".csv")))

png(file.path(ROOT, "out", "figs", "nmbca2", paste0(fo, ".png")), 
    width = 2000, height = 2000)
op <- par(mfrow=c(2,1), mar=c(1,1,1,1)+0.1)

#x <- plam[,"Mean"]
x <- plam[,"Median"]
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))
br <- Lc_quantile(x, probs=probs, type="L")
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(x)) else as.integer(cut(x, breaks=br))
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)
legend("topright", bty = "n", legend=rev(TEXT), 
    fill=rev(Col), border=1, cex=3, 
    #title=paste(spp, "mean abundance"))
    title=paste(spp, "median abundance"))

br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
TEXT <- paste0(100*br[-length(br)], "-", 100*br[-1], "%")
TEXT[length(TEXT)] <- paste0(">", 100*br[length(br)-1], "%") 
#CoV <- plam[,"SD"] / plam[,"Mean"]
CoV <- plam[,"IQR"] / plam[,"Median"]
zval <- cut(CoV, breaks=br)
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)
legend("topright", bty = "n", legend=rev(TEXT), 
    fill=rev(Col), border=1, cex=3, 
    #title=paste(spp, "SD / mean"))
    title=paste(spp, "IQR / median"))

par(op)
dev.off()

#}
}

pe <- data.frame(ids[1:6,1:3],
    Year=BASE_YEAR,
    t(sapply(ttt, function(z) 
        c(mean=mean(colSums(z)/10^6), median=median(colSums(z)/10^6)))))
save(ttt, file=file.path(ROOT, "out", "figs", "nmbca2", 
    paste0(paste0("popsize-", spp, "-", Stage, "-", BASE_YEAR, "-", Date), ".Rdata")))

## 2003

pe03 <- structure(list(TEXT = structure(c(1L, 2L, 1L, 2L, 1L, 2L), .Label = c("gfw", 
    "fre"), class = "factor"), SEXT = structure(c(1L, 1L, 1L, 1L, 
    1L, 1L), .Label = c("can", "nam"), class = "factor"), LCTU = structure(c(1L, 
    1L, 2L, 2L, 3L, 3L), .Label = c("nlc", "lcc", "eos"), class = "factor"), 
        Year = c(2003, 2003, 2003, 2003, 2003, 2003), mean = c(7.04675889420143, 
        8.1948604617475, 7.30885952161788, 8.12954353917993, 7.32127463952691, 
        8.04625186575813), median = c(7.03908748814014, 7.90960556092234, 
        7.28424257344771, 7.97434394721777, 7.29854788856675, 7.70893150271355
        )), .Names = c("TEXT", "SEXT", "LCTU", "Year", "mean", "median"
    ), row.names = c(NA, 6L), class = "data.frame")

pe13 <- structure(list(TEXT = structure(c(1L, 2L, 1L, 2L, 1L, 2L), .Label = c("gfw", 
    "fre"), class = "factor"), SEXT = structure(c(1L, 1L, 1L, 1L, 
    1L, 1L), .Label = c("can", "nam"), class = "factor"), LCTU = structure(c(1L, 
    1L, 2L, 2L, 3L, 3L), .Label = c("nlc", "lcc", "eos"), class = "factor"), 
        Year = c(2013, 2013, 2013, 2013, 2013, 2013), mean = c(7.03823934905092, 
        8.2201686110189, 7.29573582539286, 8.15710339698905, 7.2890423740622, 
        8.0574300369338), median = c(7.02512525548913, 7.93361330635701, 
        7.27357416838527, 7.98576510700306, 7.2667175762312, 7.70826606796275
        )), .Names = c("TEXT", "SEXT", "LCTU", "Year", "mean", "median"
    ), row.names = c(NA, 6L), class = "data.frame")






## -- old

spp <- "CAWA"
wid <- 1
Stage <- 6
BASE_YEAR <- 2015
regs <- gsub(".Rdata", "",
    gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))

xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]

fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
load(file.path(ROOT2, "species", paste0(spp, "-ver3"), fl[1]))
plam <- lam
tlam <- attr(lam, "total")
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT2, "species", paste0(spp, "-ver3"), fn))
    plam <- rbind(plam, lam)
    tlam <- rbind(tlam, attr(lam, "total"))
}
dim(plam)
sum(duplicated(rownames(plam)))

hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)
mean(colSums(tlam)/10^6)


plotfun <- function(plam, what="Median") {
    z <- plam[,what]
    #z <- plam[,"IQR"]/plam[,"Median"]
    Col <- rev(brewer.pal(10, "RdBu"))
    cz <- cut(z, breaks = c(min(z)-1, quantile(z, seq(0.1, 1, 0.1))))
    plot(XYeosd[rownames(plam),], col = Col[cz], pch=".")
    invisible(NULL)
}

plotfun(z)

z <- dat1$HAB
Col <- brewer.pal(9, "Set3")
plot(dat[,2:3], col = Col[z], pch=".")

hs <- rowSums(pg4x4$nalc[rownames(dat0), c("Decid", "Mixed")])
plot(hs, aa)
table(unlist(lapply(sres$nalc, "[[", "hi")))







png(file.path(ROOT, "out", "figs", "CAWA-ver3-2015.png"), width = 2000, height = 2000)
op <- par(mfrow=c(2,1), mar=c(1,1,1,1)+0.1)

library(RColorBrewer)
x <- plam[,"Mean"]
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))
br <- Lc_quantile(x, probs=probs, type="L")
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(x)) else as.integer(cut(x, breaks=br))
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)


br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br)
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)

par(op)
dev.off()

yrs <- c(2001, 2005, 2009, 2013)
gfw <- list()
for (BASE_YEAR in yrs) {
    fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
    load(file.path(ROOT2, "species", spp, fl[1]))
    plam <- lam
    tlam <- attr(lam, "total")
    for (fn in fl[-1]) {
        cat("loading", fn, "\n");flush.console()
        load(file.path(ROOT2, "species", spp, fn))
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, attr(lam, "total"))
    }
    gfw[[as.character(BASE_YEAR)]] <- tlam
}

hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)

ps <- sapply(gfw, function(z) {
    tmp <- colSums(z)/10^6
    c(Mean=mean(tmp), quantile(tmp, c(0.025, 0.975)))
    })

ps <- rbind(ps, Percent=(100 * ps[1,] / ps[1,1]) - 100)
ps <- rbind(ps, Decadal=10*ps[4,]/c(yrs-yrs[1]))


#spp <- "CAWA"

getDistPred <- function(spp) {
    xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]
    fd <- function(i) {
        xy <- XYeosd[i, ]
        min(sqrt((xy[1] - xy1[,1])^2 + (xy[2] - xy1[,2])^2)) / 1000
    }
    pbsapply(seq_len(nrow(XYeosd)), fd)
}

dp_cawa <- getDistPred("CAWA")
save(dp_cawa, file=file.path(ROOT2, "dp-cawa.Rdata"))


png(file.path(ROOT, "out", "figs", "CAWA-d-2015.png"), width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)

col <- rev(brewer.pal(9, "Reds"))
z <- cut(dp_cawa, c(-1, 0, 5, 10, 50, 100, 200, 500, 1000, Inf))
plot(XYeosd, pch=".", col=col[z],
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)

par(op)
dev.off()


spp <- "CAWA"
Stage <- 7
BASE_YEAR <- 2015
est <- list()
map <- list()
for (wid in c(0,2,3)) {
    fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
    load(file.path(ROOT2, "species", spp, fl[1]))
    plam <- lam
    tlam <- attr(lam, "total")
    for (fn in fl[-1]) {
        cat("loading", fn, "\n");flush.console()
        load(file.path(ROOT2, "species", spp, fn))
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, attr(lam, "total"))
    }
    est[[as.character(wid)]] <- tlam
    map[[as.character(wid)]] <- plam
}

ps <- sapply(est, function(z) {
    tmp <- colSums(z)/10^6
    c(Mean=mean(tmp), quantile(tmp, c(0.025, 0.975)))
    })
colnames(ps) <- c("NALC","EOSD","LCC")
ps

png(file.path(ROOT, "out", "figs", "CAWA-all-7-2015.png"), width = 2000, height = 3000)
op <- par(mfrow=c(3,1), mar=c(1,1,1,1)+0.1)

for (i in 1:3) {
    x <- map[[i]][,"Mean"]
    probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
    TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
    Col <- rev(brewer.pal(5, "RdYlBu"))
    br <- Lc_quantile(x, probs=probs, type="L")
    zval <- if (length(unique(round(br,10))) < 5)
        rep(1, length(x)) else as.integer(cut(x, breaks=br))
    plot(XYeosd[rownames(map[[i]]),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
    points(xy1, pch=19, cex=1.2)
}

par(op)
dev.off()


library(RColorBrewer)
ROOT <- "c:/bam/May2015"
load("e:/peter/bam/pred-2015/pg-clim.Rdata")
load("e:/peter/bam/pred-2015/pg-loss.Rdata")
clim$YearFire <- loss$YearFire[match(clim$pointid, loss$pointid)]
clim$YearLoss <- loss$YearLoss[match(clim$pointid, loss$pointid)]

for (fn in c("CMIJJA", "CMI", "TD", "DD0", 
    "DD5", "EMT", "MSP", "CTI", "SLP")) {

    png(file.path(ROOT, "out", "figs", paste0("x-", fn,".png")), width = 2000, height = 1000)
    op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)

    Col <- brewer.pal(5, "OrRd")
    z <- cut(clim[,fn], breaks=quantile(clim[,fn], seq(0,1,0.2), na.rm=TRUE))
    plot(clim[,c("POINT_X","POINT_Y")], col = Col[z], pch=".",
        ann=FALSE, axes=FALSE)
    title(main=fn)
    legend("bottomleft", bty = "n", legend=rev(levels(z)), fill=rev(Col))
    

    par(op)
    dev.off()

}


## mapping, fid = 4

library(RColorBrewer)
library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")
load(file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))
load(file.path(ROOT2, "XYfull.Rdata"))

spp <- "CAWA"
wid <- 1
Stage <- 6
BASE_YEAR <- 2015
regs <- gsub(".Rdata", "",
    gsub("pgdat-full-", "", list.files(file.path(ROOT2, "chunks2"))))

xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]

fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
load(file.path(ROOT2, "species", paste0(spp, "-4"), fl[1]))
plam <- lam
tlam <- attr(lam, "total")
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT2, "species", paste0(spp, "-4"), fn))
    plam <- rbind(plam, lam)
    tlam <- rbind(tlam, attr(lam, "total"))
}
dim(plam)
sum(duplicated(rownames(plam)))


hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)
mean(colSums(tlam)/10^6)

z <- ifelse(rownames(XYfull) %in% names(Brandt), 2, 1)
plot(XYfull, col = z, pch=".", ann=FALSE, axes=FALSE)

z <- plam[,"Median"]
z[z >= 100] <- max(z[z < 100])
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))
br <- Lc_quantile(z, probs=probs, type="L")

png(file.path(ROOT, "out", "figs", "CAWAfid4-2015.png"), width = 2000, height = 2000)
op <- par(mfrow=c(2,1), mar=c(1,1,1,1)+0.1)

zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(z)) else as.integer(cut(z, breaks=br))
plot(XYfull[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)


br2 <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col2 <- rev(brewer.pal(10, "RdYlGn"))
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br2)
plot(XYfull[rownames(plam),], col = Col2[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)

par(op)
dev.off()

## marginal plots

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)

load(file.path(ROOT, "pg-main.Rdata"))
x <- x[x$EOSD_COVER == 1,]
rownames(x) <- x$pointid

x <- x[rownames(plam),]
load(file.path(ROOT, "pg-loss.Rdata"))
ii <- loss$YearFire >= 9000 & !is.na(loss$YearFire)
loss$YearFire[ii] <- loss$YearFire[ii] - 8000
x$YearFire <- loss$YearFire[match(x$pointid, loss$pointid)]
x$YearLoss <- loss$YearLoss[match(x$pointid, loss$pointid)]

i <- sample.int(nrow(x), 5000)
tc <- c(rgb(1,0,0,alpha=0.2), rgb(0,0,1,alpha=0.2), rgb(0,1,0,alpha=0.2), rgb(0,0,0,alpha=0.2))
j <- as.integer(x$HAB_NALC2[i])
j[] <- 4
j[x$HAB_NALC2[i] == "Decid"] <- 1
j[x$HAB_NALC2[i] == "Mixed"] <- 2
j[x$HAB_NALC2[i] == "Conif"] <- 3

boxplot(plam[i,"Mean"] ~ x$HAB_NALC2[i], col="gold", range=0, main="Habitat", ylab="D")
boxplot(plam[i,"Mean"] ~ x$TR3[i], col="gold", range=0, main="Tree", ylab="D")
plot(plam[i,"Mean"] ~ jitter(x$HGT[i]), pch=19, cex=1, col=tc[j], main="Height", ylab="D")
legend("topleft", pch=19, col=tc, legend=c("Dec","Mix","Con","Else"))

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(x$LIN[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(x$POL[i]), pch=19, cex=1, col=tc[j], ylab="D")


load(file.path(ROOT, "pg-clim.Rdata"))
rownames(clim) <- clim$pointid
clim <- clim[match(x$pointid, clim$pointid),4:14]
clim <- clim[rownames(plam),]

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(clim$CTI[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(clim$SLP[i]), pch=19, cex=1, col=tc[j], ylab="D")

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(2015-x$YearFire[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(2015-x$YearLoss[i]), pch=19, cex=1, col=tc[j], ylab="D")
