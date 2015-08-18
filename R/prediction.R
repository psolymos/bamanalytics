library(RColorBrewer)
library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")

## this is for fid 1:3 (not 4 !)
#mods <- mods_gfw

#fid <- 1
Xnlist <- list()
for (fid in 1:3) {
    fl <- c("analysis_package_gfwfire-nalc-2015-08-17.Rdata",
        "analysis_package_gfwfire-eosd-2015-08-17.Rdata",
        "analysis_package_gfwfire-lcc-2015-08-17.Rdata",
        "analysis_package_fire-nalc-2015-08-17.Rdata")
    e <- new.env()
    load(file.path(ROOT, "out", "data", fl[fid]), envir=e)
    mods <- e$mods

    Terms <- getTerms(mods, "list")
    setdiff(Terms, colnames(e$DAT))
    xn <- e$DAT[1:1000,Terms]
    Xn <- model.matrix(getTerms(mods, "formula"), xn)
    colnames(Xn) <- fixNames(colnames(Xn))
    Xn <- Xn[,colnames(Xn) != "REGSouth:YR"]
    rm(e)
    if (fid == 1)
        Xnlist[["nalc"]] <- Xn
    if (fid == 2)
        Xnlist[["eosd"]] <- Xn
    if (fid == 3)
        Xnlist[["lcc"]] <- Xn
}
rm(xn, Xn, fid)

spp <- "CAWA"
wid <- 1
fid <- 4
#Stage <- which(names(mods) == "HS")
Stage <- which(names(mods) == "Dist")
# 2001, 2005, 2009, 2013
BASE_YEAR <- 2015


## fid 1,2,3
if (fid < 4) {
    regs <- gsub(".Rdata", "",
        gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))

    fn <- paste0("bam-", 1, "_", spp, ".Rdata", sep="")
    load(file.path(ROOT, "out", "results", fn))
    sres <- list(nalc=res)
    fn <- paste0("bam-", 2, "_", spp, ".Rdata", sep="")
    load(file.path(ROOT, "out", "results", fn))
    sres$eosd <- res
    fn <- paste0("bam-", 3, "_", spp, ".Rdata", sep="")
    load(file.path(ROOT, "out", "results", fn))
    sres$lcc <- res
    rm(res)

    if (FALSE) {
    Aic <- list()
    for (Stage in 1:8) {
        Aic[[names(mods)[Stage]]] <- cbind(nalc=getCaic(sres$nalc, stage = Stage),
            eosd=getCaic(sres$eosd, stage = Stage), lcc=getCaic(sres$lcc, stage = Stage))
    }
    lapply(Aic, function(aic) table(colnames(aic)[apply(aic, 1, which.min)]))
    }

    est <- list(nalc=getEst(sres$nalc, stage = Stage, X=Xnlist$nalc),
        eosd=getEst(sres$eosd, stage = Stage, X=Xnlist$eosd),
        lcc=getEst(sres$lcc, stage = Stage, X=Xnlist$lcc))
    aic <- cbind(nalc=getCaic(sres$nalc, stage = Stage),
        eosd=getCaic(sres$eosd, stage = Stage),
        lcc=getCaic(sres$lcc, stage = Stage))
    table(colnames(aic)[apply(aic, 1, which.min)])
}
## fid 4
if (fid == 4) {
    wid <- 1 # NALC only
    regs <- gsub(".Rdata", "",
        gsub("pgdat-full-", "", list.files(file.path(ROOT2, "chunks2"))))

    fn <- paste0("bam-", 4, "_", spp, ".Rdata", sep="")
    load(file.path(ROOT, "out", "results", fn))
    res4 <- res
    rm(res)

    est <- getEst(res4, stage = Stage, X=Xnlist$nalc)
    aic <- getCaic(res4, stage = Stage)
}
######


#for (BASE_YEAR in c(2001, 2005, 2009, 2013)) {
#for (wid in 3:2) {

#regi <- "6_AB"
for (regi in regs) {

cat(spp, "--- Stage:", Stage, "--- setup:", 
    wid, "--- base yr:", BASE_YEAR, "--- region:", regi, "\n")
flush.console()

if (fid < 4) {
    load(file.path(ROOT2, "chunks", paste0("pgdat-", regi, ".Rdata")))
    if (Stage <= 6)
        rm(pg4x4)
}
if (fid == 4)
    load(file.path(ROOT2, "chunks2", paste0("pgdat-full-", regi, ".Rdata")))
gc()

## placeholders: HSH, HSH2, isDM, isNF
dat$HSH <- 0
dat$HSH2 <- 0
dat$HAB <- 0
dat$isDM <- 0
dat$isNF <- 0
dat$ROAD <- 0L
dat$ARU <- 0L

## YR
dat$YR <- BASE_YEAR - 1997

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

dat$BRN <- ifelse(dat$YSF <= 10, 1L, 0L)
dat$LSS <- ifelse(dat$YSL <= 10, 1L, 0L)
dat$LSS[dat$YEAR < 2000] <- NA
dat$DTB <- ifelse(dat$YSD <= 10, 1L, 0L)
dat$DTB[dat$YEAR < 2000] <- NA

AA <- nrow(dat)
lamfun <- function(mu, tr=0.99) {
    lam <- exp(mu)
    q <- quantile(lam, tr)
    lam[lam > q] <- q
    Mean <- rowMeans(lam)
    SD <- apply(lam, 1, sd)
    qq <- apply(lam, 1, quantile, c(0.25, 0.5, 0.75))
#    qq <- apply(lam, 1, quantile, c(0, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 1))
    IQR <- qq[3,] - qq[1,]
    out <- cbind(Mean=Mean, SD=SD, Median=qq[2,], IQR=IQR, One=lam[,1])
    ## km^2 vs ha diff is 100
    ## division by aa is to bump it back to original extent
    attr(out, "total") <- colSums(lam) * 100 / aa
    out
}

## fid != 4 because that is not a combo
if (fid < 4) {

    dat0 <- dat[rowSums(is.na(dat)) == 0,]
    #aa <- nrow(dat0) / AA
    aa <- 1

    if (wid %in% c(0,1)) {
    dat1 <- dat0
    dat1$HAB <- dat0$HAB_NALC
    dat1$isDM <- dat0$isDM_NALC
    dat1$isNF <- dat0$isNF_NALC
    Xn1 <- model.matrix(getTerms(mods, "formula"), dat1)
    colnames(Xn1) <- fixNames(colnames(Xn1))
    rm(dat1)
    NR <- nrow(Xn1)
    }

    if (wid %in% c(0,2)) {
    dat2 <- dat0
    dat2$HAB <- dat0$HAB_EOSD
    dat2$isDM <- dat0$isDM_EOSD
    dat2$isNF <- dat0$isNF_EOSD
    Xn2 <- model.matrix(getTerms(mods, "formula"), dat2)
    colnames(Xn2) <- fixNames(colnames(Xn2))
    rm(dat2)
    NR <- nrow(Xn2)
    }

    if (wid %in% c(0,3)) {
    dat3 <- dat0
    dat3$HAB <- dat0$HAB_LCC
    dat3$isDM <- dat0$isDM_LCC
    dat3$isNF <- dat0$isNF_LCC
    Xn3 <- model.matrix(getTerms(mods, "formula"), dat3)
    colnames(Xn3) <- fixNames(colnames(Xn3))
    rm(dat3)
    NR <- nrow(Xn3)
    }

    aicv <- aic
    if (wid == 1)
        aicv[,1] <- -Inf
    if (wid == 2)
        aicv[,2] <- -Inf
    if (wid == 3)
        aicv[,3] <- -Inf

    mu4 <- matrix(0, NR, 240)
    for (j in 1:240) {
        best <- which.min(aicv[j,])
        if (best == 1) {
            if (Stage > 6) {
                Xn1[,"HSH"] <- rowSums(pg4x4$nalc[rownames(dat0), sres$nalc[[j]]$hi])
                Xn1[,"HSH2"] <- Xn1[,"HSH2"]^2
            }
            mu4[,j] <- drop(Xn1 %*% est$nalc[j,])
        }
        if (best == 2) {
            if (Stage > 6) {
                Xn2[,"HSH"] <- rowSums(pg4x4$eosd[rownames(dat0), sres$eosd[[j]]$hi])
                Xn2[,"HSH2"] <- Xn2[,"HSH2"]^2
            }
            mu4[,j] <- drop(Xn2 %*% est$eosd[j,])
        }
        if (best == 3) {
            if (Stage > 6) {
                Xn3[,"HSH"] <- rowSums(pg4x4$lcc[rownames(dat0), sres$lcc[[j]]$hi])
                Xn3[,"HSH2"] <- Xn3[,"HSH2"]^2
            }
            mu4[,j] <- drop(Xn3 %*% est$lcc[j,])
        }
    }
    #rm(Xn1, Xn2, Xn3)
    lam <- lamfun(mu4)
    rownames(lam) <- rownames(dat0)
    rm(mu4, dat0)

    attr(lam, "spp") <- spp
    attr(lam, "stage") <- Stage
    attr(lam, "base-year") <- BASE_YEAR
    attr(lam, "bcr-jurs") <- regi
    gc()

    save(lam, file=file.path(ROOT2, "species", paste0(spp, "-ver3"), 
        paste0(paste(spp, wid, Stage, BASE_YEAR, regi, sep="-"), ".Rdata")))
    rm(lam)
}
if (fid == 4) {
    dat1 <- dat
    dat1$HAB_LCC2 <- NULL
    dat1$HAB_EOSD2 <- NULL
    dat1$isDM_LCC <- NULL
    dat1$isDM_EOSD <- NULL
    dat1$isNF_LCC <- NULL
    dat1$isNF_EOSD <- NULL
    dat1 <- dat1[rowSums(is.na(dat1)) == 0,]
    if (nrow(dat1) > 0) {
        #aa <- nrow(dat0) / AA
        aa <- 1

        dat1$HAB <- dat1$HAB_NALC
        dat1$isDM <- dat1$isDM_NALC
        dat1$isNF <- dat1$isNF_NALC
        Xn1 <- model.matrix(getTerms(mods, "formula"), dat1)
        colnames(Xn1) <- fixNames(colnames(Xn1))
        #rm(dat1)
        if (Stage > 6)
            stop("no 4x4 processed for full extent")
        mu4 <- matrix(0, nrow(Xn1), nrow(est))
        for (j in 1:nrow(est)) {
            mu4[,j] <- drop(Xn1[,colnames(est)] %*% est[j,])
        }
        rm(Xn1)
        lam <- lamfun(mu4)
        rownames(lam) <- rownames(dat1)
        rm(mu4, dat1)

    } else {
        #lam <- matrix(0, 0, 0)
        lam <- structure(numeric(0), .Dim = c(0L, 5L), .Dimnames = list(NULL, 
            c("Mean", "SD", "Median", "IQR", "One")))
    }
    attr(lam, "spp") <- spp
    attr(lam, "stage") <- Stage
    attr(lam, "base-year") <- BASE_YEAR
    attr(lam, "bcr-jurs") <- regi
    gc()
    save(lam, file=file.path(ROOT2, "species", paste0(spp, "-4"), 
        paste0(paste(spp, wid, Stage, BASE_YEAR, regi, sep="-"), ".Rdata")))
    rm(lam)
}
}
#}

as.matrix(rev(sort(sapply(ls(), function(x) object.size(get(x)))))[1:10])
sum(sapply(ls(), function(x) object.size(get(x))))


## mapping

library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")
load(file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))
load(file.path(ROOT2, "XYeosd.Rdata"))

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






Lc_quantile <- function (xx, probs=seq(0, 1, 0.1), type=c("L","p")) {
    xx <- xx[!is.na(xx)]
    o <- order(xx)
    x <- cumsum(xx[o]) / sum(xx)
    if (type=="L")
        q <- probs
    if (type=="p")
        q <- quantile(x, probs=probs, na.rm=TRUE)
    xxo <- xx[o]
    i <- sapply(q, function(z) min(xxo[x >= z]))
    i
}

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


