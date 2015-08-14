library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")

## this is for fid 1:3 (not 4 !)
mods <- mods_gfw

#fid <- 1
Xnlist <- list()
for (fid in 1:3) {
    fl <- c("analysis_package_gfwfire-nalc-2015-07-24.Rdata",
        "analysis_package_gfwfire-eosd-2015-07-24.Rdata",
        "analysis_package_gfwfire-lcc-2015-07-24.Rdata",
        "analysis_package_fire-nalc-2015-07-24.Rdata")
    e <- new.env()
    load(file.path(ROOT, "out", "data", fl[fid]), envir=e)

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
rm(xn, Xn)

spp <- "CAWA"
wid <- 0
#Stage <- which(names(mods) == "HS")
Stage <- which(names(mods) == "Dist")
BASE_YEAR <- 2015
regs <- gsub(".Rdata", "",
    gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))

######

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

est <- list(nalc=getEst(sres$nalc, stage = Stage, X=Xnlist$nalc),
    eosd=getEst(sres$eosd, stage = Stage, X=Xnlist$eosd),
    lcc=getEst(sres$lcc, stage = Stage, X=Xnlist$lcc))
aic <- cbind(nalc=getCaic(sres$nalc, stage = Stage),
    eosd=getCaic(sres$eosd, stage = Stage),
    lcc=getCaic(sres$lcc, stage = Stage))
table(colnames(aic)[apply(aic, 1, which.min)])


#regi <- "6_AB"
for (regi in regs) {

cat(spp, "--- Stage:", Stage, "--- setup:", 
    wid, "--- base yr:", BASE_YEAR, "--- region:", regi, "\n")
flush.console()

load(file.path(ROOT2, "chunks", paste0("pgdat-", regi, ".Rdata")))

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

if (wid == 1) {
    dat1 <- dat
    dat1$HAB <- dat$HAB_NALC
    dat1$isDM <- dat$isDM_NALC
    dat1$isNF <- dat$isNF_NALC
    dat1 <- dat1[!is.na(dat1$HAB),]
    aa <- nrow(dat1) / AA
    Xn1 <- model.matrix(getTerms(mods, "formula"), dat1)
    colnames(Xn1) <- fixNames(colnames(Xn1))
    mu1 <- pbapply(est$nalc, 1, function(z) Xn1 %*% z)
    lam <- lamfun(mu1)
    rownames(lam) <- rownames(dat1)
    rn(dat1, Xn1, mu1)
}
if (wid == 3) {
    dat2 <- dat
    dat2$HAB <- dat$HAB_EOSD
    dat2$isDM <- dat$isDM_EOSD
    dat2$isNF <- dat$isNF_EOSD
    dat2 <- dat2[!is.na(dat2$HAB),]
    aa <- nrow(dat2) / AA
    Xn2 <- model.matrix(getTerms(mods, "formula"), dat2)
    colnames(Xn2) <- fixNames(colnames(Xn2))
    mu2 <- pbapply(est$eosd, 1, function(z) Xn2 %*% z)
    lam <- lamfun(mu2)
    rownames(lam) <- rownames(dat2)
    rn(dat2, Xn2, mu2)
}
if (wid == 2) {
    dat3 <- dat
    dat3$HAB <- dat$HAB_LCC
    dat3$isDM <- dat$isDM_LCC
    dat3$isNF <- dat$isNF_LCC
    dat3 <- dat3[!is.na(dat3$HAB),]
    aa <- nrow(dat3) / AA
    Xn3 <- model.matrix(getTerms(mods, "formula"), dat3)
    colnames(Xn3) <- fixNames(colnames(Xn3))
    mu3 <- pbapply(est$lcc, 1, function(z) Xn3 %*% z)
    lam <- lamfun(mu3)
    rownames(lam) <- rownames(dat3)
    rn(dat3, Xn3, mu3)
}

## fid != 4 because that is not a combo
if (wid == 0) {

    dat0 <- dat[rowSums(is.na(dat)) == 0,]
    aa <- nrow(dat0) / AA

    dat1 <- dat0
    dat1$HAB <- dat0$HAB_NALC
    dat1$isDM <- dat0$isDM_NALC
    dat1$isNF <- dat0$isNF_NALC
    Xn1 <- model.matrix(getTerms(mods, "formula"), dat1)
    colnames(Xn1) <- fixNames(colnames(Xn1))
    rm(dat1)

    dat2 <- dat0
    dat2$HAB <- dat0$HAB_EOSD
    dat2$isDM <- dat0$isDM_EOSD
    dat2$isNF <- dat0$isNF_EOSD
    Xn2 <- model.matrix(getTerms(mods, "formula"), dat2)
    colnames(Xn2) <- fixNames(colnames(Xn2))
    rm(dat2)

    dat3 <- dat0
    dat3$HAB <- dat0$HAB_LCC
    dat3$isDM <- dat0$isDM_LCC
    dat3$isNF <- dat0$isNF_LCC
    Xn3 <- model.matrix(getTerms(mods, "formula"), dat3)
    colnames(Xn3) <- fixNames(colnames(Xn3))
    rm(dat3)

    mu4 <- matrix(0, nrow(Xn3), 240)
    for (j in 1:240) {
        best <- which.min(aic[j,])
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
    rm(Xn1, Xn2, Xn3)
    lam <- lamfun(mu4)
    rownames(lam) <- rownames(dat0)
    rm(mu4, dat0)
}

attr(lam, "spp") <- spp
attr(lam, "stage") <- Stage
attr(lam, "base-year") <- BASE_YEAR
attr(lam, "bcr-jurs") <- regi
gc()

save(lam, file=file.path(ROOT2, "species", spp, 
    paste0(paste(spp, wid, Stage, BASE_YEAR, regi, sep="-"), ".Rdata")))
}

as.matrix(rev(sort(sapply(ls(), function(x) object.size(get(x)))))[1:10])
sum(sapply(ls(), function(x) object.size(get(x))))

library(RColorBrewer)
z <- lam[,"Median"]
#z <- dat$CTI

hist(attr(lam, "total"))

plotfun <- function(lam, what="Median") {
    z <- lam[,what]
    Col <- rev(brewer.pal(10, "RdBu"))
    cz <- cut(z, breaks = c(min(z)-1, quantile(z, seq(0.1, 1, 0.1))))
    plot(dat[rownames(lam),2:3], col = Col[cz], pch=".")
    invisible(NULL)
}

plotfun(z)

z <- dat1$HAB
Col <- brewer.pal(9, "Set3")
plot(dat[,2:3], col = Col[z], pch=".")

hs <- rowSums(pg4x4$nalc[rownames(dat0), c("Decid", "Mixed")])
plot(hs, aa)
table(unlist(lapply(sres$nalc, "[[", "hi")))
