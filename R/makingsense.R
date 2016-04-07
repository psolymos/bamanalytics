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
ids$prefix <- with(ids, paste0( 
    TEXT, "_", SEXT, "_", LCTU, "_", Date))
rownames(ids) <- 1:8


Stage <- 6
fid <- 1

## height histograms
if (FALSE) {
Height <- xn$HGT*50
HabClass <- xn$HAB
library(lattice)
densityplot(~ Height | HabClass)

}

for (fid in 1:6) {

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
off <- e$OFF
bb <- e$BB
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)
est <- getEst(res, stage = Stage, X=Xn)
#est_hab <- getEst(res, stage = 4, X=Xn)

## output for Sam

su <- getSummary(res)
write.csv(su, file.path(ROOT, "out", "nmbca-samuel", 
    paste0("CoefSE_", as.character(ids$prefix[fid]), ".csv")))
mt <- getFancyMidTab(res, mods)
write.csv(mt, file.path(ROOT, "out", "nmbca-samuel", 
    paste0("MIDtab_", as.character(ids$prefix[fid]), ".csv")))

pdf(file.path(ROOT, "out", "nmbca-samuel", 
    paste0("MIDplot_", as.character(ids$prefix[fid]), ".pdf")))
plotMid(res, mods, web=TRUE)
dev.off()

}
write.csv(ids, file.path(ROOT, "out", "nmbca-samuel", "model-ids.csv"))

for (fid in 1:2) {
e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
off <- e$OFF
bb <- e$BB
rm(e)
aa <- getFancyModsTab(mods)
write.csv(aa, file.path(ROOT, "out", "nmbca-samuel", 
    paste0("model-table-", as.character(ids$TEXT[fid]), ".csv")))
}

## habitat association graphs for report
for (fid in 1:6) {

cat(fid, "\n");flush.console()

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
off <- e$OFF
bb <- e$BB
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)
est <- getEst(res, stage = 2, X=Xn) # Hab + Road

## road
xn1 <- expand.grid(HAB=levels(xn$HAB), TR3=levels(xn$TR3))
xn1$ROAD <- 1
Xn1 <- model.matrix(getTerms(mods[2], "formula"), xn1)
colnames(Xn1) <- fixNames(colnames(Xn1))
est1 <- est[,colnames(Xn1)]
Xn1 <- rbind(0, Xn1)
Xn1[1,1] <- 1
rownames(Xn1) <- c("Reference", paste0(xn1$TR3, ".", xn1$HAB))

pr <- t(apply(est1, 1, function(z) Xn1 %*% z))
colnames(pr) <- rownames(Xn1)
pr <- exp(pr - pr[,1])
pr[pr>2] <- 2

png(file.path(ROOT, "out", "nmbca-samuel", 
    paste0("CAWAfig_Road_", as.character(ids$prefix[fid]), ".png")))
op <- par(mar=c(5,8,2,2), las=1)
boxplot(pr[,-1], horizontal=TRUE, range=0, 
    xlab="Expected abundance: On-road / Off-road",
    col=rep(terrain.colors(3), each=nlevels(xn$HAB)),
    main=as.character(ids$prefix[fid]))
abline(v=1, col=2)
par(op)
dev.off()

## habitat (HGT taken as mean)

xn2 <- expand.grid(HAB=levels(xn$HAB), TR3=levels(xn$TR3))
xn2$ROAD <- 0
xn2$isNF <- ifelse(xn2$HAB %in% 
    c("Agr", "Barren", "Burn", "Devel", "Grass", "Wet"), 1L, 0L)
xn2$isDM <- ifelse(xn2$HAB %in% c("Decid", "Mixed"), 1, 0)
xn2$HGT <- ifelse(xn2$isNF == 0, mean(xn$HGT[xn$isNF==0]), mean(xn$HGT[xn$isNF==1]))
xn2$HGT2 <- xn2$HGT^2
Xn2 <- model.matrix(getTerms(mods[1:2], "formula"), xn2)
colnames(Xn2) <- fixNames(colnames(Xn2))
est2 <- est[,colnames(Xn2)]
rownames(Xn2) <- paste0(xn2$TR3, ".", xn2$HAB)

pr <- exp(t(apply(est2, 1, function(z) Xn2 %*% z)))
colnames(pr) <- rownames(Xn2)

#pr <- pr[,1:nlevels(xn$HAB)]
#pr <- pr[,(nlevels(xn$HAB)+1):(2*nlevels(xn$HAB))]
pr <- pr[,(2*nlevels(xn$HAB)+1):(3*nlevels(xn$HAB))]
colnames(pr) <- levels(xn$HAB)
pr <- pr[,order(colMeans(pr))]

png(file.path(ROOT, "out", "nmbca-samuel", 
    paste0("CAWAfig_Lcover_", as.character(ids$prefix[fid]), ".png")))
op <- par(mar=c(5,8,2,2), las=1)
boxplot(pr, horizontal=TRUE, range=0, 
    xlab="Expected density (males / ha)",
    col=rev(terrain.colors(nlevels(xn$HAB))),
    main=as.character(ids$prefix[fid]))
par(op)
dev.off()


HGT <- seq(0,0.5,by=0.01)
xn2 <- expand.grid(HAB=factor(c("Conif", "Decid", "Mixed", "Wet"), levels(xn$HAB)), 
    HGT=HGT)
xn2$TR3 <- factor("Dense", levels(xn$TR3))
xn2$ROAD <- 0
xn2$isNF <- ifelse(xn2$HAB %in% 
    c("Agr", "Barren", "Burn", "Devel", "Grass", "Wet"), 1L, 0L)
xn2$isDM <- ifelse(xn2$HAB %in% c("Decid", "Mixed"), 1, 0)
xn2$HGT2 <- xn2$HGT^2
xn2$Height <- 50*xn2$HGT
Xn2 <- model.matrix(getTerms(mods[1:2], "formula"), xn2)
colnames(Xn2) <- fixNames(colnames(Xn2))
est2 <- est[,colnames(Xn2)]


pr <- exp(t(apply(est2, 1, function(z) Xn2 %*% z)))
xn2$Density <- colMeans(pr)
xn2$lcl <- apply(pr, 2, quantile, 0.05)
xn2$ucl <- apply(pr, 2, quantile, 0.95)

png(file.path(ROOT, "out", "nmbca-samuel", 
    paste0("CAWAfig_Height_", as.character(ids$prefix[fid]), ".png")))
op <- par(las=1)
lam <- t(matrix(xn2$Density, nrow=4))
matplot(HGT*50, lam, lty=1, type="l", lwd=2, ylim=c(0, 1.2*max(lam)),
    ylab="Density (males/ha)", xlab="Height (m)", main=as.character(ids$prefix[fid]))
legend("topright",
    lty=1, lwd=2, bty="n", col=1:4, legend=c("Conif", "Decid", "Mixed", "Wet"))
par(op)
dev.off()

}


#mods <- if (fid == 4)
#    mods_fire else mods_gfw
fn <- paste0("bam-", fid, "_", spp, ".Rdata", sep="")
load(file.path(ROOT, "out", folder, fn))

sum(getOK(res)) / length(res)

## need to load data for xn, Xn
#est <- getEst(res, stage=NULL)

getCaic(res)
printCoefmat(getSummary(res))

round(cbind(getSummary(res, show0=TRUE)[,c(1,3)], 
    getConfint(res, 0.9, "quantile", show0=TRUE))[c("CTI","CTI2",
    "SLP","SLP2"),], 3)

round(cbind(getSummary(res, show0=TRUE)[,c(1,3)], 
    getConfint(res, 0.9, "quantile", show0=TRUE))[c("LIN","POL",
    "BRN","LSS","DTB",  "YSD","YSL","YSF"),], 3)

getMidPure(res, mods)
getFancyMid(res, mods)
getFancyMidTab(res, mods)

getFancyModsTab(mods)

plotMid(res, mods, web=TRUE)

## ARU effect
est6 <- getEst(res, stage=6)
aru <- est6[,"ARU"]
## exclude very small values (bimodal!)
## those came from (?) 0 obs per bootstrap run ???
summary(sqrt(exp(aru[aru > -10])))
#summary(sqrt(exp(aru)))

## year effect

trend_res <- list()

for (fid in 1:6) {

fo <- paste0(spp, "-", ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date)
cat(fo, "\n");flush.console()

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
off <- e$OFF
bb <- e$BB
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)

est <- getEst(res, stage=6)
est_yr <- getEst(res, stage=7)
## decadal
#tw <- 100 * (exp(est_yr[,"YR"]) - 1)
#te <- 100 * (exp(est_yr[,"YR"] + est_yr[,"EWE:YR"]) - 1)
## yearly
tw <- 100 * (exp(0.1*est_yr[,"YR"]) - 1)
te <- 100 * (exp(0.1*(est_yr[,"YR"] + est_yr[,"EWE:YR"])) - 1)

summary(tw)
summary(te)
fstat <- round(cbind(West=c(Median=median(tw), quantile(tw, c(0.05, 0.95))),
    East=c(Median=median(te), quantile(te, c(0.05, 0.95)))),2)

stopifnot(all(colnames(Xn) == colnames(est)))
lam <- exp(pbapply(est, 1, function(z) Xn %*% z))
lam_hat <- pbapply(lam, 1, median)
pr <- log(lam_hat) + off[,spp]
y <- yy[,spp]
if (spp == "CAWA")
    y[y>5] <- 6
r <- (y - exp(pr)) / sqrt(y + 0.5)

NN <- aggregate(y, list(YR=10*xn$YR+2001,EW=xn$EW), mean)
RR <- aggregate(r, list(YR=10*xn$YR+2001,EW=xn$EW), mean)
NN <- NN[NN$YR >= 2002 & NN$YR <= 2012,]
RR <- RR[RR$YR >= 2002 & RR$YR <= 2012,]

dd <- data.frame(y=y, pr=pr, YR=xn$YR, EW=xn$EW)

resid_yr <- function(i) {
    m <- glm(y ~ offset(pr) + YR + EW:YR, dd[bb[,i],], family=poisson)
    c(mtw = 100 * (exp(0.1*coef(m)["YR"]) - 1),
        mte = 100 * (exp(0.1*(coef(m)["YR"] + coef(m)["YR:EWE"])) - 1))
}
#resid_yr(1)
yr_res <- t(pbsapply(1:ncol(bb), resid_yr))
fstatm <- round(cbind(West=c(Median=median(yr_res[,"mtw.YR"]), 
    quantile(yr_res[,"mtw.YR"], c(0.05, 0.95))),
    East=c(Median=median(yr_res[,"mte.YR"]), 
    quantile(yr_res[,"mte.YR"], c(0.05, 0.95)))),2)

pdf(file.path(ROOT2, "species", "cawa-nmbca-trend", paste0("trend-", fo, ".pdf")),
    width=7, height=10)
op <- par(mfrow=c(4,2))

plot(x ~ YR, NN[NN$EW=="W",], ylab="Annual Mean Abundance Index", xlab="Year", 
    type="b", col=1, pch=19, main="CAWA, West", ylim=c(0, max(NN$x)))
abline(lm(x ~ YR, NN[NN$EW=="W",]), col="red4")
plot(x ~ YR, NN[NN$EW=="E",], ylab="Annual Mean Abundance Index", xlab="Year", 
    type="b", col=1, pch=19, main="CAWA, East", ylim=c(0, max(NN$x)))
abline(lm(x ~ YR, NN[NN$EW=="E",]), col="red4")

plot(x ~ YR, RR[RR$EW=="W",], ylab="Std. Residuals", xlab="Year", 
    type="b", col=1, pch=19, main="", ylim=range(RR$x))
abline(lm(x ~ YR, RR[RR$EW=="W",]), col="red4")
plot(x ~ YR, RR[RR$EW=="E",], ylab="Std. Residuals", xlab="Year", 
    type="b", col=1, pch=19, main="", ylim=range(RR$x))
abline(lm(x ~ YR, RR[RR$EW=="E",]), col="red4")

hist(tw, col="gold", xlab="Modeled Annual Trend (%)", main="", 
    xlim=range(c(yr_res,tw,te)), freq=FALSE)
abline(v=fstat[1,"West"], col="red4", lty=1, lwd=2)
abline(v=fstat[2:3,"West"], col="red4", lty=2, lwd=1)
hist(te, col="gold", xlab="Modeled Annual Trend (%)", main="", 
    xlim=range(c(yr_res,tw,te)), freq=FALSE)
abline(v=fstat[1,"East"], col="red4", lty=1, lwd=2)
abline(v=fstat[2:3,"East"], col="red4", lty=2, lwd=1)

hist(yr_res[,"mtw.YR"], col="gold", xlab="Residual Annual Trend (%)", 
    main="", xlim=range(c(yr_res,tw,te)), freq=FALSE)
abline(v=fstatm[1,"West"], col="red4", lty=1, lwd=2)
abline(v=fstatm[2:3,"West"], col="red4", lty=2, lwd=1)
hist(yr_res[,"mte.YR"], col="gold", xlab="Residual Annual Trend (%)", 
    main="", xlim=range(c(yr_res,tw,te)), freq=FALSE)
abline(v=fstatm[1,"East"], col="red4", lty=1, lwd=2)
abline(v=fstatm[2:3,"East"], col="red4", lty=2, lwd=1)

par(op)
dev.off()

trend_res[[fo]] <- list(NN=NN, RR=RR, rtr=yr_res, mtr=cbind(W=tw, E=te))
}

save(trend_res, file=file.path(ROOT2, "species", "cawa-nmbca-trend", "trend-CAWA.Rdata"))


## Province level residual trend
load("c:/bam/May2015/out/data_package_2015-08-26.Rdata")
PKEY$JURS <- SS$JURS[match(PKEY$SS, SS$SS)]
trend_res2 <- list()

for (fid in c(1,3,5)) {

fo <- paste0(spp, "-", ids$TEXT[fid], "_", ids$SEXT[fid], "_", ids$LCTU[fid], "_", Date)
cat(fo, "\n");flush.console()

e <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[fid])), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
off <- e$OFF
bb <- e$BB
bb <- bb[,1:100]
rm(e)

load(file.path(ROOT, "out", "results", as.character(ids$fn[fid])))
100 * sum(getOK(res)) / length(res)

est <- getEst(res, stage=6)
stopifnot(all(colnames(Xn) == colnames(est)))
lam <- exp(pbapply(est, 1, function(z) Xn %*% z))
lam_hat <- pbapply(lam, 1, median)
pr <- log(lam_hat) + off[,spp]
y <- yy[,spp]
if (spp == "CAWA")
    y[y>5] <- 6

xn$JURS <- droplevels(PKEY$JURS[match(rownames(xn), PKEY$PKEY)])
xn$SS <- PKEY$SS[match(rownames(xn), PKEY$PKEY)]
table(xn$JUSR)
summary(xn$JUSR)
xn$isBBS <- ifelse(substr(rownames(xn), 1, 3) == "BBS", "BBS", "BAM")
xn$YEAR <- 10*xn$YR+2001
#table(xn$YEAR, xn$isBBS)
#xn$X <- SS$X[match(xn$SS, SS$SS)]
#xn$Y <- SS$Y[match(xn$SS, SS$SS)]
#tmp <- aggregate(xn$Y, list(BBS=xn$isBBS, Jurs=xn$JURS), mean)
#barplot(matrix(tmp$x, nrow=2), beside=TRUE, ylim=c(40,60))
#summary(lm(Y ~ isBBS*JURS, xn))

dd <- data.frame(y=y, pr=pr, YR=xn$YR, YEAR=xn$YEAR, JURS=xn$JURS, isBBS=xn$isBBS)
resid_yr_jurs <- function(i, BBSonly=FALSE) {
    dat <- dd[bb[,i],]
    dat <- dat[dat$YEAR >= 2002 & dat$YEAR <= 2012, ]
    if (BBSonly)
        dat <- dat[dat$isBBS=="BBS",]
    dat <- dat[dat$JURS %in% c("AB","ON","QC"),]
    dat$JURS <- droplevels(dat$JURS)
    m <- glm(y ~ offset(pr) + YR + JURS:YR, dat, family=poisson)
    cf <- coef(m)[-1]
    cf[-1] <- cf[1] + cf[-1]
    names(cf) <- levels(dat$JURS)
    100 * (exp(0.1*cf) - 1)
}

all <- t(pbsapply(1:ncol(bb), resid_yr_jurs, BBSonly=FALSE))
bbs <- t(pbsapply(1:ncol(bb), resid_yr_jurs, BBSonly=TRUE))

trend_res2[[fo]] <- list(all=all, bbs=bbs)

}

save(trend_res2, file=file.path(ROOT2, "species", "cawa-nmbca-trend", "trend2-CAWA.Rdata"))

level <- 0.9

Tr <- rbind(t(apply(data.frame(Fid1.All=trend_res2[[1]]$all, Fid1.BBS=trend_res2[[1]]$bbs), 
    2, quantile, c(0.5, (1-level)/2, 1-(1-level)/2))),
    t(apply(data.frame(Fid3.All=trend_res2[[2]]$all, Fid3.BBS=trend_res2[[2]]$bbs), 
    2, quantile, c(0.5, (1-level)/2, 1-(1-level)/2))),
    t(apply(data.frame(Fid5.All=trend_res2[[3]]$all, Fid5.BBS=trend_res2[[3]]$bbs), 
    2, quantile, c(0.5, (1-level)/2, 1-(1-level)/2))))
write.csv(Tr, file=file.path(ROOT2, "species", "cawa-nmbca-trend", "trend2-CAWA.csv"))

#x0 <- 1:6
par(mfrow=c(3,1))
x0 <- c(1,2, 5,6, 9,10)
plot(x0, Tr[(1:6)[c(1,4,2,5,3,6)],1], col=rep(c(2, 4), 3), 
    pch=19, ylim=range(Tr),
    axes=FALSE, ann=FALSE, cex=1.5)
abline(h=0, col="grey")
segments(x0=x0, y0=Tr[(1:6)[c(1,4,2,5,3,6)],2], y1=Tr[(1:6)[c(1,4,2,5,3,6)],3],
    col=rep(c(2, 4), 3), lwd=3)
axis(2)
axis(1, at = c(1.5, 5.5, 9.5), labels = c("AB","ON","QC"), tick=FALSE)
box()
title(main="CAWA, GFW+Fire, Canada, NALCMS")
legend("topright", pch=19, lty=1, lwd=2, col=c(2,4), 
    legend=c("All", "BBS only"), bty="n")

x0 <- c(1,2, 5,6, 9,10)
plot(x0, Tr[(7:12)[c(1,4,2,5,3,6)],1], col=rep(c(2, 4), 3), 
    pch=19, ylim=range(Tr),
    axes=FALSE, ann=FALSE, cex=1.5)
abline(h=0, col="grey")
segments(x0=x0, y0=Tr[(7:12)[c(1,4,2,5,3,6)],2], y1=Tr[(7:12)[c(1,4,2,5,3,6)],3],
    col=rep(c(2, 4), 3), lwd=3)
axis(2)
axis(1, at = c(1.5, 5.5, 9.5), labels = c("AB","ON","QC"), tick=FALSE)
box()
title(main="CAWA, GFW+Fire, Canada, LCC")
legend("topright", pch=19, lty=1, lwd=2, col=c(2,4), 
    legend=c("All", "BBS only"), bty="n")

x0 <- c(1,2, 5,6, 9,10)
plot(x0, Tr[(13:18)[c(1,4,2,5,3,6)],1], col=rep(c(2, 4), 3), 
    pch=19, ylim=range(Tr),
    axes=FALSE, ann=FALSE, cex=1.5)
abline(h=0, col="grey")
segments(x0=x0, y0=Tr[(13:18)[c(1,4,2,5,3,6)],2], y1=Tr[(13:18)[c(1,4,2,5,3,6)],3],
    col=rep(c(2, 4), 3), lwd=3)
axis(2)
axis(1, at = c(1.5, 5.5, 9.5), labels = c("AB","ON","QC"), tick=FALSE)
box()
title(main="CAWA, GFW+Fire, Canada, EOSD")
legend("topright", pch=19, lty=1, lwd=2, col=c(2,4), 
    legend=c("All", "BBS only"), bty="n")

## legend for all/bbs

e1 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[1])), envir=e1)
e2 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[2])), envir=e2)
e3 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[3])), envir=e3)
e4 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[4])), envir=e4)
e5 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[5])), envir=e5)
e6 <- new.env()
load(file.path(ROOT, "out", "data", as.character(ids$data[6])), envir=e6)

dim(e1$YY)
dim(e2$YY)
dim(e3$YY)
dim(e4$YY)
dim(e5$YY)
dim(e6$YY)
