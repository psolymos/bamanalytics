library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/Apr2016/out"
source("~/repos/bamanalytics/R/makingsense_functions.R")
fstat <- function(x, level=0.95) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
}
fstatv <- function(x, level=0.95) {
    out <- data.frame(rowMeans(x),
        apply(x, 1, sd),
        x[,1],
        t(pbapply(x, 1, quantile, c(0.5, (1-level)/2, 1 - (1-level)/2))))
    colnames(out) <- c("Mean","SD","Run1","Median","LCL","UCL")
    out
}

PROJECT <- "bam"
Date <- "2016-08-16"
level <- 0.9

e <- new.env()
load(file.path(ROOT2, "data", "pack_2016-08-16.Rdata"), envir=e)

mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
yy <- e$YY
xn <- e$DAT[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
xn <- xn[rownames(Xn),]
off <- e$OFF[rownames(xn),]
#bb <- e$BB
bbb <- unique(e$BB)
rm(e)

modTab <- getFancyModsTab(mods)
xnh <- nonDuplicated(xn, HABTR, TRUE)[,c("HAB","HABTR","isNF","isDev",
    "isWet","isOpn","isDM","isDec","isMix")]
xnh <- xnh[c("ConifDense", "ConifSparse","ConifOpen",
    "DecidDense", "DecidSparse", "DecidOpen",
    "MixedDense", "MixedSparse", "MixedOpen",
    "WetDense", "WetSparse", "WetOpen",
    "Shrub", "Grass", "Barren", "Agr", "Devel"),]

spp <- "CAWA"

fn <- file.path(ROOT2, "results", "cawa",
    paste0(PROJECT, "_", spp, "_", Date, ".Rdata"))
load(fn)
100 * sum(getOK(res)) / length(res)
est_hab <- getEst(res, stage = 2, X=Xn)
est_habhgt <- getEst(res, stage = 3, X=Xn)
est_dtb <- getEst(res, stage = 4, X=Xn)
est_wet <- getEst(res, stage = 5, X=Xn)
est <- getEst(res, stage = length(mods)-1, X=Xn)
est_yr <- getEst(res, stage = length(mods), X=Xn)


printCoefmat(su <- getSummary(res))
(mt <- getFancyMidTab(res, mods))
plotMid(res, mods, web=TRUE)



## road
xn1 <- xnh
xn1$ROAD <- 1
Xn1 <- model.matrix(getTerms(mods[2], "formula", intercept=FALSE), xn1)
colnames(Xn1) <- fixNames(colnames(Xn1))
est1 <- est[,colnames(Xn1)]

pr <- t(apply(est1, 1, function(z) Xn1 %*% z))
colnames(pr) <- rownames(Xn1)
pr <- exp(pr)
pr[pr>2] <- 2

op <- par(mar=c(5,8,2,2), las=1)
boxplot(pr[,rev(colnames(pr))], horizontal=TRUE, range=0,
    xlab="Expected abundance: On-road / Off-road",
    col=terrain.colors(nlevels(xn$HABTR)),
    main=spp)
abline(v=1, col=2)
par(op)

## habitat (HGT taken as mean)

xn2 <- xnh
xn2$ROAD <- 0
Xn2 <- model.matrix(getTerms(mods[1:2], "formula"), xn2)
colnames(Xn2) <- fixNames(colnames(Xn2))
est2 <- est_hab[,colnames(Xn2)]

pr <- exp(t(apply(est2, 1, function(z) Xn2 %*% z)))
colnames(pr) <- rownames(xn2)

pr <- pr[,order(colMeans(pr))]

op <- par(mar=c(5,8,2,2), las=1)
boxplot(pr, horizontal=TRUE, range=0,
    xlab="Expected density (males / ha)",
    col=rev(terrain.colors(ncol(pr))),
    main=spp)
par(op)


HGT <- seq(0,1,by=0.01)
xn2 <- expand.grid(HABTR=factor(c("ConifDense", #"ConifSparse","ConifOpen",
    "DecidDense", #"DecidSparse", "DecidOpen",
    "MixedDense", #"MixedSparse", "MixedOpen",
    "WetDense"), #"WetSparse", "WetOpen"),
    levels(xn$HABTR)), HGT=HGT)
xn2 <- data.frame(xnh[match(xn2$HABTR, rownames(xnh)),],
    ROAD=0, HGT=xn2$HGT, HGT2=xn2$HGT^2, HGT05=sqrt(xn2$HGT))
Xn2 <- model.matrix(getTerms(mods[1:3], "formula"), xn2)
colnames(Xn2) <- fixNames(colnames(Xn2))
est2 <- est_habhgt[,colnames(Xn2)]

pr <- exp(t(apply(est2, 1, function(z) Xn2 %*% z)))
xn2$Density <- colMeans(pr)
xn2$lcl <- apply(pr, 2, quantile, 0.05)
xn2$ucl <- apply(pr, 2, quantile, 0.95)

lam <- t(matrix(xn2$Density, nrow=4))
op <- par(las=1)
matplot(HGT*25, lam, type="l", lwd=2, ylim=c(0, 1.2*max(lam)),
    ylab="Density (males/ha)", xlab="Height (m)", main=spp,
    col=1:4, lty=1)
legend("topright",
    lty=1, lwd=2, bty="n", col=1:4, legend=c("Conif", "Decid", "Mixed", "Wet"))
par(op)


printCoefmat(getSummary(res)[c("SLP","SLP2"),])
printCoefmat(getSummary(res)[c("YSD","YSF","YSL","YR"),])
summary(100 * (exp(est_yr[,"YR"]) - 1))

## Marginal plots

#pr <- exp(apply(est, 1, function(z) Xn %*% z))
pr <- exp(apply(est_wet, 1, function(z) Xn %*% z))
xn$lam_hat <- rowMeans(pr)

COL <- rgb(65/255, 105/255, 225/255, alpha=0.1)

plot(xn$SLP^2*90, xn$lam_hat, col=COL, pch=21,
    main=spp, xlab="Slope (degrees)", ylab="Density (males/ha)")
lines(lowess(xn$SLP^2*90, xn$lam_hat), col=2, lwd=3)

plot(xn$CTI*4+8, xn$lam_hat, col=COL, pch=21,
    main=spp, xlab="CTI", ylab="Density (males/ha)")
lines(lowess(xn$CTI*4+8, xn$lam_hat), col=2, lwd=3)

#plot(xn$SLP^2*90, xn$CTI*4+8)


boxplot(lam_hat ~ HAB, xn, range=0)

boxplot(lam_hat ~ ROAD, xn, range=0)

## exclude YS=0 ???
par(mfrow=c(1,3))
plot(50*(1-xn$YSL), xn$lam_hat, col=COL, pch=21,
    main=spp, xlab="Years since last disturbance: loss", ylab="Density (males/ha)")
lines(lowess(50*(1-xn$YSL), xn$lam_hat), col=2, lwd=3)
plot(50*(1-xn$YSF), xn$lam_hat, col=COL, pch=21,
    main=spp, xlab="Years since last disturbance: fire", ylab="Density (males/ha)")
lines(lowess(50*(1-xn$YSF), xn$lam_hat), col=2, lwd=3)
plot(50*(1-xn$YSD), xn$lam_hat, col=COL, pch=21,
    main=spp, xlab="Years since last disturbance: both", ylab="Density (males/ha)")
lines(lowess(50*(1-xn$YSD), xn$lam_hat), col=2, lwd=3)


par(mfrow=c(4,1))
plot(lam_hat ~ I(jitter(HGT*25)), xn[xn$HAB == "Decid",], col=COL, xlim=range(xn$HGT*25), ylim=c(0,0.2))
plot(lam_hat ~ I(jitter(HGT*25)), xn[xn$HAB == "Mixed",], col=COL, xlim=range(xn$HGT*25), ylim=c(0,0.2))
plot(lam_hat ~ I(jitter(HGT*25)), xn[xn$HAB == "Conif",], col=COL, xlim=range(xn$HGT*25), ylim=c(0,0.2))
plot(lam_hat ~ I(jitter(HGT*25)), xn[xn$HAB == "Wet",], col=COL, xlim=range(xn$HGT*25), ylim=c(0,0.2))

## GoF

mn <- matrix(0, nrow(xn), length(mods)+1)
colnames(mn) <- c("NULL", names(mods))
rownames(mn) <- rownames(xn)
for (i in 0:length(mods)) {
    gc()
    cat(i, "/", length(mods), "\n");flush.console()
    est_i <- getEst(res, stage = i, X=Xn)
    col_keep <- colSums(abs(est_i) > 0) != 0
    pr <- exp(pbsapply(1:nrow(est_i), function(j)
        Xn[,colnames(est_i[,col_keep,drop=FALSE]),drop=FALSE] %*%
        est_i[j,col_keep]))
    #st <- fstatv(pr, level)
    mn[,i+1] <- rowMeans(pr)
}

ss1 <- which(!(1:nrow(xn) %in% bbb)) # test portion
#ss1 <- which((1:nrow(xn) %in% bbb)) # other portion
library(pROC)
Y <- yy[,spp]
Y1 <- ifelse(yy[,spp]>0, 1, 0)
off1 <- off[,spp]

## ROC/AUC
rocAll1 <- pblapply(1:ncol(mn), function(i) {
      pp <- mn[ss1,i] * exp(off1[ss1])
      roc(Y1[ss1], pp)
  })
names(rocAll1) <- c("NULL", names(mods))
auc <- sapply(rocAll1, function(z) as.numeric(z$auc))
barplot(auc, ylim=c(0,1), space=0, ylab="AUC", xlab="Stages")
lines(0:length(mods)+0.5, auc, col=2, lwd=2)

## QQ
vals <- seq(min(Y[ss1]), max(Y[ss1]), 1)
pobs <- numeric(length(vals))
names(pobs) <- vals
tab <- table(Y[ss1])
for (k in vals) {
    pobs[as.character(k)] <- if (as.character(k) %in% names(tab))
        tab[as.character(k)]/length(ss1) else 0
}
pobs <- cumsum(pobs)

pexp <- matrix(0, length(vals), ncol(mn))
rownames(pexp) <- vals
colnames(pexp) <- colnames(mn)
for (i in 1:ncol(mn)) {
    m_pexp <- matrix(0, length(ss1), length(vals))
    colnames(m_pexp) <- vals
    for (k in vals) {
        m_pexp[,as.character(k)] <- dpois(x=rep(k, length(ss1)),
            lambda=mn[ss1,i] * exp(off1[ss1]))
    }
    m_pexp[,ncol(m_pexp)] <- 1 - rowSums(m_pexp[,-ncol(m_pexp)])
    pexp[,i] <- cumsum(colMeans(m_pexp))
}

## plot
pro <- table(Y[ss1])/length(ss1)
p_min <- min(cbind(pobs,pexp))
par(mfrow=c(3,2))
boxplot(mn[ss1,"NULL"] * exp(off1[ss1]) ~ Y[ss1], range=0,
    at=cumsum(2*pro^0.2), width=2*pro^0.2, main="NULL",
    xlab="Observed count", ylab="Corrected density")
boxplot(mn[ss1,"Clim"] * exp(off1[ss1]) ~ Y[ss1], range=0,
    at=cumsum(2*pro^0.2), width=2*pro^0.2, main="Clim",
    xlab="Observed count", ylab="Corrected density")
plot(pobs, pexp[,"NULL"], type="b", pch=rownames(pexp),
    xlab="Observed", ylab="Expected", xlim=c(p_min, 1), ylim=c(p_min, 1))
abline(0, 1, col="grey")
plot(pobs, pexp[,"Clim"], type="b", pch=rownames(pexp),
    xlab="Observed", ylab="Expected", xlim=c(p_min, 1), ylim=c(p_min, 1))
abline(0, 1, col="grey")
plot(rocAll1[["NULL"]])
plot(rocAll1[["Clim"]])

## Std Resid Trend


## -- old

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

## maps

library(mefa4)
library(mapview)
library(sp)
ROOT2 <- "e:/peter/bam/Apr2016/out"
e <- new.env()
load(file.path(ROOT2, "data", "pack_2016-08-16.Rdata"), envir=e)

xsp <- e$DAT[,c("X","Y")]
coordinates(xsp) <- ~X+Y
proj4string(xsp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

m <- mapview(xsp)



