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


Stage <- 6
fid <- 1

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


#mods <- if (fid == 4)
#    mods_fire else mods_gfw
fn <- paste0("bam-", fid, "_", spp, ".Rdata", sep="")
load(file.path(ROOT, "out", folder, fn))

sum(getOK(res)) / length(res)

## need to load data for xn, Xn

est <- getEst(res, stage=NULL)

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
resid_yr(1)
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

