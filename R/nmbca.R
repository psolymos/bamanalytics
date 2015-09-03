library(mefa4)
ROOT <- "c:/bam/May2015"
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

#fid <- 1
Est <- list()
Aic <- list()
Mid <- list()

for (fid in 1:6) {

cat(fid, "\n");flush.console()

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

## need to load data for xn, Xn
Est[[fid]] <- getEst(res, stage=NULL)

Aic[[fid]] <- sapply(1:7, function(z) getCaic(res, stage=z))
colnames(Aic[[fid]]) <- names(mods)

Mid[[fid]] <- getFancyMid(res, mods)

}


Aic7 <- sapply(Aic, function(z) z[,"Year"])
t(apply(Aic7, 1, function(z) z - min(z)))
table(apply(Aic7, 1, which.min))

Aic6 <- sapply(Aic, function(z) z[,"Clim"])
t(apply(Aic6, 1, function(z) z - min(z)))
table(apply(Aic6, 1, which.min))

## percent annual change from year effect
int <- 1/10
lapply(Est, function(est) {
    tw <- 100 * (exp(int*est[,"YR"]) - 1)
    te <- 100 * (exp(int*(est[,"YR"] + est[,"EWE:YR"])) - 1)
    round(cbind(West=c(Mean=mean(tw), quantile(tw, c(0.025, 0.5, 0.975))),
        East=c(Mean=mean(te), quantile(te, c(0.025, 0.5, 0.975)))),2)
})

load("c:/bam/May2015/out/figs/nmbca2/popsize-CAWA-6-2013-2015-09-02.Rdata")
ttt13 <- ttt
load("c:/bam/May2015/out/figs/nmbca2/popsize-CAWA-6-2003-2015-09-02.Rdata")
ttt03 <- ttt

#i <- 1
fun <- function(i) {
int <- 1/10
x <- colSums(ttt13[[i]]) / colSums(ttt03[[i]])
yw <- exp(Est[[i]][,"YR"])
ye <- exp(Est[[i]][,"YR"] + Est[[i]][,"EWE:YR"])
t <- 100 * (exp(int*log(x)) - 1)
tw <- 100 * (exp(int*log(yw)) - 1)
te <- 100 * (exp(int*log(ye)) - 1)
ttw <- 100 * (exp(int*log(x * yw)) - 1)
tte <- 100 * (exp(int*log(x * ye)) - 1)
colMeans(cbind(t=t, tw=tw, ttw=ttw, te=te, tte=tte))
}

## pop size & trend
data.frame(ids[1:6,1:3],
    round(cbind(N=sapply(ttt13, function(z) mean(colSums(z)/10^6)),
    t(sapply(1:6, fun))), 3))


summary(cbind(t, tw, ttw))
summary(cbind(t, te, tte))s
boxplot(cbind(t, tw, ttw, NA, t, te, tte))
abline(h=0,col=2)




lapply(Mid, "[[", "Clim")

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
est <- getEst(res, stage=NULL)
## decadal
tw <- 100 * (exp(est[,"YR"]) - 1)
te <- 100 * (exp(est[,"YR"] + est[,"EWE:YR"]) - 1)
## yearly
tw <- 100 * (exp(0.1*est[,"YR"]) - 1)
te <- 100 * (exp(0.1*(est[,"YR"] + est[,"EWE:YR"])) - 1)

summary(tw)
summary(te)
round(cbind(West=c(Mean=mean(tw), quantile(tw, c(0.025, 0.975))),
    East=c(Mean=mean(te), quantile(te, c(0.025, 0.975)))),2)

