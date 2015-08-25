library(mefa4)
ROOT <- "c:/bam/May2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")

spp <- "CAWA"
folder <- "results3"
fid <- 1

fl <- c("analysis_package_gfwfire-nalc-2015-08-17.Rdata",
    "analysis_package_gfwfire-eosd-2015-08-17.Rdata",
    "analysis_package_gfwfire-lcc-2015-08-17.Rdata",
    "analysis_package_fire-nalc-2015-08-17.Rdata")
e <- new.env()
load(file.path(ROOT, "out", "data", fl[fid]), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:1000,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

#mods <- if (fid == 4)
#    mods_fire else mods_gfw
fn <- paste0("bam-", fid, "_", spp, ".Rdata", sep="")
load(file.path(ROOT, "out", folder, fn))

sum(getOK(res)) / length(res)

## need to load data for xn, Xn

est <- getEst(res, stage=NULL)

getCaic(res)
printCoefmat(getSummary(res))

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
