library(mefa4)
ROOT <- "c:/bam/May2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")

spp <- "CAWA"
fid <- 3

fl <- c("analysis_package_gfwfire-nalc-2015-07-24.Rdata",
    "analysis_package_gfwfire-eosd-2015-07-24.Rdata",
    "analysis_package_gfwfire-lcc-2015-07-24.Rdata",
    "analysis_package_fire-nalc-2015-07-24.Rdata")
e <- new.env()
load(file.path(ROOT, "out", "data", fl[fid]), envir=e)

Terms <- getTerms(mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:5000,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

mods <- if (fid == 4)
    mods_fire else mods_gfw
fn <- paste0("bam-", fid, "_", spp, ".Rdata", sep="")
load(file.path(ROOT, "out", "results", fn))

sum(getOK(res)) / length(res)

## need to load data for xn, Xn

est <- getEst(res)

getCaic(res)
printCoefmat(getSummary(res))

getMidPure(res, mods)
getFancyMid(res, mods)
getFancyMidTab(res, mods)

getFancyModsTab(mods)

plotMid(res, mods, web=TRUE)


## --

aic <- cbind(aic1, aic2, aic3)
table(apply(aic, 1, which.min))
