library(mefa4)

source("~/repos/mep/R/diagnostics-functions.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-models.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")
load("e:/peter/bam/bcr4/bcr4-data.RData")

Terms <- getTerms(mods, "list")
setdiff(Terms, colnames(DAT))
xn <- DAT[BB[,1],Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))

spp <- "RCKI"
load(paste0("e:/peter/bam/bcr4/results/", spp, ".RData"))
100 * sum(getOK(out)) / length(out)
est1 <- getEst(out, stage = 2, X=Xn)
est2 <- getEst(out, stage = 4, X=Xn)
m <- out[[1]]$best

## explore the 1st run
summary(m)
#map(m)
boxplot(fitted(m) ~ hab, xn, range=0, ylab="Density")

## checking results
#getMid(out, mods)
getFancyMidTab(out, mods)

par(mfrow=c(1,3))
plotMid(out, mods)
visualize_road()
visualize_ysd()

