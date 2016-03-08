library(mefa4)
ROOT <- "c:/bam/May2015"
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/bamanalytics/R/analysis_functions.R")

load(file.path(ROOT, "out", "data", "pack_gfw_nam_nlc_2015-09-02.Rdata"))


table(DAT$ROAD, DAT$HAB)
SPP <- colnames(YY)[colSums(YY>0) > 100 & TAX$Order == "Passeriformes"]

res <- list()
for (spp in SPP) {
    cat(spp, date(), "\n");flush.console()
    m0 <- glm_skeleton(glm(YY[,spp] ~ 1, data=DAT, offset=OFF[,spp],family=poisson), CAICalpha=1)
    m1 <- glm_skeleton(glm(YY[,spp] ~ HAB, data=DAT, offset=OFF[,spp],family=poisson), CAICalpha=1)
    m2 <- glm_skeleton(glm(YY[,spp] ~ HAB + ROAD, data=DAT, offset=OFF[,spp],family=poisson), CAICalpha=1)
    m3 <- glm_skeleton(glm(YY[,spp] ~ HAB * ROAD, data=DAT, offset=OFF[,spp],family=poisson), CAICalpha=1)
    res[[spp]] <- list(m0=m0, m1=m1, m2=m2, m3=m3)
}


fid <- 1
fl <- c("analysis_package_gfwfire-nalc-2015-08-14.Rdata",
    "analysis_package_gfwfire-eosd-2015-08-14.Rdata",
    "analysis_package_gfwfire-lcc-2015-08-14.Rdata",
    "analysis_package_fire-nalc-2015-08-14.Rdata")
load(file.path(ROOT, "out", "data", fl[fid]))
load(file.path(ROOT, "out", "analysis_package_distances.Rdata"))

Ymean <- as.matrix(t(groupMeans(YY[,SPP], 1, DAT$ROAD)))
rownames(Ymean) <- SPP
colnames(Ymean) <- c("Off-road","Roadside")
Beta <- sort(log(Ymean[,2] / Ymean[,1]))
plot(Beta, 1:length(Beta), type="n", axes=FALSE, ann=FALSE)
abline(v=0)
text(Beta, 1:length(Beta), names(Beta), cex=0.7, col=ifelse(Beta<0, 4,2))
axis(1)
title(xlab="Roadside bias")


Ymean <- as.matrix(t(groupMeans(YY[,SPP], 1, interaction(DAT$ROAD, DAT$HAB))))
Ymean <- Ymean[,sort(colnames(Ymean))]
BetaH <- log(Ymean[,10:18]/Ymean[,1:9])
BetaH <- BetaH[names(Beta),]
BetaH[!is.finite(BetaH)] <- NA
BetaH <- BetaH[,order(apply(BetaH, 2, median, na.rm=TRUE))]
colnames(BetaH) <- gsub("1.", "", colnames(BetaH), fixed=TRUE)

library(plotrix)
boxplot(BetaH)
abline(h=0)

par(las=1)
ladderplot(BetaH[rowSums(is.na(BetaH))==0,], vertical=FALSE, col=3)
abline(v=0)
