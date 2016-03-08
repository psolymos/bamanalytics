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

save(res, file=file.path(ROOT, "out", "data", "roadside-bias-estimates.Rdata"))

aic <- t(sapply(res, function(z) sapply(z, "[[", "aic")))
best <- find_min(aic)
table(best$index)

beta <- sapply(res, function(z) z$m2$coef["ROAD"])
names(beta) <- names(res)

betaH0 <- t(sapply(res, function(z) z$m3$coef))
betaH <- cbind(betaH0[,"ROAD"], betaH0[,"ROAD"]+betaH0[,grepl(":", colnames(betaH0))])
colnames(betaH) <- levels(DAT$HAB)
betaH <- betaH[,order(apply(betaH, 2, mean))]

Ymean <- as.matrix(t(groupMeans(YY[,SPP], 1, DAT$ROAD)))
rownames(Ymean) <- SPP
colnames(Ymean) <- c("Off-road","Roadside")
Beta <- sort(log(Ymean[,2] / Ymean[,1]))
beta <- beta[names(Beta)]
Col <- rep("green", length(beta))
Col[Beta < 0 & beta < 0] <- "blue"
Col[Beta > 0 & beta > 0] <- "red"

plot(Beta, 1:length(Beta), type="n", axes=FALSE, ann=FALSE, xlim=range(Beta, beta))
abline(v=0)
text(Beta, 1:length(Beta), names(Beta), cex=0.7, col=Col)
axis(1)
title(xlab="Roadside bias")
points(beta, 1:length(beta), pch=4, col=Col)

plot(beta[SPP], Beta[SPP])


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

betaH <- betaH[rownames(BetaH),]
betaH <- betaH[,order(apply(betaH, 2, mean))]

par(las=1, mfrow=c(1,2))
ladderplot(BetaH[rowSums(is.na(BetaH))==0,], vertical=FALSE, col=3, xlab="Roadside bias",
    ylim=c(-5,5))
abline(v=0)
points(apply(BetaH[rowSums(is.na(BetaH))==0,], 2, median), 1:ncol(BetaH), pch="|", cex=2.5, col=2)

ladderplot(betaH, vertical=FALSE, col=3, xlab="Roadside bias", 
    ylim=c(-5,5))
abline(v=0)
points(apply(betaH, 2, mean), 1:ncol(betaH), pch="|", cex=2.5, col=2)
