## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2016"

## Load required packages
library(MASS)
library(mefa4)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2016-03-02.Rdata"))

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","YEAR","PCODE","X","Y")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])

## models to consider
NAMES <- list(
    "0"="INTERCEPT",
    "1"=c("INTERCEPT", "JDAY"),
    "2"=c("INTERCEPT", "TSSR"),
    "3"=c("INTERCEPT", "JDAY", "JDAY2"),
    "4"=c("INTERCEPT", "TSSR", "TSSR2"),
    "5"=c("INTERCEPT", "JDAY", "TSSR"),
    "6"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR"),
    "7"=c("INTERCEPT", "JDAY", "TSSR", "TSSR2"),
    "8"=c("INTERCEPT", "JDAY", "JDAY2", "TSSR", "TSSR2"),
    "9"=c("INTERCEPT", "TSLS"),
    "10"=c("INTERCEPT", "TSLS", "TSLS2"),
    "11"=c("INTERCEPT", "TSLS", "TSSR"),
    "12"=c("INTERCEPT", "TSLS", "TSLS2", "TSSR"),
    "13"=c("INTERCEPT", "TSLS", "TSSR", "TSSR2"),
    "14"=c("INTERCEPT", "TSLS", "TSLS2", "TSSR", "TSSR2"))
ff <- list(
    ~ 1,
    ~ JDAY,
    ~ TSSR,
    ~ JDAY + I(JDAY^2),
    ~ TSSR + I(TSSR^2),
    ~ JDAY + TSSR,
    ~ JDAY + I(JDAY^2) + TSSR,
    ~ JDAY + TSSR + I(TSSR^2),
    ~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2),
    ~ TSLS,
    ~ TSLS + I(TSLS^2),
    ~ TSLS + TSSR,
    ~ TSLS + I(TSLS^2) + TSSR,
    ~ TSLS + TSSR + I(TSSR^2),
    ~ TSLS + I(TSLS^2) + TSSR + I(TSSR^2))
names(ff) <- names(NAMES)

## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL

## map
library(rworldmap)
rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
X0 <- pkDur[rn,]
D <- ltdur$end[match(X0$DURMETH, rownames(ltdur$end)),]
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(pkDur[, c("X","Y")], pch=".",
    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255))
points(X0[rowSums(!is.na(D)) > 1, c("X","Y")], pch=19,
    col="red", cex=0.3)


e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_rem.rda"), envir=e)
.BAMCOEFSrem <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS

## species where rem model sample size is at least 25

compare_sets(names(.BAMCOEFSrem$sra_n), names(.BAMCOEFSmix$sra_n))

SPPfull <- sort(names(.BAMCOEFSrem$sra_n)[.BAMCOEFSrem$sra_n >= 25])
SPPfull <- SPPfull[!(SPPfull %in% c("CBCH","CORE","PSFL","RBSA"))]

sptab <- .BAMCOEFSrem$spp_table[SPPfull,]
sptab$nfull <- .BAMCOEFSrem$sra_n[SPPfull]

SPPmix <- sort(names(.BAMCOEFSmix$sra_n)[.BAMCOEFSmix$sra_n >= 25])
SPPmix <- SPPmix[!(SPPmix %in% c("BOBO","CLSW","DUFL","SAPH","STGR","VGSW"))]

sptab$model <- factor("rem", c("rem","mix","both"))
sptab[SPPmix, "model"] <- "both"

SPP <- rownames(sptab)[sptab$model=="both"]

cfall0 <- sapply(SPPfull, function(spp) {
    exp(unname(.BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients))
    })
dim(cfall0) <- c(length(SPPfull), 1)
dimnames(cfall0) <- list(SPPfull, "phi_0")

cfallb <- t(sapply(SPP, function(spp) {
    tmp <- unname(.BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients)
    c(phi_b=exp(tmp[1]), c=plogis(tmp[2]))
    }))

sptab$M0_phi <- cfall0[,1]
sptab$Mb_phi <- cfallb[match(SPPfull, SPP),"phi_b"]
sptab$Mb_c <- cfallb[match(SPPfull, SPP),"c"]

## 

plot(X0$JDAY, X0$TSSR, pch=19, cex=1.2, col=rgb(0,0,0,0.1))
contour(kde2d(X0$JDAY, X0$TSSR), add=TRUE, col=2)

Projects <- table(droplevels(X0$PCODE))
SppOcc <- sapply(SPPfull, function(z) sum(rowSums(xtDur[[z]][rn,], na.rm=TRUE)>0)/length(rn))
SppYmean <- sapply(SPPfull, function(z) mean(rowSums(xtDur[[z]][rn,], na.rm=TRUE)))


## sra 0 vs b models

## aic is really AICc to account for small sample sizes
aic0 <- .BAMCOEFSrem$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
df0 <- matrix(.BAMCOEFSrem$sra_df, nrow(aic0), 15, byrow=TRUE)
dfb <- matrix(.BAMCOEFSmix$sra_df, nrow(aicb), 15, byrow=TRUE)
n0 <- .BAMCOEFSrem$sra_n
nb <- .BAMCOEFSmix$sra_n

aicc0 <- aic0 + (2*df0*(df0+1)) / (n0-df0-1)
aiccb <- aicb + (2*dfb*(dfb+1)) / (nb-dfb-1)
colnames(aicc0) <- paste0("m0_", colnames(aic0))
colnames(aiccb) <- paste0("mb_", colnames(aicb))

aic0 <- aicc0[SPPfull,]
aicb <- aiccb[SPP,]
aic <- cbind(aicc0[SPP,], aiccb[SPP,])

waic0 <- t(apply(aic0, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicb <- t(apply(aicb, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waic <- t(apply(aic, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
best <- colnames(aic)[apply(aic, 1, which.min)]
names(best0) <- SPPfull
names(bestb) <- names(best) <- SPP

par(mfrow=c(1,3))
plot(sptab[SPPfull, "nfull"], waic0[,1], log="x", ylim=c(0,1), pch=ifelse(best0=="0", "o", "+"))
plot(sptab[SPP, "nfull"], waicb[,1], log="x", ylim=c(0,1), pch=ifelse(bestb=="0", "o", "+"))
plot(sptab[SPP, "nfull"], waic[,"m0_0"] + waic[,"mb_0"], log="x", ylim=c(0,1),
    pch=ifelse(best %in% c("m0_0", "mb_0"), "o", "+"))

MAX <- 5000
nn <- 25:MAX
ww <- sapply(nn, function(z) mean((rowSums(waic[,grepl("_0", colnames(waic))]))[np >= z]))
ww2 <- sapply(nn, function(z) mean((rowSums(waic[,grepl("mb_", colnames(waic))]))[np >= z]))

par(mfrow=c(1,2))
plot(nn, 100*(1-ww), type="l", ylim=100*c(0.75, 1), xlab="Number of detections",
    ylab="% time varying", xlim=c(0,MAX))
rug(np)
plot(nn, 100*ww2, type="l", ylim=100*c(0.75, 1), xlab="Number of detections",
    ylab="% mixture", xlim=c(0,MAX))
rug(np)

table(Timevar=!grepl("_0", best), Mixture=grepl("mb_", best))

hasJD <- paste0(rep(c("m0_", "mb_"), each=6), rep(c(1, 3, 5:8), 2))
hasSR <- paste0(rep(c("m0_", "mb_"), each=10), rep(c(2, 4:8, 11:14), 2))
hasLS <- paste0(rep(c("m0_", "mb_"), each=6), rep(c(9:14), 2))

waicx <- waic[SPP,]
bestx <- best[SPP]
best0x <- best0[SPP]
bestbx <- bestb[SPP]
best0bx <- sapply(strsplit(bestx, "_"), "[[", 2)
Support <- t(sapply(as.character(0:14), function(z) {
    c(m0=sum(best0x == z), mb=sum(bestbx == z), Combined=sum(best0bx == z))
}))
Support

## m0 and mb combined here, only spp >=25 det
## constant
sum(grepl("_0", bestx)) * 100 / length(SPP)
## JDAY
sum(bestx %in% hasJD) * 100 / length(SPP)
## TSLS
sum(bestx %in% hasLS) * 100 / length(SPP)
## TSSR
sum(bestx %in% hasSR) * 100 / length(SPP)
## JDAY & TSSR
sum(bestx %in% intersect(hasJD, hasSR)) * 100 / length(SPP)
## TSLS & TSSR
sum(bestx %in% intersect(hasLS, hasSR)) * 100 / length(SPP)
## JDAY or TSLS
sum(bestx %in% c(hasJD, hasLS)) * 100 / length(SPP)
## JDAY or TSLS & TSSR
sum(bestx %in% intersect(c(hasJD, hasLS), hasSR)) * 100 / length(SPP)


sptab$M0_best <- best0
sptab$Mb_best <- bestb[match(SPPfull, SPP)]
sptab$Both_best <- best[match(SPPfull, SPP)]
sptab$Occ <- SppOcc
sptab$Ymean <- SppYmean
write.csv(sptab, row.names=FALSE, file=file.path(ROOT2, "spptab.csv"))

tt <- seq(0, 11, len=1000)
mat0 <- sapply(SPPfull, function(z) 1-exp(-tt*cfall0[z, 1]))
matb <- sapply(SPP, function(z) 1-cfallb[z, "c"]*exp(-tt*cfallb[z, "phi_b"]))

par(mfrow=c(1,2))
matplot(tt, mat0, type="l", lty=1, main="A",
    xlab="Point count duration (minutes)", ylab="P(availability)",
    col=rgb(50,50,50,50,maxColorValue=255), lwd=2)
abline(v=c(3,5,10))
matplot(tt, matb, type="l", lty=1, main="B",
    xlab="Point count duration (minutes)", ylab="P(availability)",
    col=rgb(50,50,50,50,maxColorValue=255), lwd=2)
abline(v=c(3,5,10))

save(sptab, Projects, Support,
    file=file.path(ROOT2, "this-and-that.Rdata"))











load(file.path(ROOT2, "var-bias-res.Rdata"))
res2 <- res[!sapply(res, inherits, "try-error")]

aaa <- data.frame(Var=as.numeric(sapply(res2, "[[", "Var")),
    MSE=as.numeric(sapply(res2, "[[", "MSE")),
    Bias=as.numeric(sapply(res2, "[[", "Bias")),
    Model=rep(c("0","b","0t","bt"), each=2),
    Duration=c("3","5"),
    Species=as.factor(rep(names(res2), each=8)))
aaa$n <- nob[match(aaa$Species, names(nob))]
aaa$logn <- log(aaa$n)
dim(aaa)
aaa <- aaa[!is.na(aaa$Var) & aaa$n >= 2,]
summary(aaa)
dim(aaa)
aaa$Mixture <- ifelse(aaa$Model %in% c("b","bt"), 1, 0)
aaa$Timevar <- ifelse(aaa$Model %in% c("0t","bt"), 1, 0)

library(lme4)
m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1m <- lmer(Var ~ (Mixture + Timevar + Duration + logn)^2 + (1 | Species), aaa)
m2m <- lmer(Bias ~ (Mixture + Timevar + Duration + logn)^2 + (1 | Species), aaa)
cbind(fixef(m1m), coef(m1))
cbind(fixef(m2m), coef(m2))
## minor diffs: use lm

#m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
#m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn)^2, aaa)
m1 <- step(m1)
m2 <- step(m2)

#m1 <- lm(Var ~ Model + Duration + logn, aaa)
#m2 <- lm(Bias ~ Model + Duration + logn, aaa)
summary(m1)
summary(m2)

a1 <- anova(update(m1, . ~ . + Species))
a1$Perc <- 100 * a1[,"Sum Sq"] / sum(a1[,"Sum Sq"])
a2 <- anova(update(m2, . ~ . + Species))
a2$Perc <- 100 * a2[,"Sum Sq"] / sum(a2[,"Sum Sq"])
a1
a2

ng <- 1:200
sf <- function(x) quantile(x, 0.9, na.rm=TRUE)
maxVar <- cbind(max3_0 = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf(aaa$Var[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxBias <- cbind(max3_0 = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf(aaa$Bias[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))

par(mfrow=c(4,2))
for (i in c(1,3,5,7)+1) {
plot(ng, maxVar[,i], main=paste("Var", colnames(maxVar)[i]),
    type="l", lwd=2, col=2, ylim=c(0, max(maxVar)))
lines(ng, maxVar[,i-1], lty=2, col=2, lwd=2)
plot(ng, maxBias[,i], main=paste("Bias", colnames(maxVar)[i]),
    type="l", lwd=2, col=2, ylim=c(min(maxBias), max(maxBias)))
lines(ng, maxBias[,i-1], lty=2, col=2, lwd=2)
}

## =============================================================================
## species specific predictions plots



## =============================================================================
## when does model fail?

load(file.path(ROOT2, "problem.Rdata"))
ss <- c("0 observation with multiple duration (1)", "1 observation with multiple duration (2)")

dat0 <- problem0[!(problem0$msg %in% ss),]
dat0$okfit <- ifelse(is.na(dat0$logphi), 0, 1)
dat0$okse <- ifelse(is.na(dat0$se_logphi), 0, 1)
dat0$msg <- NULL

datb <- problemb[!(problemb$msg %in% ss),]
datb$okfit <- ifelse(is.na(datb$logphi), 0, 1)
datb$okse <- ifelse(is.na(datb$se_logphi), 0, 1)
datb$msg <- NULL

dat0$occ <- sptab$Occ[match(dat0$spp, sptab$spp)]
dat0$ym <- sptab$Ymean[match(dat0$spp, sptab$spp)]
datb$occ <- sptab$Occ[match(datb$spp, sptab$spp)]
datb$ym <- sptab$Ymean[match(datb$spp, sptab$spp)]

summary(m0x <- glm(okfit ~ log(ndet+1)*log(ym), dat0, family=binomial))
summary(mbx <- glm(okfit ~ log(ndet+1)*log(ym), datb, family=binomial))
summary(m0xs <- glm(okse ~ log(ndet+1)*log(ym), dat0, family=binomial))
summary(mbxs <- glm(okse ~ log(ndet+1)*log(ym), datb, family=binomial))


mod00 <- glm(okfit ~ log(ndet+1), dat0, family=binomial)
mod0se0 <- glm(okse ~ log(ndet+1), dat0, family=binomial)
modb0 <- glm(okfit ~ log(ndet+1), datb, family=binomial)
modbse0 <- glm(okse ~ log(ndet+1), datb, family=binomial)

mod01 <- glm(okfit ~ log(n1+1), dat0, family=binomial)
mod0se1 <- glm(okse ~ log(n1+1), dat0, family=binomial)
mod02 <- glm(okfit ~ log(n2+1), dat0, family=binomial)
mod0se2 <- glm(okse ~ log(n2+1), dat0, family=binomial)

modb1 <- glm(okfit ~ log(n1+1), datb, family=binomial)
modbse1 <- glm(okse ~ log(n1+1), datb, family=binomial)
modb2 <- glm(okfit ~ log(n2+1), datb, family=binomial)
modbse2 <- glm(okse ~ log(n2+1), datb, family=binomial)


op <- par(mfrow=c(1,3))
for (i in 1:3) {
if (i == 1) {
    mod0 <- mod00
    modb <- modb0
    mod0se <- mod0se0
    modbse <- modbse0
    xlab <- "Number of >0 survey counts"
}
if (i == 2) {
    mod0 <- mod01
    modb <- modb1
    mod0se <- mod0se1
    modbse <- modbse1
    xlab <- "Number of >1 survey counts"
}
if (i == 3) {
    mod0 <- mod02
    modb <- modb2
    mod0se <- mod0se2
    modbse <- modbse2
    xlab <- "Number of >2 survey counts"
}
nmax <- 100
ndat <- data.frame(n=2:nmax,
    p0=plogis(coef(mod0)[1] + coef(mod0)[2]*log(1 + 2:nmax)),
    pb=plogis(coef(modb)[1] + coef(modb)[2]*log(1 + 2:nmax)),
    p0se=plogis(coef(mod0se)[1] + coef(mod0se)[2]*log(1 + 2:nmax)),
    pbse=plogis(coef(modbse)[1] + coef(modbse)[2]*log(1 + 2:nmax)))

plot(ndat[,1], ndat[,2], col=2, type="l", ylim=c(0,1), lwd=2,
    xlab=xlab, ylab="Probability",
    xlim=c(0,nmax))
lines(ndat[,1], ndat[,3], col=4, lwd=2)
lines(ndat[,1], ndat[,4], col=2, lwd=2, lty=2)
lines(ndat[,1], ndat[,5], col=4, lwd=2, lty=2)
abline(h=0.9, lty=1)
abline(v=ndat[which.min(abs(ndat[,2]-0.9)),1], col=2)
abline(v=ndat[which.min(abs(ndat[,3]-0.9)),1], col=4)
legend("bottomright", col=c(2,2,4,4), lty=c(1,2,1,2), lwd=2,
    legend=c("m0 fit", "m0 SE", "mb fit", "mb SE"), bty="n")
}
par(op)

dat0se <- dat0[!is.na(dat0$se_logphi),]
datbse <- datb[!is.na(datb$se_logphi),]

points(se_logphi ~ nobs, dat0se, ylim=c(0,50), xlim=c(0, 100),
    pch=19, col=rgb(50,50,50,50,maxColorValue=255))

plot(se_logphi ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(se_logitc ~ nobs, datbse, ylim=c(0,1000), xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
plot(abs(cor) ~ nobs, datbse, xlim=c(0, 500),
    pch=19, col=rgb(0,0,255,50,maxColorValue=255))
hist(datbse$cor)



## compare PIF time adjustment

sptab <- read.csv(file.path(ROOT2, "spptab.csv"))
rownames(sptab) <- sptab$spp

pif <- read.csv(file.path(ROOT2, "popContinental_v2_22-May-2013.csv"))
pif <- pif[,c("Common.Name","Scientific.Name","Time.Adjust")]

compare_sets(sptab$common_name, pif$Common.Name)
compare_sets(sptab$scientific_name, pif$Scientific.Name)

sptab[sptab$common_name %in% setdiff(sptab$common_name, pif$Common.Name),]

sptab$tadj <- pif$Time.Adjust[match(sptab$common_name, pif$Common.Name)]
sptab$p30 <- 1-exp(-3*sptab$M0_phi)
sptab$p3b <- 1-sptab$Mb_c*exp(-3*sptab$Mb_phi)

with(sptab, plot(p30, p3b, cex=0.2+0.02*sqrt(sptab$nfull)))
abline(0,1)

with(sptab, plot(p30, tadj, cex=0.2+0.02*sqrt(sptab$nfull)))
with(sptab, plot(p3b, tadj, cex=0.2+0.02*sqrt(sptab$nfull)))

with(sptab, plot(1/p30, tadj, cex=0.2+0.02*sqrt(sptab$nfull),
    xlim=c(1,5), ylim=c(1,5)))
abline(0,1)
with(sptab, plot(1/p3b, tadj, cex=0.2+0.02*sqrt(sptab$nfull),
    xlim=c(1,5), ylim=c(1,5)))
abline(0,1)

## m0t mbt prediction

vjd <- seq(0.38, 0.52, len=200)
vsr <- seq(-0.13, 0.36, len=200)
df <- expand.grid(JDAY=vjd, TSSR=vsr)
#df$JDAY2 <- df$JDAY^2
#df$TSSR2 <- df$TSSR^2
#df$Jday <- df$JDAY * 365
#df$Tssr <- df$TSSR * 24

X <- model.matrix(ff[["8"]], df)
sppCor <- list()

pdf(file.path(ROOT2, "spp-pred.pdf"), onefile=TRUE, width=10, height=5)

#spp <- "BBCU"
#fff <- ff[["8"]]
for (spp in SPP) {

cf0 <- .BAMCOEFSrem$sra_estimates[[spp]][["8"]]$coefficients
cfb <- .BAMCOEFSmix$sra_estimates[[spp]][["8"]]$coefficients

p30 <- 1-exp(-3*exp(drop(X %*% cf0)))
p3b <- 1-plogis(drop(X %*% cfb[-1]))*exp(-3*exp(cfb[1]))

sppCor[[spp]] <- cor(cbind(p30,p3b))[1,2]
cat(spp, "cor =", sppCor[[spp]], "\n");flush.console()

#plot(p30, p3b)
z0 <- matrix(p30, length(vjd), length(vsr))
zb <- matrix(p3b, length(vjd), length(vsr))
#pmax <- 1#0.2*max(p30, p3b)
Col0 <- grey(0.2+0.8*seq(1-min(z0), 1-max(z0), len=24))
Colb <- grey(0.2+0.8*seq(1-min(zb), 1-max(zb), len=24))
op <- par(mfrow=c(1,2), las=1)
image(vjd*365, vsr*24, z0,
    col = Col0,
    xlab="Julian days", 
    ylab="Hours since sunrise",
    main=paste(spp, "M0t"))
contour(vjd*365, vsr*24, z0, add=TRUE, col=1, labcex=1)
image(vjd*365, vsr*24, zb,
    col = Colb,
    xlab="Julian days", 
    ylab="Hours since sunrise",
    main=paste(spp, "Mbt"))
contour(vjd*365, vsr*24, zb, add=TRUE, col=1, labcex=1)
#plot(p30, p3b, main=round(sppCor[[spp]], 3))
par(op)
}
dev.off()
