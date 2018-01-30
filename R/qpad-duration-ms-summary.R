## Define root folder where data are stored
#ROOT <- "c:/bam/May2015"
#ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2016"
ROOT <- "e:/peter/bam/May2015"
#ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2017"
#ROOT <- "~/Dropbox/Public"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionMarch2017"

## Load required packages
library(MASS)
library(mefa4)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2017-03-01.Rdata"))

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","YEAR","PCODE","X","Y","SS")]
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
png(file.path(ROOT2, "tabfig", "Fig1_map.png"), height=600, width=800)
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(pkDur[, c("X","Y")], pch=".",
    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255))
points(X0[rowSums(!is.na(D)) > 1, c("X","Y")], pch=19,
    col="red", cex=0.3)
dev.off()

ld <- data.frame(table(pkDur$DURMETH))
rownames(ld) <- ld[,1]
ld <- data.frame(ld, ltdur$end[rownames(ld),])
ld <- ld[rowSums(!is.na(ltdur$end[rownames(ld),])) > 1,]

datx <- droplevels(dat[rownames(X0)[rowSums(!is.na(D)) > 1],])

## data summaries
pkDurOK <- droplevels(pkDur[pkDur$DURMETH %in% rownames(ld),])
str(pkDurOK)

if (FALSE) {
## export data for Fig 1
ss <- droplevels(nonDuplicated(pkDurOK, SS))
#ss <- droplevels(pkDurOK)
ss <- ss[,c("SS","PCODE","DURMETH","X","Y")]

ll <- apply(ltdur$end[rownames(ld),], 1, function(z) {
    paste(c(0, z[!is.na(z)]), collapse="-")
    })
ss$INTERVALS <- as.factor(ll[match(ss$DURMETH, names(ll))])
l <- data.frame(table(ss$INTERVALS))
l$Legend <- paste0(as.character(l$Var1), " (n=", l$Freq, ")")

write.csv(ss, row.names=FALSE, file="~/Downloads/removal-ms-points-for-fig-1.csv")
}

## species
e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_rem.rda"), envir=e)
.BAMCOEFSrem <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_duration_fmix.rda"), envir=e)
.BAMCOEFSfmix <- e$.BAMCOEFS

## species where rem model sample size is at least 25

## check here !!!
sb <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")

rownames(sb) <- sb$Species_ID

compare_sets(names(.BAMCOEFSrem$sra_n), names(.BAMCOEFSmix$sra_n))

NMIN <- 75

SPPfull <- sort(names(.BAMCOEFSrem$sra_n)[.BAMCOEFSrem$sra_n >= NMIN])
table(sb[SPPfull, "Singing_birds"])
## this comes from checking M0/Mb estimates (all BAM scale)
EXC1 <- c("CBCH", "CORE", "PAWR", "PSFL", "RBSA")
SPPfull <- SPPfull[!(SPPfull %in% EXC1)]


sptab <- .BAMCOEFSrem$spp_table[SPPfull,]
sptab$nfull <- .BAMCOEFSrem$sra_n[SPPfull]

SPPmix <- sort(names(.BAMCOEFSmix$sra_n)[.BAMCOEFSmix$sra_n >= NMIN])
EXC2 <- c("BEKI", "BOBO",
    "BRBL", "CLSW", "DUFL", "PAWR", "PIGR", "RBWO",
    "RNPH", "ROPI", "SAPH", "STGR", "VGSW")
SPPmix <- SPPmix[!(SPPmix %in% EXC2)]

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

#plot(pkDurOK$JDAY, pkDurOK$TSSR, pch=19, cex=1.2, col=rgb(0,0,0,0.1))
#contour(kde2d(pkDurOK$JDAY, pkDurOK$TSSR), add=TRUE, col=2)

Projects <- table(droplevels(X0$PCODE))
SppOcc <- sapply(SPPfull, function(z) sum(rowSums(xtDur[[z]][rn,], na.rm=TRUE)>0)/length(rn))
SppYmean <- sapply(SPPfull, function(z) mean(rowSums(xtDur[[z]][rn,], na.rm=TRUE)))


## sra 0 vs b (or f) models

## aic is really AICc to account for small sample sizes
aic0 <- .BAMCOEFSrem$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
aicf <- .BAMCOEFSfmix$sra_aic
df0 <- matrix(.BAMCOEFSrem$sra_df, nrow(aic0), 15, byrow=TRUE)
dfb <- matrix(.BAMCOEFSmix$sra_df, nrow(aicb), 15, byrow=TRUE)
dff <- matrix(.BAMCOEFSfmix$sra_df, nrow(aicb), 15, byrow=TRUE)
n0 <- .BAMCOEFSrem$sra_n
nb <- .BAMCOEFSmix$sra_n
nf <- .BAMCOEFSfmix$sra_n

aicc0 <- aic0 + (2*df0*(df0+1)) / (n0-df0-1)
aiccb <- aicb + (2*dfb*(dfb+1)) / (nb-dfb-1)
aiccf <- aicf + (2*dff*(dff+1)) / (nf-dff-1)
colnames(aicc0) <- paste0("m0_", colnames(aic0))
colnames(aiccb) <- paste0("mb_", colnames(aicb))
colnames(aiccf) <- paste0("mf_", colnames(aicf))

aic0 <- aicc0[SPPfull,]
aicb <- aiccb[SPP,]
aicf <- aiccf[SPP,]
aic <- cbind(aicc0[SPP,], aiccb[SPP,])
aicx <- cbind(aicc0[SPP,], aiccb[SPP,], aiccf[SPP,])
aicbf <- cbind(aiccb[SPP,], aiccf[SPP,])

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
waicf <- t(apply(aicf, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waic <- t(apply(aic, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicx <- t(apply(aicx, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicbf <- t(apply(aicbf, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waic0b <- waic[,1:15] + waic[,16:30]
colnames(waic0b) <- colnames(waic0)

H <- apply(waic, 1, function(z) sum(z^2))
H0 <- apply(waic0, 1, function(z) sum(z^2))
Hb <- apply(waicb, 1, function(z) sum(z^2))
H0b <- apply(waic0b, 1, function(z) sum(z^2))

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
bestf <- as.character(0:14)[apply(aicf, 1, which.min)]
best <- colnames(aic)[apply(aic, 1, which.min)]
bestx <- as.factor(colnames(aicx))[apply(aicx, 1, which.min)]
bestbf <- as.factor(colnames(aicbf))[apply(aicbf, 1, which.min)]
names(best0) <- SPPfull
names(bestb) <- names(best) <- names(bestf) <- names(bestx) <- names(bestbf) <- SPP

cbind(Mb=table(bestbf)[1:15], Mf=table(bestbf)[16:30])

par(mfrow=c(1,3))
plot(sptab[SPPfull, "nfull"], waic0[,1], log="x", ylim=c(0,1), pch=ifelse(best0=="0", "o", "+"))
plot(sptab[SPP, "nfull"], waicb[,1], log="x", ylim=c(0,1), pch=ifelse(bestb=="0", "o", "+"))
plot(sptab[SPP, "nfull"], waic[,"m0_0"] + waic[,"mb_0"], log="x", ylim=c(0,1),
    pch=ifelse(best %in% c("m0_0", "mb_0"), "o", "+"))

MAX <- 5000
nn <- 25:MAX
ww1 <- sapply(nn, function(z) mean((rowSums(waic[,grepl("_0",
    colnames(waic))]))[sptab[rownames(waic), "nfull"] >= z]))
ww2 <- sapply(nn, function(z) mean((rowSums(waic[,grepl("mb_",
    colnames(waic))]))[sptab[rownames(waic), "nfull"] >= z]))
ww3 <- sapply(nn, function(z) mean((rowSums(waic[,!grepl("_0",
    colnames(waic)) & grepl("mb_", colnames(waic))]))[sptab[rownames(waic), "nfull"] >= z]))

#png(file.path(ROOT2, "tabfig", "Fig4_weights.png"), width=800, height=450)
pdf(file.path(ROOT2, "tabfig", "Fig4_weights.pdf"), width=6, height=6)
op <- par(mfrow=c(1,1), las=1)
plot(nn, 100*(1-ww1), type="l", ylim=100*c(0.8, 1), xlab="Number of detections",
    ylab="Percent of species", xlim=c(0,MAX), lwd=2, col=1)
lines(nn, 100*ww2, lwd=2, col=1, lty=2)
lines(nn, 100*ww3, lwd=2, col=1, lty=3)
rug(sptab[rownames(waic), "nfull"])
#abline(v=nn[which(ww2 > 0.9)[1]], lty=2)
#abline(h=90, lty=3)
legend("bottomright", bty="n", lwd=2, lty=1:3,
    legend=c("Time varying", "Mixture", "Time varying mixture"))
par(op)
dev.off()
nn[which(ww1 > 0.9)[1]]
nn[which(ww2 > 0.9)[1]]
nn[which(ww3 > 0.9)[1]]

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
Support <- cbind(Support, w0=colMeans(waic0),
    wb=colMeans(waicb), w0b=colMeans(waic0b))
Support
#write.csv(Support, file=file.path(ROOT2, "tabfig", "support.csv"))

(best0bxx <- table(sapply(strsplit(bestx, "_"), "[[", 2),
    sapply(strsplit(bestx, "_"), "[[", 1)))
colSums(best0bxx)

t0b <- table(M0=best0x,Mb=bestbx)
addmargins(t0b[as.character(0:14),as.character(0:14)])

t0b2 <- table(Mb=best0bx,M0b=bestbx)
addmargins(t0b2[as.character(0:14),as.character(0:14)])


zzz <- table(bestx)
sum(zzz[grep("mb_", names(zzz))])
sum(zzz[grep("m0_", names(zzz))])

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

mean(rowSums(waic[,hasJD]))
mean(rowSums(waic[,hasSR]))
mean(rowSums(waic[,hasLS]))
mean(rowSums(waic[,union(hasJD, hasLS)]))
mean(rowSums(waic[,union(hasSR, hasJD)]))
mean(rowSums(waic[,union(hasSR, hasLS)]))
mean(rowSums(waic[,union(union(hasSR, hasJD), hasLS)]))

sptab$M0_best <- best0
sptab$Mb_best <- bestb[match(SPPfull, SPP)]
sptab$Both_best <- best[match(SPPfull, SPP)]
sptab$Occ <- SppOcc
sptab$Ymean <- SppYmean

## check if values are sensible
rownames(sptab)[sptab$M0_phi < 0.01]
rownames(sptab)[!is.na(sptab$Mb_phi) & sptab$Mb_phi < 0.01]
rownames(sptab)[!is.na(sptab$Mb_phi) & sptab$Mb_c < 0.01]


tt <- seq(0, 11, len=1000)
mat0 <- sapply(SPPfull, function(z) 1-exp(-tt*cfall0[z, 1]))
matb <- sapply(SPP, function(z) 1-cfallb[z, "c"]*exp(-tt*cfallb[z, "phi_b"]))

#png(file.path(ROOT2, "tabfig", "Fig2_probs.png"), width=800, height=450)
pdf(file.path(ROOT2, "tabfig", "Fig2_probs.pdf"), width=10, height=5)
op <- par(mfrow=c(1,2), las=1)
matplot(tt, mat0, type="l", lty=1, main="",
    xlab="Point count duration (minutes)", ylab="P(availability)",
    col=rgb(50,50,50,50,maxColorValue=255), lwd=2)
abline(v=c(3,5,10),lty=3)
text(0.5,0.95,expression(M[0]))
matplot(tt, matb, type="l", lty=1, main="",
    xlab="Point count duration (minutes)", ylab="P(availability)",
    col=rgb(50,50,50,50,maxColorValue=255), lwd=2)
abline(v=c(3,5,10),lty=3)
text(0.5,0.95,expression(M[f]))
par(op)
dev.off()

tt <- c(3,5,10)
mat0 <- t(sapply(SPPfull, function(z) 1-exp(-tt*cfall0[z, 1])))
matb <- t(sapply(SPP, function(z) 1-cfallb[z, "c"]*exp(-tt*cfallb[z, "phi_b"])))

summary(mat0)
summary(matb)

round(apply(mat0, 2, mean), 2)
round(apply(mat0, 2, sd), 3)
round(apply(matb, 2, mean), 2)
round(apply(matb, 2, sd), 3)



## model support
data.frame(w_avg=sort(colMeans(waic0b)))
best0b <- find_max(waic0b)$index
names(best0b) <- rownames(waic0b)
df0b <- (1+c(0,1,1,2,2,2,3,3,4,1,2,2,3,3,4))[best0b]

save(sptab, Projects, Support, waic0, waicb, waic, waic0b, best0, bestb, best, best0b,
    file=file.path(ROOT2, "this-and-that.Rdata"))

plot(jitter(df0b), log(sptab[names(best0b),"nfull"]))
abline(lm(log(sptab[names(best0b),"nfull"]) ~ df0b))

load(file.path(ROOT2, "this-and-that.Rdata"))

lht <- read.csv("~/Dropbox/bam/LifeHistoryToUpdate.csv")
compare_sets(lht$SPECIES, SPPfull)
setdiff(SPPfull, lht$SPECIES)
rownames(lht) <- lht$SPECIES
lht <- droplevels(lht[SPPfull,1:6])
lht[is.na(lht$DATABASE_MIG_TYPE),]

load(file.path(ROOT2, "var-bias-res.Rdata"))
load(file.path(ROOT2, "var-bias-res-2.Rdata"))
load(file.path(ROOT2, "var-bias-res-3.Rdata"))

resx <- res[!sapply(res, inherits, "try-error")] # ymin=1, y3min=0
resx <- res2[!sapply(res2, inherits, "try-error")] # y3min=1
resx <- res3[!sapply(res3, inherits, "try-error")] # aggregate
resx <- res4[!sapply(res4, inherits, "try-error")] # aggregate + estimate

aaa <- data.frame(Var=as.numeric(sapply(resx, "[[", "Var")),
    MSE=as.numeric(sapply(resx, "[[", "MSE")),
    Bias=as.numeric(sapply(resx, "[[", "Bias")),
    Model=rep(c("0","b","0t","bt"), each=2),
    Duration=c("3","5"),
    Species=as.factor(rep(names(resx), each=8)))
aaa$n <- nob[match(aaa$Species, names(nob))]
aaa$logn <- log(aaa$n)
dim(aaa)
aaa <- aaa[!is.na(aaa$Var) & aaa$n >= 2,]
summary(aaa)
dim(aaa)
aaa$Mixture <- ifelse(aaa$Model %in% c("b","bt"), 1, 0)
aaa$Timevar <- ifelse(aaa$Model %in% c("0t","bt"), 1, 0)
#aaa$logBias <- log(aaa$Bias + 0.5)
#aaa$logVar <- log(aaa$Var + 0.5)
aaa$MigT <- lht$DATABASE_MIG_TYPE[match(aaa$Species, lht$SPECIES)]

aaa$phi0 <- sptab$M0_phi[match(aaa$Species, rownames(sptab))]
#aaa$phi0 <- as.integer(cut(aaa$phi0, c(0, 0.2, 0.3, 0.5, 1))) - 1
aaa$phi0 <- (as.integer(cut(aaa$phi0, c(0, 0.2, 0.4, Inf))) - 1) / 2
table(aaa$phi0)

## retain species where both models worked fine
aaa <- droplevels(aaa[aaa$Species %in% SPP,])

m1 <- lm(Var ~ (Mixture + Timevar + Duration)^2 + (logn + phi0)^2, aaa)
m2 <- lm(Bias ~ (Mixture + Timevar + Duration)^2 + (logn + phi0)^2, aaa)

#m1 <- lm(Var ~ Mixture * Timevar + Duration + logn * phi0, aaa)
#m2 <- lm(Bias ~ Mixture * Timevar + Duration + logn * phi0, aaa)
#m1 <- lm(Var ~ (Mixture + Timevar + Duration + logn + phi0)^2, aaa)
#m2 <- lm(Bias ~ (Mixture + Timevar + Duration + logn + phi0)^2, aaa)
#m1 <- step(m1)
#m2 <- step(m2)
summary(m2)
summary(m1)

a1 <- anova(m1)
#a1 <- anova(update(m1, . ~ . + Species))
a1$Perc <- 100 * a1[,"Sum Sq"] / sum(a1[,"Sum Sq"])
a2 <- anova(m2)
#a2 <- anova(update(m2, . ~ . + Species))
a2$Perc <- 100 * a2[,"Sum Sq"] / sum(a2[,"Sum Sq"])
rownames(a1) <- gsub("Duration", "Duration5", rownames(a1))
rownames(a2) <- gsub("Duration", "Duration5", rownames(a2))

#summary(m1)
#a1
#summary(m2)
#a2

RN <- union(rownames(coef(summary(m2))), rownames(a2))
#RN1 <- union(rownames(coef(summary(m1))), rownames(a1))
#RN <- union(RN1, RN2)
tmp <- coef(summary(m2))[match(RN, rownames(coef(summary(m2)))),c(1,2,4)]
rownames(tmp)[nrow(tmp)] <- "Residuals"
bbb <- round(data.frame(
    Bias=tmp,
    Bias.Perc=a2[match(RN, rownames(a2)),"Perc"],
    Var=coef(summary(m1))[match(RN, rownames(coef(summary(m1)))),c(1,2,4)],
    Var.Perc=a1[match(RN, rownames(a1)),"Perc"]), 4)
rownames(bbb) <- RN
bbb
write.csv(bbb, file=file.path(ROOT2, "tabfig", "var-bias-tab.csv"))


## looking at residuals to figure out what species have high/low var/bias
aaa$var_hat <- fitted(m1)
aaa$var_res <- aaa$Var - aaa$var_hat
aaa$bias_hat <- fitted(m2)
aaa$bias_res <- aaa$Bias - aaa$bias_hat

zz <- aggregate(aaa[,c("var_res", "bias_res")], list(Spp=aaa$Species), mean)
rownames(zz) <- zz[,1]
zz <- zz[,-1]
round(zz, 3)[order(abs(zz[,1]), decreasing=TRUE),]

#sptab$Resid_bias <- zz$bias_res[match(rownames(sptab), rownames(zz))]
#sptab$Resid_var <- zz$var_res[match(rownames(sptab), rownames(zz))]


#ng <- 2:400
ng <- c(20, 50, 100, 200, 400)
sf0 <- function(x) quantile(x, 0.5, na.rm=TRUE)
sf1 <- function(x) quantile(x, 0.9, na.rm=TRUE)
sf2 <- function(x) quantile(x, 0.95, na.rm=TRUE)
sf3 <- function(x) quantile(x, 0.05, na.rm=TRUE)
maxVar1 <- cbind(max3_0 = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf2(aaa$Var[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxVar2 <- cbind(max3_0 = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf3(aaa$Var[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxVar0 <- cbind(max3_0 = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf0(aaa$Var[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxBias0 <- cbind(max3_0 = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf0(aaa$Bias[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxBias1 <- cbind(max3_0 = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf2(aaa$Bias[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))
maxBias2 <- cbind(max3_0 = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0" & aaa$n >= z])),
    max5_0 = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0" & aaa$n >= z])),
    max3_b = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "3" & aaa$Model == "b" & aaa$n >= z])),
    max5_b = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "5" & aaa$Model == "b" & aaa$n >= z])),
    max3_0t = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "3" & aaa$Model == "0t" & aaa$n >= z])),
    max5_0t = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "5" & aaa$Model == "0t" & aaa$n >= z])),
    max3_bt = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "3" & aaa$Model == "bt" & aaa$n >= z])),
    max5_bt = sapply(ng, function(z)
        sf3(aaa$Bias[aaa$Duration == "5" & aaa$Model == "bt" & aaa$n >= z])))

par(mfrow=c(4,2))
for (i in c(1,3,5,7)+1) {
plot(ng, maxVar1[,i], main=paste("Var", colnames(maxVar1)[i]),
    type="l", lwd=2, col=2, ylim=c(0, max(maxVar1)))
lines(ng, maxVar1[,i-1], lty=2, col=2, lwd=2)
lines(ng, maxVar0[,i], lty=1, col=4, lwd=2)
lines(ng, maxVar0[,i-1], lty=2, col=4, lwd=2)
plot(ng, maxBias1[,i], main=paste("Bias", colnames(maxBias1)[i]),
    type="l", lwd=2, col=2, ylim=range(maxBias1,maxBias2))
lines(ng, maxBias1[,i-1], lty=2, col=2, lwd=2)
lines(ng, maxBias2[,i], lty=1, col=2, lwd=2)
lines(ng, maxBias2[,i-1], lty=2, col=2, lwd=2)
lines(ng, maxBias0[,i], lty=1, col=4, lwd=2)
lines(ng, maxBias0[,i-1], lty=2, col=4, lwd=2)
abline(h=0)
}

MOD <- c("M0", "M0", "Mb", "Mb", "M0t", "M0t", "Mbt", "Mbt")
#png(file.path(ROOT2, "tabfig", "Fig5_var-bias.png"), width=600, height=4*250)
pdf(file.path(ROOT2, "tabfig", "Fig5_var-bias.pdf"), width=10, height=4*3.5)
par(mfrow=c(4,2), las=1)
for (i in c(1,3,5,7)+1) {
plot(ng, maxBias1[,i], main=paste("Bias", MOD[i]),
    xlab="Sample size", ylab="Bias",
    type="n", lwd=2, col=2,
    ylim=max(abs(maxBias1),abs(maxBias2))*c(-1,1), xlim=c(2,400))
polygon(c(ng, rev(ng)), c(maxBias1[,i-1], rev(maxBias2[,i-1])),
    col="grey", border="grey")
polygon(c(ng, rev(ng)), c(maxBias1[,i], rev(maxBias2[,i])),
    col="black")
lines(ng, maxBias0[,i], lty=1, col="white", lwd=2)
lines(ng, maxBias0[,i-1], lty=2, col="white", lwd=2)
box()

plot(ng, maxVar1[,i], main=paste("Variance", MOD[i]),
    xlab="Sample size", ylab="Variance",
    type="n", lwd=2, col=2, ylim=c(0,max(maxVar1, maxVar2)), xlim=c(2,400))
polygon(c(ng, rev(ng)), c(maxVar1[,i-1], rev(maxVar2[,i-1])),
    col="grey", border="grey")
polygon(c(ng, rev(ng)), c(maxVar1[,i], rev(maxVar2[,i])),
    col="black")
#lines(ng, maxVar0[,i], lty=1, col="white", lwd=2)
#lines(ng, maxVar0[,i-1], lty=2, col="white", lwd=2)
box()
}
dev.off()


pdf(file.path(ROOT2, "tabfig", "Fig5_var-bias2.pdf"), width=10, height=2*5)
op <- par(mfrow=c(2,1), las=1, mar=c(0,4,3,3))

ct <- 6*1:length(ng)-6

w <- 0.3 # width
plot(ct, rep(0, length(ng)), xlab="Sample size", ylab="Bias",
    type="n", axes=FALSE,
    ylim=max(abs(maxBias1),abs(maxBias2))*c(-1,1), xlim=c(ct[1]-2,ct[length(ct)]+2))
box()
axis(2)
#axis(1, ct, ng, tick=FALSE)
for (i in 1:length(ng)) {
    for (j in 1:4) {
        ## 3 min
        c3 <- c("max3_0", "max3_b", "max3_0t", "max3_bt")[j]
        xj <- ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j] - w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxBias1[i,c3], maxBias1[i,c3], maxBias2[i,c3], maxBias2[i,c3]),
            border="grey", col="grey")
        lines(xj+c(-w,w)/2, rep(maxBias0[i,c3], 2), lwd=2, col=1, lend=2)
        ## 5 min
        c5 <- c("max5_0", "max5_b", "max5_0t", "max5_bt")[j]
        xj <- ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j] + w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxBias1[i,c5], maxBias1[i,c5], maxBias2[i,c5], maxBias2[i,c5]),
            border=1, col=1)
        lines(xj+c(-w,w)/2, rep(maxBias0[i,c5], 2), lwd=2, col="grey", lend=2)
        ## text
        if (i == 1) {
            ex <- list(
                expression(M[0]),
                expression(M[f]),
                expression(M[0]^varphi),
                expression(M[f]^c))
            text(ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j],
                maxBias1[i,c("max3_0", "max3_b", "max3_0t", "max3_bt")[j]] + 0.02,
                ex[[j]], cex=1)
        }
    }
}
abline(h=0, lty=2)
legend("topright", fill=c("grey", "black"), legend=c("3-min", "5-min"), bty="n")

par(mar=c(5,4,0,3))
plot(ct, rep(0, length(ng)), xlab="Sample size", ylab="Variance",
    type="n", axes=FALSE,
    ylim=c(0, 1.1*max(maxVar1,maxVar2)), xlim=c(ct[1]-2,ct[length(ct)]+2))
box()
axis(2)
axis(1, ct, ng, tick=FALSE)
for (i in 1:length(ng)) {
    for (j in 1:4) {
        ## 3 min
        c3 <- c("max3_0", "max3_b", "max3_0t", "max3_bt")[j]
        xj <- ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j] - w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxVar1[i,c3], maxVar1[i,c3], maxVar2[i,c3], maxVar2[i,c3]),
            border="grey", col="grey")
        #lines(xj+c(-w,w)/2, rep(maxVar0[i,c3], 2), lwd=3, col=1, lend=2)
        ## 5 min
        c5 <- c("max5_0", "max5_b", "max5_0t", "max5_bt")[j]
        xj <- ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j] + w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxVar1[i,c5], maxVar1[i,c5], maxVar2[i,c5], maxVar2[i,c5]),
            border=1, col=1)
        #lines(xj+c(-w,w)/2, rep(maxVar0[i,c5], 2), lwd=3, col="grey", lend=2)
        ## text
        if (i == 1) {
            ex <- list(
                expression(M[0]),
                expression(M[b]),
                expression(M[0]^varphi),
                expression(M[f]^c))
            text(ct[i] + c(-1.5, -0.5, 0.5, 1.5)[j],
                maxVar1[i,c("max3_0", "max3_b", "max3_0t", "max3_bt")[j]]+0.003,
                ex[[j]], cex=1)
        }
    }
}
#legend("topright", fill=c("grey", "black"), legend=c("3-min", "5-min"), bty="n")
par(op)
dev.off()

## use res4 for aaa
maxBias0[maxBias0 == -1] <- 0
maxBias1[maxBias1 == -1] <- 0
maxBias2[maxBias2 == -1] <- 0

pdf(file.path(ROOT2, "tabfig", "Fig5_var-bias2-est.pdf"), width=10, height=2*5)
op <- par(mfrow=c(2,1), las=1, mar=c(0,4,3,3))

ct <- 3*1:length(ng)-3

w <- 0.4 # width
plot(ct, rep(0, length(ng)), xlab="Sample size", ylab="Bias",
    type="n", axes=FALSE,
    ylim=max(abs(maxBias1),abs(maxBias2))*c(-1,1), xlim=c(ct[1]-2,ct[length(ct)]+2))
box()
axis(2)
for (i in 1:length(ng)) {
    for (j in 1:2) {
        ## 3 min
        c3 <- c("max3_0", "max3_b")[j]
        xj <- ct[i] + c(-0.5, 0.5)[j] - w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxBias1[i,c3], maxBias1[i,c3], maxBias2[i,c3], maxBias2[i,c3]),
            border="grey", col="grey")
        lines(xj+c(-w,w)/2, rep(maxBias0[i,c3], 2), lwd=2, col=1, lend=2)
        ## 5 min
        c5 <- c("max5_0", "max5_b")[j]
        xj <- ct[i] + c(-0.5, 0.5)[j] + w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxBias1[i,c5], maxBias1[i,c5], maxBias2[i,c5], maxBias2[i,c5]),
            border=1, col=1)
        lines(xj+c(-w,w)/2, rep(maxBias0[i,c5], 2), lwd=2, col="grey", lend=2)
        ## text
        if (i == 1) {
            ex <- list(
                expression(M[0]),
                expression(M[f]))
            text(ct[i] + c(-0.5, 0.5)[j],
                maxBias1[i,c("max3_0", "max3_b")[j]] + 0.02,
                ex[[j]], cex=1)
        }
    }
}
abline(h=0, lty=2)
legend("topright", fill=c("grey", "black"), legend=c("3-min", "5-min"), bty="n")

par(mar=c(5,4,0,3))
plot(ct, rep(0, length(ng)), xlab="Sample size", ylab="Variance",
    type="n", axes=FALSE,
    ylim=c(0, 1.1*max(maxVar1,maxVar2)), xlim=c(ct[1]-2,ct[length(ct)]+2))
box()
axis(2)
axis(1, ct, ng, tick=FALSE)
for (i in 1:length(ng)) {
    for (j in 1:2) {
        ## 3 min
        c3 <- c("max3_0", "max3_b")[j]
        xj <- ct[i] + c(-0.5, 0.5)[j] - w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxVar1[i,c3], maxVar1[i,c3], maxVar2[i,c3], maxVar2[i,c3]),
            border="grey", col="grey")
        ## 5 min
        c5 <- c("max5_0", "max5_b")[j]
        xj <- ct[i] + c(-0.5, 0.5)[j] + w/2
        polygon(xj+c(-w, w, w, -w)/2,
            c(maxVar1[i,c5], maxVar1[i,c5], maxVar2[i,c5], maxVar2[i,c5]),
            border=1, col=1)
        ## text
        if (i == 1) {
            ex <- list(
                expression(M[0]),
                expression(M[f]))
            text(ct[i] + c(-0.5, 0.5)[j],
                maxVar1[i,c("max3_0", "max3_b")[j]]+0.003,
                ex[[j]], cex=1)
        }
    }
}
par(op)
dev.off()


mbv <- data.frame(MedB=colMeans(maxBias0),MedV=colMeans(maxVar0))
round(mbv,6)

## =============================================================================
## when does model fail?

load(file.path(ROOT2, "problem.Rdata"))
ss <- c("0 observation with multiple duration (1)", "1 observation with multiple duration (2)")

dat0 <- problem0[!(problem0$msg %in% ss),]
dat0$okfit <- ifelse(is.na(dat0$logphi) |
    exp(dat0$logphi) < 0.001 | exp(dat0$logphi) > 5, 0, 1)
dat0$okse <- ifelse(is.na(dat0$se_logphi) | dat0$se_logphi > 10^4, 0, 1)
dat0$okse[dat0$okfit == 0] <- 0
dat0$msg <- NULL
dat0 <- dat0[dat0$ndet > 1,]

datb <- problemb[!(problemb$msg %in% ss),]
datb$okfit <- ifelse(is.na(datb$logphi) |
    exp(datb$logphi) < 0.001 | exp(datb$logphi) > 5, 0, 1)
datb$okse <- ifelse(is.na(datb$se_logphi) |
    datb$se_logphi > 10^4 | datb$se_logitc > 10^4, 0, 1)
datb$okse[datb$okfit == 0] <- 0
datb$msg <- NULL
datb <- datb[datb$ndet > 1,]

dat0$phi <- sptab$M0_phi[match(dat0$spp, sptab$spp)]
dat0$occ <- sptab$Occ[match(dat0$spp, sptab$spp)]
dat0$ym <- sptab$Ymean[match(dat0$spp, sptab$spp)]
datb$phi <- sptab$M0_phi[match(datb$spp, sptab$spp)]
datb$occ <- sptab$Occ[match(datb$spp, sptab$spp)]
datb$ym <- sptab$Ymean[match(datb$spp, sptab$spp)]

dat0 <- droplevels(dat0[dat0$spp %in% SPPfull,])
datb <- droplevels(datb[datb$spp %in% SPPfull,])

#m1 <- glm(okfit ~ ndet, dat0, family=binomial)
#m2 <- glm(okfit ~ log(ndet), dat0, family=binomial)
#m3 <- glm(okfit ~ sqrt(ndet), dat0, family=binomial)
#AIC(m1,m2,m3) # log(n) is best

summary(m0x <- glm(okfit ~ log(ndet), dat0, family=binomial))
summary(m0xs <- glm(okse ~ log(ndet), dat0, family=binomial))
summary(mbx <- glm(okfit ~ log(ndet), datb, family=binomial))
summary(mbxs <- glm(okse ~ log(ndet), datb, family=binomial))

vnd <- seq(2, 1400, len=100)
df <- data.frame(ndet=vnd)
X <- model.matrix(delete.response(terms(m0x)), df)

fit0 <- plogis(drop(X %*% coef(m0x)))
fit0s <- plogis(drop(X %*% coef(m0xs)))
fitb <- plogis(drop(X %*% coef(mbx)))
fitbs <- plogis(drop(X %*% coef(mbxs)))

#png(file.path(ROOT2, "tabfig", "Fig3_failure.png"), width=450, height=450)
pdf(file.path(ROOT2, "tabfig", "Fig3_failure.pdf"), width=6, height=6)
par(las=1)
plot(vnd, fit0, col=1, ylim=c(0,1), type="n", lwd=2,
    xlab="Sample size", ylab="Probability of success")
abline(h=0.9, lty=1, col="grey")
lines(vnd, fit0, col=1, lty=1, lwd=2)
lines(vnd, fitb, col=1, lty=2, lwd=2)
lines(vnd, fitbs, col=1, lty=3, lwd=2)
#abline(v=vnd[which.min(abs(fit0-0.9))], col=2, lty=2)
#abline(v=vnd[which.min(abs(fitb-0.9))], col=4, lty=2)
#abline(v=vnd[which.min(abs(fit0s-0.9))], col=2, lty=2)
#abline(v=vnd[which.min(abs(fitbs-0.9))], col=4, lty=2)
legend("bottomright", col=1, lty=c(1,2,3), lwd=2, bty="n",
    legend=c(
        paste0("M0 fit & SE (90% = ", ceiling(vnd[which.min(abs(fit0-0.9))]), ")"),
        paste0("Mb fit (90% = ", ceiling(vnd[which.min(abs(fitb-0.9))]), ")"),
        paste0("Mb SE (90% = ", ceiling(vnd[which.min(abs(fitbs-0.9))]), ")")))
dev.off()

ceiling(vnd[which.min(abs(fit0-0.9))])
ceiling(vnd[which.min(abs(fitb-0.9))])
ceiling(vnd[which.min(abs(fit0s-0.9))])
ceiling(vnd[which.min(abs(fitbs-0.9))])

## SE for prediction
se_fun <- function(mod) {
    rn <- mvrnorm(4999, coef(mod), vcov(mod))
    fit <- cbind(plogis(drop(X %*% coef(mod))),
        apply(rn, 1, function(z) plogis(drop(X %*% z))))
    apply(fit, 2, function(z) ceiling(vnd[which.min(abs(z-0.9))]))
}
rbind(m0=c(est=se_fun(m0x)[1], SE=round(sd(se_fun(m0x)), 1)),
    m0s=c(se_fun(m0xs)[1], round(sd(se_fun(m0xs)), 1)),
    mb=c(se_fun(mbx)[1], round(sd(se_fun(mbx)), 1)),
    mbs=c(se_fun(mbxs)[1], round(sd(se_fun(mbxs)), 1)))


with(dat0[dat0$okse==1,], plot(logphi, se_logphi))
with(dat0[dat0$okse==1,], plot(ym, se_logphi))

summary(m0x <- glm(okfit ~ log(ndet)*log(ymedian), dat0, family=binomial))
summary(m0xs <- glm(okse ~ log(ndet)*log(ymedian), dat0, family=binomial))
summary(mbx <- glm(okfit ~ log(ndet)*log(ymedian), datb, family=binomial))
summary(mbxs <- glm(okse ~ log(ndet)*log(ymedian), datb, family=binomial))

vnd <- seq(2, 500, len=100)
vym <- seq(1, 4, len=100)
df <- expand.grid(ndet=vnd, ymedian=vym)
X <- model.matrix(delete.response(terms(m0x)), df)

fit0 <- matrix(plogis(drop(X %*% coef(m0x))), length(vnd), length(vym))
fit0s <- matrix(plogis(drop(X %*% coef(m0xs))), length(vnd), length(vym))
fitb <- matrix(plogis(drop(X %*% coef(mbx))), length(vnd), length(vym))
fitbs <- matrix(plogis(drop(X %*% coef(mbxs))), length(vnd), length(vym))


filled.contour(vnd, vym, fit0)

par(mfrow=c(2,2))
image(vnd, vym, -fit0, main="m0 fit")
contour(vnd, vym, fit0, add=TRUE)
image(vnd, vym, -fit0s, main="m0 se")
contour(vnd, vym, fit0s, add=TRUE)
image(vnd, vym, -fitb, main="mb fit")
contour(vnd, vym, fitb, add=TRUE)
image(vnd, vym, -fitbs, main="mb se")
contour(vnd, vym, fitbs, add=TRUE)


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

## average responses across species

library(MASS)
library(KernSmooth)

pf <- function(var, mod, n=10^4, resol=0.1) {
    if (mod == "0")
        sppPred <- sppPred0
    if (mod == "b")
        sppPred <- sppPredb
    var <- switch(var,
        "TSSR"=pkDur$TSSR[iii]*24,
        "JDAY"=pkDur$JDAY[iii]*365,
        "TSLS"=pkDur$TSLS[iii]*365)
    ix <- rep(var, ncol(sppPred))
    iy <- as.numeric(sppPred)
    if (is.null(n))
        n <- length(ix)
    is <- sample.int(length(ix), n)

    d <- KernSmooth::bkde2D(cbind(ix[is], iy[is]), bandwidth=c(diff(range(ix))/50, resol))
    #d <- kde2d(ix[is], iy[is], n=c(50, round(1/resol)))
    names(d) <- c("x", "y", "z")
    d1 <- d2 <- d
    d1$z <- d$z / rowSums(d$z)
    for (i in 1:nrow(d$z))
        d2$z[i,] <- cumsum(d1$z[i,])
    list(kde=d, std=d1, cumul=d2)
}
pf_simple <- function(var, mod, nx=21, ny=11) {
    if (mod == "0")
        sppPred <- sppPred0
    if (mod == "b")
        sppPred <- sppPredb
    xval <- switch(var,
        "TSSR"=pkDur$TSSR[iii]*24,
        "JDAY"=pkDur$JDAY[iii]*365,
        "TSLS"=pkDur$TSLS[iii]*365)
    xbr <- seq(min(xval), max(xval), length.out=nx)
    ybr <- seq(0, 1, length.out=ny)
    ymid <- ybr[-1]-diff(ybr)/2
    xmid <- xbr[-1]-diff(xbr)/2
    xc <- cut(xval, breaks=xbr, include.lowest = TRUE)
    Means <- groupMeans(sppPred, 1, xc)
    Qs <- t(apply(Means, 1, function(z) table(cut(z, breaks=ybr, include.lowest = TRUE))))
    out <- list(x=xmid, y=ymid, z=Qs)
}

iii <- rep(TRUE, nrow(pkDur))
iii[pkDur$TSSR < quantile(pkDur$TSSR, 0.025)] <- FALSE
iii[pkDur$TSSR > quantile(pkDur$TSSR, 0.975)] <- FALSE
iii[pkDur$JDAY < quantile(pkDur$JDAY, 0.025)] <- FALSE
iii[pkDur$JDAY > quantile(pkDur$JDAY, 0.975)] <- FALSE
iii[pkDur$TSLS < quantile(pkDur$TSLS, 0.025)] <- FALSE
iii[pkDur$TSLS > quantile(pkDur$TSLS, 0.975)] <- FALSE

TT <- 3
RESOL <- 0.1
mSPPfull <- SPPfull
mSPP <- SPP
OUT <- list()

#migtype <- lht$DATABASE_MIG_TYPE
#names(migtype) <- rownames(lht)
#mSPPfull <- SPPfull[migtype[SPPfull]=="W"]
#mSPP <- SPP[migtype[SPP]=="W"]
#mSPPfull <- SPPfull[migtype[SPPfull]!="W"]
#mSPP <- SPP[migtype[SPP]!="W"]


for (WHAT in c("TSSR","JDAY","TSLS")) {
    cat(WHAT, "\n");flush.console()

    Xpk2 <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), pkDur[iii,])
    sppPred0 <- matrix(NA, nrow(Xpk2), length(SPPfull))
    colnames(sppPred0) <- SPPfull
    sppPredb <- matrix(NA, nrow(Xpk2), length(SPP))
    colnames(sppPredb) <- SPP

    mc <- which(!grepl(WHAT, colnames(Xpk2)))[-1]
    for (cc in mc)
        Xpk2[,cc] <- mean(Xpk2[,cc])
    COOL <- names(ff)[sapply(NAMES, function(z) any(grepl(WHAT, z)))]

    for (spp in mSPPfull) {
        best0 <- as.character((0:14)[which.max(waic0[spp,])])
        if (best0 %in% COOL) {
            cf02 <- .BAMCOEFSrem$sra_estimates[[spp]][[best0]]$coefficients
            sppPred0[,spp] <- 1-exp(-TT*exp(drop(Xpk2[,gsub("log.phi_", "",
                names(cf02)),drop=FALSE] %*% cf02)))
        }
    }
    sppPred0 <- sppPred0[,colSums(is.na(sppPred0))==0]

    for (spp in mSPP) {
        bestb <- as.character((0:14)[which.max(waicb[spp,])])
        if (bestb %in% COOL) {
            cfb2 <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$coefficients
            sppPredb[,spp] <- 1-plogis(drop(Xpk2[,gsub("logit.c_", "",
                names(cfb2)[-1]),drop=FALSE] %*%
                cfb2[-1])) * exp(-TT*exp(cfb2[1]))
        }
    }
    sppPredb <- sppPredb[,colSums(is.na(sppPredb))==0]

    gc()
    b0 <- pf(WHAT, "0", n=100000, resol=RESOL)
#    b0 <- pf_simple(WHAT, "0", nx=11+2, ny=11)
    gc()
    bb <- pf(WHAT, "b", n=100000, resol=RESOL)
#    bb <- pf_simple(WHAT, "b", nx=11+2, ny=11)

OUT[[WHAT]] <- list("0"=b0, "b"=bb)
}

#sort(apply(sppPred0, 2, max))
#sort(apply(sppPredb, 2, max))

col <- colorRampPalette(c("white", "black"))(100)[c(1,1,1,1,1:66)]
nl <- 5

plf2 <- function(b, ...) {
    b1 <- b$std
    b2 <- b$cumul
    image(b1, col=col, ...)

    o <- order(b1$z)
    i <- order(o)
    f <- b1$z
    f <- f/sum(f)
    cs <- cumsum(f[o])[i]
    dim(cs) <- dim(b1$z)
    b3 <- b1
    b3$z <- cs

    contour(b3, add=TRUE, lwd=1, levels=0.75, drawlabels=FALSE)
    contour(b2, add=TRUE, lty=2, levels=c(0.05, 0.25, 0.5, 0.75, 0.95))

    box()
}

#png(file.path(ROOT2, "tabfig", paste0("Fig4_responses", TT, "min.png")),
#    height=1200, width=1600, res=150)
pdf(file.path(ROOT2, "tabfig", paste0("Fig4_responses", TT, "min-bkde2d.pdf")),
    height=7, width=8)
op <- par(mfrow=c(2,3), las=1, mar=c(5, 4, 4, 2)+0.1, oma=c(2,6,1,2), xpd=NA)

par(mar=c(5, 0, 4, 0))
plf2(OUT[["TSSR"]][["0"]], axes=FALSE, ann=FALSE)
axis(1)
axis(2)
title(ylab="P(availability)", xlab="Time since sunrise (h)")
text(-6.5, 1.2, expression(M[0]^varphi), cex=1.25)
par(mar=c(5, 0, 4, 0))
plf2(OUT[["JDAY"]][["0"]], axes=FALSE, ann=FALSE)
axis(1)
title(xlab="Ordinal day")
par(mar=c(5, 0, 4, 0))
plf2(OUT[["TSLS"]][["0"]], axes=FALSE, ann=FALSE)
axis(1)
title(xlab="Days since spring")

par(mar=c(5, 0, 4, 0))
plf2(OUT[["TSSR"]][["b"]], axes=FALSE, ann=FALSE)
axis(1)
axis(2)
title(ylab="P(availability)", xlab="Time since sunrise (h)")
text(-6.5, 1.2, expression(M[f]^c), cex=1.25)
par(mar=c(5, 0, 4, 0))
plf2(OUT[["JDAY"]][["b"]], axes=FALSE, ann=FALSE)
axis(1)
title(xlab="Ordinal day")
par(mar=c(5, 0, 4, 0))
plf2(OUT[["TSLS"]][["b"]], axes=FALSE, ann=FALSE)
axis(1)
title(xlab="Days since spring")

par(op)
dev.off()



## compare PIF time adjustment

sptab <- read.csv(file.path(ROOT2, "tabfig", "spptab.csv"))
rownames(sptab) <- sptab$spp

library(lhreg)
lhrd <- lhreg_data[match(rownames(sptab), lhreg_data$spp),]
table(Migr=lhrd$Mig2, Has_TSSR=sptab$Both_best %in% hasSR, useNA="a")
table(Migr=lhrd$Mig2, Has_JDAY=sptab$Both_best %in% hasJD, useNA="a")
table(Migr=lhrd$Mig2, Has_DSLS=sptab$Both_best %in% hasLS, useNA="a")
table(Migr=lhrd$Mig2, Has_LSorJD=sptab$Both_best %in% hasLS | sptab$Both_best %in% hasJD, useNA="a")
#lht2 <- lht[match(rownames(sptab), lht$SPECIES),]

pif <- read.csv(file.path(ROOT2, "popContinental_v2_22-May-2013.csv"))
pif <- pif[,c("Common.Name","Scientific.Name","Time.Adjust")]

compare_sets(sptab$common_name, pif$Common.Name)
compare_sets(sptab$scientific_name, pif$Scientific.Name)

sptab[sptab$common_name %in% setdiff(sptab$common_name, pif$Common.Name),]

sptab$tadj <- pif$Time.Adjust[match(sptab$common_name, pif$Common.Name)]
sptab$Inv_p30 <- 1/(1-exp(-3*sptab$M0_phi))
sptab$Inv_p3b <- 1/(1-sptab$Mb_c*exp(-3*sptab$Mb_phi))
#sptab$bias0 <- (1/sptab$p30) / sptab$tadj
#sptab$biasb <- (1/sptab$p3b) / sptab$tadj
write.csv(sptab, row.names=FALSE, file=file.path(ROOT2, "tabfig", "spptab.csv"))

boxplot(sptab[,c("tadj","Inv_p30","Inv_p3b")])


## response curves

nnn <- 200
TT <- 3
spt <- .BAMCOEFSrem$spp_table[SPPfull,]

ii <- intersect(rownames(xtDur[[1]]), rownames(pkDur))
pkd <- pkDur[ii,]
Range <- list()
for (spp in SPPfull) {
    xt <- xtDur[[spp]][ii,]
    Range[[spp]] <- apply(pkd[rowSums(xt, na.rm=TRUE) > 0,c("TSSR","JDAY","TSLS")], 2, range, na.rm=TRUE)
}

vsr <- seq(-0.2, 0.5, len=nnn)
vjd <- seq(0.35, 0.5, len=nnn)
vls <- seq(0, 0.3, len=nnn)
mjd <- 0.45
msr <- 0.13
mls <- 0.12
dsr <- data.frame(TSSR=vsr, JDAY=mjd, TSLS=mls)
djd <- data.frame(TSSR=msr, JDAY=vjd, TSLS=mls)
dls <- data.frame(TSSR=msr, JDAY=mjd, TSLS=vls)
Xsr <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), dsr)
Xjd <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), djd)
Xls <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), dls)

OKsr <- names(ff)[sapply(NAMES, function(z) any(grepl("TSSR", z)))]
OKjd <- names(ff)[sapply(NAMES, function(z) any(grepl("JDAY", z)))]
OKls <- names(ff)[sapply(NAMES, function(z) any(grepl("TSLS", z)))]

r0sr <- matrix(NA, nnn, length(SPPfull))
colnames(r0sr) <- SPPfull
rbjd <- rbls <- rbsr <- r0jd <- r0ls <- r0sr

for (spp in SPPfull) {
cat(spp, "\n");flush.console()

rr <- Range[[spp]]
Rsr <- dsr[,"TSSR"] >= rr[1,"TSSR"] & dsr[,"TSSR"] <= rr[2,"TSSR"]
Rjd <- djd[,"JDAY"] >= rr[1,"JDAY"] & djd[,"JDAY"] <= rr[2,"JDAY"]
Rls <- dls[,"TSLS"] >= rr[1,"TSLS"] & dls[,"TSLS"] <= rr[2,"TSLS"]

best0 <- as.character((0:14)[which.max(waic0[spp,])])
cf02 <- .BAMCOEFSrem$sra_estimates[[spp]][[best0]]$coefficients

#if (best0 == "0")
#r0sr[,spp] <- 1-exp(-TT*exp(cf02))
if (best0 %in% OKsr)
r0sr[Rsr,spp] <- 1-exp(-TT*exp(drop(Xsr[Rsr,gsub("log.phi_", "", names(cf02)),drop=FALSE] %*% cf02)))
if (best0 %in% OKjd)
r0jd[Rjd,spp] <- 1-exp(-TT*exp(drop(Xjd[Rjd,gsub("log.phi_", "", names(cf02)),drop=FALSE] %*% cf02)))
if (best0 %in% OKls)
r0ls[Rls,spp] <- 1-exp(-TT*exp(drop(Xls[Rls,gsub("log.phi_", "", names(cf02)),drop=FALSE] %*% cf02)))

if (spp %in% SPP) {
bestb <- as.character((0:14)[which.max(waicb[spp,])])
cfb2 <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$coefficients

#if (bestb == "0")
#r0sr[,spp] <- 1-plogis(cfb2[-1])*exp(-TT*exp(cfb2[1]))
if (bestb %in% OKsr)
rbsr[Rsr,spp] <- 1-plogis(drop(Xsr[Rsr,gsub("logit.c_", "", names(cfb2)[-1]),drop=FALSE] %*%
    cfb2[-1])) * exp(-TT*exp(cfb2[1]))
if (bestb %in% OKjd)
rbjd[Rjd,spp] <- 1-plogis(drop(Xjd[Rjd,gsub("logit.c_", "", names(cfb2)[-1]),drop=FALSE] %*%
    cfb2[-1])) * exp(-TT*exp(cfb2[1]))
if (bestb %in% OKls)
rbls[Rls,spp] <- 1-plogis(drop(Xls[Rls,gsub("logit.c_", "", names(cfb2)[-1]),drop=FALSE] %*%
    cfb2[-1])) * exp(-TT*exp(cfb2[1]))
}
}


col <- rgb(0.75,0.75,0.75,0.1) # "grey"
pdf(file.path(ROOT2, "tabfig", "responses.pdf"), width=7, height=5, onefile=TRUE)
for (spp in SPPfull) {
op <- par(mfrow=c(2,3))
matplot(vsr*24, r0sr, type=if (all(is.na(r0sr[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Time since sunrise (h)", main=spt[spp,"common_name"])
lines(vsr*24, r0sr[,spp], lwd=3, col=1)
matplot(vjd*365, r0jd, type=if (all(is.na(r0jd[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Julian day")
lines(vjd*365, r0jd[,spp], lwd=3, col=1)
matplot(vls*365, r0ls, type=if (all(is.na(r0ls[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Days since spring")
lines(vls*365, r0ls[,spp], lwd=3, col=1)
matplot(vsr*24, rbsr, type=if (all(is.na(rbsr[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Time since sunrise (h)")
lines(vsr*24, rbsr[,spp], lwd=3, col=1)
matplot(vjd*365, rbjd, type=if (all(is.na(rbjd[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Julian day")
lines(vjd*365, rbjd[,spp], lwd=3, col=1)
matplot(vls*365, rbls, type=if (all(is.na(rbls[,spp]))) "n" else "l",
    lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Days since spring")
lines(vls*365, rbls[,spp], lwd=3, col=1)
par(op)
}
dev.off()




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

pdf(file.path(ROOT2, "tabfig", "spp-pred.pdf"), onefile=TRUE, width=10, height=5)

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

## time of day vs Tadj
#vsr <- seq(-0.13, 0.36, len=200)
#df <- data.frame(TSSR=vsr)
#X <- model.matrix(ff[["4"]], df)
Xpk <- model.matrix(ff[["4"]], pkDur)
Xpk2 <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), pkDur)
sppTadj <- list()

for (spp in rownames(sptab)[sptab$model == "both" & !is.na(sptab$tadj)]) {
cf0 <- .BAMCOEFSrem$sra_estimates[[spp]][["4"]]$coefficients
cfb <- .BAMCOEFSmix$sra_estimates[[spp]][["4"]]$coefficients
best0 <- as.character((0:14)[which.max(waic0[spp,])])
bestb <- as.character((0:14)[which.max(waicb[spp,])])
cf02 <- .BAMCOEFSrem$sra_estimates[[spp]][[best0]]$coefficients
cfb2 <- .BAMCOEFSmix$sra_estimates[[spp]][[bestb]]$coefficients

#p30 <- 1-exp(-3*exp(drop(X %*% cf0)))
#p3b <- 1-plogis(drop(X %*% cfb[-1]))*exp(-3*exp(cfb[1]))
p30pk <- 1-exp(-3*exp(drop(Xpk %*% cf0)))
p3bpk <- 1-plogis(drop(Xpk %*% cfb[-1]))*exp(-3*exp(cfb[1]))
p30pk2 <- 1-exp(-3*exp(drop(Xpk2[,gsub("log.phi_", "", names(cf02)),drop=FALSE] %*% cf02)))
p3bpk2 <- 1-plogis(drop(Xpk2[,gsub("logit.c_", "", names(cfb2)[-1]),drop=FALSE] %*% cfb2[-1])) *
    exp(-3*exp(cfb2[1]))

tadj0 <- max(p30pk) / mean(p30pk)
tadjb <- max(p3bpk) / mean(p3bpk)
tadj02 <- max(p30pk2) / mean(p30pk2)
tadjb2 <- max(p3bpk2) / mean(p3bpk2)
sppTadj[[spp]] <- c(PIF=sptab[spp, "tadj"],
    M0t_sr=tadj0, Mbt_sr=tadjb,
    M0t_jdsr=tadj02, Mbt_jdsr=tadjb2,
    M0t_sr_mean=1/mean(p30pk), Mbt_sr_mean=1/mean(p3bpk),
    M0t_jdsr_mean=1/mean(p30pk2), Mbt_jdsr_mean=1/mean(p3bpk2))
}

sppTadj0 <- data.frame(do.call(rbind, sppTadj))

with(sppTadj0, plot(PIF, M0t_sr_mean, ylim=c(0,5), xlim=c(0,5)))

cn <- c("Mbt_jdsr_mean", "M0t_jdsr_mean", "Mbt_sr", "M0t_sr", "PIF")

sppTadj <- sppTadj0[!(rownames(sppTadj0) %in% c("EVGR", "MAWR", "WIPT")),rev(cn)] # MAWR WIPT EVGR

cor.test(sppTadj[,"M0t_sr"], sppTadj[,"PIF"])
cor.test(sppTadj[,"Mbt_sr"], sppTadj[,"PIF"])

cor.test(sppTadj[,"M0t_jdsr_mean"], sppTadj[,"PIF"], method="spearman")
cor.test(sppTadj[,"Mbt_jdsr_mean"], sppTadj[,"PIF"], method="spearman")
cor.test(sppTadj[,"Mbt_jdsr_mean"], sppTadj[,"M0t_jdsr_mean"], method="spearman")

summary(sppTadj[,"M0t_jdsr_mean"]/ sppTadj[,"PIF"])
summary(sppTadj[,"Mbt_jdsr_mean"]/ sppTadj[,"PIF"])
summary(sppTadj[,"Mbt_jdsr_mean"]/ sppTadj[,"M0t_jdsr_mean"])

tmp <- sppTadj[,cn]
tmp[tmp > 5] <- 5
pairs(tmp)
round(cor(tmp, method="spearman"), 3)

#png(file.path(ROOT2, "tabfig", "FigX_tadj.png"), height=450, width=450)
pdf(file.path(ROOT2, "tabfig", "Fig6_tadj.pdf"), height=6, width=6)
par(las=1, mar=c(5, 6, 4, 2) + 0.1)
boxplot(sppTadj[,cn],
    xlab="Time adjustment",
    ylim=c(1,5), horizontal=TRUE, col="grey",
    names=rev(c("Tadj (PIF)", "Tadj (M0t)", "Tadj (Mbt)", "1/p (M0t)", "1/p (Mbt)")))
abline(v=median(sppTadj[,1], na.rm=TRUE))
dev.off()

op <- par(mfrow=c(3,3))
plot(density(sppTadj[sppTadj[,"PIF"]<5,"PIF"]), xlim=c(0,5))
plot.new()
plot.new()

plot(sppTadj[,"PIF"], sppTadj[,"M0t_jdsr_mean"], xlim=c(0,5), ylim=c(0,5))
plot(density(sppTadj[sppTadj[,"M0t_jdsr_mean"]<5,"M0t_jdsr_mean"]), xlim=c(0,5))
plot.new()

plot(sppTadj[,"PIF"], sppTadj[,"Mbt_jdsr_mean"], xlim=c(0,5), ylim=c(0,5))
plot(sppTadj[,"M0t_jdsr_mean"], sppTadj[,"Mbt_jdsr_mean"], xlim=c(0,5), ylim=c(0,5))
plot(density(sppTadj[sppTadj[,"Mbt_jdsr_mean"]<5,"Mbt_jdsr_mean"]), xlim=c(0,5))
par(op)

##





## sra m0 vs mb models


R <- 1000
t <- 0:100/10

pdf(file.path(ROOT2, "tabfig", "pred.pdf"), onefile=TRUE, width=10, height=12)
for (spp in SPP) {

cat(spp, "\n");flush.console()

## model weights
wp <- waic0[spp,]
wq <- waicb[spp,]
names(wp) <- colnames(waic0)
names(wq) <- colnames(waicb)
wph <- waic[spp,1:15]
wqh <- waic[spp,16:30]
names(wp) <- 0:14
names(wq) <- 0:14


## covariate effects
mi0 <- as.character(0:14)[which.min(aicc0[spp,])]
cfi0 <- .BAMCOEFSrem$sra_estimates[[spp]][[mi0]]$coefficients
vci0 <- .BAMCOEFSrem$sra_estimates[[spp]][[mi0]]$vcov

mib <- as.character(0:14)[which.min(aiccb[spp,])]
cfib <- .BAMCOEFSmix$sra_estimates[[spp]][[mib]]$coefficients
vcib <- .BAMCOEFSmix$sra_estimates[[spp]][[mib]]$vcov

#     TSSR             JDAY            TSLS
# Min.   :-0.315   Min.   :0.351   Min.   :-0.101
# 1st Qu.: 0.063   1st Qu.:0.433   1st Qu.: 0.103
# Median : 0.149   Median :0.455   Median : 0.131
# Mean   : 0.141   Mean   :0.455   Mean   : 0.133
# 3rd Qu.: 0.234   3rd Qu.:0.479   3rd Qu.: 0.164
# Max.   : 0.520   Max.   :0.641   Max.   : 0.442
# NA's   :8455     NA's   :5804    NA's   :17255
jd <- seq(0.35, 0.55, 0.01)
ts <- seq(-0.25, 0.5, 0.01)
ls <- seq(0, 0.3, len=length(jd))

xp1 <- expand.grid(JDAY=jd, # ---------- CHECK !!!
    TSSR=ts)
xp1$JDAY2 <- xp1$JDAY^2
xp1$TSSR2 <- xp1$TSSR^2
xp1$Jday <- xp1$JDAY * 365
xp1$Tssr <- xp1$TSSR * 24

xp2 <- expand.grid(TSLS=ls, # ---------- CHECK !!!
    TSSR=ts)
xp2$TSLS2 <- xp2$TSLS^2
xp2$TSSR2 <- xp2$TSSR^2
xp2$Tsls <- xp2$TSLS * 365
xp2$Tssr <- xp2$TSSR * 24

Xp1 <- model.matrix(~., xp1)
colnames(Xp1)[1] <- "INTERCEPT"
Xp2 <- model.matrix(~., xp2)
colnames(Xp2)[1] <- "INTERCEPT"

names(cfi0) <- sapply(strsplit(names(cfi0), "_"), "[[", 2)
names(cfi0)[1] <- "INTERCEPT"
names(cfi0)[names(cfi0) == "I(TSSR^2)"] <- "TSSR2"
names(cfi0)[names(cfi0) == "I(TSLS^2)"] <- "TSLS2"
names(cfi0)[names(cfi0) == "I(JDAY^2)"] <- "JDAY2"
Xp <- if (mi0 %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
Xp <- Xp[,names(cfi0),drop=FALSE]
lphi1 <- drop(Xp %*% cfi0)
pmat <- matrix(exp(lphi1), length(jd), length(ts))
pmat <- 1-exp(-3 * pmat)
pmax <- 1

Xpb <- if (mib %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
cfib1 <- exp(cfib[1])
cfib2 <- cfib[-1]
names(cfib2) <- sapply(strsplit(names(cfib2), "_"), "[[", 2)
names(cfib2)[1] <- "INTERCEPT"
names(cfib2)[names(cfib2) == "I(TSSR^2)"] <- "TSSR2"
names(cfib2)[names(cfib2) == "I(TSLS^2)"] <- "TSLS2"
names(cfib2)[names(cfib2) == "I(JDAY^2)"] <- "JDAY2"

Xpb <- Xpb[,names(cfib2),drop=FALSE]
lphi1b <- 1-plogis(drop(Xpb %*% cfib2))*exp(-3*cfib1)
pmatb <- matrix(lphi1b, length(jd), length(ts))

## CI for m0, m0
cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
phi00 <- exp(rnorm(R, cfi00, sqrt(vci00)))
ci00 <- sapply(phi00, function(z) 1-exp(-t*z))

cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
pcf1b <- mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat)
ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))

np <- .BAMCOEFSrem$sra_n

op <- par(las=1, mfrow=c(3,2))

barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
    main=paste0(as.character(sptab[spp,"common_name"]), " (n=", np[spp], ")",
    " w0=", round(sum(wph),2)),
    ylab="Model weight", xlab="Model ID")
barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", np[spp], ")",
    " wb=", round(sum(wqh),2)),
    ylab="Model weight", xlab="Model ID")

plot(t, 1-exp(-t*exp(cfi00)), type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, pp0, col="grey", lwd=1, lty=1)
matlines(t, ci00, col="grey", lwd=1, lty=1)
lines(t, 1-exp(-t*exp(cfi00)), col=1, lwd=2)

plot(t, 1-plogis(cfi0b[2])*exp(-t*exp(cfi0b[1])), type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, ppb, col="grey", lwd=1, lty=1)
matlines(t, ci0b, col="grey", lwd=1, lty=1)
lines(t, 1-plogis(cfi0b[2])*exp(-t*exp(cfi0b[1])), col=1, lwd=2)

xval <- if (mi0 %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mi0 %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mi0))
box()

xval <- if (mib %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmatb,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mib %in% c("9","10","11","12","13","14"),
        "Days since local springs", "Julian days"),
    ylab="Hours since sunrise",
    main=paste("Best model:", mib))
box()

par(op)
}
dev.off()


## model averaged marginal responses

#spp <- "BOCH"
spt <- .BAMCOEFSrem$spp_table[SPPfull,]
TT <- 3
Nout <- 200
Xpk3 <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), pkDur)
Ctrl <- loess.control(surface="interpolate", statistics="none", trace.hat="approximate")
pdat <- data.frame(
    JDAY = seq(quantile(pkDur$JDAY, 0.025), quantile(pkDur$JDAY, 0.975), length.out=Nout),
    TSSR = seq(quantile(pkDur$TSSR, 0.025), quantile(pkDur$TSSR, 0.975), length.out=Nout),
    TSLS = seq(quantile(pkDur$TSLS, 0.025), quantile(pkDur$TSLS, 0.975), length.out=Nout))

pred_all <- list()
for (spp in SPP) {
    cat(spp, "\n");flush.console()
    pr0 <- prb <- matrix(0, 15, nrow(Xpk3))
    for (mid in 1:15) {
        cf0 <- .BAMCOEFSrem$sra_estimates[[spp]][[mid]]$coefficients
        pr0[mid,] <- 1 - exp(-TT*exp(drop(Xpk3[,gsub("log.phi_", "",
            names(cf0)),drop=FALSE] %*% cf0)))
        cfb <- .BAMCOEFSmix$sra_estimates[[spp]][[mid]]$coefficients
        prb[mid,] <- 1 - plogis(drop(Xpk3[,gsub("logit.c_", "",
            names(cfb)[-1]),drop=FALSE] %*% cfb[-1])) * exp(-TT*exp(cfb[1]))
    }
    w0 <- waic0[spp,]
    wb <- waicb[spp,]
    pr0 <- colSums(w0 * pr0)
    prb <- colSums(wb * prb)

    pred_all[[spp]] <- data.frame(
        ls0_JDAY = predict(loess(pr0 ~ JDAY, pkDur, control=Ctrl), pdat),
        ls0_TSSR = predict(loess(pr0 ~ TSSR, pkDur, control=Ctrl), pdat),
        ls0_TSLS = predict(loess(pr0 ~ TSLS, pkDur, control=Ctrl), pdat),
        lsb_JDAY = predict(loess(prb ~ JDAY, pkDur, control=Ctrl), pdat),
        lsb_TSSR = predict(loess(prb ~ TSSR, pkDur, control=Ctrl), pdat),
        lsb_TSLS = predict(loess(prb ~ TSLS, pkDur, control=Ctrl), pdat))
}
save(pred_all, pdat, file=file.path(ROOT2, "responses-averaged.Rdata"))

col <- rgb(0.5,0.5,1,1) # "grey"
pdf(file.path(ROOT2, "tabfig", "responses-averaged.pdf"), width=7, height=5, onefile=TRUE)
for (spp in SPP) {
op <- par(mfrow=c(2,3))
plot(pdat$TSSR*24, pred_all[[spp]]$ls0_TSSR,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Time since sunrise (h)", main=spt[spp,"common_name"])
plot(pdat$JDAY*365, pred_all[[spp]]$ls0_JDAY,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Julian day")
plot(pdat$TSLS*365, pred_all[[spp]]$ls0_TSLS,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (M0)", xlab="Days since spring", main=spp)
plot(pdat$TSSR*24, pred_all[[spp]]$lsb_TSSR,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Time since sunrise (h)")
plot(pdat$JDAY*365, pred_all[[spp]]$lsb_JDAY,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Julian day")
plot(pdat$TSLS*365, pred_all[[spp]]$lsb_TSLS,
    type="l", lty=1, col=col, ylim=c(0,1), lwd=3,
    ylab="Probability (Mb)", xlab="Days since spring")
par(op)
}
dev.off()

pf <- function(var, mod, resol=0.02) {
    if (mod == "0")
        sppPred <- sapply(pred_all, function(z) {
            switch(var,
                "TSSR"=z$ls0_TSSR,
                "JDAY"=z$ls0_JDAY,
                "TSLS"=z$ls0_TSLS)
        })
    if (mod == "b")
        sppPred <- sapply(pred_all, function(z) {
            switch(var,
                "TSSR"=z$lsb_TSSR,
                "JDAY"=z$lsb_JDAY,
                "TSLS"=z$lsb_TSLS)
        })
    var <- switch(var,
        "TSSR"=pdat$TSSR*24,
        "JDAY"=pdat$JDAY*365,
        "TSLS"=pdat$TSLS*365)
    ix <- rep(var, ncol(sppPred))
    iy <- as.numeric(sppPred)
    d <- kde2d(ix, iy, n=c(100, round(1/resol)))
    d1 <- d2 <- d
    d1$z <- d$z / rowSums(d$z)
    for (i in 1:nrow(d$z))
        d2$z[i,] <- cumsum(d1$z[i,])
    list(kde=d, std=d1, cumul=d2)
}
pf1 <- function(var, mod) {
    if (mod == "0")
        sppPred <- sapply(pred_all, function(z) {
            switch(var,
                "TSSR"=z$ls0_TSSR,
                "JDAY"=z$ls0_JDAY,
                "TSLS"=z$ls0_TSLS)
        })
    if (mod == "b")
        sppPred <- sapply(pred_all, function(z) {
            switch(var,
                "TSSR"=z$lsb_TSSR,
                "JDAY"=z$lsb_JDAY,
                "TSLS"=z$lsb_TSLS)
        })
    var <- switch(var,
        "TSSR"=pdat$TSSR*24,
        "JDAY"=pdat$JDAY*365,
        "TSLS"=pdat$TSLS*365)
    list(x=var, y=sppPred)
}

plf2 <- function(b, ...) {
#    contour(b1, levels=levels, lty=1, ...)
#    contour(b2, add=TRUE, levels=levels, lty=2)
    #col <- colorRampPalette(c("white", "black"))(30)[c(rep(1, 9),1:20)]
    col <- colorRampPalette(c("white", "black"))(30)[1:25]
    #b1 <- b$std
    b1 <- b$kde
    b2 <- b$cumul
    image(b1, col=col, ...)
    contour(b2, add=TRUE, levels=seq(0.1, 0.9, by=0.1))
    #contour(b$std, add=TRUE)
    box()
}


png(file.path(ROOT2, "tabfig", paste0("Fig4_responses-avg-", TT, "min.png")),
    height=1200, width=1600, res=150)
op <- par(mfrow=c(2,3), las=1)
plf2(pf("TSSR", "0"),
    ylab="Probability (M0)", xlab="Time since sunrise (h)")
plf2(pf("JDAY", "0"),
    ylab="Probability (M0)", xlab="Julian day")
plf2(pf("TSLS", "0"),
    ylab="Probability (M0)", xlab="Days since spring")
plf2(pf("TSSR", "b"),
    ylab="Probability (Mb)", xlab="Time since sunrise (h)")
plf2(pf("JDAY", "b"),
    ylab="Probability (Mb)", xlab="Julian day")
plf2(pf("TSLS", "b"),
    ylab="Probability (Mb)", xlab="Days since spring")
par(op)
dev.off()

col <- rgb(0,0,0,0.1) # "grey"
wd <- 10
png(file.path(ROOT2, "tabfig", paste0("Fig4_responses-spaghetti-", TT, "min.png")),
    height=1200, width=1600, res=150)
op <- par(mfrow=c(2,3), las=1)
tmp <- pf1("TSSR", "0")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (M0)", xlab="Time since sunrise (h)")
tmp <- pf1("JDAY", "0")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (M0)", xlab="Julian day")
tmp <- pf1("TSLS", "0")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (M0)", xlab="Days since spring")
tmp <- pf1("TSSR", "b")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (Mb)", xlab="Time since sunrise (h)")
tmp <- pf1("JDAY", "b")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (Mb)", xlab="Julian day")
tmp <- pf1("TSLS", "b")
matplot(tmp$x, tmp$y, lty=1, type="l", col=col, lwd=wd,
    ylab="Probability (Mb)", xlab="Days since spring")
par(op)
dev.off()

## corrected counts and CIs

load(file.path(ROOT2, "var-bias-countsummary.Rdata"))

t <- 0:100/10
R <- 10^4
ii <- c('3'=which(t==3), '5'=which(t==5), '10'=which(t==10))
spp <- "OVEN"

SPPs <- names(res_sum)
SPPss <- SPPs[colSums(sapply(res_sum, function(z) colSums(z)[1:3]) > 0) == 3]
SPPss <- SPPss[SPPss %in% SPP]

pdf(file.path(ROOT2, "tabfig", "corrected-counts-0.pdf"), onefile=TRUE, width=10, height=5)
for (spp in SPPss) {
    cat(spp, "\n");flush.console()
    YYsum <- colSums(res_sum[[spp]])
    YYmean <- YYsum[1:3] / YYsum["n"]
    nn <- sum(rowSums(res_sum[[spp]][,-4]) > 0)

    ## CI for m0, m0
    cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
    phi00 <- exp(c(cfi00, rnorm(R, cfi00, sqrt(vci00))))
    ci00 <- sapply(phi00, function(z) 1-exp(-t*z))
    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))


    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    pcf1b <- rbind(cfi0b, mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))

    p00 <- CI00[ii,]
    p0b <- CI0b[ii,]

    yc00 <- YYmean / p00
    yc0b <- YYmean / p0b

    op <- par(mfrow=c(1,2))
    plot(0, type="n", ylim=c(0,1), xlim=c(0,10), main=spp,
        xlab="Duration (min)", ylab="Probability")
    polygon(c(t, rev(t)), c(CI00[,2], rev(CI00[,3])), border="#0000FF40", col="#0000FF40")
    lines(t, CI00[,1], col="#0000FF")
    polygon(c(t, rev(t)), c(CI0b[,2], rev(CI0b[,3])), border="#FF000040", col="#FF000040")
    lines(t, CI0b[,1], col="#FF0000")
    legend("bottomright", bty="n", fill=c("#0000FF40", "#FF000040"), border=NA,
        lty=1, col=c("#0000FF", "#FF0000"), legend=c("M0", "Mb"))

    plot(0, type="n", ylim=c(0, max(yc00,yc0b)), xlim=c(0,10.2), main=paste("n =", nn),
        xlab="Duration (min)", ylab="Corrected Mean Count")
    polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
        border="#0000FF40", col="#0000FF40")
    polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
        border="#FF000040", col="#FF000040")
    segments(x0=c(3,5,10)-0.1, y0=yc00[,2], y1=yc00[,3], col="#0000FF")
    segments(x0=c(3,5,10)+0.1, y0=yc0b[,2], y1=yc0b[,3], col="#FF0000")
    points(c(3,5,10)-0.1, yc00[,1], col="#0000FF", cex=1.2)
    points(c(3,5,10)+0.1, yc0b[,1], col="#FF0000", cex=1.2)
    par(op)

}
dev.off()

## todo: bootstrap for M0 and Mb to better capture uncertainty

## resDur_mix, resDur_rem
load("e:/peter/bam/May2015/out/estimates_SRA_QPAD_v2016-ABMIfix-boot0rem.Rdata")
load("e:/peter/bam/May2015/out/estimates_SRA_QPAD_v2016-ABMIfix-boot0mix.Rdata")

SPPss <- sort(SPPss[SPPss %in% names(resDur_mix)])
pdf(file.path(ROOT2, "tabfig", "corrected-counts-0boot.pdf"), onefile=TRUE, width=10, height=10)
for (spp in SPPss) {
    cat(spp, "\n");flush.console()
    YYsum <- colSums(res_sum[[spp]])
    YYmean <- YYsum[1:3] / YYsum["n"]
    nn <- sum(rowSums(res_sum[[spp]][,-4]) > 0)
    op <- par(mfrow=c(2,2))

    ## CI for m0, m0: asymptotics
    cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
    phi00 <- exp(c(cfi00, rnorm(R, cfi00, sqrt(vci00))))
    ci00 <- sapply(phi00, function(z) 1-exp(-t*z))
    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))


    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    pcf1b <- rbind(cfi0b, mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))

    p00 <- CI00[ii,]
    p0b <- CI0b[ii,]

    yc00 <- YYmean / p00
    yc0b <- YYmean / p0b

    plot(0, type="n", ylim=c(0,1), xlim=c(0,10), main=paste(spp, "Wald-type"),
        xlab="Duration (min)", ylab="Probability")
    polygon(c(t, rev(t)), c(CI00[,2], rev(CI00[,3])), border="#0000FF40", col="#0000FF40")
    lines(t, CI00[,1], col="#0000FF")
    polygon(c(t, rev(t)), c(CI0b[,2], rev(CI0b[,3])), border="#FF000040", col="#FF000040")
    lines(t, CI0b[,1], col="#FF0000")
    legend("bottomright", bty="n", fill=c("#0000FF40", "#FF000040"), border=NA,
        lty=1, col=c("#0000FF", "#FF0000"), legend=c("M0", "Mb"))

    plot(0, type="n", ylim=c(0, max(yc00,yc0b)), xlim=c(0,10.2), main=paste("n =", nn),
        xlab="Duration (min)", ylab="Corrected Mean Count")
    polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
        border="#0000FF40", col="#0000FF40")
    polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
        border="#FF000040", col="#FF000040")
    segments(x0=c(3,5,10)-0.1, y0=yc00[,2], y1=yc00[,3], col="#0000FF")
    segments(x0=c(3,5,10)+0.1, y0=yc0b[,2], y1=yc0b[,3], col="#FF0000")
    points(c(3,5,10)-0.1, yc00[,1], col="#0000FF", cex=1.2)
    points(c(3,5,10)+0.1, yc0b[,1], col="#FF0000", cex=1.2)

    ## CI for m0, m0: bootstrap
    tmp0 <- resDur_rem[[spp]]
    tmp0 <- tmp0[!sapply(tmp0, inherits, "try-error")]
    phi00 <- exp(sapply(tmp0, "[[", "coefficients"))
    ci00 <- sapply(phi00, function(z) 1-exp(-t*z))
    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))

    tmpb <- resDur_mix[[spp]]
    tmpb <- tmpb[!sapply(tmpb, inherits, "try-error")]
    pcf1b <- t(sapply(tmpb, "[[", "coefficients"))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))

    p00 <- CI00[ii,]
    p0b <- CI0b[ii,]

    yc00 <- YYmean / p00
    yc0b <- YYmean / p0b

    plot(0, type="n", ylim=c(0,1), xlim=c(0,10), main=paste(spp, "bootstrap"),
        xlab="Duration (min)", ylab="Probability")
    polygon(c(t, rev(t)), c(CI00[,2], rev(CI00[,3])), border="#0000FF40", col="#0000FF40")
    lines(t, CI00[,1], col="#0000FF")
    polygon(c(t, rev(t)), c(CI0b[,2], rev(CI0b[,3])), border="#FF000040", col="#FF000040")
    lines(t, CI0b[,1], col="#FF0000")
    legend("bottomright", bty="n", fill=c("#0000FF40", "#FF000040"), border=NA,
        lty=1, col=c("#0000FF", "#FF0000"), legend=c("M0", "Mb"))

    plot(0, type="n", ylim=c(0, max(yc00,yc0b)), xlim=c(0,10.2), main=paste("n =", nn),
        xlab="Duration (min)", ylab="Corrected Mean Count")
    polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
        border="#0000FF40", col="#0000FF40")
    polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
        border="#FF000040", col="#FF000040")
    segments(x0=c(3,5,10)-0.1, y0=yc00[,2], y1=yc00[,3], col="#0000FF")
    segments(x0=c(3,5,10)+0.1, y0=yc0b[,2], y1=yc0b[,3], col="#FF0000")
    points(c(3,5,10)-0.1, yc00[,1], col="#0000FF", cex=1.2)
    points(c(3,5,10)+0.1, yc0b[,1], col="#FF0000", cex=1.2)

    par(op)

}
dev.off()

## composite figure for ms

load(file.path(ROOT2, "var-bias-countsummary.Rdata"))

t <- 0:100/10
R <- 10^4
ii <- c('3'=which(t==3), '5'=which(t==5), '10'=which(t==10))


pdf(file.path(ROOT2, "tabfig", "FigZ_corrected-counts.pdf"), onefile=TRUE, width=7, height=9)
op <- par(mfrow=c(3,2), las=1, mar=c(5,5,1,1))
for (spp in c("CONW", "WEWP", "RUBL")) {

    YYsum <- colSums(res_sum[[spp]])
    YYmean <- YYsum[1:3] / YYsum["n"]
    nn <- sum(rowSums(res_sum[[spp]][,-4]) > 0)

    ## CI for m0, m0: asymptotics
    cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
    phi00 <- exp(c(cfi00, rnorm(R, cfi00, sqrt(vci00))))
    ci00 <- sapply(phi00, function(z) 1-exp(-t*z))
    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))

    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    pcf1b <- rbind(cfi0b, mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))

    p00 <- CI00[ii,]
    p0b <- CI0b[ii,]

    yc00 <- YYmean / p00
    yc0b <- YYmean / p0b

    Ref <- yc00[3,1]
    yc00 <- yc00 / Ref
    yc0b <- yc0b / Ref


    col <- "#80808080"
    plot(0, type="n", ylim=c(0,1), xlim=c(0,10),
        xlab="Duration (min)", ylab="Probability")
    polygon(c(t, rev(t)), c(CI00[,2], rev(CI00[,3])), border=col, col=col)
    lines(t, CI00[,1], col=1, lty=1, lwd=2)
    polygon(c(t, rev(t)), c(CI0b[,2], rev(CI0b[,3])), border=col, col=col)
    lines(t, CI0b[,1], col=1, lty=2, lwd=2)
    abline(v=c(3,5),lty=3)
    legend("bottomright", bty="n", lty=c(1,2), col=1,
        title=paste(spp, "\nn =", nn), lwd=2,
        legend=c(expression(M[0]), expression(M[f])))

    plot(0, type="n", ylim=c(0, 2), xlim=c(0,10.2),
        xlab="Duration (min)", ylab="Corrected Relative Mean Count")
    polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
        border=col, col=col)
    polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
        border=col, col=col)
    segments(x0=c(3,5,10)-0.1, y0=yc00[,2], y1=yc00[,3], col=1, lwd=2)
    segments(x0=c(3,5,10)+0.1, y0=yc0b[,2], y1=yc0b[,3], col=1, lwd=2, lty=1)
    abline(v=c(3,5),lty=3)
    points(c(3,5,10)-0.1, yc00[,1], col=1, cex=1.2, pch=19)
    points(c(3,5,10)+0.1, yc0b[,1], col=1, cex=1.2, pch=21)
    legend("bottomright", bty="n", pch=c(19,21), col=1, title=spp, lty=c(1,1),
        legend=c(expression(M[0]), expression(M[f])))
}
par(op)
dev.off()

## figure based on res4

pdf(file.path(ROOT2, "tabfig", "FigZZ_corrected-counts.pdf"), onefile=TRUE, width=7, height=9)
op <- par(mfrow=c(3,2), las=1, mar=c(5,5,1,1))
for (spp in c("CONW", "WEWP", "RUBL")) {

    tmp <- res4[[spp]]

    YYsum <- colSums(res_sum[[spp]])
    YYmean <- YYsum[1:3] / YYsum["n"]
    nn <- sum(rowSums(res_sum[[spp]][,-4]) > 0)

    ## CI for m0, m0: asymptotics
    #cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    #vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
    cfi00 <- attr(tmp, "est")$m0cf
    vci00 <- attr(tmp, "est")$m0vc
    phi00 <- exp(c(cfi00, rnorm(R, cfi00, sqrt(vci00))))
    ci00 <- sapply(phi00, function(z) 1-exp(-t*z))
    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))

    #cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    #vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    cfi0b <- attr(tmp, "est")$mbcf
    vci0b <- attr(tmp, "est")$mbvc
    pcf1b <- rbind(cfi0b, mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))

    p00 <- CI00[ii,]
    p0b <- CI0b[ii,]

    yc00 <- YYmean / p00
    yc0b <- YYmean / p0b

    Ref <- yc00[3,1]
    yc00 <- yc00 / Ref
    yc0b <- yc0b / Ref


    col <- "#80808080"
    plot(0, type="n", ylim=c(0,1), xlim=c(0,10),
        xlab="Duration (min)", ylab="Probability")
    polygon(c(t, rev(t)), c(CI00[,2], rev(CI00[,3])), border=col, col=col)
    polygon(c(t, rev(t)), c(CI0b[,2], rev(CI0b[,3])), border=col, col=col)
    lines(t, CI00[,1], col=1, lty=1, lwd=1.5)
    lines(t, CI0b[,1], col=1, lty=2, lwd=1.5)
    abline(v=c(3,5),lty=3)
    legend("topleft", bty="n", lty=c(1,2), col=1,
        title=paste0(spp, " (n =", nn, ")"), lwd=1.5,
        legend=c(expression(M[0]), expression(M[f])))

    plot(0, type="n", ylim=c(0, max(2, yc00, yc0b)), xlim=c(0,10.2),
        xlab="Duration (min)", ylab="Corrected Relative Mean Count")
    polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
        border=col, col=col)
    polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
        border=col, col=col)
    segments(x0=c(3,5,10)-0.15, y0=yc00[,2], y1=yc00[,3], col=1, lwd=1.5)
    segments(x0=c(3,5,10)+0.15, y0=yc0b[,2], y1=yc0b[,3], col=1, lwd=1.5, lty=1)
    abline(v=c(3,5),lty=3)
    points(c(3,5,10)-0.15, yc00[,1], col=1, cex=1.2, pch=19)
    points(c(3,5,10)+0.15, yc0b[,1], col=1, cex=1.2, pch=21)
    legend("topleft", bty="n", pch=c(19,21), col=1, title=spp, lty=c(1,1), lwd=1.5,
        legend=c(expression(M[0]), expression(M[f])))
}
par(op)
dev.off()

