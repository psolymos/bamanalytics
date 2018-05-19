## this version incorporates 3 model types
## M0, Mf, Mb
## Define root folder where data are stored
ROOT <- "e:/peter/bam/May2015"
ROOT2 <- "~/GoogleWork/bam/duration_ms/revisionMarch2017"

## Load required packages
library(MASS)
library(mefa4)
library(detect)

## Load functions kept in separate file
#source("~/repos/bamanalytics/R/dataprocessing_functions.R")

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
#png(file.path(ROOT2, "tabfig", "Fig1_map.png"), height=600, width=800)
#plot(getMap(resolution = "low"),
#    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
#points(pkDur[, c("X","Y")], pch=".",
#    col=rgb(70, 130, 180, alpha=255*0.15, maxColorValue=255))
#points(X0[rowSums(!is.na(D)) > 1, c("X","Y")], pch=19,
#    col="red", cex=0.3)
#dev.off()

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

## species where rem model sample size is at least NMIN

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

Projects <- table(droplevels(X0$PCODE))
SppOcc <- sapply(SPPfull, function(z) sum(rowSums(xtDur[[z]][rn,], na.rm=TRUE)>0)/length(rn))
SppYmean <- sapply(SPPfull, function(z) mean(rowSums(xtDur[[z]][rn,], na.rm=TRUE)))


## sra 0 vs f or b models

## aic is really AICc to account for small sample sizes
aic0 <- .BAMCOEFSrem$sra_aic
aicf <- .BAMCOEFSfmix$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
df0 <- matrix(.BAMCOEFSrem$sra_df, nrow(aic0), 15, byrow=TRUE)
dff <- matrix(.BAMCOEFSfmix$sra_df, nrow(aicb), 15, byrow=TRUE)
dfb <- matrix(.BAMCOEFSmix$sra_df, nrow(aicb), 15, byrow=TRUE)
n0 <- .BAMCOEFSrem$sra_n
nf <- .BAMCOEFSfmix$sra_n
nb <- .BAMCOEFSmix$sra_n

aicc0 <- aic0 + (2*df0*(df0+1)) / (n0-df0-1)
aiccf <- aicf + (2*dff*(dff+1)) / (nf-dff-1)
aiccb <- aicb + (2*dfb*(dfb+1)) / (nb-dfb-1)
colnames(aicc0) <- paste0("M0_", colnames(aic0))
colnames(aiccf) <- paste0("Mf_", colnames(aicf))
colnames(aiccb) <- paste0("Mb_", colnames(aicb))

aic0 <- aicc0[SPPfull,]
aicf <- aiccf[SPP,]
aicb <- aiccb[SPP,]
aicx <- cbind(aicc0[SPP,], aiccb[SPP,], aiccf[SPP,]) # all 3 combined

waic0 <- t(apply(aic0, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicf <- t(apply(aicf, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicb <- t(apply(aicb, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))
waicx <- t(apply(aicx, 1, function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2)
    w/sum(w)
}))

H0 <- apply(waic0, 1, function(z) sum(z^2))
Hf <- apply(waicf, 1, function(z) sum(z^2))
Hb <- apply(waicb, 1, function(z) sum(z^2))
Hx <- apply(waicx, 1, function(z) sum(z^2))

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestf <- as.character(0:14)[apply(aicf, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
bestx <- as.factor(colnames(aicx))[apply(aicx, 1, which.min)]
## these 2 are identical, shouldnt have both
bestx[bestx %in% c("Mf_0", "Mb_0")]
names(best0) <- SPPfull
names(bestb) <- names(bestf) <- names(bestx) <- SPP

tmp <- strsplit(as.character(bestx), "_")
dfb <- data.frame(spp=SPP, best=bestx, type=sapply(tmp, "[[", 1),
    model=as.integer(sapply(tmp, "[[", 2)))
btab <- table(dfb$model, dfb$type)
rownames(btab) <- paste(rownames(btab), as.character(ff))
btab <- addmargins(btab)

## appendix table ---------------------------------

tb <- read.csv("~/GoogleWork/bam/duration_ms/revisionMarch2017/tabfig/spptab.csv")
rownames(tb) <- tb$spp
tb$Mf_best <- as.integer(bestf)[match(rownames(tb), names(bestf))]
tb$Best3 <- bestx[match(rownames(tb), names(bestx))]
#write.csv(tb, row.names=FALSE, file="~/GoogleWork/bam/duration_ms/revisionMarch2017/internal/Appendix-table.csv")

## constant model figure (M0/Mf) ---------------------------------

tt <- seq(0, 11, len=1000)
mat0 <- sapply(SPPfull, function(z) 1-exp(-tt*cfall0[z, 1]))
matb <- sapply(SPP, function(z) 1-cfallb[z, "c"]*exp(-tt*cfallb[z, "phi_b"]))

pdf(file.path(ROOT2, "internal", "Fig2-M0-mf-probs.pdf"), width=10, height=5)
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

## p(avail) stats for M0/Mf ---------------------------------

sf <- function(y)
    t(apply(y, 2, function(x) c(Mean=mean(x), SD=sd(x), min=min(x), max=max(x))))
round(sf(t(sapply(SPPfull, function(z)
    1-exp(-c(t3=3,t5=5,t10=10)*cfall0[z, 1])))), 3)
round(sf(t(sapply(SPP, function(z)
    1-cfallb[z, "c"]*exp(-c(t3=3,t5=5,t10=10)*cfallb[z, "phi_b"])))), 3)

## model failure (est and SE) for M0/Mf ---------------------------------

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

#pdf(file.path(ROOT2, "tabfig", "Fig3_failure.pdf"), width=6, height=6)
par(las=1)
plot(vnd, fit0, col=1, ylim=c(0,1), type="n", lwd=2,
    xlab="Sample size", ylab="Probability of success")
abline(h=0.9, lty=1, col="grey")
lines(vnd, fit0, col=1, lty=1, lwd=2)
lines(vnd, fitb, col=1, lty=2, lwd=2)
lines(vnd, fitbs, col=1, lty=3, lwd=2)
legend("bottomright", col=1, lty=c(1,2,3), lwd=2, bty="n",
    legend=c(
        paste0("M0 fit & SE (90% = ", ceiling(vnd[which.min(abs(fit0-0.9))]), ")"),
        paste0("Mb fit (90% = ", ceiling(vnd[which.min(abs(fitb-0.9))]), ")"),
        paste0("Mb SE (90% = ", ceiling(vnd[which.min(abs(fitbs-0.9))]), ")")))
#dev.off()

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

## bias and variance for M0/Mf ---------------------------------

load(file.path(ROOT2, "var-bias-res-3.Rdata"))
resx <- res4[!sapply(res4, inherits, "try-error")] # aggregate + estimate

aaa <- data.frame(Var=as.numeric(sapply(resx, "[[", "Var")),
    MSE=as.numeric(sapply(resx, "[[", "MSE")),
    Bias=as.numeric(sapply(resx, "[[", "Bias")),
    Model=rep(c("0","b","0t","bt"), each=2),
    Duration=c("3","5"),
    Species=as.factor(rep(names(resx), each=8)))
ndet <- sapply(resx, function(z) unname(sum(attr(z, "n"))-attr(z, "n")["0"]))
aaa$n <- ndet[match(aaa$Species, names(ndet))]

ng <- c(20, 50, 100, 200, 400, 1000, 2000)
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

## use res4 for aaa
maxBias0[maxBias0 == -1] <- 0
maxBias1[maxBias1 == -1] <- 0
maxBias2[maxBias2 == -1] <- 0

pdf(file.path(ROOT2, "internal", "Fig3-var-bias-est.pdf"), width=10, height=2*5)
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

rownames(maxBias0) <- rownames(maxVar0) <- paste0("ndet=",ng)
round(100*maxBias2[,1:4],2) # 5%
round(100*maxBias0[,1:4],2) # 50%
round(100*maxBias1[,1:4],2) # 95%
round(100*maxVar2[,1:4],2) # 5%
round(100*maxVar0[,1:4],2) # 50%
round(100*maxVar1[,1:4],2) # 95%
round(100*maxVar0[,1:4],2)



## count correction, example species (m0/Mf) ---------------------------------

## figure is based on res4
tt <- 0:1000/100
R <- 10^4
ii <- c('3'=which(tt==3), '5'=which(tt==5), '10'=which(tt==10))

pdf(file.path(ROOT2, "internal", "Fig4-corrected-counts.pdf"), onefile=TRUE, width=7, height=9)
op <- par(mfrow=c(3,2), las=1, mar=c(5,5,1,1))
for (spp in rev(c("CONW", "WEWP", "RUBL"))) {

    tmp <- res4[[spp]]
    YYmean <- attr(tmp, "y")
#    YYsum <- colSums(res_sum[[spp]])
#    YYmean <- YYsum[1:3] / YYsum["n"]
#    nn <- sum(rowSums(res_sum[[spp]][,-4]) > 0)
    nn <- sum(attr(tmp, "n")) - attr(tmp, "n")["0"]

    ## CI for m0, m0: asymptotics
    #cfi00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
    #vci00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$vcov
    cfi00 <- attr(tmp, "est")$m0cf
    vci00 <- attr(tmp, "est")$m0vc
    phi00 <- exp(c(cfi00, rnorm(R, cfi00, sqrt(vci00))))
    ci00 <- sapply(phi00, function(z) 1-exp(-tt*z))

    #cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    #vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    cfi0b <- attr(tmp, "est")$mbcf
    vci0b <- attr(tmp, "est")$mbvc
    pcf1b <- rbind(cfi0b, mvrnorm(R, cfi0b, Matrix::nearPD(vci0b)$mat))
    ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-tt*exp(z[1])))

    CI00 <- cbind(Est=ci00[,1], t(apply(ci00, 1, quantile, c(0.025, 0.975))))
    CI0b <- cbind(Est=ci0b[,1], t(apply(ci0b, 1, quantile, c(0.025, 0.975))))
    ref1 <- FALSE
    if (ref1) {

        p00 <- CI00[ii,]
        p0b <- CI0b[ii,]

        yc00 <- YYmean / p00
        yc0b <- YYmean / p0b

        Ref <- yc00[3,1]
        yc00 <- yc00 / Ref
        yc0b <- yc0b / Ref
    } else {
        dp00 <- ci00[ii,] / YYmean
        dp00 <- t(dp00) / dp00[3,]
        dp0b <- ci0b[ii,] / YYmean
        dp0b <- t(dp0b) / dp0b[3,]
        yc00 <- cbind(Est=dp00[1,], t(apply(dp00, 2, quantile, c(0.025, 0.975))))
        yc0b <- cbind(Est=dp0b[1,], t(apply(dp0b, 2, quantile, c(0.025, 0.975))))
    }


    col <- "#80808080"
    plot(0, type="n", ylim=c(0,1), xlim=c(0,10),
        xlab="Duration (min)", ylab="P(availability)")
    polygon(c(tt, rev(tt)), c(CI00[,2], rev(CI00[,3])), border=col, col=col)
    polygon(c(tt, rev(tt)), c(CI0b[,2], rev(CI0b[,3])), border=col, col=col)
    lines(tt, CI00[,1], col=1, lty=1, lwd=1.5)
    lines(tt, CI0b[,1], col=1, lty=2, lwd=1.5)
    abline(v=c(3,5),lty=3)
    legend("topleft", bty="n", lty=c(1,2), col=1,
        title=paste0(spp, " (n=", nn, ")"), lwd=1.5,
        legend=c(expression(M[0]), expression(M[f])))

    ylim <- range(yc00, yc0b)
    if (ylim[1] >= 0.8)
        ylim[1] <- 0.8
    if (ylim[2] <= 1.2)
        ylim[2] <- 1.2
    plot(0, type="n", ylim=ylim, xlim=c(0,10.2),
        xlab="Duration (min)", ylab="Corrected relative count")
    if (ref1) {
        polygon(c(-1, 11, 11, -1), c(yc00[3,2], yc00[3,2], yc00[3,3], yc00[3,3]),
            border=col, col=col)
        polygon(c(-1, 11, 11, -1), c(yc0b[3,2], yc0b[3,2], yc0b[3,3], yc0b[3,3]),
            border=col, col=col)
        ss <- 1:3
    } else {
        ss <- 1:2
    }
    segments(x0=c(3,5,10)[ss]-0.15, y0=yc00[ss,2], y1=yc00[ss,3], col=1, lwd=1.5)
    segments(x0=c(3,5,10)[ss]+0.15, y0=yc0b[ss,2], y1=yc0b[ss,3], col=1, lwd=1.5, lty=1)
    abline(v=c(3,5),lty=3)
    abline(h=1)
    points(c(3,5,10)[ss]-0.15, yc00[ss,1], col=1, cex=1.2, pch=19)
    points(c(3,5,10)[ss]+0.15, yc0b[ss,1], col=1, cex=1.2, pch=21)
    legend("topleft", bty="n", pch=c(19,21), col=1, title=spp, lty=c(1,1), lwd=1.5,
        legend=c(expression(M[0]), expression(M[f])))
}
par(op)
dev.off()

## stats for the figures
cf <- lapply(c("CONW", "WEWP", "RUBL"), function(z) attr(resx[[z]], "est"))
names(cf) <- c("CONW", "WEWP", "RUBL")

ff <- function(x, R=10^4) {
    z1 <- exp(rnorm(R, x$m0cf, x$m0vc[1,1]))
    z2 <- MASS::mvrnorm(R, x$mbcf, x$mbvc)
    z2[,1] <- exp(z2[,1])
    z2[,2] <- plogis(z2[,2])
    c(m0_phi=exp(x$m0cf), M0_phi_SE=sd(z1),
      mb_phi=exp(x$mbcf[1]), mb_phi_SE=sd(z2[,1]),
      mb_c=plogis(x$mbcf[2]), mb_phi_SE=sd(z2[,2]),
      mv_cor=cor(z2)[2,1])
}
round(sapply(cf, ff), 3)


## time varying responses ---------------------------------

library(MASS)
library(KernSmooth)

pf <- function(var, mod, n=10^4, resol=0.1) {
    if (mod == "0")
        sppPred <- sppPred0
    if (mod == "b")
        sppPred <- sppPredb
    if (mod == "f")
        sppPred <- sppPredf
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

for (WHAT in c("TSSR","JDAY","TSLS")) {
    cat(WHAT, "\n");flush.console()

    Xpk2 <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), pkDur[iii,])
    sppPred0 <- matrix(NA, nrow(Xpk2), length(SPPfull))
    colnames(sppPred0) <- SPPfull
    sppPredb <- matrix(NA, nrow(Xpk2), length(SPP))
    colnames(sppPredb) <- SPP
    sppPredf <- sppPredb

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

    for (spp in mSPP) {
        bestf <- as.character((0:14)[which.max(waicf[spp,])])
        if (bestf %in% COOL) {
            cff2 <- .BAMCOEFSfmix$sra_estimates[[spp]][[bestf]]$coefficients
            c2 <- plogis(cff2["logit.c"])
            logphi2 <- cff2[-length(cff2)]
            sppPredf[,spp] <- 1-c2*exp(-TT*exp(drop(Xpk2[,gsub("log.phi_", "",
                names(logphi2)),drop=FALSE] %*% logphi2)))
        }
    }
    sppPredf <- sppPredf[,colSums(is.na(sppPredf))==0]

    gc()
    b0 <- pf(WHAT, "0", n=100000, resol=RESOL)
    gc()
    bb <- pf(WHAT, "b", n=100000, resol=RESOL)
    gc()
    bf <- pf(WHAT, "f", n=100000, resol=RESOL)

OUT[[WHAT]] <- list("0"=b0, "b"=bb, "f"=bf)
}

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

pdf(file.path(ROOT2, "internal", paste0("Fig5-responses", TT, "min-bkde2d.pdf")),
    height=10, width=8)
op <- par(mfrow=c(3,3), las=1, mar=c(5, 4, 4, 2)+0.1, oma=c(2,6,1,2), xpd=NA)

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
plf2(OUT[["TSSR"]][["f"]], axes=FALSE, ann=FALSE)
axis(1)
axis(2)
title(ylab="P(availability)", xlab="Time since sunrise (h)")
text(-6.5, 1.2, expression(M[f]^varphi), cex=1.25)
par(mar=c(5, 0, 4, 0))
plf2(OUT[["JDAY"]][["f"]], axes=FALSE, ann=FALSE)
axis(1)
title(xlab="Ordinal day")
par(mar=c(5, 0, 4, 0))
plf2(OUT[["TSLS"]][["f"]], axes=FALSE, ann=FALSE)
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



## compare PIF time adjustment -----------------------------

sptab <- read.csv("~/GoogleWork/bam/duration_ms/revisionMarch2017/internal/Appendix-table.csv")
rownames(sptab) <- sptab$spp

spp <- "OVEN"
TT <- 3
psumm <- list()
for (spp in SPPfull) {
    cat(spp, "\n");flush.console()
    bmod <- strsplit(as.character(bestx[spp]), "_")[[1]]
    if (all(is.na(bmod)))
        bmod <- c("M0", "0")
    if (bmod[1] == "M0" && bmod[2] == "0") {
        cf00 <- .BAMCOEFSrem$sra_estimates[[spp]][["0"]]$coefficients
        p <- 1-exp(-TT*exp(cf00))
    } else {
        Xpk2 <- model.matrix(~JDAY + I(JDAY^2) + TSSR + I(TSSR^2) + TSLS + I(TSLS^2), pkDur[iii,])
        if (bmod[1] == "M0") {
            cf02 <- .BAMCOEFSrem$sra_estimates[[spp]][[bmod[2]]]$coefficients
            p <- 1-exp(-TT*exp(drop(Xpk2[,gsub("log.phi_", "",
                names(cf02)),drop=FALSE] %*% cf02)))
        }
        if (bmod[1] == "Mb") {
            cfb2 <- .BAMCOEFSmix$sra_estimates[[spp]][[bmod[2]]]$coefficients
            p <- 1-plogis(drop(Xpk2[,gsub("logit.c_", "",
                names(cfb2)[-1]),drop=FALSE] %*%
                cfb2[-1])) * exp(-TT*exp(cfb2[1]))
        }
        if (bmod[1] == "Mf") {
            cff2 <- .BAMCOEFSfmix$sra_estimates[[spp]][[bmod[2]]]$coefficients
            c2 <- plogis(cff2["logit.c"])
            logphi2 <- cff2[-length(cff2)]
            p <- 1-c2*exp(-TT*exp(drop(Xpk2[,gsub("log.phi_", "",
                names(logphi2)),drop=FALSE] %*% logphi2)))
        }
    }
    psumm[[spp]] <- summary(p)
}
pp <- do.call(rbind, psumm)

sptab$Inv_p3best <- 1/pp[rownames(sptab),"Mean"]
sptab$Inv_p3bestadj <- 1/(pp[rownames(sptab),"Mean"]/pp[rownames(sptab),"Max."])
ta <- sptab[!is.na(sptab$tadj),c("tadj","Inv_p3best","Inv_p3bestadj")]
ta <- ta[ta$Inv_p3bestadj < 5 & ta$tadj < 5,]

write.csv(sptab, row.names=FALSE,
    file="~/GoogleWork/bam/duration_ms/revisionMarch2017/internal/Appendix-table-pif.csv")

spt <- read.csv("~/GoogleWork/bam/duration_ms/revisionMarch2017/internal/Appendix-table-pif.csv")
rownames(spt) <- spt$spp

x <- spt[!is.na(spt$tadj),c("tadj","Inv_p3best","Inv_p3bestadj")]
colnames(x) <- c("Tadj", "Uinv", "Uadj")
summary(x)
max(x$Tadj)
x <- x[x$Uadj <= max(x$Tadj),]

cor.test(x$Tadj, x$Uinv, method="spearman")
cor.test(x$Tadj, x$Uadj, method="spearman")

plot(x$Tadj, x$Uinv)

summary(x$Uinv/x$Tadj)
summary(x$Uadj/x$Tadj)

## Migratory status

if (FALSE) {
tb <- read.csv("~/GoogleWork/bam/duration_ms/revisionMarch2017/internal/Appendix-table.csv")
rownames(tb) <- tb$spp
library(lhreg)
data(lhreg_data)
tb$Mig <- lhreg_data$Mig[match(rownames(tb), lhreg_data$spp)]
tb <- write.csv(tb, row.names=FALSE,
    "~/GoogleWork/bam/duration_ms/revisionMarch2018/Appendix-table.csv")
tb[is.na(tb$Mig),]
}

tb <- read.csv("~/GoogleWork/bam/duration_ms/revisionMarch2018/Appendix-table.csv")
rownames(tb) <- tb$spp
tb$b3m <- sapply(strsplit(as.character(tb$Best3), "_"),
    function(z) if (length(z) < 2 && is.na(z)) -1 else as.integer(z[2]))
tb$m3ts <- ifelse(tb$b3m %in% (grep("TSSR", sapply(NAMES, paste, collapse=" "))-1), 1, 0)
tb$m3jd <- ifelse(tb$b3m %in% (grep("JDAY", sapply(NAMES, paste, collapse=" "))-1), 1, 0)
tb$m3ls <- ifelse(tb$b3m %in% (grep("TSLS", sapply(NAMES, paste, collapse=" "))-1), 1, 0)
tb$Mig2 <- tb$Mig
levels(tb$Mig2)[levels(tb$Mig2) %in% c("SD","LD")] <- "MI"

with(tb[!is.na(tb$Mb_phi),], table(m3jd, m3ls))
addmargins(with(tb[!is.na(tb$Mb_phi),], table(b3m, Mig2)))

addmargins(with(tb[!is.na(tb$Mb_phi),], table(m3jd, Mig2)))
with(tb[!is.na(tb$Mb_phi),], chisq.test(m3jd, Mig2))

addmargins(with(tb[!is.na(tb$Mb_phi),], table(m3ls, Mig2)))
with(tb[!is.na(tb$Mb_phi),], chisq.test(m3ls, Mig2))

## JDAY
prop.test(c(62, 13), c(62+62, 10+13))

## DSLS
prop.test(c(86, 13), c(86+38, 10+13))


xtfun2 <- function(spp, ymin=1) {
    Y <- groupSums(xtDur[[spp]][rn,], 2,
        c("0-3", "xxx", "0-3", "0-3", "xxx", "xxx", "0-3", "0-3",
        "xxx", "3-5", "3-5", "xxx", "xxx", "3-5", "5-10", "5-10",
        "5-10", "5-10", "xxx", "5-10", "5-10", "5-10", "5-10", "5-10"))
    Y <- as.matrix(Y)[,c("0-3","3-5","5-10")]
    YY <- cbind("0-3"=Y[,1], "0-5"=Y[,1]+Y[,2], "0-10"=rowSums(Y))
    YY <- YY[YY[,3] >= ymin,,drop=FALSE]
    colMeans(YY)
}
yyy <- t(pbsapply(SPP, xtfun2, ymin=1))
yyyy <- yyy/yyy[,3]

## addressing non-independence of revisits

library(detect)
load("~/GoogleWork/bam/duration_ms/pkResDur_data.Rdata")
pkDur$SSYR <- interaction(pkDur$SS, pkDur$YEAR, drop=TRUE)
f <- function(fa) {
    fa <- as.integer(droplevels(as.factor(fa)))
    o <- seq_along(fa)
    r <- sample(o)
    keep <- !duplicated(fa[r])
    r[keep]
}
g <- function(m00, m0f) {
    V <- vcov(m0f)
    unname(c(coef(m00), coef(m0f),
        sqrt(vcov(m00)[1]), sqrt(diag(V)), cov2cor(V)[2,1]))
}

SPP <- names(resDurOK)
B <- 200
nmax <- 500

OUT <- list()
for (spp in SPP) {
    z <- resDurOK[["BTNW"]]
    x <- droplevels(pkDur[z$pkey,])
    n <- min(nlevels(x$SS), nmax)
    if (n >= 200) {

        out0 <- matrix(0, B, 7)
        colnames(out0) <- c("logphi0", "logphi", "logitc",
                            "SE_logphi0", "SE_logphi", "SE_logitc", "cor")
        out2 <- out1 <- out0

        for (j in 1:B) {
            cat(spp, j, "/", B, "\n");flush.console()
            i0 <- sample.int(z$n, n, replace=FALSE) # ss/yr/visit
            i1 <- f(x$SSYR)[1:n]                 # ss/yr
            i2 <- f(x$SS)[1:n]                   # ss

            m00 <- cmulti(z$Y[i0,] | z$D[i0,] ~ 1, type="rem")
            m0f <- cmulti(z$Y[i0,] | z$D[i0,] ~ 1, type="mix")
            m10 <- cmulti(z$Y[i1,] | z$D[i1,] ~ 1, type="rem")
            m1f <- cmulti(z$Y[i1,] | z$D[i1,] ~ 1, type="mix")
            m20 <- cmulti(z$Y[i2,] | z$D[i2,] ~ 1, type="rem")
            m2f <- cmulti(z$Y[i2,] | z$D[i2,] ~ 1, type="mix")

            out0[j,] <- g(m00, m0f)
            out1[j,] <- g(m10, m1f)
            out2[j,] <- g(m20, m2f)
        }

        OUT <- array(0, c(7, 3, 3))
        dimnames(OUT) <- list(colnames(out0), c("50%", "2.5%", "97.5"),
            c("ss/yr/visit", "ss/yr", "ss"))
        OUT[,,1] <- t(apply(out0, 2, quantile, c(0.5, 0.025, 0.975)))
        OUT[,,2] <- t(apply(out1, 2, quantile, c(0.5, 0.025, 0.975)))
        OUT[,,3] <- t(apply(out2, 2, quantile, c(0.5, 0.025, 0.975)))
    }
}

par(mfrow=c(2,3))
for (k in colnames(out0)[1:6]) {
    out <- cbind(out0[,k], out1[,k], out2[,k])
    #summary(out)

    d0 <- density(out[,1])
    d1 <- density(out[,2])
    d2 <- density(out[,3])
    ylim <- range(d0$y, d1$y, d2$y)
    xlim <- range(d0$x, d1$x, d2$x)
    plot(d0, xlim=xlim, ylim=ylim, col=1, main=k)
    lines(d1, xlim=xlim, ylim=ylim, col=2)
    lines(d2, xlim=xlim, ylim=ylim, col=4)
}
