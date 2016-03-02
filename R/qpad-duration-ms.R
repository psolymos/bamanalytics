
## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
ROOT2 <- "~/Dropbox/bam/duration_ms/revisionOct2015"

B <- 1000

## Load required packages
library(mefa4)
library(pbapply)
library(detect)
library(MASS)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2015-10-08.Rdata"))

### Removal sampling

## non NA subset for duration related estimates
pkDur <- dat[,c("PKEY","JDAY","TSSR","TSLS","DURMETH","PCODE")]
pkDur <- droplevels(pkDur[rowSums(is.na(pkDur)) == 0,])
sort(table(pkDur$PCODE))
## only 1 obs in project (all others are >20)
pkDur <- droplevels(pkDur[pkDur$PCODE != "GLDRPCCL01",])
## strange methodology where all counts have been filtered
## thus this only leads to 0 total count and exclusion
pkDur <- droplevels(pkDur[pkDur$DURMETH != "J",])


## crosstab for species
xtDur <- Xtab(ABUND ~ PKEY + dur + SPECIES, pc)
xtDur[["NONE"]] <- NULL

load_BAM_QPAD(1)
.BAMCOEFS$version
load(file.path(ROOT2, "BAMCOEFS_QPAD_v3.rda"))
.BAMCOEFS$version

e <- new.env()
load(file.path(ROOT2, "BAMCOEFS_QPAD_v3_mix.rda"), envir=e)
.BAMCOEFSmix <- e$.BAMCOEFS
.BAMCOEFSmix$version

aic0 <- .BAMCOEFS$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))

aic0 <- aic0[SPP,]
aicb <- aicb[SPP,]
aic <- cbind(aic0[SPP,], aicb[SPP,])

best0 <- as.character(0:14)[apply(aic0, 1, which.min)]
bestb <- as.character(0:14)[apply(aicb, 1, which.min)]
best <- colnames(aic)[apply(aic, 1, which.min)]
names(best) <- names(best0) <- names(bestb) <- SPP

rn <- intersect(rownames(pkDur), rownames(xtDur[[1]]))
pkDur <- droplevels(pkDur[rn,])

## random PCODE
set.seed(1000)
pkDur$PCODErnd <- pkDur$PCODE[sample.int(nrow(pkDur))]

D <- ltdur$end[match(pkDur$DURMETH, rownames(ltdur$end)),]
## exclude 0 sum and <1 interval rows
nOK <- table(PC=pkDur$PCODE, n_int=rowSums(!is.na(D)))

nPC <- rowSums(nOK) - nOK[,"1"]
nPC <- sort(nPC[nPC > 0])
PC <- names(nPC)

## read all the suff, read in lines 1-86 (RND defs there!)

## random PCODE IDs or not (projects as they are)
isRND <- FALSE

## rem -- subsets are for sets of species to be combined
e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-1.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-1rnd.Rdata"), envir=e)
resDurBAMless1 <- e$resDurBAMless1
resDurPcode1 <- e$resDurPcode1

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-2.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-2rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-3.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-3rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-rem-4.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-rem-4rnd.Rdata"), envir=e)
resDurBAMless1 <- c(resDurBAMless1, e$resDurBAMless1)
resDurPcode1 <- c(resDurPcode1, e$resDurPcode1)

## mix -- subsets are for sets of species to be combined
e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-1.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-1rnd.Rdata"), envir=e)
resDurBAMless1_mix <- e$resDurBAMless1_mix
resDurPcode1_mix <- e$resDurPcode1_mix

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-2.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-2rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-3.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-3rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

e <- new.env()
if (!isRND)
    load(file.path(ROOT2, "xval-mix-4.Rdata"), envir=e)
if (isRND)
    load(file.path(ROOT2, "xval-mix-4rnd.Rdata"), envir=e)
resDurBAMless1_mix <- c(resDurBAMless1_mix, e$resDurBAMless1_mix)
resDurPcode1_mix <- c(resDurPcode1_mix, e$resDurPcode1_mix)

ff <- list(
        "0"=~ 1,
        "1"=~ JDAY,
        "2"=~ TSSR,
        "3"=~ JDAY + I(JDAY^2),
        "4"=~ TSSR + I(TSSR^2),
        "5"=~ JDAY + TSSR,
        "6"=~ JDAY + I(JDAY^2) + TSSR,
        "7"=~ JDAY + TSSR + I(TSSR^2),
        "8"=~ JDAY + I(JDAY^2) + TSSR + I(TSSR^2),
        "9"=~ TSLS,
        "10"=~ TSLS + I(TSLS^2),
        "11"=~ TSLS + TSSR,
        "12"=~ TSLS + I(TSLS^2) + TSSR,
        "13"=~ TSLS + TSSR + I(TSSR^2),
        "14"=~ TSLS + I(TSLS^2) + TSSR + I(TSSR^2))
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

aic_fun <- function(x) {
#    if (inherits(x, "try-error"))
    if (is.character(x))
        Inf else -2*x$loglik + 2*x$p
}
coef_fun <- function(x, id) {
#    if (inherits(x, "try-error"))
    if (is.character(x[[id]]))
        NA else x[[id]]$coefficients
}
vcov_fun <- function(x, id) {
#    if (inherits(x, "try-error"))
    if (is.character(x[[id]]))
        NA else x[[id]]$vcov
}
waic_fun <- function(z) {
    dAIC <- z - min(z)
    w <- exp(-dAIC/2) 
    w/sum(w)
}

rownames(TAX) <- TAX$Species_ID

## when does model fail?

problem <- list()
for (spp in names(resDurPcode1)) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        yy <- rowSums(xtDur[[spp]][ii, ])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                logphi=NA, se_logphi=NA,
                    nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=NA, se_logphi=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=unname(tmp3$coefficients), 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problem0 <- do.call(rbind, problem)

problem <- list()
for (spp in names(resDurPcode1_mix)) {
    cat(spp, "\n")
    flush.console()
    tmp <- resDurPcode1_mix[[spp]]
    for (PCi in names(tmp)) {
        tmp2 <- tmp[[PCi]]
        ii <- rownames(pkDur)[pkDur$PCODE==PCi]
        yy <- rowSums(xtDur[[spp]][ii, ])
        if (inherits(tmp2, "try-error")) {
            out <- data.frame(spp=spp, pc=PCi, 
                ntot=length(yy),
                ndet=sum(yy > 0), 
                n1=sum(yy > 1), 
                n2=sum(yy > 2), 
                logphi=NA, 
                logitc=NA, 
                se_logphi=NA,
                se_logitc=NA,
                cor=NA,
                nobs=NA,
                msg=as.character(tmp2))
        } else {
            tmp3 <- tmp2[["0"]]
            if (inherits(tmp3, "try-error") || is.character(tmp3)) {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=NA, 
                    logitc=NA, 
                    se_logphi=NA,
                    se_logitc=NA,
                    cor=NA,
                    nobs=NA,
                    msg=as.character(tmp3))
            } else {
                out <- data.frame(spp=spp, pc=PCi, 
                    ntot=length(yy),
                    ndet=sum(yy > 0), 
                    n1=sum(yy > 1), 
                    n2=sum(yy > 2), 
                    logphi=unname(tmp3$coefficients)[1], 
                    logitc=unname(tmp3$coefficients)[2], 
                    se_logphi=sqrt(tmp3$vcov[1,1]),
                    se_logitc=sqrt(tmp3$vcov[2,2]),
                    cor=cov2cor(tmp3$vcov)[1,2],
                    nobs=tmp3$nobs,
                    msg="")
            }
        }
        problem[[paste(spp, PCi, sep=":")]] <- out
    }
}
problemb <- do.call(rbind, problem)

save(problem0, problemb, file=file.path(ROOT2, "problem.Rdata"))

ss <- c("0 observation with multiple duration (1)", "1 observation with multiple duration (2)")

dat0 <- problem0[!(problem0$msg %in% ss),]
dat0$okfit <- ifelse(is.na(dat0$logphi), 0, 1)
dat0$okse <- ifelse(is.na(dat0$se_logphi), 0, 1)
dat0$msg <- NULL
mod0 <- glm(okfit ~ log(ndet+1), dat0, family=binomial)
mod0se <- glm(okse ~ log(ndet+1), dat0, family=binomial)

datb <- problemb[!(problemb$msg %in% ss),]
datb$okfit <- ifelse(is.na(datb$logphi), 0, 1)
datb$okse <- ifelse(is.na(datb$se_logphi), 0, 1)
#datb$okcor <- ifelse(abs(datb$cor) > 0.99, 0, 1)
datb$msg <- NULL
modb <- glm(okfit ~ log(ndet+1), datb, family=binomial)
modbse <- glm(okse ~ log(ndet+1), datb, family=binomial)
#modbcor <- glm(okcor ~ log(ndet+1), datb, family=binomial)

#modb <- glm(okfit ~ log(n1+1), datb, family=binomial)
#modbse <- glm(okse ~ log(n1+1), datb, family=binomial)
#modb <- glm(okfit ~ log(n2+1), datb, family=binomial)
#modbse <- glm(okse ~ log(n2+1), datb, family=binomial)

nmax <- 100
ndat <- data.frame(n=2:nmax, 
    p0=plogis(coef(mod0)[1] + coef(mod0)[2]*log(1 + 2:nmax)),
    pb=plogis(coef(modb)[1] + coef(modb)[2]*log(1 + 2:nmax)),
    p0se=plogis(coef(mod0se)[1] + coef(mod0se)[2]*log(1 + 2:nmax)),
    pbse=plogis(coef(modbse)[1] + coef(modbse)[2]*log(1 + 2:nmax)),
    pbcor=plogis(coef(modbcor)[1] + coef(modbcor)[2]*log(1 + 2:nmax)))

plot(ndat[,1], ndat[,2], col=2, type="l", ylim=c(0,1), lwd=2,
    xlab="Number of >0 survey counts", ylab="Probability",
    xlim=c(0,nmax))
lines(ndat[,1], ndat[,3], col=4, lwd=2)
lines(ndat[,1], ndat[,4], col=2, lwd=2, lty=2)
lines(ndat[,1], ndat[,5], col=4, lwd=2, lty=2)
#lines(ndat[,1], ndat[,6], col=4, lwd=2, lty=3)
abline(h=0.9, lty=1)
abline(v=ndat[which.min(abs(ndat[,2]-0.9)),1], col=2)
abline(v=ndat[which.min(abs(ndat[,3]-0.9)),1], col=4)
#abline(v=ndat[which.min(abs(ndat[,4]-0.9)),1], col=2, lty=2)
#abline(v=ndat[which.min(abs(ndat[,5]-0.9)),1], col=4, lty=2)
#legend("bottomright", col=c(2,2,4,4,4), lty=c(1,2,1,2,3), lwd=2, 
#    legend=c("m0 fit", "m0 SE", "mb fit", "mb SE", "mb cor"), bty="n")
legend("bottomright", col=c(2,2,4,4), lty=c(1,2,1,2), lwd=2, 
    legend=c("m0 fit", "m0 SE", "mb fit", "mb SE"), bty="n")

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


## sra m0 vs mb models


aic0 <- .BAMCOEFS$sra_aic
aicb <- .BAMCOEFSmix$sra_aic
colnames(aic0) <- paste0("m0_", colnames(aic0))
colnames(aicb) <- paste0("mb_", colnames(aicb))
SPP <- sort(intersect(rownames(aic0), rownames(aicb)))
aic <- cbind(aic0[SPP,], aicb[SPP,])
np <- sapply(SPP, function(z) selectmodelBAMspecies(z)$sra$nobs[1])

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

best <- colnames(aic)[apply(aic, 1, which.min)]
data.frame(table(best))
by(np, best %in% c("m0_0", "mb_0"), summary)


R <- 1000

cfall0 <- sapply(SPP, function(spp) exp(coefBAMspecies(spp, 0, 0)$sra))
dim(cfall0) <- c(length(SPP), 1)
dimnames(cfall0) <- list(SPP, "phi_0")
cfallb <- t(sapply(SPP, function(spp) {
    tmp <- unname(.BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients)
    c(phi_b=exp(tmp[1]), c=plogis(tmp[2]))
    }))

t <- seq(0, 10, 0.1)
pp0 <- sapply(SPP, function(spp) 1-exp(-t*cfall0[spp,"phi_0"]))
ppb <- sapply(SPP, function(spp) 1-cfallb[spp,"c"]*exp(-t*cfallb[spp,"phi_b"]))

f <- function(spp) {
    cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
    vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
    !inherits(try(mvrnorm(R, cfi0b, vci0b)), "try-error")
}
OK <- sapply(SPP, f)
table(OK)
SPP <- SPP[OK]

w00 <- waic0[SPP,"m0_0"]
w0b <- waicb[SPP,"mb_0"]
w0 <- pmax(w00, w0b)
H0 <- apply(waic0[SPP,], 1, function(z) sum(z^2))
Hb <- apply(waicb[SPP,], 1, function(z) sum(z^2))
H <- apply(waic[SPP,], 1, function(z) sum(z^2))
n <- np[SPP]

par(mfrow=c(1,3))
plot(n, w00, pch=21, cex=1+2*H0, log="x")
lines(n[order(n)], fitted(glm(w00 ~ n, family=binomial, weights=H0))[order(n)], col=2)
plot(n, w0b, pch=21, cex=1+2*Hb, log="x")
lines(n[order(n)], fitted(glm(w0b ~ n, family=binomial, weights=Hb))[order(n)], col=2)
plot(n, w0, pch=21, cex=1+2*H, log="x")
lines(n[order(n)], fitted(glm(w0 ~ n, family=binomial, weights=H))[order(n)], col=2)

## best model is 0
table((1:ncol(waic))[apply(waic[SPP,], 1, which.max)] <= 15)


pdf(paste0(ROOT2, "/m0-vs-mb.pdf"), onefile=TRUE, width=10, height=12)
for (spp in SPP) {

cat(spp, "\n");flush.console()

## model weights
wp <- waic0[spp,]
wq <- waicb[spp,]
names(wp) <- colnames(waic0)
names(wq) <- colnames(waicb)
wph <- waic[spp,names(wp)]
wqh <- waic[spp,names(wq)]
names(wp) <- 0:14
names(wq) <- 0:14


## covariate effects
mi <- bestmodelBAMspecies(spp, type="AIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

mib <- as.character(0:14)[which.min(aicb[spp,])]
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
ts <- seq(-0.3, 0.5, 0.01)
ls <- seq(-0.1, 0.4, len=length(jd))

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

Xp <- if (mi$sra %in% c("9","10","11","12","13","14"))
    Xp2 else Xp1
Xp <- Xp[,names(cfi$sra),drop=FALSE]
lphi1 <- drop(Xp %*% cfi$sra)
pmat <- matrix(exp(lphi1), length(jd), length(ts))
pmax <- sra_fun(10, max(exp(lphi1)))
pmat <- sra_fun(3, pmat)
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
cfi00 <- coefBAMspecies(spp, "0", mi$edr)$sra
vci00 <- drop(vcovBAMspecies(spp, "0", mi$edr)$sra)
phi00 <- exp(rnorm(R, cfi00, sqrt(vci00)))
ci00 <- sapply(phi00, function(z) 1-exp(-t*z))

cfi0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$coefficients
vci0b <- .BAMCOEFSmix$sra_estimates[[spp]][["0"]]$vcov
pcf1b <- mvrnorm(R, cfi0b, vci0b)
ci0b <- apply(pcf1b, 1, function(z) 1-plogis(z[2])*exp(-t*exp(z[1])))


op <- par(las=1, mfrow=c(3,2))

barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", np[spp], ") v", getBAMversion(),
    " w0=", round(sum(wph),2)),
    ylab="Model weight", xlab="Model ID")
barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
    main=paste0(spp, " (n=", np[spp], ") v", getBAMversion(),
    " wb=", round(sum(wqh),2)),
    ylab="Model weight", xlab="Model ID")

plot(t, pp0[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, pp0, col="grey", lwd=1, lty=1)
matlines(t, ci00, col="grey", lwd=1, lty=1)
lines(t, pp0[,spp], col=1, lwd=2)

plot(t, ppb[,spp], type="n", ylim=c(0,1),
     xlab="Point count duration (min)",
     ylab="Probability of singing")
#polygon(c(t, rev(t)), c(p[,2], rev(p[,3])),
#        col="grey", border="grey")
#matlines(t, ppb, col="grey", lwd=1, lty=1)
matlines(t, ci0b, col="grey", lwd=1, lty=1)
lines(t, ppb[,spp], col=1, lwd=2)

xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
    ls*365 else jd*365
image(xval, ts*24, pmat,
    col = rev(grey(seq(0, pmax, len=12))),
    xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"), 
        "Days since local springs", "Julian days"), 
    ylab="Hours since sunrise",
    main=paste("Best model:", mi$sra))
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





