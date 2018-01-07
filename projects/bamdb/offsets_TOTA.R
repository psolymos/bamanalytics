library(QPAD) # latest estimates
library(maptools) # for sunrise calculations
library(mefa4) # data manipulation

dat <- read.csv("~/Downloads/new.dat.for.offsets.csv")

## max duration
dat$MAXDUR <- dat$MaxDuration
table(dat$MAXDUR, useNA="always")

## max distance
dat$MAXDIS <- droplevels(dat$Maxdist)
levels(dat$MAXDIS) <- toupper(levels(dat$MAXDIS))
levels(dat$MAXDIS)[levels(dat$MAXDIS) == "UNLIMITED"] <- "Inf"
dat$MAXDIS <- as.numeric(as.character(dat$MAXDIS)) / 100
table(dat$MAXDIS, useNA="always")

## Date/time components
TM <- strptime(as.character(dat$StartTime.x), "%I:%M:%S %p")
chr <- paste0(dat$YYYY, "-", dat$MM.x, "-", dat$DD.x, " ", dat$StartTime.x)
DD <- strptime(chr, "%Y-%m-%e %I:%M:%S %p")
dat$DATE <- DD
summary(dat$DATE)
class(DD)
class(dat$DATE)

dat$JULIAN <- dat$DATE$yday
dat$JDAY <- DD$yday / 365
summary(dat$JDAY)
## prevent too far extrapolation
#dat$JDAY[dat$JDAY < 0.35 | dat$JDAY > 0.55] <- NA
hist(dat$JDAY)

Coor <- as.matrix(cbind(as.numeric(dat$POINT_X), as.numeric(dat$POINT_Y)))
JL <- as.POSIXct(DD)
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dat$srise <- NA
dat$srise[subset] <- sr
dat$start_time <- dat$DATE$hour + dat$DATE$min/60

dat$MDT_offset <- -6
dat$MDT_offset[dat$TZ == "America/Dawson_Creek"] <- -7
dat$MDT_offset[dat$TZ == "America/Vancouver"] <- -8
dat$MDT_offset <- dat$MDT_offset + 6
table(dat$MDT_offset, useNA="always")

dat$TSSR <- (dat$start_time - dat$srise + dat$MDT_offset) / 24
dat$TSSR_orig <- dat$TSSR # keep a full copy
dat$TSSR[dat$start_time > 12] <- NA # after noon
summary(dat$TSSR)
summary(dat$start_time)
hist(dat$TSSR)

dat$TREE <- dat$tree
summary(dat$TREE)
dat$TREE[dat$TREE > 100] <- NA
dat$TREE[dat$TREE < 0] <- NA
dat$TREE <- dat$TREE / 100
summary(dat$TREE)
hist(dat$TREE)

(ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv"))
table(dat$NALCMS05, useNA="always")
dat$NALCMS05[dat$NALCMS05 < 0] <- 0
compare_sets(dat$NALCMS05, ltnalc$Value)
dat$LCC2 <- reclass(dat$NALCMS05, ltnalc[,c("Value", "LCC2")], allow_NA=TRUE)
table(dat$NALCMS05, dat$LCC2, useNA="always")
dat$LCC4 <- reclass(dat$NALCMS05, ltnalc[,c("Value", "LCC4")], allow_NA=TRUE)
table(dat$NALCMS05, dat$LCC4, useNA="always")
boxplot(TREE ~ LCC4, dat)

load_BAM_QPAD(3)
SPP <- getBAMspecieslist()
dat$JDAY2 <- dat$JDAY^2
dat$TSSR2 <- dat$TSSR^2

Xp <- cbind("(Intercept)"=1, as.matrix(dat[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=dat$TREE,
    LCC2OpenWet=ifelse(dat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(dat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(dat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(dat$LCC4=="Wet", 1, 0))
OFF <- matrix(NA, nrow(dat), length(SPP))
rownames(OFF) <- dat$PKEY
colnames(OFF) <- SPP

## load TOTA results
load("~/Dropbox/Public/TOTA_BAMCOEFS_QPAD_v3.rda")
(mods <- getBAMmodellist())
(sra_mods <- names(mods$sra)[!grepl("DSLS", mods$sra)])
getBAMspecieslist()

spp <- "TOTA"
p <- rep(NA, nrow(dat))
A <- q <- p
## constant for NA cases
(cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0))))
## best model
(mi <- bestmodelBAMspecies(spp, model.sra=sra_mods, type="BIC"))
(cfi <- coefBAMspecies(spp, mi$sra, mi$edr))

Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0

p[!OKp] <- sra_fun(dat$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(dat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(dat$MAXDIS[!OKq], cf0[2]))

phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(dat$MAXDUR[OKp], phi1)
unlim <- ifelse(dat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(dat$MAXDIS[OKq], tau1))

ii <- which(p == 0)
p[ii] <- sra_fun(dat$MAXDUR[ii], cf0[1])

CORRECTION <- data.frame(p=p, A=A, q=q)
summary(CORRECTION)
OFFSET <- rowSums(log(CORRECTION))
summary(OFFSET)
dat[is.na(OFFSET),] # methodology is unknown --> drop!

out <- data.frame(dat$PKEY, Offset=OFFSET)
write.csv(out, row.names=FALSE, "~/Dropbox/Public/offsets_TOTA_20180106.csv")
