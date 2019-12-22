#' ---
#' title: "Analysis data NWT CONI models"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "`r as.Date(Sys.time())`"
#' output: pdf_document
#' ---
#+ echo=FALSE
#knitr::opts_chunk$set(eval=FALSE)
#ROOT <- "~/GoogleWork/collaborations/coni-nwt"
par(las=1)
#'
#' # Load packages
#'
library(mgcv)
library(gbm)
library(raster)
library(dismo)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
#'
#' # Load data
#'
#' Observations at `key` level: station+date+time
x <- read.csv("~/GoogleWork/collaborations/coni-nwt/CONI-AB-NWT-withOffsets.csv")
x$key <- gsub(" ", "+", paste0(x$station, "+", x$datetime))
#' station level data turned into a spatial object: more efficient to extract
#' at station level and repeat the values
s <- x[!duplicated(x$station),]
rownames(s) <- s$station
coordinates(s) <- ~long+lat
proj4string(s) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#'
#' # Extract GIS info
#'
#' This is the stack we will use for prediction
ROOT <- "d:/bam/BAM_data_v2019/gnm"
pr60 <- stack(file.path(ROOT, "data", "subunits", paste0("bcr60all_1km.grd")))
pr61 <- stack(file.path(ROOT, "data", "subunits", paste0("bcr61all_1km.grd")))
pr <- mosaic(pr60, pr61, fun=mean)

#' Reproject station object
s <- spTransform(s, proj4string(pr))
pt <- extract(pr, s, method="simple")
rownames(pt) <- rownames(s@data)
sum(is.na(pt))

#' Check data for layers to drop: constant varoables (SD=0)
SD <- apply(pt, 2, sd)
ncol(pt) # 222
sum(SD > 0) # down to 111
#' Check data for layers to drop: highly correlated variables (r_max=0.9)
get_cn <- function(z, rmax=0.9) {
    SD <- apply(z, 2, sd)
    COR <- cor(z[,SD > 0])
    cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
    cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
    cr[,1] <- as.character(cr[,1])
    cr[,2] <- as.character(cr[,2])
    cr$Landsc1 <- startsWith(cr[,1], "Landsc750_")
    cr$Landsc2 <- startsWith(cr[,2], "Landsc750_")
    cr1 <- cr[cr$Landsc1 | cr$Landsc2,]
    cr2 <- cr[!(cr$Landsc1 | cr$Landsc2),]
    while(any(abs(cr1$dist) > rmax)) {
        i <- if (cr1$Landsc1[1])
            cr1[1,1] else cr1[1,2]
        j <- cr1[,1] == i | cr1[,2] == i
        cr1 <- cr1[!j,]
    }
    cr3 <- rbind(cr1, cr2)
    cr3 <- cr3[order(abs(cr3$dist), decreasing = TRUE),]
    while(any(abs(cr3$dist) > rmax)) {
        i <- if (cr3$Landsc1[1])
            cr3[1,1] else cr3[1,2]
        j <- cr3[,1] == i | cr3[,2] == i
        cr3 <- cr3[!j,]
    }
    union(as.character(cr3[,1]), as.character(cr3[,2]))
}
#' Keep NALC and lf, these are discrete
CN <- unique(c("nalc", "lf", get_cn(pt[,SD > 0], rmax=0.9)))
length(CN) # down to 50

#' Organizing into a data frame, set NALC and lf as factors
dat <- data.frame(
    x[,c("detection", "Offset", "station", "datetime","long","lat")],
    pt[match(x$station, rownames(pt)),CN])
dat$nalc <- as.factor(dat$nalc) # !!! fixme -- strange values
nlevels(dat$nalc)
dat$lf <- as.factor(dat$lf)
nlevels(dat$lf)
str(dat)
save(dat, x, file="~/GoogleWork/collaborations/coni-nwt/CONI-AB-NWT-datawithpredictors.RData")

#' Variable importance based on cross-validated BRT.
#' Note: this is a crude model with logit link, just for exploration
b <- gbm.step(dat, gbm.y = 1, gbm.x = 5:ncol(dat),
    family = "bernoulli", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.5,
    offset=dat$Offset)

vi <- b$contributions
rownames(vi) <- NULL
vi$cumsum <- cumsum(vi$rel.inf)
vi

## analysis

#    Randomly selected 10% of the recordings scanned
#    Filtered out recordings that had booms detected in them
#    Reviewed all recognizer detections for those recordings
#    Assigned individual ID to each boom detection in those recordings
# There's definitely more than one individual in some of these recordings.
# Here's the distribution of the results:
#    1 CONI: 29 sites
#    2 CONI: 11 sites
#    3 CONI: 2 sites
#    4 CONI: 1 site
# The one thing to keep in mind is that these are almost entirely
# from the Alberta fire site. I suspect the densities will be lower at the
# other province/treatment combos.

library(mefa4)
library(lme4)
library(pbapply)

load("~/GoogleWork/collaborations/coni-nwt/CONI-AB-NWT-datawithpredictors.RData")
x <- read.csv("~/GoogleWork/collaborations/coni-nwt/CONI-AB-NWT-withOffsets.csv")
x <- nonDuplicated(x, station, TRUE)
dat$lat <- x$lat[match(dat$station, x$station)]
dat$long <- x$long[match(dat$station, x$station)]
rm(x)

str(dat)

tmp <- strsplit(as.character(dat$station), "-")
n <- sapply(tmp, length)
table(n)
head(tmp[n==3])
head(tmp[n==4])
tmp[n==5]
ss1 <- as.factor(sapply(tmp, function(z) paste0(z[1])))
table(table(ss1))
ss2 <- as.factor(sapply(tmp, function(z) paste0(z[1], "-", z[2])))
table(table(ss2))
ss3 <- as.factor(sapply(tmp, function(z) paste0(z[1], "-", z[2], "-", z[3])))
table(table(ss3))

plot(dat$long, dat$lat, col=ss3)

yy <- rep(1:4, c(29, 11, 2, 1)) # mean 1.418605

levs <- c("1"="ConTmp",
    "2"="ConTai",
    "5"="Dec",
    "6"="Mix",
    "8"="Shr",
    "10"="Grs",
    "14"="Wet",
    "16"="Bar",
    "18"="Wat")
dat$alc <- dat$nalc
levels(dat$alc) <- levs[levels(dat$nalc)]

levs2 <- c("1"="ConTmp",
    "2"="ConTai",
    "5"="DM",
    "6"="DM",
    "8"="OP",
    "10"="OP",
    "14"="OP",
    "16"="OP",
    "18"="OP")
dat$lc <- dat$nalc
levels(dat$lc) <- levs2[levels(dat$lc)]

addmargins(table(dat$alc, dat$lc)/30)

m0 <- glmer(detection ~ (1 | station), data=dat,
    offset=dat$Offset, family=binomial("cloglog"))
m1 <- glmer(detection ~ lc + (1 | station), data=dat,
    offset=dat$Offset, family=binomial("cloglog"))
g0 <- glm(detection ~ 1, data=dat,
    offset=dat$Offset, family=binomial("cloglog"))
g1 <- glm(detection ~ lc, data=dat,
    offset=dat$Offset, family=binomial("cloglog"))
AIC(m0, m1, g0, g1)

round(sort(exp(fixef(m1))), 3)
round(sort(exp(coef(g1))), 3)

#d <- mefa4::sum_by(dat$detection, dat$station)

f <- function() {
    d <- nonDuplicated(dat[sample.int(nrow(dat)),], station, TRUE)
    m <- glm(detection ~ lc-1, data=d,
        offset=d$Offset, family=binomial("cloglog"))
    coef(m)
}

B <- 1000
res <- t(pbreplicate(B, f()))

cf <- fixef(m1)
est <- exp(c(cf[1], cf[1]+cf[-1]))
names(est) <- levels(dat$lc)
est

cf <- coef(g1)
est2 <- exp(c(cf[1], cf[1]+cf[-1]))
names(est2) <- levels(dat$lc)
est2

boxplot(exp(res))
points(1:4,est, pch=4, col=2, cex=2)
points(1:4,est2, pch=4, col=4, cex=2)
