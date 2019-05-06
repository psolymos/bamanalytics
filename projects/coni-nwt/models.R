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
i <- 6 # BCR 6
ROOT <- "d:/bam/BAM_data_v2019/gnm"
pr <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "_1km.grd")))
#prc <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "cat_1km.grd")))
prc <- stack(file.path(ROOT, "data", "_new", "cat2", paste0("bcr", i, "cat_1km.grd")))
pr <- pr[[which(names(pr) != "bcr")]]
prc <- prc[[which(names(prc) != "bcr")]]
pr <- addLayer(pr, prc)
rm(prc)
#' Reproject station object
s <- spTransform(s, proj4string(pr))
pt <- extract(pr, s)
rownames(pt) <- rownames(s@data)
#' Check data for layers to drop: constant varoables (SD=0)
SD <- apply(pt, 2, sd)
ncol(pt) # 147
sum(SD > 0) # down to 73
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
length(CN) # down to 46

#' Organizing into a data frame, set NALC and lf as factors
dat <- data.frame(
    x[,c("detection", "Offset", "station", "datetime")],
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


