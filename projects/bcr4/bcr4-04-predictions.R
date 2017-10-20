library(mefa4)
library(pbapply)
library(sp)
library(raster)
library(rgdal)
library(rgeos)

## predictions

source("~/repos/mep/R/diagnostics-functions.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-models.R")
#source("~/repos/bamanalytics/projects/bcr4/bcr4-models2.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")
load("e:/peter/bam/bcr4/bcr4-data.RData")

Terms <- getTerms(mods, "list")
setdiff(Terms, colnames(DAT))
#xn <- DAT[BB[,1],Terms]
#Xn <- model.matrix(getTerms(mods, "formula"), xn)
#colnames(Xn) <- fixNames(colnames(Xn))

slp <- raster("e:/peter/bam/bcr4/alfresco_data/iem_prism_slope_1km.tiff")
load("e:/peter/bam/bcr4/bam_data/pred_grid_subset.RData")

## this is how you can take a set of x y coordinates with known projection
## and reproject according to the raster template's projection
xy <- Subset[,c("POINT_X", "POINT_Y")]
coordinates(xy) <- c("POINT_X", "POINT_Y")
proj4string(xy) <- "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
xy <- spTransform(xy, proj4string(slp))

## this is how you can set up a grid for intersections
## using a raster as template
pts <- as(slp, "SpatialPointsDataFrame")
PredData <- data.frame(coordinates(pts), Slope=pts@data[,1])
## if raster do not align, use the x y coordinates from pts
## to extract values for the same points

## this is how you can extract values from a raster given
## a set of xy coordinates
SLP <- extract(slp, xy)
Subset$slp <- sqrt(abs(SLP))
Subset$slp2 <- Subset$slp^2
Subset$slp3 <- Subset$slp^3
## set road to be 0
Subset$ROAD <- 0
## habitat reclass
lt <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
lt <- nonDuplicated(lt[!is.na(lt$Label),], Label)
Subset$hab <- reclass(Subset$HAB_NALC2, lt[,c("Label", "BCR4")])
table(Subset$HAB_NALC2, Subset$hab)
Subset$isForest <- ifelse(Subset$hab %in% c("Decid", "Conif", "Wet"), 1, 0)
## ysd
BASE_YEAR <- 2017
Subset$ysd <- ifelse(is.na(Subset$YearFire), 100, BASE_YEAR - Subset$YearFire)
Subset$ysd[Subset$ysd < 0] <- 100
Subset$ysd10 <- ifelse(Subset$ysd <= 10, 1, 0)
Subset$ysd60 <- pmax(0, 1 - Subset$ysd / 60)
Subset$ysd100 <- pmax(0, 1 - Subset$ysd / 100)
Subset$ysdcat <- cut(Subset$ysd, c(0, 10, 99, 100), include.lowest=TRUE)
levels(Subset$ysdcat) <- c("early", "mid", "old")
Subset$ysdcat <- relevel(Subset$ysdcat, "old")
Subset$ysdmid <- ifelse(Subset$ysdcat == "mid", 1, 0)
Subset$ysd2 <- Subset$ysd^2
#Subset$ysd3 <- Subset$ysd^3
Subset$ysd05 <- sqrt(Subset$ysd)

setdiff(Terms, colnames(Subset))
xn <- Subset[,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
xn <- xn[rownames(Xn),]
xy <- Subset[rownames(xn),c("POINT_X", "POINT_Y")]
xy$pred <- 0
coordinates(xy) <- c("POINT_X", "POINT_Y")
proj4string(xy) <- "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
xy <- spTransform(xy, proj4string(slp))

od <- setwd("e:/peter/AB_data_v2017/data/raw/xy/bcr/")
BCR <- readOGR(".", "BCR_Terrestrial_master_International") # rgdal
setwd(od)
BCR <- spTransform(BCR, proj4string(xy))
BCR <- gSimplify(BCR, tol=500, topologyPreserve=TRUE)


## load species data
Col <- colorRampPalette(c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506'))(10)
#spp <- "BCCH"
SPP <- gsub("\\.RData", "", list.files("e:/peter/bam/bcr4/results/"))
for (spp in SPP) {
    gc()
    cat(spp, "\n");flush.console()
    load(paste0("e:/peter/bam/bcr4/results/", spp, ".RData"))
    est <- getEst(out, stage = 4, X=Xn)
    est <- est[1:25,]

    ## making prediction based on 1st run
    pr <- exp(drop(Xn %*% est[1,]))

    ## making prediction based on all the bootstrap runs
    mat <- pbapply(est, 1, function(z) Xn %*% z)
    ## median and 90% CI
    #qpr <- exp(t(apply(mat, 1, quantile, c(0.5, 0.05, 0.95))))
    ## bootstrap mean is much quicker -- use this for exploring the patterns
    mpr <- exp(rowMeans(mat))

    ## usually truncating at 99% helps get rid of outliers
    q <- quantile(pr, 0.99)
    pr[pr > q] <- q
    q <- quantile(mpr, 0.99)
    mpr[mpr > q] <- q

    ## See 1st and averaged prediction is different
    #i <- sample(nrow(qpr), 5000)
    #plot(pr[i], qpr[i,1])

    ## save raster and corresponding map
    ## rasterize predictions
    #xy@data$pred <- mpr
    #rpr <- rasterize(xy, slp, field="pred")
    #writeRaster(rpr, paste0("e:/peter/bam/bcr4/maps/", spp, ".tif"))
    #png(paste0("e:/peter/bam/bcr4/maps/", spp, ".png"))
    #plot(rpr, main=spp)
    #dev.off()

    ## save a map without creating a raster (save time)
    q <- quantile(mpr, seq(0, 1, by=0.1))
    ii <- cut(mpr, q, include.lowest=TRUE)
    png(paste0("e:/peter/bam/bcr4/maps/", spp, "-fancy-map.png"))
    plot(xy, pch=".", col=Col[ii], main=spp)
    plot(BCR, add=TRUE, border="grey")
    legend("topright", bty="n", fill=rev(Col), border=NA,
        legend=paste0(">", round(rev(q[1:9]), 6)))
    dev.off()
}

## ALFRESCO

library(sp)
library(raster)

slp <- raster("e:/peter/bam/bcr4/alfresco_data/iem_prism_slope_1km.tiff")
veg <- raster("e:/peter/bam/bcr4/alfresco_data/iem_vegetation_model_input_v0_4.tiff")
lcc <- raster("e:/peter/bam/bcr4/alfresco_data/LandCover_iem_ALFRESCO_2005.tiff")
lcc[lcc==255] <- NA

veg1 <- raster("e:/peter/bam/bcr4/alfresco_data/Veg_171_2050.tif")
age1 <- raster("e:/peter/bam/bcr4/alfresco_data/Age_171_2050.tif")
age1[age1 < 0] <- 100


V4.0
0 | Not Modeled
1 | Black Spruce Forest
2 | White Spruce Forest
3 | Deciduous Forest
4 | Shrub Tundra
5 | Graminoid Tundra
6 | Wetland Tundra
7 | Barren lichen-moss
8 | Heath
9 | Maritime Upland Forest
10 | Maritime Forested Wetland
11 | Maritime Fen
12 | Maritime Alder Shrubland

v2.0
Final Legend:
0 - Not Modeled
1 - Black Spruce
2 - White Spruce
3 - Deciduous
4 - Shrub Tundra
5 - Graminoid  Tundra
6 - Wetland Tundra
7 - Barren lichen-moss
8 - Temperate Rainforest
255 - out of bounds

