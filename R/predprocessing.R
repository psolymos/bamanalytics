## Data processing for nationam BAM analyses: prediction

ROOT <- "e:/peter/bam/pred-2015"
ROOT2 <- "e:/peter/bam/Apr2016"

library(mefa4)
#library(pbapply)
library(data.table)
getOption("datatable.fread.datatable")
options(datatable.fread.datatable=FALSE)
getOption("datatable.fread.datatable")

################# PREPROCESSING ###############################

## Climate preprocessing -------------------------
z <- fread(file.path(ROOT, "cti", "ECGrid_climMarch2014.csv"))
#z <- read.csv(file.path(ROOT, "cti", "ECGrid_climMarch2014.csv"))
cn1 <- c("pointid", "POINT_X", "POINT_Y",
    "NORM_6190_CMDpm", "NORM_6190_CMIJJApm", "NORM_6190_CMIpm",
    "NORM_6190_DD0", "NORM_6190_DD5", "NORM_6190_MWMT", "NORM_6190_MCMT",
    "NORM_6190_EMT", "NORM_6190_TD", "NORM_6190_MSP")

write.csv(z[,cn1], file=file.path(ROOT, "cti", "clim-only.csv"))

## climate, cti, slope
clim <- fread(file.path(ROOT, "cti", "clim-only.csv"))
clim$V1 <- NULL

## climate
clim$CMIJJA <- (clim$NORM_6190_CMIJJApm - 0) / 50
clim$NORM_6190_CMIJJApm <- NULL
clim$CMI <- (clim$NORM_6190_CMIpm - 0) / 50
clim$NORM_6190_CMIpm <- NULL
clim$TD <- (clim$NORM_6190_TD - 300) / 100
clim$NORM_6190_TD <- NULL
clim$DD0 <- (clim$NORM_6190_DD0 - 1000) / 1000
clim$NORM_6190_DD0 <- NULL
clim$DD5 <- (clim$NORM_6190_DD5 - 1600) / 1000
clim$NORM_6190_DD5 <- NULL
clim$EMT <- (clim$NORM_6190_EMT + 400) / 100
clim$NORM_6190_EMT <- NULL
clim$MSP <- (clim$NORM_6190_MSP - 400) / 200
clim$NORM_6190_MSP <- NULL
#clim$DD02 <- clim$DD0^2
#clim$DD52 <- clim$DD5^2
#clim$CMI2 <- clim$CMI^2
#clim$CMIJJA2 <- clim$CMIJJA^2
clim$NORM_6190_CMDpm <- NULL
clim$NORM_6190_MWMT <- NULL
clim$NORM_6190_MCMT <- NULL

climN <- read.csv("e:/peter/bam/Apr2016/north/climlu_NorthXY_2016.csv")
climN <- climN[,c("pointid", "east", "north", "CMD",
    "CMI", "CMIJJA", "dd01", "dd51", "EMT", "mcmt", "msp", "mwmt", "td")]
climN$POINT_X <- climN$east
climN$east <- NULL
climN$POINT_Y <- climN$north
climN$north <- NULL
climN$CMIJJA <- (climN$CMIJJA - 0) / 50
climN$CMI <- (climN$CMI - 0) / 50
climN$TD <- (climN$td - 300) / 100
climN$td <- NULL
climN$DD0 <- (climN$dd01 - 1000) / 1000
climN$dd01 <- NULL
climN$DD5 <- (climN$dd51 - 1600) / 1000
climN$dd51 <- NULL
climN$EMT <- (climN$EMT + 400) / 100
climN$MSP <- (climN$msp - 400) / 200
climN$msp <- NULL
climN <- climN[,colnames(clim)]

clim <- rbind(clim, climN)

## loading CTI
pg2 <- fread(file.path(ROOT, "cti", "EC1kclip_cti90.csv"))
pg2$X <- NULL
pg2$Y <- NULL
pg2$grid_code <- NULL
gc()
## additional CTI values
pg2x <- fread(file.path(ROOT, "cti", "EC1kclip_cti90_extra_NB.csv"))
pg2x$X <- NULL
pg2x$Y <- NULL
pg2x$grid_code <- NULL
pg2 <- rbind(pg2,pg2x)
rm(pg2x)
head(pg2)

pg2$lat <- NULL
pg2$lon <- NULL

pgN <- read.csv("e:/peter/bam/Apr2016/north/NorthernXY_2016_terrain90.csv")
pg2 <- rbind(pg2,pgN[,c("pointid","cti90")])
rownames(pg2) <- pg2$pointid

pg2$CTI <- (pg2$cti90 - 8) / 4
clim$CTI <- pg2$CTI[match(clim$pointid, pg2$pointid)]
clim$CTI2 <- clim$CTI^2
rm(pg2)
gc()

## loading slope values
pg3 <- fread(file.path(ROOT, "cti", "EC1kclip_slope90.csv"))
pg3$X <- NULL
pg3$Y <- NULL
pg3$grid_code <- NULL
pg3$lon <- NULL
pg3$lat <- NULL
gc()
## additional slope values
pg3a <- fread(file.path(ROOT, "cti", "EC1kclip_slope90_extra.csv"))[,c("pointid","slope90")]
pg3b <- fread(file.path(ROOT, "cti", "EC1kclip_slope90_extra_NB.csv"))[,c("pointid","slope90")]
pg3ab <- rbind(pg3a,pg3b)
pg3 <- rbind(pg3,pg3ab)
pg3 <- rbind(pg3,pgN[,c("pointid","slope90")])


## slope, cti, elev
pg3$SLP <- sqrt(pg3$slope90)
clim$SLP <- pg3$SLP[match(clim$pointid, pg3$pointid)]
clim$SLP2 <- clim$SLP^2
rm(pg3,pg3a,pg3ab,pg3b,pgN)
gc()

save(clim, file=file.path(ROOT, "pg-clim.Rdata"))

## adding in missing Northern belt points for CTI and Slope
load(file.path(ROOT, "pg-clim.Rdata"))
sum(is.na(clim$CTI))
rownames(clim) <- clim$pointid
#ctin <- fread("e:/peter/bam/Apr2016/north2/NorthernXYgap_2016_terrain90.csv")
#compare_sets(clim$pointid, ctin$pointid)
#compare_sets(clim$pointid[is.na(clim$CTI)], ctin$pointid)

ctin <- fread("e:/peter/bam/Apr2016/north3/cracksXY2_2016_cti90.csv")
slpin <- fread("e:/peter/bam/Apr2016/north3/cracksXY2_2016_slope90.csv")
ctin$slope90 <- slpin$slope90[match(ctin$pointid, slpin$pointid)]
rm(slpin)

## Nov 18 2016 update from DS
#ctin2 <- fread("e:/peter/bam/Apr2016/north3/missingXY2_2016_cti90.csv")
#slpin <- fread("e:/peter/bam/Apr2016/north3/missingXY2_2016_slope90.csv")
ctin2 <- fread("e:/peter/bam/Apr2016/north4/MissingXY_2016_cti90.csv")
slpin <- fread("e:/peter/bam/Apr2016/north4/missingXY_2016_slope90.csv")
ctin2$slope90 <- slpin$slope90[match(ctin2$pointid, slpin$pointid)]
rm(slpin)

ct1 <- fread("e:/peter/bam/Apr2016/north2/NorthernXYgap_2016_terrain90.csv")

ct2 <- fread("e:/peter/bam/Apr2016/north2/NorthernXYgap2_2016_cti90.csv")
ct3 <- fread("e:/peter/bam/Apr2016/north2/NorthernXYgap2_2016_slope90.csv")
ct2$slope90 <- ct3$slope90[match(ct2$pointid, ct3$pointid)]
rm(ct3)

ct4 <- fread("e:/peter/bam/Apr2016/north/NorthernXY_2016_terrain90.csv")

ct5 <- fread("e:/peter/bam/Apr2016/north4/missingXY4_2016_cti90.csv")
ct6 <- fread("e:/peter/bam/Apr2016/north4/missingXY4_2016_slope90.csv")
ct5$slope90 <- ct6$slope90[match(ct5$pointid, ct6$pointid)]
rm(ct6)

if (FALSE) {
xy <-  fread("e:/peter/bam/Apr2016/north3/missing-cti.csv")
with(xy,plot(POINT_X,POINT_Y,pch=".",col=1))
with(ctin,points(x,y,pch=".",col=3))
with(ctin2,points(x,y,pch=".",col=4))
with(ct1,points(x,y,pch=".",col=2))
with(ct4,points(x,y,pch=".",col="gold"))
with(ct2,points(x,y,pch=".",col="lightblue"))
with(ct5,points(x,y,pch=".",col="pink"))
}

#compare_sets(ctin$pointid, ctin2$pointid)
ctin0 <- ctin
cn <- c("pointid", "cti90", "slope90")
ctin <- rbind(ctin0[,cn], ctin2[,cn], ct1[,cn], ct2[,cn], ct4[,cn], ct5[,cn])
dim(ctin);sum(duplicated(ctin$pointid))
ctin <- nonDuplicated(ctin, pointid, TRUE)
#rownames(clim) <- clim$pointid


#pid <- as.character(intersect(clim$pointid[is.na(clim$CTI)], ctin$pointid))
#compare_sets(clim$pointid, ctin$pointid)
pid <- as.character(intersect(clim$pointid, ctin$pointid))
length(pid)

#with(clim, plot(POINT_X, POINT_Y, pch="."))
#with(clim[is.na(clim$CTI),], points(POINT_X, POINT_Y, pch=".", col=2))
#with(clim[pid,], points(POINT_X, POINT_Y, pch=".", col=3))

sum(is.na(clim$CTI))
sum(is.na(clim$SLP))
ctin2 <- ctin[pid,]
ctin2$CTI <- (ctin2$cti90 - 8) / 4
ctin2$CTI2 <- ctin2$CTI^2
ctin2$SLP <- sqrt(ctin2$slope90)
ctin2$SLP2 <- ctin2$SLP^2
clim[pid,"CTI"] <- ctin2[pid,"CTI"]
clim[pid,"CTI2"] <- ctin2[pid,"CTI2"]
clim[pid,"SLP"] <- ctin2[pid,"SLP"]
clim[pid,"SLP2"] <- ctin2[pid,"SLP2"]
sum(is.na(clim$CTI))
sum(is.na(clim$SLP))

sum(duplicated(clim$pointid))
save(clim, file=file.path(ROOT, "pg-clim-gap.Rdata"))

## -----------------------------------

## CEC
x <- fread(file.path(ROOT2, "north2", "AllPredPointsCECFrameworkOct23_2016.csv"))
x$LEVEL3[x$LEVEL3 == ""] <- NA
x$LEVEL3[x$LEVEL3 == "Water"] <- NA
x$LEVEL3 <- as.factor(x$LEVEL3)
cec <- x[,c("pointid", "LEVEL3", "POINT_X", "POINT_Y")]
rownames(cec) <- cec$pointid
rm(x)
gc()

save(cec, file=file.path(ROOT, "cec.Rdata"))

## -------------------------------------

x <- fread(file.path(ROOT, "PredictionIntersections",
    "CovariatesPredGridOneforpeter.txt"))
rownames(x) <- x$pointid

x$BCRNAME <- NULL
x$COUNTRY <- as.factor(x$COUNTRY)
x$PROVINCE_S <- as.factor(x$PROVINCE_S)

## REG (based on BCR and Jurs)
x$JURS <- x$PROVINCE_S

x$eosdbam <- NULL
x$PROVINCE_S <- NULL
x$COUNTRY <- NULL
x$LCCV1_3Can <- NULL
#x$JURS <- NULL
#x$BCR_JURS <- NULL
#x$BCR_JURS0 <- NULL
#x$BCR <- NULL
gc()

## HAB: NALC
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
x$HAB_NALC2 <- ltnalc$Label[match(x$nalcms_05, ltnalc$Value)]
x$HAB_NALC2 <- relevel(x$HAB_NALC2, "Conif")
table(x$nalcms_05, x$HAB_NALC2,useNA="a")
x$nalcms_05 <- NULL
#levels(x$HAB_NALC2) <- c(levels(x$HAB_NALC2), "0-density")
#x$HAB_NALC2[is.na(x$HAB_NALC2)] <- "0-density"
table(x$HAB_NALC2)

## HGT, HGT2
x$HGT <- x$SimardG / 25
x$HGT2 <- x$HGT^2
x$SimardG <- NULL

## TREE fix comes here -------------------
tf <- fread(file.path(ROOT2, "tree_update",
    "TreePredNov29.csv"))
rownames(tf) <- tf$pointid
x$tree_wrong <- x$tree
x$tree <- tf$Tree[match(x$pointid, tf$pointid)]
sum(is.na(x$tree_wrong))
sum(is.na(x$tree))
#with(x[sample(nrow(x), 10^5),], plot(jitter(tree_wrong), jitter(tree),
#    xlim=c(0,100), ylim=c(0,100), col=rgb(0,0,1,0.05), pch=19, cex=1))

## TR3
#x$TREE <- x$tree
x$tree[x$tree > 100] <- NA
x$tree[x$tree < 0] <- NA
x$tree <- x$tree / 100
x$TR3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
x$TR3[x$tree < 0.25] <- "Open"
x$TR3[x$tree >= 0.25 & x$tree < 0.60] <- "Sparse"
x$TR3[x$tree >= 0.60] <- "Dense"
table(x$TR3, useNA="a")
x$tree <- NULL
x$tree_wrong <- NULL

## NALC-TREE combo
x$HAB_NALC1 <- as.character(x$HAB_NALC2)
ii <- x$HAB_NALC1 %in% c("Conif", "Decid", "Mixed", "Wet")
tmp <- as.character(interaction(x$HAB_NALC2, x$TR3, sep="", drop=TRUE))
x$HAB_NALC1[ii] <- tmp[ii]
x$HAB_NALC1 <- as.factor(x$HAB_NALC1)
x$HAB_NALC1 <- relevel(x$HAB_NALC1, "ConifDense")
table(x$HAB_NALC2, x$HAB_NALC1)

## decid + mixed
x$isDM <- ifelse(x$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
x$isNF <- ifelse(x$HAB_NALC1 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub",
    "WetOpen", "DecidOpen", "ConifOpen", "MixedOpen"), 1L, 0L)
x$isDev <- ifelse(x$HAB_NALC2 %in% c("Agr", "Devel"), 1L, 0L)
x$isOpn <- x$isNF
x$isOpn[x$isDev == 1] <- 0
x$isWet <- ifelse(x$HAB_NALC2 %in% c("Wet"), 1L, 0L)
x$isDec <- ifelse(x$HAB_NALC2 %in% c("Decid"), 1L, 0L)
x$isMix <- ifelse(x$HAB_NALC2 %in% c("Mixed"), 1L, 0L)

## LIN, POL
## linear features
x$LIN <- log(x$BEAD_linear + 1)
x$BEAD_linear <- NULL
## Polygon is fine as it is (0-1)
x$POL <- x$BEAD_poly
x$BEAD_poly <- NULL

all(rownames(x) == as.character(x$pointid))

## Brandt boreal
#br <- read.csv(file.path(ROOT, "brandt", "Pred_BrandtBoreal.csv"))
br <- fread(file.path(ROOT, "brandt", "Pred_BrandtBoreal.csv"))
br$pointid <- as.character(br$pointid)
br$pointid <- gsub(",", "", br$pointid)
br <- br[!duplicated(br$pointid),]
rownames(br) <- br$pointid
table(br$TYPE)
x$Brandt <- br$TYPE[match(rownames(x), br$pointid)]
table(x$Brandt, useNA="a")

xN <- read.csv("e:/peter/bam/Apr2016/north/MissingPredValuesJune8.csv")

## "pointid", "POINT_X", "POINT_Y", "BCR", "JURS","Brandt"
x2 <- with(xN, data.frame(
    pointid=pointid,
    POINT_X=POINT_X,
    POINT_Y=POINT_Y,
    BCR=BCR,
    JURS=NAME,
    Brandt=NA
))
rownames(x2) <- x2$pointid

br2 <- fread(file.path(ROOT2, "north2", "MissingPredValuesBoreal_Oct23_2016.csv"))
rownames(br2) <- br2$pointid
compare_sets(x2$pointid, br2$pointid)
x2$Brandt <- factor(br2$TYPE[match(x2$pointid, br2$pointid)],
    c("B_ALPINE","BOREAL","H_ALPINE","HEMIBOREAL"))

## "HAB_NALC2", "HGT", "HGT2", "TR3", "HAB_NALC1", "isDM", "isNF",
## "isDev", "isOpn", "isWet", "isDec", "isMix", "LIN", "POL"

## HAB: NALC
x2$HAB_NALC2 <- ltnalc$Label[match(xN$NALCMS_05, ltnalc$Value)]
x2$HAB_NALC2 <- relevel(x2$HAB_NALC2, "Conif")
table(xN$NALCMS_05, x2$HAB_NALC2,useNA="a")
table(x2$HAB_NALC2)

## HGT, HGT2
x2$HGT <- xN$Simard / 25
x2$HGT2 <- x2$HGT^2

## TREE fix comes here -------------------
xN$tree_wrong <- xN$Tree
xN$Tree <- tf$Tree[match(xN$pointid, tf$pointid)]
sum(is.na(xN$tree_wrong))
sum(is.na(xN$Tree))

## TR3
xN$Tree[xN$Tree > 100] <- NA
xN$Tree[xN$Tree < 0] <- NA
xN$Tree <- xN$Tree / 100
x2$TR3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
x2$TR3[xN$Tree < 0.25] <- "Open"
x2$TR3[xN$Tree >= 0.25 & xN$Tree < 0.60] <- "Sparse"
x2$TR3[xN$Tree >= 0.60] <- "Dense"
table(x2$TR3, useNA="a")

## NALC-TREE combo
x2$HAB_NALC1 <- as.character(x2$HAB_NALC2)
ii <- x2$HAB_NALC1 %in% c("Conif", "Decid", "Mixed", "Wet")
tmp <- as.character(interaction(x2$HAB_NALC2, x2$TR3, sep="", drop=TRUE))
x2$HAB_NALC1[ii] <- tmp[ii]
x2$HAB_NALC1 <- as.factor(x2$HAB_NALC1)
x2$HAB_NALC1 <- relevel(x2$HAB_NALC1, "ConifDense")
table(x2$HAB_NALC2, x2$HAB_NALC1)

## decid + mixed
x2$isDM <- ifelse(x2$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
x2$isNF <- ifelse(x2$HAB_NALC1 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub",
    "WetOpen", "DecidOpen", "ConifOpen", "MixedOpen"), 1L, 0L)
x2$isDev <- ifelse(x2$HAB_NALC2 %in% c("Agr", "Devel"), 1L, 0L)
x2$isOpn <- x2$isNF
x2$isOpn[x2$isDev == 1] <- 0
x2$isWet <- ifelse(x2$HAB_NALC2 %in% c("Wet"), 1L, 0L)
x2$isDec <- ifelse(x2$HAB_NALC2 %in% c("Decid"), 1L, 0L)
x2$isMix <- ifelse(x2$HAB_NALC2 %in% c("Mixed"), 1L, 0L)

## LIN, POL
## linear features
x2$LIN <- log(xN$BEAD_L + 1)
## Polygon is fine as it is (0-1)
x2$POL <- xN$BEAD_P
x2$POL[x2$POL < 0] <- NA

setdiff(colnames(x), colnames(x2))
#compare_sets(x$pointid, x2$pointid)

x <- rbind(x, x2[,colnames(x)])
rownames(x) <- x$pointid
gc()

## new tree cover -- same as the old ...
#xtree <- fread(file.path(ROOT2, "north2", "AllPredPointsTreecover_oct23.csv"))
#compare_sets(x$pointid, xtree$pointid)
#x$tree_new <- xtree$treecov[match(x$pointid, xtree$pointid)]
#summary(x$tree_new - x$tree)


load(file.path(ROOT, "cec.Rdata"))
compare_sets(x$pointid, cec$pointid)
x$LEVEL3 <- cec$LEVEL3[match(x$pointid, cec$pointid)]
#x$LEVEL3 <- as.factor(x$LEVEL3)
#x$LEVEL3[!is.na(x$LEVEL3) & x$LEVEL3 == ""] <- NA
#x$LEVEL3[!is.na(x$LEVEL3) & x$LEVEL3 == "Water"] <- NA
#x$LEVEL3 <- droplevels(x$LEVEL3)
table(x$LEVEL3,useNA="a")

x$BCR[x$BCR < 1 | x$BCR >= 100] <- NA
levels(x$JURS) <- toupper(levels(x$JURS))
x$JURS[x$JURS == ""] <- NA
x$JURS <- droplevels(x$JURS)

table(x$Brandt, useNA="a")

save(x, file=file.path(ROOT, "pg-main-NALConly.Rdata"))
rm(cec)
gc()

#i <- "Brandt"
ii <- c("BCR", "JURS", "HAB_NALC2",
    "HGT", "TR3", "HAB_NALC1", "LIN", "POL", "Brandt", "LEVEL3")
for (i in ii) {
gc()
cat(i, "\n");flush.console()
png(paste0("e:/peter/bam/Apr2016/maps/", i, ".png"), width=2000, height=1500)
op <- par(mar=c(1,1,1,1))
plot(x$POINT_X, x$POINT_Y,pch=".",col="lightgrey", axes=FALSE, ann=FALSE, main=i)
#plot(x$POINT_X, x$POINT_Y,pch=".",col=1, axes=FALSE, ann=FALSE, main=i)
#points(x2$POINT_X, x2$POINT_Y,pch=".",col="darkgrey")
with(x[!is.na(x[,i]),], points(POINT_X, POINT_Y,pch=".",col="tomato"))
#with(x[x$TREE > 100,], points(POINT_X, POINT_Y,pch=".",col=2))
par(op)
dev.off()
}
png(paste0("e:/peter/bam/Apr2016/maps/LEVEL3-allVals.png"), width=2000, height=1500)
op <- par(mar=c(1,1,1,1))
plot(x$POINT_X, x$POINT_Y,pch=".", axes=FALSE, ann=FALSE, main="LEVEL3",
    col=as.integer(x$LEVEL3))
par(op)
dev.off()

st <- read.csv(file.path(ROOT2, "BAMCECStudyAreaEcoregionLevel2.csv"))
levs <- levels(st$LEVEL3)
iii <- x$LEVEL3 %in% levs
png(paste0("e:/peter/bam/Apr2016/maps/LEVEL3-usedVals.png"), width=2000, height=1500)
op <- par(mar=c(1,1,1,1))
plot(x$POINT_X, x$POINT_Y,pch=".", axes=FALSE, ann=FALSE, main="LEVEL3", col="lightgrey")
with(x[iii,], points(POINT_X, POINT_Y,pch=".",col=LEVEL3))
par(op)
dev.off()

st <- read.csv(file.path(ROOT2, "BAMCECStudyAreaEcoregionLevel2.csv"))
levs <- levels(st$LEVEL3)
iii <- x$LEVEL3 %in% levs
png(paste0("e:/peter/bam/Apr2016/maps/LEVEL3-usedVals2.png"), width=2000, height=1500)
op <- par(mar=c(1,1,1,1))
plot(x$POINT_X, x$POINT_Y,pch=".", axes=FALSE, ann=FALSE, main="LEVEL3", col="lightgrey")
with(x[iii,], points(POINT_X, POINT_Y,pch=".",col="tomato"))
par(op)
dev.off()


## Fire, GWF

z <- fread(file.path(ROOT, "fire", "ecgrid_fire.csv"))
z$POINT_X <- NULL
z$POINT_Y <- NULL
z$SIZE_HA <- NULL

zz <- fread(file.path(ROOT, "gfw", "GFW_Predictiongrid.csv"))
zz$OBJECTID <- NULL
zz$FID_Pred1kmIdentGrid4k <- NULL
zz$FID_Ec_1k_clip <- NULL
zz$grid_code <- NULL
zz$POINT_X <- NULL
zz$POINT_Y <- NULL
zz$FID_Grid4x4GFWall <- NULL
zz$Id <- NULL
zz$gridcode <- NULL
zz$Grid4x4 <- NULL
zz$FID_MergeLossLCC <- NULL
zz$ID_1 <- NULL
zz$pointid <- gsub(",", "", zz$pointid)
zz$pointid <- as.integer(zz$pointid)
zz$YearLoss <- as.integer(zz$GRIDCODE_1 + 2000)
zz$YearLoss[zz$GRIDCODE_1 < 1] <- NA

length(intersect(z$pointid, zz$pointid))
colnames(z)[colnames(z) == "YEAR"] <- "YearFire"
z$YearLoss <- zz$YearLoss[match(z$pointid, zz$pointid)]

## N here
lN <- read.csv("e:/peter/bam/Apr2016/north/MissingPred_GFWC_Loss.csv")
lN$YearLoss <- as.integer(lN$GFWC_noFire + 2000)
lN$YearLoss[lN$GFWC_noFire < 1] <- NA

fN <- fread(file.path(ROOT2, "north2", "fire_NorthXY_2016.csv"))
fN$YEAR[fN$YEAR==0] <- NA
lN$YearFire <- fN$YEAR[match(lN$pointid,fN$pointid)]

loss <- rbind(z[,c("pointid", "YearFire", "YearLoss")],
    lN[,c("pointid", "YearFire", "YearLoss")])
rownames(loss) <- loss$pointid
ii <- !is.na(loss$YearFire) & loss$YearFire > 9000
loss$YearFire[ii] <- loss$YearFire[ii] - 8000
save(loss, file=file.path(ROOT, "pg-loss.Rdata"))

################# PACKAGING ###############################


## save Level3 chunks

ROOT <- "e:/peter/bam/pred-2015"
ROOT2 <- "e:/peter/bam/Apr2016"
library(mefa4)

load(file.path(ROOT, "pg-main-NALConly.Rdata"))
#x <- x[x$EOSD_COVER == 1,]
rownames(x) <- x$pointid

load(file.path(ROOT, "pg-loss.Rdata"))
ii <- loss$YearFire >= 9000 & !is.na(loss$YearFire)
loss$YearFire[ii] <- loss$YearFire[ii] - 8000
x$YearFire <- loss$YearFire[match(x$pointid, loss$pointid)]
x$YearLoss <- loss$YearLoss[match(x$pointid, loss$pointid)]
rm(loss)

#load(file.path(ROOT, "pg-clim.Rdata"))
load(file.path(ROOT, "pg-clim-gap.Rdata"))
#with(clim[!is.na(clim$CTI),],plot(POINT_X,POINT_Y,pch=".",col=1))
#with(clim[is.na(clim$CTI),],points(POINT_X,POINT_Y,pch=".",col=2))
sum(duplicated(clim$pointid))
rownames(clim) <- clim$pointid
sum(duplicated(x$pointid))
clim <- clim[match(x$pointid, clim$pointid),]
x <- data.frame(x, clim)
rm(clim)
gc()
#cti_missing <- x[is.na(x$CTI),c("pointid","POINT_X","POINT_Y")]
#write.csv(cti_missing, row.names=FALSE, file="e:/peter/bam/Apr2016/cti-missing.csv")

if (FALSE) {

    #tmp <- is.na(x$CTI) & !is.na(clim$CTI)
    table(clim=is.na(x$CTI), clim_gap=!is.na(clim$CTI))
    aa <- table(x$LEVEL3, gap=is.na(x$CTI) & !is.na(clim$CTI))
    gap3 <- rownames(aa[aa[,2]>0,])
    ## gap3
    c("2.1.5", "2.1.9", "2.2.1", "2.2.2", "2.3.1", "2.4.1", "2.4.2",
    "2.4.4", "3.1.1", "3.1.3", "3.2.1", "3.2.3", "3.3.1", "3.4.5", "6.1.1")
    ## gap3 %in% levs
    c("2.1.9", "2.2.1", "2.2.2", "2.3.1", "2.4.1", "2.4.2", "2.4.4",
    "3.1.1", "3.1.3", "3.2.1", "3.2.3", "3.3.1", "3.4.5", "6.1.1")


    ii <- colnames(x)
    jj <- sample.int(nrow(x), 10^6)
    xtmp <- x[jj,]

    for (i in ii) {
    gc()
    cat(i, "\n");flush.console()
    png(paste0("e:/peter/bam/Apr2016/maps-sub/", i, ".png"), width=2000, height=1500)
    op <- par(mar=c(1,1,1,1))
    plot(xtmp$POINT_X, xtmp$POINT_Y,pch=".",col="lightgrey", axes=FALSE, ann=FALSE, main=i)
    with(xtmp[!is.na(xtmp[,i]),], points(POINT_X, POINT_Y,pch=".",col="tomato"))
    par(op)
    dev.off()
    }

nn <- c("BCR", "JURS",
    "HAB_NALC2", "HGT", "HGT2", "TR3", "HAB_NALC1",
    "isDM", "isNF", "isDev", "isOpn",
    "isWet", "isDec", "isMix",
    "CMIJJA", "CMI", "TD", "DD0", "DD5", "EMT", "MSP", "CTI", "CTI2",
    "SLP", "SLP2", "DD02", "DD52", "CMI2", "CMIJJA2")
IS_OK <- rowSums(is.na(x[,nn])) == 0
table(IS_NA)
(aa <- data.frame(colSums(is.na(x[,nn]))))

png(paste0("e:/peter/bam/Apr2016/maps-sub/all_NA.png"), width=2000, height=1500)
op <- par(mar=c(1,1,1,1))
with(x[IS_OK,], plot(POINT_X, POINT_Y, pch=".",col="grey", axes=FALSE, ann=FALSE))
with(x[!IS_OK,], points(POINT_X, POINT_Y,pch=".",col="tomato"))
par(op)
dev.off()

nnn <- c("pointid", "POINT_X", "POINT_Y", "JURS", "HAB_NALC2", "HAB_NALC1", "CMIJJA", "CTI")
x <- x[,nnn]
#i <- "CTI"
for (i in colnames(x)[-(1:3)]) {
    gc()
    cat(i, "\n");flush.console()
    IS_OK <- !is.na(x[,i])
    png(paste0("e:/peter/bam/Apr2016/maps-sub/fix", i, ".png"), width=2000, height=1500)
    op <- par(mar=c(1,1,1,1))
    #plot(x$POINT_X, x$POINT_Y,pch=".",col="lightgrey", axes=FALSE, ann=FALSE, main=i)
    #with(x[!is.na(x[,i]),], points(POINT_X, POINT_Y,pch=".",col="tomato"))
    with(x[IS_OK,], plot(POINT_X, POINT_Y, pch=".",col="grey", axes=FALSE, ann=FALSE))
    with(x[!IS_OK,], points(POINT_X, POINT_Y,pch=".",col="tomato"))

#    with(x, plot(POINT_X, POINT_Y, pch=".",col="tomato", axes=FALSE, ann=FALSE))
#    with(x[IS_OK,], points(POINT_X, POINT_Y,pch=".",col="grey"))
    par(op)
    dev.off()
}

}

st <- read.csv(file.path(ROOT2, "BAMCECStudyAreaEcoregionLevel2.csv"))
levs <- levels(st$LEVEL3)
iii <- x$LEVEL3 %in% levs
table(iii)
x <- x[iii,]
gc()

#x <- x[!is.na(x$LEVEL3),]

x$DD02 <- x$DD0^2
x$DD52 <- x$DD5^2
x$CMI2 <- x$CMI^2
x$CMIJJA2 <- x$CMIJJA^2
x$TR3[is.na(x$TR3)] <- "Open" # this is global
x$pointid.1 <- NULL
x$POINT_X.1 <- NULL
x$POINT_Y.1 <- NULL

aa <- data.frame(colSums(is.na(x)))

reg <- levels(droplevels(x$LEVEL3))
for (i in reg) {
    cat(i, "\n");flush.console()
    ii <- x$LEVEL3 == i
    dat <- x[ii,]
    save(dat, file=file.path("e:/peter/bam/pred-2016", "chunks", paste0("pgdat-", i, ".Rdata")))
    gc()
}

## fill-in NAs
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
reg <- levels(droplevels(st$LEVEL3))
for (i in reg) {
    load(file.path("e:/peter/bam/pred-2016", "chunks", paste0("pgdat-", i, ".Rdata")))
    cat(i, sum(is.na(dat)));flush.console()
    dat$LIN <- NULL
    dat$POL <- NULL
    datCheck <- dat
    datCheck$Brandt <- NULL
    datCheck$YearFire <- NULL
    datCheck$YearLoss <- NULL
    jj <- which(rowSums(is.na(datCheck)) > 0)
    for (j in jj) {
        d <- sqrt((dat$POINT_X - dat$POINT_X[j])^2 + (dat$POINT_Y - dat$POINT_Y[j])^2)
        d[jj] <- Inf
        jfill <- which.min(d)
        wc <- is.na(dat[j,,drop=TRUE])
        wc[c("Brandt","YearFire","YearLoss")] <- FALSE
        dat[j, wc] <- dat[jfill, wc]
    }
    cat(" -->", sum(is.na(dat)), "\n");flush.console()
    save(dat, file=file.path("e:/peter/bam/pred-2016", "chunks", paste0("pgdat-", i, ".Rdata")))
    gc()
}

## check NAs
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
reg <- levels(droplevels(st$LEVEL3))
for (i in reg) {
    load(file.path("e:/peter/bam/pred-2016", "chunks", paste0("pgdat-", i, ".Rdata")))
    datCheck <- dat
    datCheck$Brandt <- NULL
    datCheck$YearFire <- NULL
    datCheck$YearLoss <- NULL
    jj <- which(rowSums(is.na(datCheck)) > 0)
    cat(i, " -->", sum(is.na(datCheck)), "\n");flush.console()
}


################# POSTPROCESSING ###############################


## placeholders: HSH, HSH2, isDM, isNF
x$HAB <- 0
x$HABTR <- 0
x$isDM <- 0
x$isDec <- 0
x$isMix <- 0
x$isWet <- 0
x$isNF <- 0
x$isDev <- 0
x$isOpn <- 0

BASE_YEAR <- 2012

## YR
x$YR <- BASE_YEAR - 1997

## disturbance
SS$YearFire[is.na(SS$YearFire)] <- BASE_YEAR - 200
SS$YearLoss[is.na(SS$YearLoss)] <- BASE_YEAR - 200

## years since fire
DAT$YSF <- DAT$YEAR - DAT$YearFire
DAT$YSF[DAT$YSF < 0] <- 100
## years since loss
DAT$YSL <- DAT$YEAR - DAT$YearLoss
DAT$YSL[DAT$YSL < 0] <- 100
## years since most recent burn or loss
DAT$YSD <- pmin(DAT$YSF, DAT$YSL)

DAT$BRN <- ifelse(DAT$YSF <= 10, 1L, 0L)
DAT$LSS <- ifelse(DAT$YSL <= 10, 1L, 0L)
DAT$LSS[DAT$YEAR < 2000] <- NA
DAT$DTB <- ifelse(DAT$YSD <= 10, 1L, 0L)
DAT$DTB[DAT$YEAR < 2000] <- NA

## check NA occurrences, especially in HAB_*
plot(x[,2:3],pch=".",col=ifelse(is.na(x$HAB_EOSD2),2,1))

plot(x[,2:3],pch=".",col=ifelse(is.na(x$EMT),2,1)) # coastal
plot(x[,2:3],pch=".",col=ifelse(is.na(x$SLP),2,1)) # coastal + N YK


## fill in unknown BCRs with closest / do this by CEC level3
if (FALSE) {
x$BCR[x$BCR == 100] <- 0
x$BCR[is.na(x$BCR)] <- 0
x$BCRo <- x$BCR
table(x$BCR==0, x$JURS=="")

ii <- which(x$BCR == 0)
Prc <- 0
for (j in 1:length(ii)) {
    i <- ii[j]
    d <- sqrt((x$POINT_X - x$POINT_X[i])^2 + (x$POINT_Y - x$POINT_Y[i])^2)
    d[ii] <- Inf # keep ALL unknowns out of possible options
    x$BCR[i] <- x$BCR[which.min(d)]
    x$JURS[i] <- x$JURS[which.min(d)]
    Prc2 <- round(100*j/length(ii), 1)
    if (Prc2 > Prc) {
        Prc <- Prc2
        cat(Prc, "%\n")
        flush.console()
    }
}

## subset to include study region only
# c(2:14, 23)
save(x, file=file.path(ROOT, "pg-main-raw-NALConly.Rdata"))
}
