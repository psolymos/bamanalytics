## Data processing for nationam BAM analyses: prediction

ROOT <- "e:/peter/bam/pred-2015"
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

## NALC 4x4
if (FALSE) { # no surrounding stuff
cn2 <- c("pointid", "POINT_X", "POINT_Y",
    "lc1", "lc10", "lc11", "lc12", "lc13", "lc14", "lc15", "lc16", "lc17",
    "lc18", "lc19", "lc2", "lc5", "lc6", "lc8")
z <- z[,cn2]
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
rownames(ltnalc) <- paste0("lc", ltnalc$Value)
ltnalc <- ltnalc[c("lc1", "lc10", "lc11", "lc12", "lc13", "lc14", "lc15", "lc16", "lc17",
    "lc18", "lc19", "lc2", "lc5", "lc6", "lc8"),]
ltnalc$Label <- as.character(ltnalc$Label)
ltnalc$Label[is.na(ltnalc$Label)] <- "NONE"
m <- as.matrix(z[,c("lc1", "lc10", "lc11", "lc12", "lc13", "lc14", "lc15", "lc16", "lc17",
    "lc18", "lc19", "lc2", "lc5", "lc6", "lc8")])
m <- groupSums(m, 2, ltnalc$Label)
rownames(m) <- z$pointid
write.csv(m, file=file.path(ROOT, "cti", "nalc4x4-processed.csv"))
}

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

## -----------------------------------

x <- fread(file.path(ROOT, "PredictionIntersections",
    "CovariatesPredGridOneforpeter.txt"))
rownames(x) <- x$pointid

x$BCRNAME <- NULL
x$COUNTRY <- as.factor(x$COUNTRY)
x$PROVINCE_S <- as.factor(x$PROVINCE_S)

## REG (based on BCR and Jurs)
x$JURS <- x$PROVINCE_S
#x$JURS <- droplevels(x$JURS)

#levels(x$JURS) <- c("", "AK", "AB", "BC", "CA",
#    "CT", "DE", "ID", "IL", "IN", "IA",
#    "ME", "MB", "MD", "MA", "MI",
#    "MN", "MO", "MT", "NE", "NV", "NB",
#    "NH", "NJ", "NY", "NL", "ND",
#    "NT", "NS", "NU", "OH", "ON",
#    "OR", "PA", "PE", "QC", "RI",
#    "SK", "SD", "UT", "VT", "VA",
#    "WA", "WV", "WI", "WY", "YK")

## fill in unknown BCRs with closest
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

#x <- droplevels(x[x$BCR %in% c(2:14, 23),])

#x$BCR_JURS <- interaction(x$BCR, x$JURS, drop=TRUE, sep="_")
#levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("0_",   "100_"  )] <- "UNKNOWN"

## EOSC/LCC/NALC skeleton: gridcode and eosd extent definition !!!

#x$LCCV1_3Can <- NULL
x$eosdbam <- NULL

#xx <- x[sample.int(nrow(tmp), 10000),c("POINT_X","POINT_Y","eosdbam")]
#plot(xx[,1:2], col=ifelse(x[,3]==0, 1, 2), pch=".")
#x$EOSD_COVER <- ifelse(x$eosdbam > 0 & x$COUNTRY == "CANADA" & x$REG != "South", 1L, 0L)
#x <- x[x$eosdbam > 0 & x$COUNTRY == "CANADA" & x$REG != "South",]
#x$BCR <- NULL
x$PROVINCE_S <- NULL
x$COUNTRY <- NULL
#x$JURS <- NULL
#x$BCR_JURS <- NULL
#x$BCR_JURS0 <- NULL
#x$BCR <- NULL
gc()

if (FALSE) { # no EOSD used
## HAB: EOSD
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
x$eosdbam[x$eosdbam < 1] <- NA
x$HAB_EOSD2 <- lteosd$Reclass_label2[match(x$eosdbam, lteosd$Value)]
x$HAB_EOSD2 <- relevel(x$HAB_EOSD2, "Conif")
table(x$eosdbam, x$HAB_EOSD2,useNA="a")
x$eosdbam <- NULL
#levels(x$HAB_EOSD2) <- c(levels(x$HAB_EOSD2), "0-density")
#x$HAB_EOSD2[is.na(x$HAB_EOSD2)] <- "0-density"
table(x$HAB_EOSD2)
}

## HAB: NALC
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
x$HAB_NALC2 <- ltnalc$Label[match(x$nalcms_05, ltnalc$Value)]
x$HAB_NALC2 <- relevel(x$HAB_NALC2, "Conif")
table(x$nalcms_05, x$HAB_NALC2,useNA="a")
x$nalcms_05 <- NULL
#levels(x$HAB_NALC2) <- c(levels(x$HAB_NALC2), "0-density")
#x$HAB_NALC2[is.na(x$HAB_NALC2)] <- "0-density"
table(x$HAB_NALC2)

if (FALSE) { # no LCC used
## HAB: LCC
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
x$LCCV1_3Can[x$LCCV1_3Can < 1 | x$LCCV1_3Can > 39] <- NA
x$LCCV1_3Can[x$COUNTRY == "USA"] <- NA
x$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(x$LCCV1_3Can, ltlcc$lcc05v1_2)]
x$HAB_LCC2 <- relevel(x$HAB_LCC2, "Conif")
table(x$LCCV1_3Can, x$HAB_LCC2,useNA="a")
x$LCCV1_3Can <- NULL
#levels(x$HAB_LCC2) <- c(levels(x$HAB_LCC2), "0-density")
#x$HAB_LCC2[is.na(x$HAB_LCC2)] <- "0-density"
table(x$HAB_LCC2)
}

## isDM, isNF
## decid + mixed
if (FALSE) {
x$isDM_LCC <- ifelse(x$HAB_LCC2 %in% c("Decid", "Mixed"), 1L, 0L)
x$isDM_EOSD <- ifelse(x$HAB_EOSD2 %in% c("Decid", "Mixed"), 1L, 0L)
x$isDM_NALC <- ifelse(x$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
x$isNF_LCC <- ifelse(x$HAB_LCC2 %in%
    c("Agr", "Barren", "Burn", "Devel", "Grass", "Wet"), 1L, 0L)
x$isNF_EOSD <- ifelse(x$HAB_EOSD2 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet"), 1L, 0L)
x$isNF_NALC <- ifelse(x$HAB_NALC2 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet"), 1L, 0L)
}

## HGT, HGT2
x$HGT <- x$SimardG / 25
x$HGT2 <- x$HGT^2
x$SimardG <- NULL

## TR3
x$tree[x$tree > 100] <- NA
x$tree <- x$tree / 100
x$TR3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
x$TR3[x$tree < 0.25] <- "Open"
x$TR3[x$tree >= 0.25 & x$tree < 0.60] <- "Sparse"
x$TR3[x$tree >= 0.60] <- "Dense"
table(x$TR3, useNA="a")
x$tree <- NULL

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
br <- read.csv(file.path(ROOT, "brandt", "Pred_BrandtBoreal.csv"))
levels(br$pointid) <- gsub(",", "", levels(br$pointid))
br <- br[!duplicated(br$pointid),]
rownames(br) <- br$pointid
table(br$TYPE)
x$Brandt <- br$TYPE[match(rownames(x), br$pointid)]

x$LEVEL3 <- NA
xN <- read.csv("e:/peter/bam/Apr2016/north/MissingPredValuesJune8.csv")
xN$LCCV1_3Can <- NA

c("pointid", "POINT_X", "POINT_Y", "BCR", "LCCV1_3Can", "JURS",
"HAB_NALC2", "HGT", "HGT2", "TR3", "HAB_NALC1", "isDM", "isNF",
"isDev", "isOpn", "isWet", "isDec", "isMix", "LIN", "POL", "Brandt", "LEVEL3")

xN2 <- xN[,c()

#x <- droplevels(x[x$BCR_JURS != "UNKNOWN",])
save(x, file=file.path(ROOT, "pg-main-NALConly.Rdata"))

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

loss <- z
save(loss, file=file.path(ROOT, "pg-loss.Rdata"))

################# PACKAGING ###############################

## save bcr/jurs chunks

ROOT <- "e:/peter/bam/pred-2015"
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

load(file.path(ROOT, "pg-clim.Rdata"))
rownames(clim) <- clim$pointid
clim <- clim[match(x$pointid, clim$pointid),4:14]
x <- data.frame(x, clim)
rm(clim)

#ii <- !is.na(x$CTI) & ! is.na(x$TD)
#x <- x[ii,]

x$DD02 <- x$DD0^2
x$DD52 <- x$DD5^2

if (FALSE) {
load(file.path(ROOT, "pg-4x4.Rdata"))
tmp <- setdiff(rownames(x), rownames(lcc4x4))
plot(x[,2:3],pch=".",col=ifelse(rownames(x) %in% tmp,2,1))
x <- x[rownames(x) %in% rownames(lcc4x4),]

length(setdiff(rownames(x), rownames(lcc4x4)))
length(setdiff(rownames(x), rownames(eosd4x4)))
length(setdiff(rownames(x), rownames(nalc4x4)))

x$REG <- droplevels(x$REG)
x$BCR_JURS0 <- droplevels(x$BCR_JURS0)
lcc4x4 <- lcc4x4[rownames(x),]
eosd4x4 <- eosd4x4[rownames(x),]
nalc4x4 <- nalc4x4[rownames(x),]
}

x$TR3[is.na(x$TR3)] <- "Open" # this is global

#XYeosd <- as.matrix(x[,2:3])
#rownames(XYeosd) <- x[,1]
#save(XYeosd, file=file.path(ROOT, "XYeosd.Rdata"))

reg <- levels(droplevels(x$BCR_JURS))
for (i in reg) {
    cat(i, "\n");flush.console()
    ii <- x$BCR_JURS == i
    dat <- x[ii,]
    save(dat, file=file.path(ROOT, "chunks3", paste0("pgdat-", i, ".Rdata")))
    gc()
}

## packaging: NALC full extent

## use eosd coverage to save bcr/jurs0 chunks

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)

load(file.path(ROOT, "pg-main.Rdata"))
#x <- x[x$EOSD_COVER == 1,]
rownames(x) <- x$pointid

load(file.path(ROOT, "pg-loss.Rdata"))
ii <- loss$YearFire >= 9000 & !is.na(loss$YearFire)
loss$YearFire[ii] <- loss$YearFire[ii] - 8000
x$YearFire <- loss$YearFire[match(x$pointid, loss$pointid)]
x$YearLoss <- loss$YearLoss[match(x$pointid, loss$pointid)]
rm(loss)

load(file.path(ROOT, "pg-clim.Rdata"))
rownames(clim) <- clim$pointid
clim <- clim[match(x$pointid, clim$pointid),4:14]
x <- data.frame(x, clim)
rm(clim)
gc()

ii <- !is.na(x$CTI) & ! is.na(x$TD)
x <- x[ii,]

x$DD02 <- x$DD0^2
x$DD52 <- x$DD5^2

x$REG <- droplevels(x$REG)
x$BCR_JURS0 <- droplevels(x$BCR_JURS0)

x$TR3[is.na(x$TR3)] <- "Open" # this is global

x$pointid <- NULL

XYfull <- as.matrix(x[,c("POINT_X","POINT_Y")])
rownames(XYfull) <- rownames(x)

br <- read.csv(file.path(ROOT2, "brandt", "Pred_BrandtBoreal.csv"))
levels(br$pointid) <- gsub(",", "", levels(br$pointid))
br <- br[!duplicated(br$pointid),]
rownames(br) <- br$pointid
table(br$TYPE)
Brandt <- br$TYPE
names(Brandt) <- rownames(br)

save(XYfull, Brandt, file=file.path(ROOT, "XYfull.Rdata"))

reg <- levels(x$BCR_JURS0)
for (i in reg) {
    gc()
    cat(i, "\n");flush.console()
    ii <- x$BCR_JURS0 == i
    dat <- x[ii,]
    save(dat, file=file.path(ROOT, "chunks2", paste0("pgdat-full-", i, ".Rdata")))
}
dat <- x
save(dat, file=file.path(ROOT, paste0("pgdat-full.Rdata")))


################# POSTPROCESSING ###############################


## placeholders: HSH, HSH2, isDM, isNF
x$HSH <- 0
x$HSH2 <- 0
x$HAB <- 0
x$isDM <- 0
x$isNF <- 0

BASE_YEAR <- 2015

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


