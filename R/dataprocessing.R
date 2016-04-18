##---
##title: "Data processing for nationam BAM analyses"
##author: "Peter Solymos"
##date: "Apr 18, 2016"
##output: 
##  pdf_document: 
##    toc: true
##    toc_depth: 2
##---

### Preliminaries

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/Apr2016"

## Load required packages
library(mefa4)
library(RODBC)
library(maptools)
library(QPAD)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

### Pulling in tables

## Define MS Access database connection
con <- odbcConnectAccess2007(file.path(ROOT, "BAM_BayneAccess_BAMBBScore.accdb"))
con2 <- odbcConnectAccess2007(file.path(ROOT2, "BBS_V4draft_2016_aprilForErin.accdb"))

#### Species lookup table

TAX <- sqlFetch(con, "dbo_DD_BAM_AOU52_SPECIES_TBL")
TAX$SSMA_TimeStamp <- NULL

#### SS level stuff for BAM+BBS combined

## This table has time zone, BCR, jurisdiction, XY
SS01 <- sqlFetch(con, "dbo_BBSBAM_V4_XYTBL_ATTRIBUTES1")
SS01 <- nonDuplicated(SS01, SS, TRUE)
SS01$COUNTRY <- ifelse(SS01$JURSALPHA %in% c("AB","BC","MB","NB",
    "NL","NS","NT","NU","ON","PE","QC","SK","YT"), "CAN", "USA")

## Tree proportions
SS02 <- sqlFetch(con, "dbo_TREE_BBSBAM_V4_tbl")
SS02 <- nonDuplicated(SS02, SS, TRUE)
SS02 <- SS02[rownames(SS01),]
SS02$TREE[SS02$TREE > 100] <- NA
SS02$TREE <- SS02$TREE / 100
SS02$TREE3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
SS02$TREE3[SS02$TREE < 0.25] <- "Open"
SS02$TREE3[SS02$TREE >= 0.25 & SS02$TREE < 0.60] <- "Sparse"
SS02$TREE3[SS02$TREE >= 0.60] <- "Dense"

## Point level land cover
SS03 <- sqlFetch(con, "dbo_BAMBBS_LANDCOVER_PTS")
SS03 <- nonDuplicated(SS03, SS, TRUE)
SS03 <- SS03[rownames(SS01),]

## comment out all LCC05 and EOSD
if (FALSE) {
## Reclass LCC05
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
SS03$LCC05_PT[SS03$LCC05_PT < 1 | SS03$LCC05_PT > 39] <- NA
SS03$LCC05_PT[SS01$COUNTRY == "USA"] <- NA
SS03$HAB_LCC1 <- ltlcc$BAMLCC05V2_label1[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC1 <- relevel(SS03$HAB_LCC1, "ConifDense")
SS03$HAB_LCC2 <- relevel(SS03$HAB_LCC2, "Conif")

## Reclass EOSD
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
levels(SS03$EOSD_PT) <- sub(",", "", levels(SS03$EOSD_PT))
SS03$EOSD_PT <- as.integer(as.character(SS03$EOSD_PT))
SS03$EOSD_PT[SS03$EOSD_PT < 1] <- NA
SS03$HAB_EOSD1 <- lteosd$Reclass_label1[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD2 <- lteosd$Reclass_label2[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD1 <- relevel(SS03$HAB_EOSD1, "ConifDense")
SS03$HAB_EOSD2 <- relevel(SS03$HAB_EOSD2, "Conif")
}

## Reclass NALCMS
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
SS03$HAB_NALC2 <- ltnalc$Label[match(SS03$NALCMS_PT, ltnalc$Value)]
tmp <- as.character(interaction(SS03$HAB_NALC2, SS02$TREE3, sep="", drop=TRUE))
SS03$HAB_NALC1 <- as.character(SS03$HAB_NALC2)
ii <- SS03$HAB_NALC1 %in% c("Conif", "Decid", "Mixed")
SS03$HAB_NALC1[ii] <- tmp[ii]
SS03$HAB_NALC1 <- as.factor(SS03$HAB_NALC1)
SS03$HAB_NALC1 <- relevel(SS03$HAB_NALC1, "ConifDense")
SS03$HAB_NALC2 <- relevel(SS03$HAB_NALC2, "Conif")

## NALC is used for QPADv3
if (FALSE) {
## LCC for offsets
SS03$LCC_OFF1 <- as.factor(ltlcc$qpad_num[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)])
SS03$LCC_OFF2 <- factor(5, 1:5)
SS03$LCC_OFF2[is.na(SS03$HAB_NALC1)] <- NA
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("DecidSparse")] <- "4"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("ConifSparse","MixedSparse")] <- "3"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("DecidDense")] <- "2"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("ConifDense","MixedDense")] <- "1"
SS03$LCC_combo <- SS03$LCC_OFF1
SS03$LCC_combo[is.na(SS03$LCC_OFF1)] <- SS03$LCC_OFF2[is.na(SS03$LCC_OFF1)]
}

## 4x4 stuff is not done
if (FALSE) {
## Grid ID 4x4 km
SS_grid <- read.csv(file.path(ROOT, "BAMBBS_Gridcode.csv"))
rownames(SS_grid) <- SS_grid$SS
compare.sets(rownames(SS01), rownames(SS_grid))
SS_grid <- SS_grid[rownames(SS01),"gridcode",drop=FALSE]
levels(SS_grid$gridcode) <- gsub(",", "", levels(SS_grid$gridcode))
}

## Road: dist to, class, #lanes, surface
SS_road <- sqlFetch(con, "dbo_BAMBBS_2015_NearDistanceRoadJoin1000M")
rownames(SS_road) <- SS_road$SS
compare.sets(rownames(SS01), rownames(SS_road))
SS_road <- SS_road[rownames(SS01),]
SS_road$d2road <- SS_road[["Distance to Road"]]
table(SS_road$ROADCLASS, SS_road$d2road > 0, useNA="a")
table(SS_road$NBRLANES, SS_road$d2road > 0, useNA="a")
table(SS_road$PAVSTATUS, SS_road$d2road > 0, useNA="a")
table(SS_road$ROADCLASS, SS_road$NBRLANES, useNA="a")
SS_road <- SS_road[,c("d2road","ROADCLASS","NBRLANES","PAVSTATUS")]
## need to exclude # lanes >2

## Fire: year, size
SS_fire <- sqlFetch(con, "dbo_BBSBAM_2015_FIRE")
rownames(SS_fire) <- SS_fire$SS
compare.sets(rownames(SS01), rownames(SS_fire))
SS_fire <- SS_fire[rownames(SS01),]
SS_fire <- SS_fire[,c("Year","SIZE_HA")]
colnames(SS_fire) <- c("YearFire","FIRE_HA")

## Terrain: slope, twi, elev
SS_terr <- sqlFetch(con, "dbo_BBSBAM_2015_TERRAIN90")
rownames(SS_terr) <- SS_terr$SS
compare.sets(rownames(SS01), rownames(SS_terr))
SS_terr <- SS_terr[rownames(SS01),]
t(table(is.na(SS_terr$cti90), SS01$PCODE)) # mostly affects BBS in some states
SS_terr <- SS_terr[,c("slp90","cti90","elv90")]

## Climate variables from DS and NALCMS 4x4 level
SS_clim <- sqlFetch(con, "dbo_BBSBAM_2015__CLIMLU")
rownames(SS_clim) <- SS_clim$SS
compare.sets(rownames(SS01), rownames(SS_clim))
SS_clim <- SS_clim[rownames(SS01),]
tmp <- as.matrix(SS_clim[,grepl("NALCMS05_", colnames(SS_clim))])
SS_clim <- SS_clim[,!grepl("NALCMS05_", colnames(SS_clim))]
colnames(tmp) <- gsub("NALCMS05_", "", colnames(tmp))
Col <- as.character(ltnalc$Label)[match(colnames(tmp), as.character(ltnalc$Value))]
Col[is.na(Col)] <- "Water"
## 4x4 stuff is not done
#SS_NALC4x4 <- data.frame(groupSums(tmp, 2, Col, na.rm=TRUE))
#colnames(SS_NALC4x4) <- paste0("GRID4_NALC_", colnames(SS_NALC4x4))
SS_clim$NALCMS05 <- NULL
SS_clim$PCODE <- NULL
SS_clim$SS <- NULL

## 4x4 stuff is not done
if (FALSE) {
## LCC05 4x4 level
SS_LCC4x4 <- sqlFetch(con, "dbo_BBSBAM_V4_LCC05CND_4X4SUMM")
SS_LCC4x4 <- nonDuplicated(SS_LCC4x4, SS, TRUE)
#rownames(SS_LCC4x4) <- SS_LCC4x4$SS
compare.sets(rownames(SS01), rownames(SS_LCC4x4))
SS_LCC4x4 <- SS_LCC4x4[rownames(SS01),]
SS_LCC4x4$SS <- NULL
SS_LCC4x4$gridcode <- NULL
SS_LCC4x4$LCCVVSUM <- NULL
SS_LCC4x4 <- as.matrix(SS_LCC4x4)
colnames(SS_LCC4x4) <- gsub("LCCVV", "", colnames(SS_LCC4x4))
#Col <- as.character(ltlcc$BAMLCC05V2_label1)[match(colnames(SS_LCC4x4), 
#    as.character(ltlcc$lcc05v1_2))]
Col <- as.character(ltlcc$BAMLCC05V2_label2)[match(colnames(SS_LCC4x4), 
    as.character(ltlcc$lcc05v1_2))]
Col[is.na(Col)] <- "BARREN"
SS_LCC4x4 <- data.frame(groupSums(SS_LCC4x4, 2, Col, na.rm=TRUE))
SS_LCC4x4[is.na(SS03$HAB_LCC2),] <- NA
colnames(SS_LCC4x4) <- paste0("GRID4_LCC_", colnames(SS_LCC4x4))

## EOSD 4x4 level
SS_EOSD4x4 <- sqlFetch(con, "dbo_BBSBAM_V4_EOSD_4X4SUMM")
SS_EOSD4x4$upsize_ts <- NULL
rownames(SS_EOSD4x4) <- SS_EOSD4x4$SS
compare.sets(rownames(SS01), rownames(SS_EOSD4x4))
SS_EOSD4x4 <- SS_EOSD4x4[match(rownames(SS01), rownames(SS_EOSD4x4)),]
rownames(SS_EOSD4x4) <- SS01$SS

SS_EOSD4x4 <- as.matrix(SS_EOSD4x4[,grepl("eosdVV", colnames(SS_EOSD4x4))])
colnames(SS_EOSD4x4) <- gsub("eosdVV", "", colnames(SS_EOSD4x4))

#Col <- as.character(lteosd$Reclass_label1)[match(colnames(SS_EOSD4x4), 
#    as.character(lteosd$Value))]
Col <- as.character(lteosd$Reclass_label2)[match(colnames(SS_EOSD4x4), 
    as.character(lteosd$Value))]
Col[is.na(Col)] <- "BARREN"
SS_EOSD4x4 <- data.frame(groupSums(SS_EOSD4x4, 2, Col, na.rm=TRUE))
SS_EOSD4x4[is.na(SS03$HAB_EOSD2),] <- NA
colnames(SS_EOSD4x4) <- paste0("GRID4_EOSD_", colnames(SS_EOSD4x4))
}

## HEIGHT (Simard)
SS_height <- sqlFetch(con, "dbo_Height")
SS_height <- nonDuplicated(SS_height, SS, TRUE)
compare.sets(rownames(SS01), rownames(SS_height))
SS_height <- SS_height[rownames(SS01),]
SS_height <- SS_height[,"HEIGHTSIMARD",drop=FALSE]

if (FALSE) {
## Nature Serve range: 3 spp (Can clipped range used 0/1)
SS_nserv <- sqlFetch(con, "dbo_BBSBAM_SARSPPLOCATIONSRange")
SS_nserv <- nonDuplicated(SS_nserv, SS, TRUE)
compare.sets(rownames(SS01), rownames(SS_nserv))
SS_nserv <- SS_nserv[rownames(SS01),]
SS_nserv <- SS_nserv[,c("CAWAINOUT","OSFLINOUT","CONIINOUT")]
}

## GFW yearly loss intersections and 1st year of loss
#SS_gfw <- sqlFetch(con, "dbo_BAMBBS_GFWLossYear")
SS_gfw <- read.csv(file.path(ROOT, "GFWLossYear.csv"))
SS_gfw <- nonDuplicated(SS_gfw, SS, TRUE)
compare.sets(rownames(SS01), rownames(SS_gfw))
levels(SS_gfw$YearLoss) <- gsub(",", "", levels(SS_gfw$YearLoss))
SS_gfw$YearLoss <- as.integer(as.character(SS_gfw$YearLoss))
SS_gfw <- SS_gfw[rownames(SS01),]
SS_gfw <- SS_gfw[,"YearLoss",drop=FALSE]

## Pasher disturbance
SS_pash <- read.csv(file.path(ROOT, "bambbs2015beadandpasher.csv"))
SS_pash <- nonDuplicated(SS_pash, SS, TRUE)
compare.sets(rownames(SS01), rownames(SS_pash))
SS_pash <- SS_pash[rownames(SS01),]
SS_pash <- SS_pash[,c("BEADTotalL","BEADtotPol")]

## Local spring

SS_sprng <- read.csv("e:/peter/bam/May2015/NRCAN_SG_001_BAMBBS2015_71_13.csv")
SS_sprng <- SS_sprng[,c("SS","RASTERVALU")]
SS_sprng$SPRNG <- SS_sprng$RASTERVALU
levels(SS_sprng$SPRNG) <- gsub(",", "", levels(SS_sprng$SPRNG))
SS_sprng$SPRNG <- as.numeric(as.character(SS_sprng$SPRNG))
SS_sprng$SPRNG[SS_sprng$SPRNG < 0] <- NA
rownames(SS_sprng) <- SS_sprng$SS
SS_sprng <- SS_sprng[rownames(SS01),]

## Put together the main SS level object
SS <- data.frame(
    PCODE=SS01$PCODE,
    SS=SS01$SS,
    X=SS01$X_GEONAD83,
    Y=SS01$Y_GEONAD83,
    Xcl=SS01$X_CLCC,
    Ycl=SS01$Y_CLCC,
    JURS=SS01$JURSALPHA,
    COUNTRY=SS01$COUNTRY,
    TZONE=SS01$TZONE_CODE,
    BOREALLOC=SS01$BOREALLOC,
    BCR=as.factor(SS01$BCR),
    TREE=SS02$TREE,
    TREE3=SS02$TREE3,
    SPRNG=SS_sprng$SPRNG,
    #LCC05_PT=SS03$LCC05_PT, # -- FOR NICOLE
    SS03[,c("HAB_NALC2", "HAB_NALC1")],
    #SS_grid,
    #SS_nserv,
    SS_road, 
    SS_terr, 
    SS_fire, 
    SS_clim, 
    SS_pash, 
    SS_gfw, 
    SS_height)
    #SS_NALC4x4, 
    #SS_LCC4x4, 
    #SS_EOSD4x4)

#### Project summary table

## This table needed local tweaks to be operable
#PCODE <- sqlFetch(con, "dbo_National_Proj_Summary_V4_2015")
#PCODE$SSMA_TimeStamp <- NULL
PCODE <- read.csv(file.path(ROOT,"proj.csv"))
levels(PCODE$Maxdist) <- tolower(levels(PCODE$Maxdist))
levels(PCODE$Maxdist)[levels(PCODE$Maxdist)=="unlimited"] <- "Inf"
PCODE$Maxdist[PCODE$Maxdist=="unknown"] <- NA
PCODE$Maxdist <- droplevels(PCODE$Maxdist)
PCODE$Maxdist <- as.numeric(as.character(PCODE$Maxdist))
PCODE$Maxdur <- pmin(PCODE$MaxDuration, 10)

#### Survey level fields

## Pull in point count tables
pkbam <- sqlFetch(con, "dbo_National_PKEY_V4_2015")
pkbam$SSMA_TimeStamp <- NULL
#pkbbs <- sqlFetch(con, "dbo_PKEY_BBS_V3_2015")
pkbbs <- sqlFetch(con2, "PKEY_BBS_V4_2016")

## What columns to retain
PKCOLS <- c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND",
    "YEAR","MONTH","DAY","HOUR","MIN","PART","MAXDUR","MAXDIS")

colnames(pkbam) <- toupper(colnames(pkbam))
pkbam$MAXDUR <- PCODE$Maxdur[match(pkbam$METHOD, PCODE$Method)]
pkbam$MAXDIS <- PCODE$Maxdist[match(pkbam$METHOD, PCODE$Method)]
pkbam$PART <- 1L
pkbam$MONTH <- pkbam$MM
pkbam$DAY <- pkbam$DD
pkbam$HOUR <- pkbam$HR
pkbam$YEAR <- pkbam$YYYY
levels(pkbam$ROUND) <- sub("[[:alpha:]]+$", "", levels(pkbam$ROUND))
pkbam$ROUND <- as.integer(as.character(pkbam$ROUND))
pkbam <- pkbam[,PKCOLS]

colnames(pkbbs) <- toupper(colnames(pkbbs))
pkbbs$MAXDUR <- 3
pkbbs$MAXDIS <- Inf
pkbbs$PART <- 2L
pkbbs$METHOD <- as.factor("BBS:9999")
pkbbs$SITE <- as.factor(pkbbs$SITE)
pkbbs$YEAR <- pkbbs$YYYY
pkbbs$MONTH <- pkbbs$MM
pkbbs$DAY <- pkbbs$DD
pkbbs$HOUR <- pkbbs$HR
pkbbs <- pkbbs[,PKCOLS]

PKEY <- rbind(pkbam, pkbbs)
#rm(pkbam, pkbbs)
#gc()

## Map `METHOD` field from project summary table onto `PKEY$METHOD`
## so that duration and distance method can be carried forward to 
## point count table
levels(PCODE$Method)[levels(PCODE$Method) == "QCAtlas:118"] <- "QCATLAS:118"
compare.sets(PCODE$Method, PKEY$METHOD)
setdiff(PKEY$METHOD, PCODE$Method)
setdiff(PCODE$Method, PKEY$METHOD)
PKEY$DURMETH <- PCODE$DURMETH[match(PKEY$METHOD, PCODE$Method)]
PKEY$DISMETH <- PCODE$DISTMETH[match(PKEY$METHOD, PCODE$Method)]

## Identifying roadside surveys
PKEY$ROAD <- 0L
treat.as.bbs <- c("HOBBS","CF","MNBBA", levels(pkbbs$PCODE))
PKEY$ROAD[PKEY$PCODE %in% treat.as.bbs] <- 1L

#### Offset specific variables

## Date/time components
PKEY$MIN[is.na(PKEY$MIN)] <- 0 # min that critical, is not -- said Yoda
MM <- ifelse(PKEY$MONTH < 10, paste0("0", PKEY$MONTH), as.character(PKEY$MONTH))
HH <- ifelse(PKEY$HOUR < 10, paste0("0", PKEY$HOUR), as.character(PKEY$HOUR))
mm <- ifelse(PKEY$MIN < 10, paste0("0", PKEY$MIN), as.character(PKEY$MIN))
#mm[is.na(mm) & !is.na(HH)] <- "00"
DD <- with(PKEY, paste0(YEAR, "-", MM, "-", DAY, " ", HH, ":", mm, ":00"))
DD <- strptime(DD, "%Y-%m-%e %H:%M:%S")
PKEY$DATE <- DD
## Julian day
PKEY$JULIAN <- DD$yday # this is kept as original
PKEY$JDAY <- DD$yday / 365
summary(PKEY$JDAY)
## prevent too far extrapolation
PKEY$JDAY[PKEY$JDAY < 0.35 | PKEY$JDAY > 0.55] <- NA
## TSSR = time since sunrise
Coor <- as.matrix(cbind(as.numeric(SS$X),as.numeric(SS$Y)))[match(PKEY$SS, rownames(SS)),]
JL <- as.POSIXct(DD)
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
PKEY$srise <- NA
PKEY$srise[subset] <- sr
PKEY$start_time <- PKEY$HOUR + PKEY$MIN/60
TZ <- SS$TZONE[match(PKEY$SS, rownames(SS))]
lttz <- read.csv("~/repos/bamanalytics//lookup/tzone.csv")
lttz <- nonDuplicated(lttz, Timezone, TRUE)
PKEY$MDT_offset <- lttz$MDT_offset[match(TZ, rownames(lttz))]
table(TZ, PKEY$MDT_offset)
PKEY$TSSR <- (PKEY$start_time - PKEY$srise + PKEY$MDT_offset) / 24
PKEY$TSSR_orig <- PKEY$TSSR # keep a full copy
PKEY$TSSR[PKEY$start_time > 12] <- NA ## after noon
summary(PKEY$TSSR)
summary(PKEY$start_time)

PKEY <- PKEY[PKEY$DURMETH != "J",] # unknown duration
PKEY <- PKEY[PKEY$DISMETH != "O",] # unknown distance
PKEY <- droplevels(PKEY)

## QC Atlas problem
#with(PKEY, table(PCODE, is.na(MAXDUR)))
#with(PKEY, table(PCODE, is.na(MAXDIS)))


#### Point count tables and methodology

## Some of these tables are also tweaked locally (see below)
BEH <- sqlFetch(con, "dbo_DD_DescripBEH")
#DISINT <- sqlFetch(con, "dbo_DD_DescripDistance")
#DISINT$SSMA_TimeStamp <- NULL
#DURINT <- sqlFetch(con, "dbo_DD_DescripPeriod")
DISMET <- sqlFetch(con, "dbo_DD_distance_codes_methodology")
DURMET <- sqlFetch(con, "dbo_DD_duration_codes_methodology")
## Trailing whitespace removed from factor levels
levels(DISMET$DISTANCECODE) <- gsub(" *$", "", levels(DISMET$DISTANCECODE))
levels(DURMET$DURATIONCODE) <- gsub(" *$", "", levels(DURMET$DURATIONCODE))

#### Point count tables
pcbam <- sqlFetch(con, "dbo_National_PtCount_V4_2015")
pcbam$SSMA_TimeStamp <- NULL

pcbbs <- sqlFetch(con2, "POINTCOUNT_BBS_V4_2016")
colnames(pcbbs) <- toupper(colnames(pcbbs))
colnames(pcbbs)[colnames(pcbbs) == "SPECIES_ID"] <- "SPECIES"
colnames(pcbbs)[colnames(pcbbs) == "PERIOD"] <- "DURATION"
pcbbs$PCODE <- "BBS"

## Columns to keep
pccols <- c("PCODE","SS","PKEY","DURATION","DISTANCE",
    "SPECIES","ABUND","BEH")
## Close the database connections
close(con)
close(con2)

## Duration and distance intervals (locally tweaked)
DURINT <- read.csv("~/repos/bamanalytics/lookup/durint.csv")
DISINT <- read.csv("~/repos/bamanalytics/lookup/disint.csv")
DURINT$dur <- paste0(DURINT$Dur_Start, "-", DURINT$DUR_end)
DISINT$dis <- paste0(DISINT$DIST_START, "-", DISINT$DIST_END)
rownames(DURINT) <- DURINT[,1]
rownames(DISINT) <- DISINT[,1]
## Combined point count table
PCTBL <- rbind(pcbam[,pccols], pcbbs[,pccols])
## Mapping duration and distance intervals 
PCTBL$dur <- as.factor(DURINT$dur[match(PCTBL$DURATION, rownames(DURINT))])
PCTBL$dis <- as.factor(DISINT$dis[match(PCTBL$DISTANCE, rownames(DISINT))])
## Methodology
PCTBL$DISMETH <- droplevels(PKEY$DISMETH[match(PCTBL$PKEY, PKEY$PKEY)])
PCTBL$DURMETH <- droplevels(PKEY$DURMETH[match(PCTBL$PKEY, PKEY$PKEY)])

## Filtering surveys (need to exclude PKEY)
keeppkey <- rep(TRUE, nrow(PCTBL))
## 11=0-20
## 8=unk
keeppkey[PCTBL$DURATION %in% c(11,8)] <- FALSE
## Excluding unknown distance bands
keeppkey[PCTBL$DISTANCE %in% c(4,5,9)] <- FALSE
## Excluding unknown duration methodology
keeppkey[PCTBL$DISMETH == "J"] <- FALSE
## Excluding unknown distance methodology
keeppkey[PCTBL$DURMETH == "O"] <- FALSE
## Actual filtering -- but dropping PKEYs
PCTBL <- droplevels(PCTBL[keeppkey,])

## Filtering within survey (do not exclude PKEY)
## Filtering behaviour
#sort(100 * table(PCTBL$BEH) / sum(table(PCTBL$BEH)))
## 1=Heard
## 11=no birds observed at station - added 2011
## 6=seen and heard
## Excluding non-aerial detections
keep <- rep(TRUE, nrow(PCTBL))
keep[!(PCTBL$BEH %in% c("1","6","11"))] <- FALSE
## this is fake, but there is no other option until a fix
keep[is.na(PCTBL$BEH)] <- TRUE # dont know what this is

if (FALSE) {
PCTBLx <- PCTBL
PCTBL$YEAR <- PKEY$YEAR[match(PCTBL$PKEY, PKEY$PKEY)]
PCTBL$JURS <- SS$JURS[match(PCTBL$SS, SS$SS)]
PCTBL <- PCTBL[PCTBL$JURS=="AB" & PCTBL$PCODE=="BBS",]
table(PCTBL$YEAR, PCTBL$BEH, useNA="a")
PCTBL$YEAR <- PKEY$YEAR[match(PCTBL$PKEY, PKEY$PKEY)]
PCTBL$JURS <- SS$JURS[match(PCTBL$SS, SS$SS)]
with(PCTBL, aggregate(ABUND, list(Yr=YEAR), sum))
}

## check species with high visual detection rates
if (FALSE) {
    xtx <- Xtab(ABUND ~ SPECIES + BEH, PCTBL)
    save(BEH, xtx, file=file.path(ROOT, "out", "spp-beh.Rdata"))
}
## Excluding >10 min intervals
## 10=10-20
## 3=before or after
## 9=10-15
keep[PCTBL$DURATION %in% c(10,3,9)] <- FALSE
## Excluding NA values
keep[is.na(PCTBL$dur)] <- FALSE
keep[is.na(PCTBL$dis)] <- FALSE
keep[is.na(PCTBL$ABUND)] <- FALSE
## Actual filtering -- but keeping PKEYs (do not drop levels)
#PCTBL$keep <- keep
PCTBL <- PCTBL[keep,]

## Excluding/dropping species

PCTBL$SPECIES <- droplevels(PCTBL$SPECIES)
levels(PCTBL$SPECIES) <- toupper(levels(PCTBL$SPECIES))
compare.sets(PCTBL$SPECIES, TAX$Species_ID)
setdiff(PCTBL$SPECIES, TAX$Species_ID)
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "YWAR"] <- "YEWA"

levels(TAX$Species_ID)[levels(TAX$Species_ID) == "SCJU"] <- "DEJU" # change SCJU to DEJU
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "MYWA"] <- "YRWA" # change MYWA to YRWA
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "COSN"] <- "WISN" # change COSN to WISN


sp <- levels(PCTBL$SPECIES)
spp.to.exclude <- rep(FALSE, length(sp))
spp.to.exclude[sp %in% setdiff(PCTBL$SPECIES, TAX$Species_ID)] <- TRUE
spp.to.exclude[sp == "NONE"] <- TRUE
spp.to.exclude[sp == "RESQ"] <- TRUE
spp.to.exclude[sp == "DUCK"] <- TRUE
spp.to.exclude[nchar(sp) > 4] <- TRUE
spp.to.exclude[substr(sp, 1, 2) == "UN"] <- TRUE
levels(PCTBL$SPECIES)[spp.to.exclude] <- "NONE"
compare.sets(PCTBL$SPECIES, TAX$Species_ID)
setdiff(PCTBL$SPECIES, TAX$Species_ID)
## Excluding columns
PCTBL$DURATION <- NULL
PCTBL$DISTANCE <- NULL
PCTBL$BEH <- NULL
PCTBL$dur <- droplevels(PCTBL$dur)
PCTBL$dis <- droplevels(PCTBL$dis)

compare.sets(SS$SS, PKEY$SS)
compare.sets(SS$SS, PCTBL$SS)
compare.sets(PKEY$PKEY, PCTBL$PKEY)

### ABMI data processing

pcabmi <- read.csv(file.path(ROOT, "abmi", "birds.csv"))
ssabmi <- read.csv(file.path(ROOT, "abmi", "sitemetadata.csv"))

## Labels etc
## SK alpac sites to exclude (xy is hard to track down)
pcabmi <- droplevels(pcabmi[!grepl("ALPAC-SK", pcabmi$SITE_LABEL),])

tmp <- do.call(rbind, sapply(levels(pcabmi$SITE_LABEL), strsplit, "_"))
colnames(tmp) <- c("Protocol", "OnOffGrid", "DataProvider", "SiteLabel", "YYYY", "Visit", "SubType", "BPC")
tmp2 <- sapply(tmp[,"SiteLabel"], strsplit, "-")
tmp3 <- sapply(tmp2, function(z) if (length(z)==1) "ABMI" else z[2])
tmp4 <- sapply(tmp2, function(z) if (length(z)==1) z[1] else z[3])
tmp <- data.frame(tmp, ClosestABMISite=tmp4)
tmp$DataProvider <- as.factor(tmp3)
tmp$Label <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, "PC", BPC, sep="_"))
tmp$Label2 <- with(tmp, paste(OnOffGrid, DataProvider, SiteLabel, YYYY, Visit, sep="_"))
tmp$ClosestABMISite <- as.integer(as.character(tmp$ClosestABMISite))
tmp$lat <- ssabmi$PUBLIC_LATTITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
tmp$long <- ssabmi$PUBLIC_LONGITUDE[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
#tmp$NatReg <- ssabmi$NATURAL_REGIONS[match(tmp$ClosestABMISite, ssabmi$SITE_ID)]
#tmp$boreal <- tmp$NatReg %in% c(c("Boreal", "Canadian Shield", "Foothills", "Rocky Mountain"))

pcabmi <- data.frame(pcabmi, tmp[match(pcabmi$SITE_LABEL, rownames(tmp)),])

## PKEY table and proper date format
PKEY_abmi <- nonDuplicated(pcabmi, pcabmi$Label, TRUE)
tmp <- PKEY_abmi$ADATE
tmp <- sapply(as.character(tmp), strsplit, split="-")
for (i in 1:length(tmp)) {
    if (length(tmp[[i]])<3) {
        tmp[[i]] <- rep("99", 3)
    }
}
table(sapply(tmp, "[[", 2))
for (i in 1:length(tmp)) {
    tmp[[i]][2] <- switch(tmp[[i]][2],
        "May"=5, "Jun"=6, "Jul"=7, "Aug"=8, "99"=99)
}
tmp <- sapply(tmp, function(z) paste("20",z[3],"-",z[2],"-",z[1], sep=""))
tmp[tmp=="2099-99-99"] <- NA
PKEY_abmi$Date <- as.POSIXct(tmp, tz="America/Edmonton")

## TSSR
Coor <- as.matrix(cbind(as.numeric(PKEY_abmi$long),as.numeric(PKEY_abmi$lat)))
JL <- as.POSIXct(PKEY_abmi$Date, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
PKEY_abmi$srise_MDT <- NA
PKEY_abmi$srise_MDT[subset] <- sr

tmp <- strsplit(as.character(PKEY_abmi$TBB_START_TIME), ":")
id <- sapply(tmp,length)==2
tmp <- tmp[id]
tmp <- as.integer(sapply(tmp,"[[",1)) + as.integer(sapply(tmp,"[[",2))/60
PKEY_abmi$start_time <- NA
PKEY_abmi$start_time[id] <- tmp
PKEY_abmi$srise <- PKEY_abmi$srise_MDT
PKEY_abmi$TSSR <- (PKEY_abmi$start_time - PKEY_abmi$srise) / 24 # MDT offset is 0

## Julian day
PKEY_abmi$jan1 <- as.Date(paste(PKEY_abmi$YEAR, "-01-01", sep=""))
PKEY_abmi$JULIAN <- as.numeric(as.Date(PKEY_abmi$Date)) - as.numeric(PKEY_abmi$jan1) + 1
PKEY_abmi$JULIAN[PKEY_abmi$JULIAN > 365] <- NA
PKEY_abmi$JDAY <- PKEY_abmi$JULIAN / 365

if (FALSE) {

aa <- PKEY_abmi
CL <- rgb(210, 180, 140, alpha=0.25*255, max=255)

op <- par(mfrow=c(3,1))
boxplot(JULIAN~YEAR,aa, ylab="Julian day")
abline(h=seq(0, 360, by=30), col="grey")
text(rep(0.5, 12), seq(0, 360, by=30)+25, c("Jan","Feb","Mar","Apr","May","Jun",
    "Jul", "Aug","Sep","Oct","Nov","Dec"))
boxplot(JULIAN~YEAR,aa, add=TRUE, col=CL)
abline(h=mean(aa$JULIAN, na.rm=TRUE), col=2, lty=2, lwd=2)

boxplot(start_time~YEAR,aa, ylab="Start time (hours)")
abline(h=seq(0, 24, by=2), col="grey")
boxplot(start_time~YEAR,aa, add=TRUE, col=CL)
abline(h=mean(aa$start_time, na.rm=TRUE), col=2, lty=2, lwd=2)

plot(start_time ~ srise, aa, pch=19, col=CL,
    xlab="Sunrise time (24 hour clock)", ylab="Start time (24 hour clock)")
abline(0, 1)
for (i in c(1:10*2))
    abline(i, 1, col="grey")
abline(lm(start_time ~ srise, aa), col=2, lwd=2, lty=2)
box()
par(op)

## Check EMCLA

rownames(PKEY) <- PKEY$PKEY
i <- grepl("EMCLA", rownames(PKEY))
x <- droplevels(PKEY[i,])
hist(x$start_time)
hist(x$JDAY*365)

plot(x$JDAY*365, x$HOUR, col=CL, pch=19, type="n", main="EMCLA",
    ylab="Start time (24 hour clock)", xlab="Julian day")
abline(v=seq(0, 360, by=30), col="grey")
text(seq(0, 360, by=30)+17, rep(13, 12), c("Jan","Feb","Mar","Apr","May","Jun",
    "Jul", "Aug","Sep","Oct","Nov","Dec"))
abline(h=seq(0, 24, by=5), col="grey")
points(x$JDAY*365, x$HOUR, col=CL, pch=19)
}

## counts

PCTBL_abmi <- pcabmi
levels(PCTBL_abmi$COMMON_NAME)[levels(PCTBL_abmi$COMMON_NAME) == "Black and White Warbler"] <- "Black-and-white Warbler"
compare.sets(TAX$English_Name, PCTBL_abmi$COMMON_NAME)
compare.sets(TAX$Scientific_Name, PCTBL_abmi$SCIENTIFIC_NAME)
setdiff(PCTBL_abmi$COMMON_NAME, TAX$English_Name)
setdiff(PCTBL_abmi$SCIENTIFIC_NAME, TAX$Scientific_Name)

PCTBL_abmi$SPECIES <- TAX$Species_ID[match(PCTBL_abmi$COMMON_NAME, TAX$English_Name)]
levels(PCTBL_abmi$SPECIES) <- c(levels(PCTBL_abmi$SPECIES), "NONE")
PCTBL_abmi$SPECIES[PCTBL_abmi$SPECIES == "TERN_UNI"] <- "NONE"
PCTBL_abmi$SPECIES <- droplevels(PCTBL_abmi$SPECIES)

PCTBL_abmi$TBB_TIME_1ST_DETECTED[PCTBL_abmi$TBB_TIME_1ST_DETECTED %in% 
    c("DNC", "NONE", "VNA")] <- NA
PCTBL_abmi$TBB_TIME_1ST_DETECTED <- as.numeric(as.character(PCTBL_abmi$TBB_TIME_1ST_DETECTED))
PCTBL_abmi$period1st <- as.numeric(cut(PCTBL_abmi$TBB_TIME_1ST_DETECTED, c(-1, 200, 400, 600)))

PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_1 %in% c("0","1"),]
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_2 %in% c("0","1"),]
PCTBL_abmi <- PCTBL_abmi[PCTBL_abmi$TBB_INTERVAL_3 %in% c("0","1"),]
PCTBL_abmi$TBB_INTERVAL_1 <- as.integer(PCTBL_abmi$TBB_INTERVAL_1) - 1L
PCTBL_abmi$TBB_INTERVAL_2 <- as.integer(PCTBL_abmi$TBB_INTERVAL_2) - 1L
PCTBL_abmi$TBB_INTERVAL_3 <- as.integer(PCTBL_abmi$TBB_INTERVAL_3) - 1L

tmp <- col(matrix(0,nrow(PCTBL_abmi),3)) * 
    PCTBL_abmi[,c("TBB_INTERVAL_1","TBB_INTERVAL_2","TBB_INTERVAL_3")]
tmp[tmp==0] <- NA
tmp <- cbind(999,tmp)
PCTBL_abmi$period123 <- apply(tmp, 1, min, na.rm=TRUE)
with(PCTBL_abmi, table(period1st, period123))
PCTBL_abmi$period1 <- pmin(PCTBL_abmi$period1st, PCTBL_abmi$period123)
with(PCTBL_abmi, table(period1st, period1))
with(PCTBL_abmi, table(period123, period1))


## Data package for new offsets
dat <- data.frame(PKEY[,c("PCODE","PKEY","SS","YEAR","TSSR","JDAY","JULIAN",
    "srise","start_time","MAXDUR","MAXDIS","METHOD","DURMETH","DISMETH","ROAD")],
    SS[match(PKEY$SS, rownames(SS)),c("TREE","TREE3","LCC_combo","HAB_NALC1","HAB_NALC2",
    "BCR","JURS","SPRNG","DD51","X","Y")], NR=NA)
dat <- dat[dat$ROAD == 0,]
rownames(dat) <- dat$PKEY
ii <- intersect(dat$PKEY, levels(PCTBL$PKEY))
dat <- droplevels(dat[ii,])
summary(dat)
colSums(is.na(dat))
## sra and edr might have different NA patterns -- it is OK to exclude them later
#dat <- dat[rowSums(is.na(dat)) == 0,]
dat <- droplevels(dat)
dat$TSLS <- (dat$JULIAN - dat$SPRNG) / 365
dat$DD5 <- (dat$DD51 - 1600) / 1000
dat$DD51 <- NULL

## nat regions to filter grasslands
luf <- read.csv("~/repos/abmianalytics/lookup/sitemetadata.csv")
PKEY_abmi$NR <- luf$NATURAL_REGIONS[match(PKEY_abmi$ClosestABMISite, luf$SITE_ID)]
table(PKEY_abmi$NR)

dat2 <- with(PKEY_abmi, data.frame(
    PCODE="ABMI",
    PKEY=as.factor(Label),
    SS=as.factor(Label2),
    YEAR=YEAR,
    TSSR=TSSR,
    JDAY=JDAY,
    JULIAN=JULIAN,
    srise=srise,
    start_time=start_time,
    MAXDUR=10,
    MAXDIS=Inf,
    METHOD="ABMI:1",
    DURMETH="X",
    DISMETH="D",
    ROAD=0,
    TREE=NA,
    TREE3=NA,
    LCC_combo=NA,
    HAB_NALC1=NA,
    HAB_NALC2=NA,
    BCR=factor("6", levels(dat$BCR)),
    JURS=factor("AB", levels=levels(dat$JURS)),
    SPRNG=NA,
    X=long,
    Y=lat,
    TSLS=NA,
    DD5=NA,
    NR=NR))
rownames(dat2) <- dat2$PKEY
#write.csv(dat2, row.names=FALSE, file="ABMI-XY.csv")
ls_abmi <- read.csv("e:/peter/bam/May2015/ABMI_XY_JDStart_DD5.csv")
rownames(ls_abmi) <- ls_abmi$PKEY
ls_abmi <- ls_abmi[rownames(dat2),]
dat2$SPRNG <- ls_abmi$JSD_00_13
dat2$TSLS <- (dat2$JULIAN - dat2$SPRNG) / 365

## besides `dat` we also need specific aoutput from `PCTBL`
## and also the methodology x interval lookup

compare.sets(DISMET$DISTANCECODE, PKEY$DISMETH)
compare.sets(DURMET$DURATIONCODE, PKEY$DURMETH)

## getting problematic ones to TF
if (FALSE) {

PCTBL$issue <- character(nrow(PCTBL))
PCTBL$issue[with(PCTBL, DURMETH=="A" & dur=="0-3")] <- "DURMETH=A: 0-3 -> 0-10"
PCTBL$issue[with(PCTBL, DURMETH=="B" & dur=="5-8")] <- "DURMETH=B: 5-8 -> 0-5"
PCTBL$issue[with(PCTBL, DURMETH=="X" & dur=="10-10")] <- "DURMETH=X: 10-10 -> 6.66-10"

PCTBL$issue[with(PCTBL, DISMETH=="B" & dis=="0-Inf")] <- "DISMETH=B: 0-Inf -> 0-50 //best guess"
PCTBL$issue[with(PCTBL, DISMETH=="C" & dis=="0-Inf")] <- "DISMETH=C: 0-Inf -> 0-50 //best guess"
PCTBL$issue[with(PCTBL, DISMETH=="F")] <- "DISMETH=F: all -> 0-100 //all kinds of weird stuff"
PCTBL$issue[with(PCTBL, DISMETH=="I" & dis=="100-125")] <- "DISMETH=I: 100-125 -> 0-25"
PCTBL$issue[with(PCTBL, DISMETH=="I" & dis=="100-Inf")] <- "DISMETH=I: 100-Inf -> 0-25 //best guess"
PCTBL$issue[with(PCTBL, DISMETH=="M" & dis=="0-Inf")] <- "DISMETH=M: 0-Inf -> 150-Inf"
PCTBL$issue[with(PCTBL, DISMETH=="U" & dis=="0-50")] <- "DISMETH=U: 0-50 -> 40-50"
PCTBL$issue[with(PCTBL, DISMETH=="U" & dis=="100-150")] <- "DISMETH=U: 100-150 -> 125-150"
PCTBL$issue[with(PCTBL, DISMETH=="U" & dis=="100-Inf")] <- "DISMETH=U: 100-Inf -> 150-Inf"
PCTBL$issue[with(PCTBL, DISMETH=="U" & dis=="50-100")] <- "DISMETH=U: 50-100 -> 90-100"
PCTBL$issue[with(PCTBL, DISMETH=="U" & dis=="50-Inf")] <- "DISMETH=U: 50-Inf -> 150-Inf"
PCTBL$issue[with(PCTBL, DISMETH=="W" & dis=="150-Inf")] <- "DISMETH=W: 150-Inf -> 100-Inf"
PCTBL$issue[with(PCTBL, DISMETH=="W" & dis=="100-125")] <- "DISMETH=W: 100-125 -> 100-Inf"

ISSUE <- PCTBL[PCTBL$issue != "",]

write.csv(ISSUE, row.names=FALSE, file=file.path(ROOT, "out",
    paste0("issues-toTF-", Sys.Date(), ".csv")))

}

## Oddities that should not happen:
PCTBL$dur <- as.character(PCTBL$dur)
PCTBL$dur[with(PCTBL, DURMETH=="A" & dur=="0-3")] <- "0-10"
# PCTBL$dur[with(PCTBL, DURMETH=="B" & dur=="5-8")] <- "0-5" -- fixed on proj summ
PCTBL$dur[with(PCTBL, DURMETH=="X" & dur=="10-10")] <- "6.66-10"
PCTBL$dur <- as.factor(PCTBL$dur)

PCTBL$dis <- as.character(PCTBL$dis)
PCTBL$dis[with(PCTBL, DISMETH=="B" & dis=="0-Inf")] <- "0-50" # best guess
PCTBL$dis[with(PCTBL, DISMETH=="C" & dis=="0-Inf")] <- "0-50" # best guess
PCTBL$dis[with(PCTBL, DISMETH=="F")] <- "0-100" # all kinds of weird stuff
PCTBL$dis[with(PCTBL, DISMETH=="I" & dis=="100-125")] <- "0-25"
PCTBL$dis[with(PCTBL, DISMETH=="I" & dis=="100-Inf")] <- "0-25" # best guess
#PCTBL$dis[with(PCTBL, DISMETH=="L" & dis=="150-Inf")] <- "100-150" # no >150
PCTBL$dis[with(PCTBL, DISMETH=="M" & dis=="0-Inf")] <- "150-Inf"
#PCTBL$dis[with(PCTBL, DISMETH=="T" & dis=="100-Inf")] <- "-" # no >100
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="0-50")] <- "40-50"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="100-150")] <- "125-150"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="100-Inf")] <- "150-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="50-100")] <- "90-100"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="50-Inf")] <- "150-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="W" & dis=="150-Inf")] <- "100-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="W" & dis=="100-125")] <- "100-Inf"
PCTBL$dis <- as.factor(PCTBL$dis)

#PCTBL$dis <- droplevels(PCTBL$dis)
#PCTBL$dur <- droplevels(PCTBL$dur)

pc <- droplevels(PCTBL[PCTBL$PKEY %in% levels(dat$PKEY),])
levels(pc$PKEY) <- c(levels(pc$PKEY), setdiff(levels(dat$PKEY), levels(pc$PKEY)))

pc2 <- with(PCTBL_abmi, data.frame(
    PCODE="ABMI",
    PKEY=as.factor(Label),
    SS=as.factor(Label2),
    SPECIES=SPECIES,
    ABUND=1,
    dur=factor(period1),
    dis="0-Inf",
    DISMETH="D",
    DURMETH="X"))
levels(pc2$dur) <- c("0-3.33","3.33-6.66","6.66-10")


## combine dat, dat2 and pc pc2
dat$BOREALLOC <- SS$BOREALLOC[match(dat$SS, SS$SS)]
dat2$BOREALLOC <- NA
dat0 <- dat
dat <- rbind(
#    dat[!(dat$BCR %in% c(11, 22)),], 
    dat0[!is.na(dat$BOREALLOC) & dat$BOREALLOC != "OUT",], 
    dat2[dat2$NR != "Grassland", colnames(dat0)])
pc <- rbind(pc, pc2[,colnames(pc)])
pc <- droplevels(pc)

durmat <- as.matrix(Xtab(~ DURMETH + dur, pc))
durmat[durmat > 0] <- 1
dismat <- as.matrix(Xtab(~ DISMETH + dis, pc))
dismat[dismat > 0] <- 1

ltdur <- arrange.intervals(durmat)
ltdis <- arrange.intervals(dismat)
## divide by 100
ltdis$end <- ltdis$end / 100

if (FALSE) {
pcc <- nonDuplicated(pc, PKEY, TRUE)
ii <- intersect(rownames(pcc), rownames(dat))
pkk <- dat[ii,]
pcc <- pcc[ii,]
table(pcc=droplevels(pcc$DISMET), pkk=droplevels(pkk$DISMET), useNA="a")
table(pcc=droplevels(pcc$DURMET), pkk=droplevels(pkk$DURMET), useNA="a")

}

chf <- function() {
PCTBL$YEAR <- PKEY$YEAR[match(PCTBL$PKEY, PKEY$PKEY)]
PCTBL$JURS <- SS$JURS[match(PCTBL$SS, SS$SS)]
with(PCTBL[PCTBL$JURS=="AB" & PCTBL$PCODE=="BBS",], aggregate(ABUND, list(Yr=YEAR), sum))
}
chf()

save(dat2, pc2, 
    file=file.path(ROOT, "out",
    paste0("abmi_data_package_", Sys.Date(), ".Rdata")))

if (FALSE) {
PCTBL$YEAR <- PKEY$YEAR[match(PCTBL$PKEY, PKEY$PKEY)]
PCTBL$JURS <- SS$JURS[match(PCTBL$SS, SS$SS)]
with(PCTBL[PCTBL$JURS=="AB" & PCTBL$PCODE=="BBS",], aggregate(ABUND, list(Yr=YEAR), sum))


pcbbs$YEAR <- PKEY$YEAR[match(pcbbs$PKEY, PKEY$PKEY)]
pcbbs$JURS <- SS$JURS[match(pcbbs$SS, SS$SS)]
with(pcbbs[pcbbs$JURS=="AB",], aggregate(ABUND, list(Yr=YEAR), sum))
}

save(SS, PKEY, PCTBL, TAX,
    file=file.path(ROOT, "out",
    paste0("data_package_", Sys.Date(), ".Rdata")))

save(dat, pc, ltdur, ltdis, TAX,
    file=file.path(ROOT, "out",
    paste0("new_offset_data_package_", Sys.Date(), ".Rdata")))


#### Calculate the offsets (optional)

if (FALSE) { # BEGIN offset calculations !!! NEW VERSION ------------------- !!!!!!!!!!!!!!

ROOT <- "c:/bam/May2015"
source("~/repos/bamanalytics/R/dataprocessing_functions.R")
load(file.path(ROOT, "out", "data_package_2015-07-24.Rdata"))

## summary of detections by species
if (FALSE) {
    library(mefa4)
    ROOT2 <- "e:/peter/bam/pred-2015"
    #spp <- "CAWA"
    xt_by_ss <- Xtab(ABUND ~ SS + SPECIES, PCTBL)
    xt_by_ss[xt_by_ss > 0] <- 1
    xt_by_pk <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)
    xt_by_pk[xt_by_ss > 0] <- 1
    SS$BJ <- interaction(SS$BCR, SS$JURS, sep="_", drop=TRUE)
    PKEY$BJ <- SS$BJ[match(PKEY$SS, SS$SS)]
    bj_by_ss <- SS$BJ[match(rownames(xt_by_ss), SS$SS)]
    bj_by_pk <- PKEY$BJ[match(rownames(xt_by_pk), PKEY$PKEY)]
    xt_by_ss <- xt_by_ss[!is.na(bj_by_ss),]
    bj_by_ss <- bj_by_ss[!is.na(bj_by_ss)]
    xt_by_pk <- xt_by_pk[!is.na(bj_by_pk),]
    bj_by_pk <- bj_by_pk[!is.na(bj_by_pk)]
    
    det_ss <- as.matrix(groupSums(xt_by_ss, 1, bj_by_ss))
    det_pk <- as.matrix(groupSums(xt_by_pk, 1, bj_by_pk))
    all(rownames(det_ss) == rownames(det_pk))

    write.csv(det_ss, row.names=TRUE,
        file=file.path(ROOT2, "species", "cawa-nmbca-tabs", 
        "detections-by-species_ss-level.csv"))
    write.csv(det_pk, row.names=TRUE,
        file=file.path(ROOT2, "species", "cawa-nmbca-tabs", 
        "detections-by-species_pkey-level.csv"))
    
}

library(QPAD)
load_BAM_QPAD(3)
getBAMversion()
sppp <- getBAMspecieslist()

rownames(TAX) <- TAX$Species_ID
tax <- droplevels(TAX[sppp,])
tax$keep <- TRUE
tax$keep[tax$Order %in% c("Gaviiformes","Pelecaniformes","Podicipediformes",
"Gruiformes","Galliformes","Anseriformes","Charadriiformes")] <- FALSE

#data.frame(n=sort(table(tax$Family_CmName_EN)))
#data.frame(n=sort(table(tax$Order)))
## OK
#Apodiformes        1 ok
#Caprimulgiformes   1 ok
#Coraciiformes      1 ok
#Columbiformes      2 ok
#Cuculiformes       2 ok
#Piciformes        10 ok
#Passeriformes    130 ok
## birds of prey
#Falconiformes      2 ???
#Strigiformes       3 ???
#Accipitriformes    7 ???
## water birds
#Gaviiformes        2 xxx
#Pelecaniformes     2 xxx
#Podicipediformes   2 xxx
#Gruiformes         3 xxx
#Galliformes        6 xxx
#Anseriformes       7 xxx
#Charadriiformes   17 xxx

load(file=file.path(ROOT, "out", "spp-beh.Rdata"))
xtx <- as.matrix(xtx[sppp,])
audi <- rowSums(xtx[,c("1","6","11")]) / rowSums(xtx)
by(audi,tax$Order,summary)

sppp <- sppp[tax$keep]

offdat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS")],
    SS[match(PKEY$SS, rownames(SS)),c("TREE","TREE3","HAB_NALC1","HAB_NALC2")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
summary(offdat)
#summary(offdat[PKEY$PCODE=="QCATLAS",]) -- problem solved: QCAtlas vs QCATLAS

#load_BAM_QPAD(version=1)
#BAMspp <- getBAMspecieslist()
#load("~/Dropbox/abmi/intactness/dataproc/BAMCOEFS25.Rdata")

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp

for (i in sppp) {
    cat(i, date(), "\n");flush.console()
    tmp <- try(offset_fun(j=1, i, offdat))
    if (!inherits(tmp, "try-error"))
        OFF[,i] <- tmp
}
## 99-100 percentile can be crazy high (~10^5), thus reset
for (i in sppp) {
    q <- quantile(OFF[,i], 0.99, na.rm=TRUE)
    OFF[!is.na(OFF[,i]) & OFF[,i] > q, i] <- q
}
colSums(is.na(OFF))/nrow(OFF)
apply(exp(OFF), 2, range, na.rm=TRUE)


save(OFF, file=file.path(ROOT, "out",
    paste0("offsets_allspp_BAMBBS_", Sys.Date(), ".Rdata")))

} # END offset calculations

if (FALSE) { # BEGIN offset calculations !!! OLD VERSION ------------------- !!!!!!!!!!!!!!

load(file.path(ROOT, "out", "data_package_2015-07-24.Rdata"))

offdat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS")],
    SS[match(PKEY$SS, rownames(SS)),c("LCC_combo","TREE")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
summary(offdat)
#summary(offdat[PKEY$PCODE=="QCATLAS",]) -- problem solved: QCAtlas vs QCATLAS

load_BAM_QPAD(version=1)
BAMspp <- getBAMspecieslist()
load("~/Dropbox/abmi/intactness/dataproc/BAMCOEFS25.Rdata")

(sppp <- union(BAMspp, BAMCOEFS25$spp))

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp
for (i in sppp) {
    cat(i, date(), "\n");flush.console()
    tmp <- try(offset_fun(j=1, i, offdat))
    if (!inherits(tmp, "try-error"))
        OFF[,i] <- tmp
}
## 99-100 percentile can be crazy high (~10^5), thus reset
for (i in sppp) {
    q <- quantile(OFF[,i], 0.99, na.rm=TRUE)
    OFF[!is.na(OFF[,i]) & OFF[,i] > q, i] <- q
}
colSums(is.na(OFF))/nrow(OFF)
apply(exp(OFF), 2, range, na.rm=TRUE)


save(OFF, file=file.path(ROOT, "out",
    paste0("offsets_allspp_BAMBBS_", Sys.Date(), ".Rdata")))
write.csv(OFF, file=file.path(ROOT, "out",
    paste0("offsets_allspp_BAMBBS_", Sys.Date(), ".csv")))
save(offdat, file=file.path(ROOT, "out",
    paste0("offset_covariates_", Sys.Date(), ".Rdata")))
write.csv(offdat, row.names=FALSE, file=file.path(ROOT, "out",
    paste0("offset_covariates_", Sys.Date(), ".csv")))

} # END offset calculations


######## These are the transformations #################

library(mefa4)
ROOT <- "c:/bam/May2015"
#load(file.path(ROOT, "out", "data_package_2015-07-24.Rdata"))
load(file.path(ROOT, "out", "data_package_2015-08-14.Rdata"))
load(file.path(ROOT, "out", "offsets_allspp_BAMBBS_2015-07-24.Rdata"))

colnames(OFF)[colnames(OFF) == "YWAR"] <- "YEWA"

rownames(TAX) <- TAX$Species_ID
TAX <- droplevels(TAX[colnames(OFF),])
#TAX$keep <- TRUE
#TAX$keep[TAX$Order %in% c("Gaviiformes","Pelecaniformes","Podicipediformes",
#"Gruiformes","Galliformes","Anseriformes","Charadriiformes")] <- FALSE
#TAX <- TAX[TAX$keep,]

SPP <- colnames(OFF)

YY <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)
YY <- YY[,SPP]
YYSS <- Xtab(ABUND ~ SS + SPECIES, PCTBL)
YYSS[YYSS > 1] <- 1
YYSS <- YYSS[,SPP]
rm(PCTBL)

PKEY <- PKEY[,c("PKEY", "SS", "PCODE", "METHOD", "SITE", "ROUND", "YEAR",
    "ROAD")]

SS$ALL_HAB_OK <- ifelse(rowSums(is.na(SS[,c("HAB_LCC1", "HAB_LCC2", "HAB_EOSD1", 
    "HAB_EOSD2", "HAB_NALC2", "HAB_NALC1")])) == 0, 0L, 1L)

with(SS, table(ROADCLASS, PAVSTATUS))
with(SS, table(NBRLANES, PAVSTATUS))
with(SS, table(ROADCLASS, NBRLANES))
SS$ROAD_OK <- ifelse(SS$NBRLANES > 2, 0L, 1L)
SS$d2road <- NULL
SS$ROADCLASS <- NULL
SS$NBRLANES <- NULL
SS$PAVSTATUS <- NULL

hist(SS$BEADTotalL)
hist(SS$BEADtotPol)
## linear features
SS$LIN <- log(SS$BEADTotalL + 1)
SS$BEADTotalL <- NULL
## Polygon is fine as it is (0-1)
SS$POL <- SS$BEADtotPol
SS$BEADtotPol <- NULL

## slope, cti, elev
SS$SLP <- sqrt(SS$slp90)
SS$SLP2 <- SS$SLP^2
SS$CTI <- (SS$cti90 - 8) / 4
SS$CTI2 <- SS$CTI^2
SS$ELV <- SS$elv90 / 1000
SS$ELV2 <- SS$ELV^2
SS$slp90 <- NULL
SS$cti90 <- NULL
SS$elv90 <- NULL

## height
SS$HGT <- SS$HEIGHTSIMARD / 50
SS$HGT2 <- SS$HGT^2
SS$HEIGHTSIMARD <- NULL

SS$TR3 <- SS$TREE3
SS$TREE3 <- NULL

## climate
SS$CMIJJA <- (SS$CMIJJA - 0) / 50
SS$CMI <- (SS$CMI - 0) / 50
SS$TD <- (SS$TD - 300) / 100
SS$DD0 <- (SS$DD01 - 1000) / 1000
SS$DD5 <- (SS$DD51 - 1600) / 1000
SS$EMT <- (SS$EMT + 400) / 100
SS$MSP <- (SS$MSP - 400) / 200

SS$DD02 <- SS$DD0^2
SS$DD52 <- SS$DD5^2

SS$CMD <- NULL
#CMI
#CMIJJA
SS$DD01 <- NULL
SS$DD51 <- NULL
#EMT
SS$FFP <- NULL
SS$ID <- NULL
SS$MAP <- NULL
SS$MAT <- NULL
SS$MCMT <- NULL
#MSP
SS$MWMT <- NULL
SS$NFFD <- NULL
SS$PAS <- NULL
SS$PET <- NULL
SS$PPT_sm <- NULL
SS$PPT_wt <- NULL
#TD

SS$FIRE_HA <- NULL
SS$YearLoss[SS$YearLoss == 0] <- NA

SS$YearFire[is.na(SS$YearFire)] <- 2015 - 200
SS$YearLoss[is.na(SS$YearLoss)] <- 2015 - 200

XYSS <- as.matrix(SS[,c("X","Y","Xcl","Ycl")])
rownames(XYSS) <- SS$SS
ii <- intersect(rownames(YYSS), rownames(XYSS))
YYSS <- YYSS[ii,]
XYSS <- XYSS[ii,]
rm(ii)

DAT <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL

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


## decid + mixed
DAT$isDM_LCC <- ifelse(DAT$HAB_LCC2 %in% c("Decid", "Mixed"), 1L, 0L)
DAT$isDM_EOSD <- ifelse(DAT$HAB_EOSD2 %in% c("Decid", "Mixed"), 1L, 0L)
DAT$isDM_NALC <- ifelse(DAT$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
DAT$isNF_LCC <- ifelse(DAT$HAB_LCC2 %in% 
    c("Agr", "Barren", "Burn", "Devel", "Grass", "Wet"), 1L, 0L)
DAT$isNF_EOSD <- ifelse(DAT$HAB_EOSD2 %in% 
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet"), 1L, 0L)
DAT$isNF_NALC <- ifelse(DAT$HAB_NALC2 %in% 
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet"), 1L, 0L)

DAT$HSH <- 0 # placeholder
DAT$HSH2 <- 0
DAT$HAB <- 0
DAT$isDM <- 0
DAT$isNF <- 0

## see above for filtering
DAT <- DAT[DAT$YEAR >= 1997,]

## year effect
#plot(table(DAT$YEAR))
DAT$YR5 <- 0
DAT$YR5[DAT$YEAR > 2000] <- 1
DAT$YR5[DAT$YEAR > 2004] <- 2
DAT$YR5[DAT$YEAR > 2009] <- 3
table(DAT$YEAR,DAT$YR5)
table(DAT$YR5)

DAT$YR5 <- as.factor(DAT$YR5)
DAT$YR <- DAT$YEAR - 1997

DAT$BCR[DAT$BCR == "0"] <- NA
DAT <- droplevels(DAT[!is.na(DAT$BCR),])
table(DAT$JURS, DAT$BCR)

DAT$BCR_JURS <- interaction(DAT$BCR, DAT$JURS, drop=TRUE, sep="_")

levels(DAT$BCR_JURS)[grepl("_AK", levels(DAT$BCR_JURS))] <- "AK"
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("9_BC","9_ID","9_WA",
    "5_BC","5_WA","10_AB","10_BC","10_ID","10_MT","10_WA")] <- "Mtn"
#levels(DAT$BCR_JURS)[grepl("11_", levels(DAT$BCR_JURS))] <- "Pra" # Prairies
#levels(DAT$BCR_JURS)[grepl("17_", levels(DAT$BCR_JURS))] <- "Pra"
#levels(DAT$BCR_JURS)[grepl("22_", levels(DAT$BCR_JURS))] <- "Pra"
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("11_AB", "11_MB", "11_MT", "11_SK",
    "17_MT", "17_ND", "17_SD", "11_ND", "11_SD")] <- "Pra_west" # Prairies
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("11_MN", 
    "22_IL", "22_IN", "22_MI", "22_MN", "22_OH", "22_WI")] <- "Pra_east"


#levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("4_BC",
#    "4_NT","4_YK","6_YK","6_NT")] <- "4+6_YK+NT"
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("4_BC","4_NT","4_YK")] <- "4_all"

levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("6_AB","6_BC","6_SK","6_MB",
    "6_YK","6_NT")] <- "6_all"

## keep northern points in
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("3_NT","7_MB","7_NT","3_NU","7_NU")] <- "3+7_west"
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("7_ON","7_QC","7_NL")] <- "3+7_east"

levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("8_MB","8_SK")] <- "8_west"
levels(DAT$BCR_JURS)[levels(DAT$BCR_JURS) %in% c("8_NL","8_ON","8_QC")] <- "8_east"

levels(DAT$BCR_JURS)[grepl("12_", levels(DAT$BCR_JURS))] <- "Grl" # Great Lakes
levels(DAT$BCR_JURS)[grepl("13_", levels(DAT$BCR_JURS))] <- "Grl"
levels(DAT$BCR_JURS)[grepl("23_", levels(DAT$BCR_JURS))] <- "Grl"

levels(DAT$BCR_JURS)[grepl("14_", levels(DAT$BCR_JURS))] <- "Mar" # Maritimes
levels(DAT$BCR_JURS)[grepl("30_", levels(DAT$BCR_JURS))] <- "Mar"

levels(DAT$BCR_JURS)[grepl("24_", levels(DAT$BCR_JURS))] <- "SEus" # SE US
levels(DAT$BCR_JURS)[grepl("26_", levels(DAT$BCR_JURS))] <- "SEus"
levels(DAT$BCR_JURS)[grepl("28_", levels(DAT$BCR_JURS))] <- "SEus"
levels(DAT$BCR_JURS)[grepl("29_", levels(DAT$BCR_JURS))] <- "SEus"

sort(levels(DAT$BCR_JURS))
table(DAT$BCR_JURS)
table(DAT$BCR_JURS, DAT$YR5)

## regions for trend
DAT$REG <- DAT$BCR_JURS
levels(DAT$REG)[levels(DAT$REG) %in% c("Mtn", "AK")] <- "Coast"
levels(DAT$REG)[levels(DAT$REG) %in% c("4_all", "3+7_west", "3+7_east")] <- "North"
levels(DAT$REG)[levels(DAT$REG) %in% c("Pra_west", "Pra_east")] <- "South"
levels(DAT$REG)[levels(DAT$REG) %in% c("8_west","6_all")] <- "West"
levels(DAT$REG)[levels(DAT$REG) %in% c("Mar","SEus","8_east","Grl")] <- "East"
table(DAT$BCR_JURS, DAT$REG)
table(DAT$REG)
DAT$REG <- relevel(DAT$REG, "West")

DAT$EW <- DAT$BCR_JURS
levels(DAT$EW)[levels(DAT$EW) %in% c("Mtn", "AK", "4_all", "3+7_west", 
    "Pra_west", "8_west","6_all")] <- "W"
levels(DAT$EW)[levels(DAT$EW) %in% c("3+7_east", "Pra_east",
    "Mar","SEus","8_east","Grl")] <- "E"
table(DAT$BCR_JURS, DAT$EW)
table(DAT$EW)
DAT$EW <- relevel(DAT$EW, "W")

library(RColorBrewer)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$REG, pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=brewer.pal(12, "Set3")[DAT$BCR_JURS], pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=ifelse(as.integer(DAT$BCR_JURS)==8,2,1), pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$EW, pch=19, cex=0.2)

## resampling blocks
DAT$bootg <- interaction(DAT$BCR_JURS, DAT$YR5, drop=TRUE)

## exclude same year revisits
#DAT$ROUND[is.na(DAT$ROUND)] <- 1
#DAT <- droplevels(DAT[DAT$ROUND == 1,])
DAT$SS_YR <- interaction(DAT$SS, DAT$YEAR, drop=TRUE)
table(table(DAT$SS_YR))
dup <- unique(DAT$SS_YR[duplicated(DAT$SS_YR)])
tmp1 <- DAT[DAT$SS_YR %in% dup, c("PKEY","SS_YR")]
rownames(tmp1) <- tmp1$PKEY
dim(tmp1)
tmp2 <- nonDuplicated(tmp1[sample.int(nrow(tmp1)),], SS_YR)
rownames(tmp2) <- tmp2$PKEY
dim(tmp2)
dim(tmp1) - dim(tmp2)
tmp1 <- tmp1[setdiff(rownames(tmp1), rownames(tmp2)),]
dim(tmp1)
dim(tmp1) + dim(tmp2)
DAT <- DAT[!(DAT$PKEY %in% rownames(tmp1)),]
DAT$SS_YR <- NULL

## unit of resampling is site x year to maintain temporal distribution
DAT$SITE <- as.character(DAT$SITE)
DAT$SITE[is.na(DAT$SITE)] <- paste0("gr.", as.character(DAT$gridcode[is.na(DAT$SITE)]))
DAT$SITE <- as.factor(DAT$SITE)

DAT$SITE_YR <- interaction(DAT$SITE, DAT$YEAR, drop=TRUE, sep="_")

DAT$ROAD_OK[is.na(DAT$ROAD_OK)] <- 1L

DAT$ARU <- ifelse(DAT$PCODE %in% levels(DAT$PCODE)[grep("EMC", levels(DAT$PCODE))],
    1L, 0L)

stopifnot(all(!is.na(DAT$PKEY)))

## placeholder for distance to nearest detection
#DAT$ND2 <- 0

summary(DAT$YSD)
summary(DAT$YSF)
summary(DAT$YSL)
AGEMAX <- 50
DAT$YSD <- pmax(0, 1 - (DAT$YSD / AGEMAX))
DAT$YSF <- pmax(0, 1 - (DAT$YSF / AGEMAX))
DAT$YSL <- pmax(0, 1 - (DAT$YSL / AGEMAX))
plot(DAT$YSD ~ DAT0$YSD)

if (FALSE) {
YS <- 0:100
YSD <- pmax(0, 1 - (YS / 50))
plot(YS, YSD, xlab="Years since disturbance", ylab="YSD covariate", type="l",
col=4, lwd=2)

}

source("~/repos/bamanalytics/R/analysis_mods.R")
source("~/repos/detect/R/hbootindex.R")
hbootindex2 <- hbootindex

nmin <- 25
B <- 239

## check max counts
apply(as.matrix(YY), 2, max)

## ---------------------------------- data packages below -----------
## SEXT: "can", "nam" # spatial extent, (canb=canadian boreal ~ eosd)
## TEXT: "gfw", "fre" # temporal extent, gfw=2001-2013, fire=1997-2014
## LCTU: "nlc", "lcc", "eos" # land cover to use

source("~/repos/bragging/R/glm_skeleton.R")
DAT0 <- DAT
TAX0 <- TAX
YY0 <- YY
OFF0 <- OFF
#Extra <- c("gridcode", "SS", "SITE_YR", "bootg", "X", "Y", "Xcl", "Ycl")
Extra <- c("SS", "SITE_YR", "X", "Y", "Xcl", "Ycl")
Save <- c("DAT", "YY", "OFF", "TAX", "mods", "BB") #, "HSH")
Date <- "2015-09-02"

if (FALSE) { # FIRE STUFF FOR NICOLE BEGIN

    Save <- c("DAT", "YY", "OFF", "TAX", "BB")
    Date <- "2016-01-14"

    rm(DAT, TAX, YY, OFF)
    gc()

    keep <- rep(TRUE, nrow(DAT0))
    keep[DAT0$ROAD_OK < 1] <- FALSE
    keep[is.na(DAT0$HAB_LCC2)] <- FALSE
    keep[DAT0$COUNTRY != "CAN"] <- FALSE

    DAT1 <- droplevels(DAT0[keep,])
    rownames(DAT1) <- DAT1$PKEY
    DAT1 <- droplevels(DAT1[intersect(rownames(DAT1), rownames(YY0)),])
    dim(DAT1)
    data.frame(n=colSums(is.na(DAT1)))
    
    DAT1 <- DAT1[,c("SS","PKEY","PCODE","METHOD","SITE","ROUND","YEAR","ROAD",
        "X","Y","JURS","BCR", "bootg")]
    ## need original SS here to be used
    DAT1$LCC05 <- SS$LCC05_PT[match(DAT1$SS, SS$SS)]
    DAT1$YearFire <- SS$YearFire[match(DAT1$SS, SS$SS)]
    DAT1$FIRE_HA <- SS$FIRE_HA[match(DAT1$SS, SS$SS)]

    set.seed(1234)
    ## make sure that time intervals are considered as blocks
    ## keep out 10% of the data for validation
    id2 <- list()
    for (l in levels(DAT1$bootg)) {
        sset <- which(DAT1$bootg == l)
        id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
    }
    KEEP_ID <- unname(unlist(id2))
    HOLDOUT_ID <- setdiff(seq_len(nrow(DAT1)), KEEP_ID)

    DAT1k <- DAT1[KEEP_ID,]
    DAT1 <- DAT1[c(KEEP_ID, HOLDOUT_ID),]
    DAT1k$SITE <- droplevels(DAT1k$SITE)
    BB1 <- hbootindex2(DAT1k$SITE, DAT1k$bootg, B=B)

    DAT <- DAT1
    pk <- rownames(DAT)
    YY <- YY0[pk,]
    YY <- YY[,colSums(YY) >= nmin]
    OFF <- OFF0[pk,colnames(YY)]
    TAX <- TAX0[colnames(YY),]
    BB <- BB1

dim(DAT)
dim(YY)
dim(OFF)
dim(TAX)
dim(BB)

    save(list = Save,
        file=file.path(ROOT, "out", "data", 
        paste0("pack_LCC05-fire-Canada_toNicole_", Date, ".Rdata")))

} # FIRE STUFF FOR NICOLE END

if (FALSE) { # ALL DATA FOR DIANA BEGIN
    DAT <- DAT0
    rownames(DAT) <- DAT$PKEY
    pk <- intersect(rownames(DAT), rownames(YY0))
    DAT <- DAT[pk,]
    YY <- YY0[pk,]
    YY <- YY[,colSums(YY) >= 0]
    OFF <- OFF0[pk,colnames(YY)]
    TAX <- TAX0[colnames(YY),]
    dim(DAT)
    dim(YY)
    dim(TAX)
    
    save(DAT, YY, TAX, OFF,
        file=file.path(ROOT, "out", "data", 
        paste0("pack_BAMBBS_toDiana_2016feb.Rdata")))

} # ALL DATA FOR DIANA END

for (TEXT in c("gfw", "fre")) {

rm(DAT, TAX, YY, OFF)
gc()
modsx <- if (TEXT == "gfw")
    mods_gfw else mods_fire

keep <- rep(TRUE, nrow(DAT0))
keep[DAT0$ROAD_OK < 1] <- FALSE
keep[is.na(DAT0$TREE)] <- FALSE
keep[is.na(DAT0$CMI)] <- FALSE
keep[is.na(DAT0$HGT)] <- FALSE
keep[is.na(DAT0$SLP)] <- FALSE
keep[is.na(DAT0$HAB_NALC2)] <- FALSE
keep2 <- keep
keep[is.na(DAT0$HAB_LCC2)] <- FALSE
keep[is.na(DAT0$HAB_EOSD2)] <- FALSE
keep[DAT0$COUNTRY != "CAN"] <- FALSE
keep[DAT0$REG == "South"] <- FALSE
#keep2[DAT0$YEAR < 1997] <- FALSE
if (TEXT == "gfw") {
    keep[DAT0$YEAR < 2001] <- FALSE # GFW years
    keep[DAT0$YEAR > 2013] <- FALSE # GFW years
    keep2[DAT0$YEAR < 2001] <- FALSE # GFW years
    keep2[DAT0$YEAR > 2013] <- FALSE # GFW years
}

#save(YYSS, XYSS,
#    file=file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))

DAT1 <- droplevels(DAT0[keep,])
rownames(DAT1) <- DAT1$PKEY
DAT1 <- droplevels(DAT1[intersect(rownames(DAT1), rownames(YY0)),])
dim(DAT1)
data.frame(n=colSums(is.na(DAT1)))
## data specific 0 year and decadal change
DAT1$YR <- (DAT1$YR - min(DAT1$YR)) / 10
table(DAT1$YR)

set.seed(1234)
## make sure that time intervals are considered as blocks
## keep out 10% of the data for validation
id2 <- list()
for (l in levels(DAT1$bootg)) {
    sset <- which(DAT1$bootg == l)
    id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
}
KEEP_ID <- unname(unlist(id2))
HOLDOUT_ID <- setdiff(seq_len(nrow(DAT1)), KEEP_ID)

DAT1k <- DAT1[KEEP_ID,]
DAT1 <- DAT1[c(KEEP_ID, HOLDOUT_ID),]
DAT1k$SITE <- droplevels(DAT1k$SITE)
BB1 <- hbootindex2(DAT1k$SITE, DAT1k$bootg, B=B)

DAT1_LCC <- DAT1
DAT1_LCC$HAB <- DAT1$HAB_LCC2
DAT1_LCC$isDM <- DAT1_LCC$isDM_LCC
DAT1_LCC$isNF <- DAT1_LCC$isNF_LCC
DAT1_LCC <- DAT1_LCC[,!grepl("_NALC", colnames(DAT1_LCC))]
DAT1_LCC <- DAT1_LCC[,!grepl("_EOSD", colnames(DAT1_LCC))]
colnames(DAT1_LCC) <- gsub("_LCC_", "_", colnames(DAT1_LCC))
DAT1_LCC <- DAT1_LCC[,!grepl("_LCC", colnames(DAT1_LCC))]

DAT1_EOSD <- DAT1
DAT1_EOSD$HAB <- DAT1$HAB_EOSD2
DAT1_EOSD$isDM <- DAT1_EOSD$isDM_EOSD
DAT1_EOSD$isNF <- DAT1_EOSD$isNF_EOSD
DAT1_EOSD <- DAT1_EOSD[,!grepl("_NALC", colnames(DAT1_EOSD))]
DAT1_EOSD <- DAT1_EOSD[,!grepl("_LCC", colnames(DAT1_EOSD))]
colnames(DAT1_EOSD) <- gsub("_EOSD_", "_", colnames(DAT1_EOSD))
DAT1_EOSD <- DAT1_EOSD[,!grepl("_EOSD", colnames(DAT1_EOSD))]

DAT1_NALC <- DAT1
DAT1_NALC$HAB <- DAT1$HAB_NALC2
DAT1_NALC$isDM <- DAT1_NALC$isDM_NALC
DAT1_NALC$isNF <- DAT1_NALC$isNF_NALC
DAT1_NALC <- DAT1_NALC[,!grepl("_LCC", colnames(DAT1_NALC))]
DAT1_NALC <- DAT1_NALC[,!grepl("_EOSD", colnames(DAT1_NALC))]
colnames(DAT1_NALC) <- gsub("_NALC_", "_", colnames(DAT1_NALC))
DAT1_NALC <- DAT1_NALC[,!grepl("_NALC", colnames(DAT1_NALC))]

DAT2 <- droplevels(DAT0[keep2,])
DAT2$HAB <- DAT2$HAB_NALC2
DAT2$isDM <- DAT2$isDM_NALC
DAT2$isNF <- DAT2$isNF_NALC
DAT2 <- DAT2[,!grepl("_LCC", colnames(DAT2))]
DAT2 <- DAT2[,!grepl("_EOSD", colnames(DAT2))]
colnames(DAT2) <- gsub("_NALC_", "_", colnames(DAT2))
DAT2 <- DAT2[,!grepl("_NALC", colnames(DAT2))]
data.frame(n=colSums(is.na(DAT2)))
rownames(DAT2) <- DAT2$PKEY
DAT2 <- droplevels(DAT2[intersect(rownames(DAT2), rownames(YY0)),])
dim(DAT2)
## data specific 0 year and decadal change
DAT2$YR <- (DAT2$YR - min(DAT2$YR)) / 10
table(DAT2$YR)

set.seed(1234)
## make sure that time intervals are considered as blocks
## keep out 10% of the data for validation
id2 <- list()
for (l in levels(DAT2$bootg)) {
    sset <- which(DAT2$bootg == l)
    id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
}
KEEP_ID <- unname(unlist(id2))
HOLDOUT_ID <- setdiff(seq_len(nrow(DAT2)), KEEP_ID)

DAT2k <- DAT2[KEEP_ID,]
DAT2 <- DAT2[c(KEEP_ID, HOLDOUT_ID),]
DAT2k$SITE <- droplevels(DAT2k$SITE)
BB2 <- hbootindex2(DAT2k$SITE, DAT2k$bootg, B=B)

## ---------------- save output -----------------------------------------

DAT <- DAT1_LCC
pk <- rownames(DAT)
YY <- YY0[pk,]
YY <- YY[,colSums(YY) >= nmin]
OFF <- OFF0[pk,colnames(YY)]
TAX <- TAX0[colnames(YY),]
mods <- modsx
HSH <- as.matrix(DAT[,grep("GRID4_", colnames(DAT))])
colnames(HSH) <- gsub("GRID4_", "", colnames(HSH))
BB <- BB1
DAT <- DAT[,c(Extra, getTerms(mods, "list"))]
cat("\n---\t", TEXT, "can", "lcc", nrow(DAT), "\n")
save(list = Save,
    file=file.path(ROOT, "out", "data", 
    paste0("pack_", TEXT, "_can_lcc_", Date, ".Rdata")))

DAT <- DAT1_EOSD
pk <- rownames(DAT)
YY <- YY0[pk,]
YY <- YY[,colSums(YY) >= nmin]
OFF <- OFF0[pk,colnames(YY)]
TAX <- TAX0[colnames(YY),]
mods <- modsx
HSH <- as.matrix(DAT[,grep("GRID4_", colnames(DAT))])
colnames(HSH) <- gsub("GRID4_", "", colnames(HSH))
BB <- BB1
DAT <- DAT[,c(Extra, getTerms(mods, "list"))]
cat("\n---\t", TEXT, "can", "eos", nrow(DAT), "\n")
save(list = Save,
    file=file.path(ROOT, "out", "data", 
    paste0("pack_", TEXT, "_can_eos_", Date, ".Rdata")))

DAT <- DAT1_NALC
pk <- rownames(DAT)
YY <- YY0[pk,]
YY <- YY[,colSums(YY) >= nmin]
OFF <- OFF0[pk,colnames(YY)]
TAX <- TAX0[colnames(YY),]
mods <- modsx
HSH <- as.matrix(DAT[,grep("GRID4_", colnames(DAT))])
colnames(HSH) <- gsub("GRID4_", "", colnames(HSH))
BB <- BB1
DAT <- DAT[,c(Extra, getTerms(mods, "list"))]
cat("\n---\t", TEXT, "can", "nlc", nrow(DAT), "\n")
save(list = Save,
    file=file.path(ROOT, "out", "data", 
    paste0("pack_", TEXT, "_can_nlc_", Date, ".Rdata")))

DAT <- DAT2
pk <- rownames(DAT)
YY <- YY0[pk,]
YY <- YY[,colSums(YY) >= nmin]
OFF <- OFF0[pk,colnames(YY)]
TAX <- TAX0[colnames(YY),]
mods <- modsx
HSH <- as.matrix(DAT[,grep("GRID4_", colnames(DAT))])
colnames(HSH) <- gsub("GRID4_", "", colnames(HSH))
BB <- BB2
DAT <- DAT[,c(Extra, getTerms(mods, "list"))]
cat("\n---\t", TEXT, "nam", "nlc", nrow(DAT), "\n")
save(list = Save,
    file=file.path(ROOT, "out", "data", 
    paste0("pack_", TEXT, "_nam_nlc_", Date, ".Rdata")))
#plot(DAT[,c("X","Y")], col=DAT$REG, pch=19, cex=0.2)

}
