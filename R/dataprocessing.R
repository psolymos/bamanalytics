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
    "NL","NS","NT","NU","ON","PEI","QC","SK","YK"), "CAN", "USA")
save(SS01, file=file.path(ROOT2, "out", "SS-regions-and-xy.Rdata"))

## Tree proportions
SS02 <- sqlFetch(con, "dbo_TREE_BBSBAM_V4_tbl")
SS02 <- nonDuplicated(SS02, SS, TRUE)
SS02 <- SS02[rownames(SS01),]

## TREE update 2016-12-01
tr <- read.csv("e:/peter/bam/Apr2016/tree_update/BAM_BBS_2015_XY_AVHRRTREE_nov30_2016.csv")
SS02$TREE_WRONG <- SS02$TREE
SS02$TREE <- tr$VCFAVHRR[match(SS02$SS, tr$SS)]

SS02$TREE[SS02$TREE > 100] <- NA
SS02$TREE[SS02$TREE < 0] <- NA
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
## Reclass LCC05
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
SS03$LCC05 <- SS03$LCC05_PT
SS03$LCC05_PT[SS03$LCC05_PT < 1 | SS03$LCC05_PT > 39] <- NA
SS03$LCC05_PT[SS01$COUNTRY == "USA"] <- NA
SS03$HAB_LCC1 <- ltlcc$BAMLCC05V2_label1[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC1 <- relevel(SS03$HAB_LCC1, "ConifDense")
SS03$HAB_LCC2 <- relevel(SS03$HAB_LCC2, "Conif")

if (FALSE) {
## for Nicole
SS_Nicole <- data.frame(
    PCODE=SS01$PCODE,
    SS=SS01$SS,
    X=SS01$X_GEONAD83,
    Y=SS01$Y_GEONAD83,
    JURS=SS01$JURSALPHA,
    COUNTRY=SS01$COUNTRY,
    BCR=as.factor(SS01$BCR),
    SS03[,c("LCC05","HAB_LCC1","HAB_LCC2")])

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
ii <- SS03$HAB_NALC1 %in% c("Conif", "Decid", "Mixed", "Wet")
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

## Grid ID 4x4 km
SS_grid <- read.csv(file.path(ROOT, "BAMBBS_Gridcode.csv"))
rownames(SS_grid) <- SS_grid$SS
compare_sets(rownames(SS01), rownames(SS_grid))
SS_grid <- SS_grid[rownames(SS01),"gridcode",drop=FALSE]
levels(SS_grid$gridcode) <- gsub(",", "", levels(SS_grid$gridcode))

## Road: dist to, class, #lanes, surface
SS_road <- sqlFetch(con, "dbo_BAMBBS_2015_NearDistanceRoadJoin1000M")
rownames(SS_road) <- SS_road$SS
compare_sets(rownames(SS01), rownames(SS_road))
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
compare_sets(rownames(SS01), rownames(SS_fire))
SS_fire <- SS_fire[rownames(SS01),]
SS_fire <- SS_fire[,c("Year","SIZE_HA")]
colnames(SS_fire) <- c("YearFire","FIRE_HA")

## Terrain: slope, twi, elev
SS_terr <- sqlFetch(con, "dbo_BBSBAM_2015_TERRAIN90")
rownames(SS_terr) <- SS_terr$SS
compare_sets(rownames(SS01), rownames(SS_terr))
SS_terr <- SS_terr[rownames(SS01),]
t(table(is.na(SS_terr$cti90), SS01$PCODE)) # mostly affects BBS in some states
SS_terr <- SS_terr[,c("slp90","cti90","elv90")]

## Climate variables from DS and NALCMS 4x4 level
SS_clim <- sqlFetch(con, "dbo_BBSBAM_2015__CLIMLU")
rownames(SS_clim) <- SS_clim$SS
compare_sets(rownames(SS01), rownames(SS_clim))
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
compare_sets(rownames(SS01), rownames(SS_LCC4x4))
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
compare_sets(rownames(SS01), rownames(SS_EOSD4x4))
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
compare_sets(rownames(SS01), rownames(SS_height))
SS_height <- SS_height[rownames(SS01),]
SS_height <- SS_height[,"HEIGHTSIMARD",drop=FALSE]

if (FALSE) {
## Nature Serve range: 3 spp (Can clipped range used 0/1)
SS_nserv <- sqlFetch(con, "dbo_BBSBAM_SARSPPLOCATIONSRange")
SS_nserv <- nonDuplicated(SS_nserv, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_nserv))
SS_nserv <- SS_nserv[rownames(SS01),]
SS_nserv <- SS_nserv[,c("CAWAINOUT","OSFLINOUT","CONIINOUT")]
}

## GFW yearly loss intersections and 1st year of loss
#SS_gfw <- sqlFetch(con, "dbo_BAMBBS_GFWLossYear")
SS_gfw <- read.csv(file.path(ROOT, "GFWLossYear.csv"))
SS_gfw <- nonDuplicated(SS_gfw, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_gfw))
levels(SS_gfw$YearLoss) <- gsub(",", "", levels(SS_gfw$YearLoss))
SS_gfw$YearLoss <- as.integer(as.character(SS_gfw$YearLoss))
SS_gfw <- SS_gfw[rownames(SS01),]
SS_gfw <- SS_gfw[,"YearLoss",drop=FALSE]

## Pasher disturbance
SS_pash <- read.csv(file.path(ROOT, "bambbs2015beadandpasher.csv"))
SS_pash <- nonDuplicated(SS_pash, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_pash))
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
    SS03[,c("LCC05","HAB_LCC1","HAB_LCC2")],
    SS03[,c("HAB_NALC2", "HAB_NALC1")],
    SS_grid,
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
compare_sets(PCODE$Method, PKEY$METHOD)
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
table(PCTBL$BEH, PCTBL$PCODE=="BBS")
keep <- rep(TRUE, nrow(PCTBL))
keep[!(PCTBL$BEH %in% c("1","6","11"))] <- FALSE
## this is fake, but there is no other option until a fix
#keep[is.na(PCTBL$BEH)] <- TRUE # dont know what this in -- FIXED in BBS_V4

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
compare_sets(PCTBL$SPECIES, TAX$Species_ID)
setdiff(PCTBL$SPECIES, TAX$Species_ID)
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "YWAR"] <- "YEWA"
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "SCJU"] <- "DEJU" # change SCJU to DEJU
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "MYWA"] <- "YRWA" # change MYWA to YRWA
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "COSN"] <- "WISN" # change COSN to WISN

levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "YWAR"] <- "YEWA"
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "SCJU"] <- "DEJU" # change SCJU to DEJU
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "MYWA"] <- "YRWA" # change MYWA to YRWA
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "COSN"] <- "WISN" # change COSN to WISN
PCTBL$SPECIES <- droplevels(PCTBL$SPECIES)
setdiff(PCTBL$SPECIES, TAX$Species_ID)


PCTBL$SPECIES_ALL <- PCTBL$SPECIES
sspp <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")
levels(PCTBL$SPECIES)[!(levels(PCTBL$SPECIES) %in% sspp$Species_ID[sspp$Singing_birds])] <- "NONE"

## Excluding columns
PCTBL$DURATION <- NULL
PCTBL$DISTANCE <- NULL
PCTBL$BEH <- NULL
PCTBL$dur <- droplevels(PCTBL$dur)
PCTBL$dis <- droplevels(PCTBL$dis)

compare_sets(SS$SS, PKEY$SS)
compare_sets(SS$SS, PCTBL$SS)
compare_sets(PKEY$PKEY, PCTBL$PKEY)

save(SS, PKEY, PCTBL, TAX,
    file=file.path(ROOT2, "out",
    #paste0("data_package_2016-04-18.Rdata")))
#    paste0("data_package_2016-07-05.Rdata")))
    paste0("data_package_2016-12-01.Rdata")))


#### Calculate the offsets (optional) -----------------------------------------
## QPADv3

## Define root folder where data are stored
ROOT <- "e:/peter/bam/Apr2016"
## Load required packages
library(mefa4)
library(QPAD)
## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

load(file.path(ROOT, "out", paste0("data_package_2016-12-01.Rdata")))

load_BAM_QPAD(3)
getBAMversion()
sppp <- getBAMspecieslist()

offdat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS","JULIAN")],
    SS[match(PKEY$SS, rownames(SS)),c("TREE","TREE3","HAB_NALC1","HAB_NALC2", "SPRNG")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
offdat$DSLS <- (offdat$JULIAN - offdat$SPRNG) / 365
summary(offdat)

offdat$JDAY2 <- offdat$JDAY^2
offdat$TSSR2 <- offdat$TSSR^2
offdat$DSLS2 <- offdat$DSLS^2
offdat$LCC4 <- as.character(offdat$HAB_NALC2)
offdat$LCC4[offdat$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
offdat$LCC4[offdat$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
offdat$LCC4 <- factor(offdat$LCC4,
        c("DecidMixed", "Conif", "Open", "Wet"))
offdat$LCC2 <- as.character(offdat$LCC4)
offdat$LCC2[offdat$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
offdat$LCC2[offdat$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
offdat$LCC2 <- factor(offdat$LCC2, c("Forest", "OpenWet"))
table(offdat$LCC4, offdat$LCC2)
offdat$MAXDIS <- offdat$MAXDIS / 100

Xp <- cbind("(Intercept)"=1, as.matrix(offdat[,c("TSSR","JDAY","DSLS","TSSR2","JDAY2","DSLS2")]))
Xq <- cbind("(Intercept)"=1, TREE=offdat$TREE,
    LCC2OpenWet=ifelse(offdat$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(offdat$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(offdat$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(offdat$LCC4=="Wet", 1, 0))
#offdat$OKp <- rowSums(is.na(offdat[,c("TSSR","JDAY","DSLS")])) == 0
#offdat$OKq <- rowSums(is.na(offdat[,c("TREE","LCC4")])) == 0
#Xp <- model.matrix(~TSSR+TSSR2+JDAY+JDAY2+DSLS+DSLS2, offdat[offdat$OKp,])
#Xq <- model.matrix(~LCC2+LCC4+TREE, offdat[offdat$OKq,])

OFF <- matrix(NA, nrow(offdat), length(sppp))
rownames(OFF) <- offdat$PKEY
colnames(OFF) <- sppp

#spp <- "OVEN"
for (spp in sppp) {
cat(spp, "\n");flush.console()
p <- rep(NA, nrow(offdat))
A <- q <- p

## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
#vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0

p[!OKp] <- sra_fun(offdat$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(offdat$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * offdat$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[!OKq], cf0[2]))

phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(offdat$MAXDUR[OKp], phi1)
unlim <- ifelse(offdat$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * offdat$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(offdat$MAXDIS[OKq], tau1))

ii <- which(p == 0)
p[ii] <- sra_fun(offdat$MAXDUR[ii], cf0[1])

OFF[,spp] <- log(p) + log(A) + log(q)

}

(Ra <- apply(OFF, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,])) # BARS GCSP
which(!is.finite(Ra[2,]))

SPP <- sppp
save(OFF, SPP,
    file=file.path(ROOT, "out", "offsets-v3_2016-12-01.Rdata"))
offdat <- offdat[,c("PKEY","TSSR","JDAY","DSLS","TREE","LCC4","MAXDUR","MAXDIS")]
save(offdat,
    file=file.path(ROOT, "out", "offsets-v3data_2016-12-01.Rdata"))

######## These are the transformations #################


## Define root folder where data are stored
ROOT <- "e:/peter/bam/Apr2016"
## Load required packages
library(mefa4)
## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

load(file.path(ROOT, "out", paste0("data_package_2016-12-01.Rdata")))
load(file.path(ROOT, "out", "offsets-v3_2016-12-01.Rdata"))

TAX <- nonDuplicated(TAX, Species_ID, TRUE)
TAX <- droplevels(TAX[SPP,])
OFF <- OFF[,SPP]

YY <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)
YY <- YY[,SPP]
YYSS <- Xtab(ABUND ~ SS + SPECIES, PCTBL)
YYSS[YYSS > 1] <- 1
YYSS <- YYSS[,SPP]
rm(PCTBL)

PKEY <- PKEY[,c("PKEY", "SS", "PCODE", "METHOD", "SITE", "ROUND", "YEAR",
    "ROAD")]

with(SS, table(ROADCLASS, PAVSTATUS))
with(SS, table(NBRLANES, PAVSTATUS))
with(SS, table(ROADCLASS, NBRLANES))

SS$ROAD_OK <- SS$NBRLANES <= 2
SS$HAB_OK <- !is.na(SS$HAB_NALC2)

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
SS$HGTorig <- SS$HEIGHTSIMARD / 25
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

## lat/lon
SS$xlon <- (SS$Xcl - 285400) / 1500000
SS$xlat <- (SS$Ycl - 1320000) / 740000
SS$xlon2 <- SS$xlon^2
SS$xlat2 <- SS$xlat^2

SS$FIRE_HA <- NULL
SS$YearLoss[SS$YearLoss == 0] <- NA

SS$YearFire[is.na(SS$YearFire)] <- 2016 - 200
SS$YearLoss[is.na(SS$YearLoss)] <- 2016 - 200

XYSS <- as.matrix(SS[,c("X","Y","Xcl","Ycl")])
rownames(XYSS) <- SS$SS
ii <- intersect(rownames(YYSS), rownames(XYSS))
YYSS <- YYSS[ii,]
XYSS <- XYSS[ii,]
rm(ii)

DAT <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL
rownames(DAT) <- DAT$PKEY

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
DAT$isDM <- ifelse(DAT$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
DAT$isNF <- ifelse(DAT$HAB_NALC1 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub",
    "WetOpen", "DecidOpen", "ConifOpen", "MixedOpen"), 1L, 0L)
DAT$isDev <- ifelse(DAT$HAB_NALC2 %in% c("Agr", "Devel"), 1L, 0L)
DAT$isOpn <- DAT$isNF
DAT$isOpn[DAT$isDev == 1] <- 0
DAT$isWet <- ifelse(DAT$HAB_NALC2 %in% c("Wet"), 1L, 0L)
DAT$isDec <- ifelse(DAT$HAB_NALC2 %in% c("Decid"), 1L, 0L)
DAT$isMix <- ifelse(DAT$HAB_NALC2 %in% c("Mixed"), 1L, 0L)

## see above for filtering
DAT <- DAT[DAT$YEAR >= 1997,]

## year effect
#plot(table(DAT$YEAR))
#DAT$YR <- DAT$YEAR - 1997

#keep <- rep(TRUE, nrow(DAT))
#keep[DAT$BCR %in% c(1,5,11,17,22,24,26,28,29,30)] <- FALSE
#keep[DAT$BCR %in% c(9, 10)] <- FALSE
#keep[DAT$BCR %in% c(9, 10) & DAT$JURS %in% c("AB","BC")] <- TRUE
#keep[DAT$BOREALLOC != "OUT"] <- TRUE
#keep[is.na(DAT$BCR) | DAT$BCR == "0"] <- FALSE

keep <- !is.na(DAT$BOREALLOC) & DAT$BOREALLOC != "OUT"
keep[!is.null(DAT$BCR) & DAT$BCR %in% c("12", "13","14")] <- TRUE

table(DAT$BCR,keep)

DAT <- droplevels(DAT[keep,])

## fill-in BCR based on closes known value
length(which(DAT$BCR == "0" | is.na(DAT$BCR)))
ssDAT <- nonDuplicated(DAT, SS)
ssDAT$xBCR <- ssDAT$BCR
ii <- which(ssDAT$BCR == "0" | is.na(ssDAT$BCR))
for (i in ii) {
    d <- sqrt((ssDAT$X - ssDAT$X[i])^2 + (ssDAT$Y - ssDAT$Y[i])^2)
    d[ii] <- Inf
    ssDAT$xBCR[i] <- ssDAT$BCR[which.min(d)]
}
DAT$xBCR <- ssDAT$xBCR[match(DAT$SS, ssDAT$SS)]
length(which(DAT$BCR == "0" | is.na(DAT$BCR)))
length(which(DAT$xBCR == "0" | is.na(DAT$xBCR)))

table(DAT$JURS, DAT$xBCR)
DAT$Units <- factor(NA, c("AK","N","W","8E","Mt","HW", "At"))
DAT$Units[DAT$xBCR %in% c("2","4")] <- "AK"
DAT$Units[DAT$xBCR %in% c("3","7")] <- "N" # Arctic/Shield
DAT$Units[DAT$xBCR %in% c("5","9","10")] <- "Mt" # Rockies

DAT$Units[DAT$xBCR %in% c("6","11")] <- "W" # West
DAT$Units[DAT$xBCR %in% c("8") & DAT$JURS %in% c("MB","SK")] <- "W"

DAT$Units[DAT$xBCR %in% c("8") & !(DAT$JURS %in% c("MB","SK"))] <- "8E" # BCR8 east

DAT$Units[DAT$xBCR %in% c("12","13","23")] <- "HW" # Hardwood
DAT$Units[DAT$xBCR %in% c("14")] <- "At" # Atlantic

table(DAT$xBCR, DAT$Units, useNA="a")
table(DAT$Units, useNA="a")

ssDAT <- nonDuplicated(DAT, SS)
table(ssDAT$Units, useNA="a")


#library(RColorBrewer)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$REG, pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=brewer.pal(12, "Set3")[DAT$BCR_JURS], pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=ifelse(as.integer(DAT$BCR_JURS)==8,2,1), pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$EW, pch=19, cex=0.2)

## resampling blocks
DAT$YR5 <- 0
DAT$YR5[DAT$YEAR > 2000] <- 1
DAT$YR5[DAT$YEAR > 2004] <- 2
DAT$YR5[DAT$YEAR > 2009] <- 3
table(DAT$YEAR,DAT$YR5)
table(DAT$YR5)
DAT$bootg <- interaction(DAT$Units, DAT$YR5, drop=TRUE)

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

DAT$HAB <- DAT$HAB_NALC2
DAT$HABTR <- DAT$HAB_NALC1
DAT$HGT <- DAT$HGTorig
DAT$HGT[DAT$HAB %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- 0
DAT$HGT2 <- DAT$HGT^2
DAT$HGT05 <- sqrt(DAT$HGT)

DAT <- DAT[DAT$ARU == 0,]

data.frame(n=colSums(is.na(DAT)))


keep <- rep(TRUE, nrow(DAT))
keep[DAT$ROAD_OK < 1] <- FALSE
keep[is.na(DAT$HABTR)] <- FALSE
keep[is.na(DAT$CMI)] <- FALSE
keep[is.na(DAT$HGT)] <- FALSE
keep[is.na(DAT$SLP)] <- FALSE
keep[DAT$YEAR < 2001] <- FALSE # GFW years
keep[DAT$YEAR > 2013] <- FALSE # GFW years

data.frame(n=colSums(is.na(DAT[keep,])))

DAT <- droplevels(DAT[keep,])
DAT$YR <- DAT$YEAR - min(DAT$YEAR)

## check NAs, exclude ARU

library(detect)
source("~/repos/bamanalytics/R/analysis_mods.R")
source("~/repos/bragging/R/glm_skeleton.R")

nmin <- 25
B <- 239
Extra <- c("SS", "SITE_YR", "X", "Y", "Xcl", "Ycl", "Units", "xBCR", "JURS")
Save <- c("DAT", "YY", "OFF", "TAX", "mods", "BB")
Date <- "2016-12-01"

pk <- intersect(rownames(DAT), rownames(YY))
DAT <- DAT[pk,]

set.seed(1234)
## make sure that time intervals are considered as blocks
## keep out 10% of the data for validation
id2 <- list()
for (l in levels(DAT$bootg)) {
    sset <- which(DAT$bootg == l)
    id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
}
KEEP_ID <- unname(unlist(id2))
HOLDOUT_ID <- setdiff(seq_len(nrow(DAT)), KEEP_ID)

DATk <- DAT[KEEP_ID,]
DAT <- DAT[c(KEEP_ID, HOLDOUT_ID),]
DATk$SITE <- droplevels(DATk$SITE)
BB <- hbootindex(DATk$SITE, DATk$bootg, B=B)

#pk <- intersect(rownames(DAT), rownames(YY))
pk <- rownames(DAT)
#DAT <- DAT[pk,]
YY <- YY[pk,]
YY <- YY[,colSums(YY) >= nmin]
apply(as.matrix(YY), 2, max)
#apply(as.matrix(YY), 2, table)

OFF <- OFF[pk,colnames(YY)]
TAX <- TAX[colnames(YY),]

DAT$CMI2 <- DAT$CMI^2
DAT$CMIJJA2 <- DAT$CMIJJA^2

DAT <- DAT[,c(Extra, getTerms(mods, "list"))]

save(list = Save,
    file=file.path(ROOT, "out", "data", paste0("pack_", Date, ".Rdata")))



ROOT <- "e:/peter/bam/Apr2016"
load(file=file.path(ROOT, "out", "data", paste0("pack_2016-04-18.Rdata")))

library(rworldmap)
X0 <- nonDuplicated(DAT, SS)
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(X0[, c("X","Y")], pch=19, cex=0.25, col=as.integer(X0$Units)+1)

library(lattice)
densityplot(~ I(HGT*25) | HABTR, DAT[DAT$HGT>0,])

