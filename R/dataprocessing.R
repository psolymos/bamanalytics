##---
##title: "Data processing for nationam BAM analyses"
##author: "Peter Solymos"
##date: "May 7, 2015"
##output: 
##  pdf_document: 
##    toc: true
##    toc_depth: 2
##---

### Preliminaries

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"

## Load required packages
library(mefa4)
library(RODBC)
library(maptools)
library(pbapply)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

### Pulling in tables

## Define MS Access database connection
con <- odbcConnectAccess2007(file.path(ROOT, "BAM_BayneAccess_BAMBBScore.accdb"))

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

## Put together the main SS level object
SS <- data.frame(
    PCODE=SS01$PCODE,
    SS=SS01$SS,
    X=SS01$X_GEONAD83,
    Y=SS01$Y_GEONAD83,
    JURS=SS01$JURSALPHA,
    COUNTRY=SS01$COUNTRY,
    TZONE=SS01$TZONE_CODE,
    BOREALLOC=SS01$BOREALLOC,
    BCR=as.factor(SS01$BCR),
    TREE=SS02$TREE,
    TREE3=SS02$TREE3,
    SS03[,c("HAB_LCC1", "HAB_LCC2", "HAB_EOSD1", "HAB_EOSD2", 
        "HAB_NALC2", "HAB_NALC1", "LCC_combo")])

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
pkbbs <- sqlFetch(con, "dbo_PKEY_BBS_V3_2015")
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
MM <- ifelse(PKEY$MONTH < 10, paste0("0", PKEY$MONTH), as.character(PKEY$MONTH))
HH <- ifelse(PKEY$HOUR < 10, paste0("0", PKEY$HOUR), as.character(PKEY$HOUR))
mm <- ifelse(PKEY$MIN < 10, paste0("0", PKEY$MIN), as.character(PKEY$MIN))
DD <- with(PKEY, paste0(YEAR, "-", MM, "-", DAY, " ", HH, ":", mm, ":00"))
DD <- strptime(DD, "%Y-%m-%e %H:%M:%S")
## Julian day
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
PKEY$TSSR[PKEY$start_time > 12] <- NA ## after noon
summary(PKEY$TSSR)
summary(PKEY$start_time)

#### Calculate the offsets (optional)
if (FALSE) { # BEGIN offset calculations

offdat <- data.frame(PKEY[,c("PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS")],
    SS[match(PKEY$SS, rownames(SS)),c("LCC_combo","TREE")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
summary(offdat)

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
pcbbs <- sqlFetch(con, "dbo_PointCount_BBS_V3_2015")
colnames(pcbbs) <- toupper(colnames(pcbbs))
colnames(pcbbs)[colnames(pcbbs) == "SPECIES_ID"] <- "SPECIES"
colnames(pcbbs)[colnames(pcbbs) == "PERIOD"] <- "DURATION"
pcbbs$PCODE <- "BBS"
## Columns to keep
pccols <- c("PCODE","SS","PKEY","DURATION","DISTANCE",
    "SPECIES","ABUND","BEH")
## Close the database connection
close(con)

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
## Filtering behaviour
#sort(100 * table(PCTBL$BEH) / sum(table(PCTBL$BEH)))
## 1=Heard
## 11=no birds observed at station - added 2011
## 6=seen and heard
## Excluding non-aerial detections
keep <- rep(TRUE, nrow(PCTBL))
keep[!(PCTBL$BEH %in% c("1","6","11"))] <- FALSE
## Excluding >10 min intervals
## 10=10-20
## 11=0-20
## 3=before or after
## 8=unk
## 9=10-15
keep[PCTBL$DURATION %in% c(10,11,3,8,9)] <- FALSE
## Excluding unknown distance bands
keep[PCTBL$DISTANCE %in% c(4,5,9)] <- FALSE
## Excluding NA values
keep[is.na(PCTBL$dur)] <- FALSE
keep[is.na(PCTBL$dis)] <- FALSE
keep[is.na(PCTBL$ABUND)] <- FALSE
## Actual filtering
#PCTBL$keep <- keep
PCTBL <- PCTBL[keep,]

## Excluding/dropping species

PCTBL$SPECIES <- droplevels(PCTBL$SPECIES)
levels(PCTBL$SPECIES) <- toupper(levels(PCTBL$SPECIES))
compare.sets(PCTBL$SPECIES, TAX$Species_ID)
setdiff(PCTBL$SPECIES, TAX$Species_ID)
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "YWAR"] <- "YEWA"

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
dat <- data.frame(PKEY[,c("PCODE","PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS","METHOD",
    "DURMETH","DISMETH","ROAD")],
    SS[match(PKEY$SS, rownames(SS)),c("TREE","TREE3","LCC_combo","HAB_NALC1","HAB_NALC2")])
dat <- dat[dat$ROAD == 0,]
rownames(dat) <- dat$PKEY
ii <- intersect(dat$PKEY, levels(PCTBL$PKEY))
dat <- droplevels(dat[ii,])
summary(dat)
colSums(is.na(dat))
## sra and edr might have different NA patterns -- it is OK to exclude them later
#dat <- dat[rowSums(is.na(dat)) == 0,]
dat <- droplevels(dat)

dat2 <- with(PKEY_abmi, data.frame(
    PCODE="ABMI",
    PKEY=as.factor(Label),
    SS=as.factor(Label2),
    TSSR=TSSR,
    JDAY=JDAY,
    MAXDIS=Inf,
    MAXDUR=10,
    METHOD="ABMI:1",
    DURMETH="X",
    DISMETH="D",
    ROAD=0,
    TREE=NA,
    TREE3=NA,
    LCC_combo=NA,
    HAB_NALC1=NA,
    HAB_NALC2=NA))
rownames(dat2) <- dat2$PKEY

## besides `dat` we also need specific aoutput from `PCTBL`
## and also the methodology x interval lookup

compare.sets(DISMET$DISTANCECODE, PKEY$DISMETH)
compare.sets(DURMET$DURATIONCODE, PKEY$DURMETH)

PCTBL$DISMETH <- droplevels(PKEY$DISMETH[match(PCTBL$PKEY, PKEY$PKEY)])
PCTBL$DURMETH <- droplevels(PKEY$DURMETH[match(PCTBL$PKEY, PKEY$PKEY)])


## Oddities that should not happen:
PCTBL$dur[with(PCTBL, DURMETH=="A" & dur=="0-3")] <- "0-10"
PCTBL$dur[with(PCTBL, DURMETH=="B" & dur=="5-8")] <- "0-5"

PCTBL$dis[with(PCTBL, DISMETH=="B" & dis=="0-Inf")] <- "0-50" # best guess
PCTBL$dis[with(PCTBL, DISMETH=="C" & dis=="0-Inf")] <- "0-50" # best guess
PCTBL$dis[with(PCTBL, DISMETH=="F")] <- "0-100" # all kinds of weird stuff
PCTBL$dis[with(PCTBL, DISMETH=="I" & dis=="100-125")] <- "0-25"
PCTBL$dis[with(PCTBL, DISMETH=="I" & dis=="100-Inf")] <- "0-25" # best guess
#PCTBL$dis[with(PCTBL, DISMETH=="L" & dis=="150-Inf")] <- "100-150" # no >150
#PCTBL$dis[with(PCTBL, DISMETH=="T" & dis=="100-Inf")] <- "-" # no >100
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="0-50")] <- "40-50"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="100-150")] <- "125-150"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="100-Inf")] <- "150-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="U" & dis=="50-Inf")] <- "150-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="W" & dis=="150-Inf")] <- "100-Inf"
PCTBL$dis[with(PCTBL, DISMETH=="W" & dis=="100-125")] <- "100-Inf"

PCTBL$dis <- droplevels(PCTBL$dis)
PCTBL$dur <- droplevels(PCTBL$dur)

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
    DURMETH="ABMI"))
levels(pc2$dur) <- c("3.33","6.66","10")


## combine dat, dat2 and pc pc2
dat <- rbind(dat, dat2)
pc <- rbind(pc, pc2)

durmat <- as.matrix(Xtab(~ DURMETH + dur, pc))
durmat[durmat > 0] <- 1
dismat <- as.matrix(Xtab(~ DISMETH + dis, pc))
dismat[dismat > 0] <- 1
ltdur <- arrange.intervals(durmat)
ltdis <- arrange.intervals(dismat)
## divide by 100
ltdis$end <- ltdis$end / 100

save(dat, pc, ltdur, ltdis, TAX,
    file=file.path(ROOT, "out",
    paste0("new_offset_data_package_", Sys.Date(), ".Rdata")))

save(SS, PKEY, PCTBL, TAX,
    file=file.path(ROOT, "out",
    paste0("data_package_", Sys.Date(), ".Rdata")))
