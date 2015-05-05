ROOT <- "c:/bam/May2015"

library(mefa4)
library(RODBC)

con <- odbcConnectAccess2007(file.path(ROOT, "BAM_BayneAccess_BAMBBScore.accdb"))

#### SS level stuff for BAM+BBS combined
## time zone, BCR, jurisdiction, XY
SS01 <- sqlFetch(con, "dbo_BBSBAM_V4_XYTBL_ATTRIBUTES1")
SS01 <- nonDuplicated(SS01, SS, TRUE)
SS01$COUNTRY <- ifelse(SS01$JURSALPHA %in% c("AB","BC","MB","NB",
    "NL","NS","NT","NU","ON","PE","QC","SK","YT"), "CAN", "USA")

## tree
SS02 <- sqlFetch(con, "dbo_TREE_BBSBAM_V4_tbl")
SS02 <- nonDuplicated(SS02, SS, TRUE)
SS02 <- SS02[rownames(SS01),]
SS02$TREE[SS02$TREE > 100] <- NA
SS02$TREE <- SS02$TREE / 100
SS02$TREE3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
SS02$TREE3[SS02$TREE < 0.25] <- "Open"
SS02$TREE3[SS02$TREE >= 0.25 & SS02$TREE < 0.60] <- "Sparse"
SS02$TREE3[SS02$TREE >= 0.60] <- "Dense"


## point level land cover
SS03 <- sqlFetch(con, "dbo_BAMBBS_LANDCOVER_PTS")
SS03 <- nonDuplicated(SS03, SS, TRUE)
SS03 <- SS03[rownames(SS01),]

## LCC05 classes
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
SS03$LCC05_PT[SS03$LCC05_PT < 1 | SS03$LCC05_PT > 39] <- NA
SS03$LCC05_PT[SS01$COUNTRY == "USA"] <- NA
SS03$HAB_LCC1 <- ltlcc$BAMLCC05V2_label1[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC1 <- relevel(SS03$HAB_LCC1, "ConifDense")
SS03$HAB_LCC2 <- relevel(SS03$HAB_LCC2, "Conif")

## reclass EOSD
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
levels(SS03$EOSD_PT) <- sub(",", "", levels(SS03$EOSD_PT))
SS03$EOSD_PT <- as.integer(as.character(SS03$EOSD_PT))
SS03$EOSD_PT[SS03$EOSD_PT < 1] <- NA
SS03$HAB_EOSD1 <- lteosd$Reclass_label1[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD2 <- lteosd$Reclass_label2[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD1 <- relevel(SS03$HAB_EOSD1, "ConifDense")
SS03$HAB_EOSD2 <- relevel(SS03$HAB_EOSD2, "Conif")

## reclass NALCMS
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

## main SS level object
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

#### project summary table -- locally tweaked
#PCODE <- sqlFetch(con, "dbo_National_Proj_Summary_V4_2015")
#PCODE$SSMA_TimeStamp <- NULL
PCODE <- read.csv(file.path(ROOT,"proj.csv"))
levels(PCODE$Maxdist) <- tolower(levels(PCODE$Maxdist))
levels(PCODE$Maxdist)[levels(PCODE$Maxdist)=="unlimited"] <- "Inf"
PCODE$Maxdist[PCODE$Maxdist=="unknown"] <- NA
PCODE$Maxdist <- droplevels(PCODE$Maxdist)
PCODE$Maxdist <- as.numeric(as.character(PCODE$Maxdist))
PCODE$Maxdur <- pmin(PCODE$MaxDuration, 10)


#### survey level stuff

pkbam <- sqlFetch(con, "dbo_National_PKEY_V4_2015")
pkbam$SSMA_TimeStamp <- NULL
pkbbs <- sqlFetch(con, "dbo_PKEY_BBS_V3_2015")

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
rm(pkbam, pkbbs)
gc()

#### offsets

MM <- ifelse(PKEY$MONTH < 10, paste0("0", PKEY$MONTH), as.character(PKEY$MONTH))
HH <- ifelse(PKEY$HOUR < 10, paste0("0", PKEY$HOUR), as.character(PKEY$HOUR))
mm <- ifelse(PKEY$MIN < 10, paste0("0", PKEY$MIN), as.character(PKEY$MIN))
DD <- with(PKEY, paste0(YEAR, "-", MM, "-", DAY, " ", HH, ":", mm, ":00"))
DD <- strptime(DD, "%Y-%m-%e %H:%M:%S")
PKEY$JDAY <- DD$yday / 365
summary(PKEY$JDAY)
## prevent too far extrapolation
PKEY$JDAY[PKEY$JDAY < 0.35 | PKEY$JDAY > 0.55] <- NA

## TSSR
library(maptools)
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

offdat <- data.frame(PKEY[,c("PKEY","SS","TSSR","JDAY","MAXDUR","MAXDIS")],
    SS[match(PKEY$SS, rownames(SS)),c("LCC_combo","TREE")])
offdat$srise <- PKEY$srise + PKEY$MDT_offset
summary(offdat)

library(pbapply)
library(detect)
load_BAM_QPAD(version=1)
BAMspp <- getBAMspecieslist()
load("~/Dropbox/abmi/intactness/dataproc/BAMCOEFS25.Rdata")

#i="OVEN";j=1
#x=DAT1[,c("QPADLCCcombo","TREE","JDAY","TSSR","MAXDIS","MAXDUR")]
offset_fun <- function(j, i, x) {
    if (i %in% BAMspp) { ## spp with offsets (>=75 obs)
        type <- "full"
    } else {
        if (i %in% BAMCOEFS25$spp) { ## SPP based on >=25 obs
            type <- "const"
        } else {
            stop("no offsets available")
        }
    }
    if (type == "full") {
        xna <- as.data.frame(lapply(x[,c("LCC_combo","TREE","JDAY","TSSR")], 
            function(z) as.integer(is.na(z))))
        xna <- with(xna, interaction(LCC_combo,TREE,JDAY,TSSR, drop=TRUE, sep=""))
        xx <- lapply(levels(xna), function(z) x[xna==z,])
        id <- lapply(levels(xna), function(z) which(xna==z))
        names(xx) <- levels(xna)
        off <- lapply(xx, function(z) offset_fun0(j, i, x=z, type=type))
        OFF <- numeric(nrow(x))
        for (k in seq_len(length(off)))
            OFF[id[[k]]] <- off[[k]]
        gc <- globalBAMcorrections(i, 
            r=x$MAXDIS[!is.finite(OFF)],
            t=x$MAXDUR[!is.finite(OFF)])
        gc <- rowSums(log(as.matrix(gc)))
        OFF[!is.finite(OFF)] <- gc
    } else {
        OFF <- offset_fun0(j, i, x=x, type=type)
    }
    OFF
}

offset_fun0 <- function(j, i, x, type=c("full","const")) {
    BOOT <- j != 1
    if (type == "full") { ## spp with offsets (>=75 obs)
        ml1 <- getBAMmodellist()$sra
        ml2 <- getBAMmodellist()$edr
        if (any(is.na(x$JDAY)))
            ml1 <- ml1[!grepl("JDAY", ml1)]
        if (any(is.na(x$TSSR)))
            ml1 <- ml1[!grepl("TSSR", ml1)]
        if (any(is.na(x$TREE)))
            ml2 <- ml2[!grepl("TREE", ml2)]
        if (any(is.na(x$LCC_combo)))
            ml2 <- ml2[!grepl("LCC", ml2)]
        best <- bestmodelBAMspecies(i, 
            model.sra=names(ml1), model.edr=names(ml2),
            type=ifelse(BOOT, "multi", "BIC"))
        out <- with(x, localBAMcorrections(i,
            r=MAXDIS,
            t=MAXDUR,
            jday=JDAY, 
            tssr=TSSR, 
            tree=TREE, 
            lcc=LCC_combo,
            model.sra=best$sra, 
            model.edr=best$edr,
            boot=BOOT, ## MVN approach of j != 1
            ver=1))
    }
    if (type == "const") { ## SPP based on >=25 obs
        if (BOOT) {
            PHI <- exp(rnorm(1, 
                mean = BAMCOEFS25$sra_estimates[[i]][["0"]]$coef[1],
                sd = BAMCOEFS25$sra_estimates[[i]][["0"]]$vcov[1,1]))
            TAU <- exp(rnorm(1, 
                mean = BAMCOEFS25$edr_estimates[[i]][["0"]]$coef[1],
                sd = BAMCOEFS25$edr_estimates[[i]][["0"]]$vcov[1,1]))
        } else {
            PHI <- exp(BAMCOEFS25$sra_estimates[[i]][["0"]]$coef[1])
            TAU <- exp(BAMCOEFS25$edr_estimates[[i]][["0"]]$coef[1])
        }
        out <- customBAMcorrections(r=x$MAXDIS,
            t=x$MAXDUR,
            phi=PHI,
            tau=TAU)
    }
    #out
    rowSums(log(as.matrix(out)))
}

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







#### methodology stuff
BEH <- sqlFetch(con, "dbo_DD_DescripBEH")
DISINT <- sqlFetch(con, "dbo_DD_DescripDistance")
DISINT$SSMA_TimeStamp <- NULL
DURINT <- sqlFetch(con, "dbo_DD_DescripPeriod")
DISMET <- sqlFetch(con, "dbo_DD_distance_codes_methodology")
DURMET <- sqlFetch(con, "dbo_DD_duration_codes_methodology")


#### counts

pc1 <- sqlFetch(con, "dbo_National_PtCount_V4_2015")
pc1$SSMA_TimeStamp <- NULL
pc2 <- sqlFetch(con, "dbo_PointCount_BBS_V3_2015")


close(con)

PCODE <- proj
SS <- xy
PKEY <- pk
PCTBL <- pc


#levels(PCODE$DISTMETH)
#levels(PCODE$DISTMETH) <- toupper(levels(PCODE$DISTMETH))
DISMET
levels(DISMET$DISTANCECODE)
levels(DISMET$DISTANCECODE) <- gsub(" *$", "", levels(DISMET$DISTANCECODE))

levels(PCODE$DURMETH)
DURMET
levels(DURMET$DURATIONCODE)
levels(DURMET$DURATIONCODE) <- gsub(" *$", "", levels(DURMET$DURATIONCODE))


## excluding non-aerial detections
BEH
sort(100 * table(PCTBL$BEH) / sum(table(PCTBL$BEH)))
## 1=Heard
## 11=no birds observed at station - added 2011
## 6=seen and heard
PCTBL <- PCTBL[PCTBL$BEH %in% c(1,11,6),]

## excluding >10 min intervals
DURINT
## 10=10-20
## 11=0-20
## 3=before or after
## 8=unk
## 9=10-15
PCTBL <- PCTBL[!(PCTBL$DURATION %in% c(10,11,3,8,9)),]

