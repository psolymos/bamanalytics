## taking pre-processed BAM+BBS data and subset it to BCR4+

library(mefa4)
load("e:/peter/bam/Apr2016/out/data/pack_2016-12-01.Rdata")
rm(BB, mods)

## visualize the location subset
library(sp)
library(raster)
library(leaflet)
tmp <- nonDuplicated(DAT, SS, TRUE)
ss <- tmp$JURS %in% c("AK", "AB", "BC", "NT", "YK") & tmp$xBCR %in% c(1:6, 10)
xsp <- tmp[ss,c("X","Y")]
coordinates(xsp) <- ~X+Y
proj4string(xsp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
m <- mapview(xsp)


## data subset
ss <- DAT$JURS %in% c("AK", "AB", "BC", "NT", "YK") & DAT$xBCR %in% c(1:6, 10)
addmargins(with(droplevels(DAT[ss, ]), table(JURS, xBCR)))

DAT <- droplevels(DAT[ss,])
YY <- YY[ss,]
OFF <- OFF[ss,]

e <- new.env()
load("e:/peter/bam/Apr2016/out/data_package_2016-12-01.Rdata", envir=e)

SS <- e$SS[,c("PCODE", "SS", "TREE", "HAB_NALC2", "HAB_NALC1",
    "slp90", "cti90", "elv90", "YearFire", "FIRE_HA", "YearLoss",
    "CMD", "CMI", "CMIJJA", "DD01", "DD51", "EMT", "FFP",
    "ID", "MAP", "MAT", "MCMT", "MSP", "MWMT", "NFFD", "PAS", "PET",
    "PPT_sm", "PPT_wt", "TD", "BEADTotalL", "BEADtotPol", "HEIGHTSIMARD")]
SS <- SS[match(DAT$SS, SS$SS),]
DAT$DATE <- e$PKEY$DATE[match(rownames(DAT), e$PKEY$PKEY)]
DAT$YEAR <- e$PKEY$YEAR[match(rownames(DAT), e$PKEY$PKEY)]

DAT <- data.frame(DAT[,c("SS", "SITE_YR", "X", "Y", "Xcl", "Ycl", "xBCR", "JURS",
    "ROAD", "YEAR", "DATE")], SS)
DAT$SS.1 <- NULL

library(RODBC)
con <- odbcConnectAccess2007("c:/bam/May2015/BAM_BayneAccess_BAMBBScore.accdb")
SS03 <- sqlFetch(con, "dbo_BAMBBS_LANDCOVER_PTS")
close(con)
#save(SS03, file="e:/peter/bam/bcr4/SS03.RData")
SS03 <- nonDuplicated(SS03, SS, TRUE)
DAT$NALC <- SS03$NALCMS_PT[match(DAT$SS, SS03$SS)]

## Resampling:
## - downsize AB (what %? -- based on Area)
## - ? within jurisdiction blocks?
## - temporal blocks (5-6 yrs)
## - use AK, BC, AB, YT+NT
## - roadside?

DAT$YR2 <- ifelse(DAT$YEAR<2007, 0, 1)
DAT$bootg <- interaction(DAT$JURS, DAT$YR2, drop=TRUE)
levels(DAT$bootg)[levels(DAT$bootg) %in% c("NT.0", "YK.0")] <- "NTYK.0"
levels(DAT$bootg)[levels(DAT$bootg) %in% c("NT.1", "YK.1")] <- "NTYK.1"
table(DAT$bootg)

bbfun2 <- function(DAT1, B, out=0.1, seed=1234) {
    set.seed(seed)
    DAT1$SS_YR <- interaction(DAT1$SS, DAT1$YEAR, drop=TRUE)
    DAT1$IDMAP <- seq_len(nrow(DAT1))
    ## randomize input
    DAT1 <- DAT1[sample.int(nrow(DAT1)),]
    kkk <- floor(nlevels(DAT1$SS_YR) * (1-out))
    DAT1k <- DAT1[as.integer(DAT1$SS_YR) <= kkk,]
    if (nlevels(droplevels(DAT1k$bootg)) != nlevels(droplevels(DAT1$bootg)))
        stop("bootg problem: pick larger blocks for validation")
    ## one run
    r1fun <- function(DAT1k, n=NULL, replace=FALSE) {
        ## get rid of resamples
        DAT1k <- DAT1k[sample.int(nrow(DAT1k)),]
        DAT1k <- nonDuplicated(DAT1k, SS_YR)
        id2 <- list()
        for (l in levels(DAT1k$bootg)) {
            sset0 <- which(DAT1k$bootg == l)
            if (is.null(n))
                n <- length(sset0)
            id2[[l]] <- if (length(sset0) < 2)
                sset0 else sample(sset0, n, replace=replace)
        }
        DAT1k$IDMAP[unname(unlist(id2))]
    }
    n <- 5700
    BB0 <- r1fun(DAT1k, n=n, replace=FALSE)
    BB1 <- pbsapply(seq_len(B), function(i) r1fun(DAT1k, n=n, replace=TRUE))
    cbind(BB0, BB1)
}

B <- 240-1
BB <- bbfun2(DAT, B)
dim(BB)
dim(DAT)

nrow(BB)/nrow(DAT)
aa <- unique(as.numeric(BB))
DAT$HeldOut <- !(seq_len(nrow(DAT)) %in% aa)
bb <- table(selected=seq_len(nrow(DAT)) %in% aa)
bb[2]/sum(bb)

YY <- YY[,colSums(YY)>0]
OFF <- OFF[,colnames(YY)]

## data transformations

## years since last disturbance
DAT$ysd <- ifelse(is.na(DAT$YearFire), 100, DAT$YEAR - DAT$YearFire)
DAT$ysd[DAT$ysd < 0] <- 100
DAT$ysd10 <- ifelse(DAT$ysd <= 10, 1, 0)
DAT$ysd60 <- pmax(0, 1 - DAT$ysd / 60)
DAT$ysd100 <- pmax(0, 1 - DAT$ysd / 100)
DAT$ysdcat <- cut(DAT$ysd, c(0, 10, 99, 100), include.lowest=TRUE)
levels(DAT$ysdcat) <- c("early", "mid", "old")
DAT$ysdcat <- relevel(DAT$ysdcat, "old")
DAT$ysdmid <- ifelse(DAT$ysdcat == "mid", 1, 0)
DAT$ysd2 <- DAT$ysd^2
DAT$ysd3 <- DAT$ysd^3
DAT$ysd05 <- sqrt(DAT$ysd)

## land cover reclass
lt <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
DAT$hab <- reclass(DAT$NALC, lt[,c("Value", "BCR4")], allow_NA=TRUE)
DAT$hab <- relevel(DAT$hab, "Conif")
with(DAT[aa,], table(xBCR, hab))
table(DAT$NALC, DAT$hab)

## indicator variables for interactions
DAT$isConif <- ifelse(DAT$hab == "Conif", 1, 0)
DAT$isDecid <- ifelse(DAT$hab == "Decid", 1, 0)
DAT$isForest <- ifelse(DAT$hab %in% c("Decid", "Conif", "Wet"), 1, 0)

## terrain
DAT$slp <- sqrt(DAT$slp90)
DAT$cti <- (DAT$cti90 - 8) / 4
DAT$elv <- sqrt(DAT$elv90) / 25 - 1
DAT$slp2 <- DAT$slp^2
DAT$slp3 <- DAT$slp^3
DAT$cti2 <- DAT$cti^2
DAT$elv2 <- DAT$elv^2
DAT$cti3 <- DAT$cti^3
tmp <- data.frame(y=YY[,"OSFL"], DAT[,c("slp", "cti", "elv")])
#siplot(y ~ ., tmp, B=0)
cor(tmp)

## climate
DAT$CMIJJA <- (DAT$CMIJJA - 0) / 50
DAT$CMI <- (DAT$CMI - 0) / 50
DAT$TD <- (DAT$TD - 300) / 100
DAT$DD0 <- (DAT$DD01 - 1000) / 1000
DAT$DD5 <- (DAT$DD51 - 1600) / 1000
DAT$EMT <- (DAT$EMT + 400) / 100
DAT$MSP <- (DAT$MSP - 400) / 200

DAT$DD02 <- DAT$DD0^2
DAT$DD52 <- DAT$DD5^2
DAT$CMI2 <- DAT$CMI^2
DAT$CMIJJA2 <- DAT$CMIJJA^2

cor(DAT[,c("CMI", "CMIJJA", "TD", "DD0", "DD5", "EMT", "MSP")])
## look up what these mean


save(YY, OFF, DAT, BB, file="e:/peter/bam/bcr4/bcr4-data.RData")
