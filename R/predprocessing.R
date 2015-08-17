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
    "NORM_6190_CMIJJApm", "NORM_6190_DD0", "NORM_6190_DD5", 
    "NORM_6190_EMT", "NORM_6190_TD", "NORM_6190_CMIpm", "NORM_6190_MSP")
write.csv(z[,cn1], file=file.path(ROOT, "cti", "clim-only.csv"))
## NALC 4x4
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


## slope, cti, elev
pg3$SLP <- sqrt(pg3$slope90)
clim$SLP <- pg3$SLP[match(clim$pointid, pg3$pointid)]
clim$SLP2 <- clim$SLP^2
rm(pg3,pg3a,pg3ab,pg3b)

save(clim, file=file.path(ROOT, "pg-clim.Rdata"))

## -----------------------------------

x <- fread(file.path(ROOT, "PredictionIntersections",
    "CovariatesPredGridOneforpeter.txt"))
x$BCRNAME <- NULL
x$COUNTRY <- as.factor(x$COUNTRY)
x$PROVINCE_S <- as.factor(x$PROVINCE_S)


## REG (based on BCR and Jurs)
x$JURS <- x$PROVINCE_S
#x$JURS <- droplevels(x$JURS)

levels(x$JURS) <- c("", "AK", "AB", "BC", "CA", 
    "CT", "DE", "ID", "IL", "IN", "IA", 
    "ME", "MB", "MD", "MA", "MI", 
    "MN", "MO", "MT", "NE", "NV", "NB", 
    "NH", "NJ", "NY", "NL", "ND", 
    "NT", "NS", "NU", "OH", "ON", 
    "OR", "PA", "PE", "QC", "RI", 
    "SK", "SD", "UT", "VT", "VA", 
    "WA", "WV", "WI", "WY", "YK")

x$BCR_JURS <- interaction(x$BCR, x$JURS, drop=TRUE, sep="_")
x$BCR_JURS0 <- x$BCR_JURS

levels(x$BCR_JURS)[grepl("_AK", levels(x$BCR_JURS))] <- "AK"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("9_BC","9_ID","9_WA",
    "5_BC","5_WA","10_AB","10_BC","10_ID","10_MT","10_WA")] <- "Mtn"
levels(x$BCR_JURS)[grepl("11_", levels(x$BCR_JURS))] <- "Pra" # Prairies
levels(x$BCR_JURS)[grepl("17_", levels(x$BCR_JURS))] <- "Pra"
levels(x$BCR_JURS)[grepl("22_", levels(x$BCR_JURS))] <- "Pra"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("4_BC","4_NT","4_YK")] <- "4_all"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("6_AB","6_BC","6_SK","6_MB",
    "6_YK","6_NT")] <- "6_all"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("6_AB","6_BC","6_SK","6_MB")] <- "6_south"

## keep northern points in
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_NT","7_MB","7_NT")] <- "3+7_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_NU","7_ON","7_QC","7_NL")] <- "3+7_east"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("8_MB","8_SK")] <- "8_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("8_NL","8_ON","8_QC")] <- "8_east"

levels(x$BCR_JURS)[grepl("12_", levels(x$BCR_JURS))] <- "Grl" # Great Lakes
levels(x$BCR_JURS)[grepl("13_", levels(x$BCR_JURS))] <- "Grl"
levels(x$BCR_JURS)[grepl("23_", levels(x$BCR_JURS))] <- "Grl"

levels(x$BCR_JURS)[grepl("14_", levels(x$BCR_JURS))] <- "Mar" # Maritimes
levels(x$BCR_JURS)[grepl("30_", levels(x$BCR_JURS))] <- "Mar"

levels(x$BCR_JURS)[grepl("24_", levels(x$BCR_JURS))] <- "Seus" # SE US
levels(x$BCR_JURS)[grepl("26_", levels(x$BCR_JURS))] <- "Seus"
levels(x$BCR_JURS)[grepl("28_", levels(x$BCR_JURS))] <- "Seus"
levels(x$BCR_JURS)[grepl("29_", levels(x$BCR_JURS))] <- "Seus"

sort(levels(x$BCR_JURS))
table(x$BCR_JURS)

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("5_YK")] <- "Mtn"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_YK")] <- "3+7_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_MB")] <- "3+7_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_NL")] <- "3+7_east"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("3_QC")] <- "3+7_east"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("6_NU")] <- "6_all"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("8_AB")] <- "8_west"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("7_SK")] <- "3+7_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("7_AB")] <- "3+7_west"
levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("7_NU")] <- "3+7_east"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("6_MN")] <- "6_all"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("5_CA",  "5_OR" ,
    "9_CA",  "9_NV", "9_OR", "9_UT",  "9_WY",
    "10_OR", "10_UT", "10_WY",
    "16_ID",  "16_UT", "16_WY")] <- "Mtn"

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("18_NE", "18_SD", "18_WY",
    "19_NE", "19_SD")] <- "Pra"      

levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("0_",   "100_"  )] <- "UNKNOWN"


## regions for trend
x$REG <- x$BCR_JURS
levels(x$REG)[levels(x$REG) %in% c("Mtn", "AK")] <- "Coast"
levels(x$REG)[levels(x$REG) %in% c("4_all", "3+7_west", "3+7_east")] <- "North"
levels(x$REG)[levels(x$REG) %in% c("Pra")] <- "South"
levels(x$REG)[levels(x$REG) %in% c("8_west","6_all")] <- "West"
levels(x$REG)[levels(x$REG) %in% c("Mar","Seus","8_east","Grl")] <- "East"

table(x$BCR_JURS, x$REG)
table(x$BCR_JURS0, x$REG)
x$REG <- relevel(x$REG, "West")
table(x$REG)
#plot(x[,2:3],pch=".",col=x$REG)

## EOSC/LCC/NALC skeleton: gridcode and eosd extent definition !!!

#xx <- x[sample.int(nrow(tmp), 10000),c("POINT_X","POINT_Y","eosdbam")]
#plot(xx[,1:2], col=ifelse(x[,3]==0, 1, 2), pch=".")
x$EOSD_COVER <- ifelse(x$eosdbam > 0 & x$COUNTRY == "CANADA" & x$REG != "South", 1L, 0L)
#x <- x[x$eosdbam > 0 & x$COUNTRY == "CANADA" & x$REG != "South",]
x$BCR <- NULL
x$PROVINCE_S <- NULL
x$COUNTRY <- NULL
x$JURS <- NULL
x$BCR_JURS <- NULL
#x$BCR_JURS0 <- NULL
x$BCR <- NULL
gc()

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

## HAB: NALC
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
x$HAB_NALC2 <- ltnalc$Label[match(x$nalcms_05, ltnalc$Value)]
x$HAB_NALC2 <- relevel(x$HAB_NALC2, "Conif")
table(x$nalcms_05, x$HAB_NALC2,useNA="a")
x$nalcms_05 <- NULL
#levels(x$HAB_NALC2) <- c(levels(x$HAB_NALC2), "0-density")
#x$HAB_NALC2[is.na(x$HAB_NALC2)] <- "0-density"
table(x$HAB_NALC2)

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

## isDM, isNF
## decid + mixed
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

## HGT, HGT2
x$HGT <- x$SimardG / 50
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

## LIN, POL
## linear features
x$LIN <- log(x$BEAD_linear + 1)
x$BEAD_linear <- NULL
## Polygon is fine as it is (0-1)
x$POL <- x$BEAD_poly
x$BEAD_poly <- NULL

x <- droplevels(x[!(x$BCR_JURS0 %in% c("0_","100_")),])
save(x, file=file.path(ROOT, "pg-main.Rdata"))

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

loss <- z
save(loss, file=file.path(ROOT, "pg-loss.Rdata"))

## LCC 4x4

fn <- file.path(ROOT, "PredictionIntersections", "Pred_LCC05_4x4CAN.csv")
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
ltlcc$vv <- paste0("LCCVV", ltlcc$lcc05v1_2)
ltlcc$BAMLCC05V2_label2 <- as.character(ltlcc$BAMLCC05V2_label2)
ltlcc$BAMLCC05V2_label2[is.na(ltlcc$BAMLCC05V2_label2)] <- "NONE"

tmp0 <- read.csv(fn, nrows=2, skip=0, header=TRUE)

f1 <- function(tmp0) {
    tmp0$LCCVVSUM <- NULL
    tmp0$LCCVV0 <- NULL
    m <- as(as.matrix(tmp0[,grepl("LCCVV", colnames(tmp0))]), "dgCMatrix")
    rownames(m) <- gsub(",", "", tmp0$pointid)
    groupSums(m, 2, ltlcc$BAMLCC05V2_label2)
}
#f1(tmp0)
#rBind(tmp0,tmp0)

nlines <- function(file) {
    ## http://r.789695.n4.nabble.com/Fast-way-to-determine-number-of-lines-in-a-file-td1472962.html
    ## needs Rtools on Windows
    if (.Platform$OS.type == "windows") { 
        nr <- as.integer(strsplit(system(paste("/RTools/bin/wc -l", 
            file), intern=TRUE), " ")[[1]][1])
    } else {
        nr <- as.integer(strsplit(system(paste("wc -l", 
            file), intern=TRUE), " ")[[1]][1])
    }
    nr
}
nlines(fn)

n <- 10^6
nr <- 7935918
cn <- colnames(tmp0)
#res <- list()
m <- floor((nr-1) / n)
mm <- (nr-1) %% n
if (mm > 0)
    m <- m+1
for (i in 1:m) {
    cat(i, "/", m, "\n");flush.console()
    tmp <- fread(fn, nrows=n, skip=(i-1)*n+1, header=FALSE)
    colnames(tmp) <- cn
    res <- f1(tmp)
    save(res, file = file.path(ROOT, paste0("lcc4x4-", i, ".Rdata")))
    rm(tmp, res)
    gc()
}

## EOSD 4x4

fn <- file.path(ROOT, "eosd4x4", "EOSD4x4Pred_eosdextent.csv")
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
lteosd$vv <- paste0("eosdVV", lteosd$Value)
lteosd$Reclass_label2 <- as.character(lteosd$Reclass_label2)
rownames(lteosd) <- lteosd$vv

tmp0 <- read.csv(fn, nrows=2, skip=0, header=TRUE)
tmp0$eosdVVSUM <- NULL
m <- as(as.matrix(tmp0[,grepl("eosdVV", colnames(tmp0))]), "dgCMatrix")
rownames(m) <- gsub(",", "", tmp0$pointid)
lteosd <- lteosd[colnames(m),]
rownames(lteosd) <- colnames(m)
lteosd$Reclass_label2[is.na(lteosd$Reclass_label2)] <- "OUT"

f1 <- function(tmp0) {
    tmp0$eosdVVSUM <- NULL
    m <- as(as.matrix(tmp0[,grepl("eosdVV", colnames(tmp0))]), "dgCMatrix")
    rownames(m) <- gsub(",", "", tmp0$pointid)
    groupSums(m, 2, lteosd$Reclass_label2)
}
#f1(tmp0)
#rBind(tmp0,tmp0)

nlines(fn)

n <- 10^6
nr <- 6357648
cn <- colnames(tmp0)
#res <- list()
m <- floor((nr-1) / n)
mm <- (nr-1) %% n
if (mm > 0)
    m <- m+1
for (i in 1:m) {
    cat(i, "/", m, "\n");flush.console()
    tmp <- fread(fn, nrows=n, skip=(i-1)*n+1, header=FALSE)
    colnames(tmp) <- cn
    res <- f1(tmp)
    save(res, file = file.path(ROOT, paste0("eosd4x4-", i, ".Rdata")))
    rm(tmp, res)
    gc()
}

nalc4x4 <- fread(file.path(ROOT, "cti", "nalc4x4-processed.csv"))
rownames(nalc4x4) <- nalc4x4[,1]
nalc4x4[,1] <- NULL
nalc4x4 <- as(as.matrix(nalc4x4), "dgCMatrix")

load(file.path(ROOT, paste0("eosd4x4-", 1, ".Rdata")))
eosd4x4 <- res
for (i in 2:7) {
    cat(i, "\t");flush.console()
    load(file.path(ROOT, paste0("eosd4x4-", i, ".Rdata")))
    eosd4x4 <- rBind(eosd4x4, res)
}
rm(res)

load(file.path(ROOT, paste0("lcc4x4-", 1, ".Rdata")))
lcc4x4 <- res
for (i in 2:8) {
    cat(i, "\t");flush.console()
    load(file.path(ROOT, paste0("lcc4x4-", i, ".Rdata")))
    lcc4x4 <- rBind(lcc4x4, res)
}
rm(res)

save(lcc4x4, eosd4x4, nalc4x4, file=file.path(ROOT, "pg-4x4.Rdata"))

################# PACKAGING ###############################

## use eosd coverage to save bcr/jurs0 chunks

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)

load(file.path(ROOT, "pg-main.Rdata"))
x <- x[x$EOSD_COVER == 1,]
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

x <- x[!is.na(x$CTI) & ! is.na(x$TD),]

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

x$TR3[is.na(x$TR3)] <- "Open" # this is global

XYeosd <- as.matrix(x[,2:3])
rownames(XYeosd) <- x[,1]
save(XYeosd, file=file.path(ROOT, "XYeosd.Rdata"))

reg <- levels(x$BCR_JURS0)
for (i in reg) {
    cat(i, "\n");flush.console()
    ii <- x$BCR_JURS0 == i
    dat <- x[ii,]
    pg4x4 <- list(lcc=lcc4x4[ii,], eosd=eosd4x4[ii,], nalc=nalc4x4[ii,])
    save(dat, pg4x4, file=file.path(ROOT, "chunks", paste0("pgdat-", i, ".Rdata")))
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

x <- x[!is.na(x$CTI) & ! is.na(x$TD),]
x$DD02 <- x$DD0^2
x$DD52 <- x$DD5^2

x$REG <- droplevels(x$REG)
x$BCR_JURS0 <- droplevels(x$BCR_JURS0)

x$TR3[is.na(x$TR3)] <- "Open" # this is global

x$pointid <- NULL

XYfull <- as.matrix(x[,c("POINT_X","POINT_Y")])
rownames(XYfull) <- x[,1]
save(XYfull, file=file.path(ROOT, "XYfull.Rdata"))

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


