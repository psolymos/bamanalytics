## Data processing for nationam BAM analyses: prediction

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)
library(data.table)
getOption("datatable.fread.datatable")
options(datatable.fread.datatable=FALSE)
getOption("datatable.fread.datatable")

## Climate preprocessing -------------------------
z <- fread(file.path(ROOT, "cti", "ECGrid_climMarch2014.csv"))
#z <- read.csv(file.path(ROOT, "cti", "ECGrid_climMarch2014.csv"))
cn1 <- c("pointid", "POINT_X", "POINT_Y", 
    "NORM_6190_CMIJJApm", "NORM_6190_DD0", "NORM_6190_DD5", 
    "NORM_6190_EMT", "NORM_6190_TD", "NORM_6190_CMIpm", "NORM_6190_MSP")
write.csv(z[,cn1], file=file.path(ROOT, "cti", "clim-only.csv"))
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
zz$pointid <- gsub(",", "", zz$pointid)
zz$pointid <- as.integer(zz$pointid)

load(file.path(ROOT, "pg-main.Rdata"))
load(file.path(ROOT, "pg-clim.Rdata"))


## placeholders: HSH, HSH2, isDM, isNF
x$HSH <- 0
x$HSH2 <- 0
x$HAB <- 0
x$isDM <- 0
x$isNF <- 0

## YR
x$YR <- 2015 - 1997

## check NA occurrences, especially in HAB_*
plot(x[,2:3],pch=".",col=ifelse(is.na(x$HAB_EOSD2),2,1))
