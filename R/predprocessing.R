## Data processing for nationam BAM analyses: prediction

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)

x <- read.csv(file.path(ROOT, "PredictionIntersections",
    "CovariatesPredGridOneforpeter.txt"))
x$BCRNAME <- NULL

## REG (based on BCR and Jurs)
#x$BCR[x$BCR == "0"] <- NA
#x$PROVINCE_S <- droplevels(x$PROVINCE_S)
#table(x$PROVINCE_S, x$BCR)

x$JURS <- x$PROVINCE_S
#x <- x[x$JURS != "",]
x$JURS <- droplevels(x$JURS)
#levels(x$JURS) <- c("AB", "BC", "MB", "NB", "NL", "NT", "NS", "NU", 
#    "ON", "PE", "QC", "SK", "YK")

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

#levels(x$BCR_JURS)[levels(x$BCR_JURS) %in% c("4_BC",
#    "4_NT","4_YK","6_YK","6_NT")] <- "4+6_YK+NT"
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

## all BCR 4 is in Coast (but 4+6_YK+NT kept separate for bootstrap blocks)
#x$REG[x$BCR == 4] <- "Coast"

table(x$BCR_JURS, x$REG)
table(x$BCR_JURS0, x$REG)
x$REG <- relevel(x$REG, "West")
#x <- x[x$REG != "South",]
#x$REG <- droplevels(x$REG)

table(x$REG)
#plot(x[,2:3],pch=".",col=x$REG)

## EOSC/LCC/NALC skeleton: gridcode and eosd extent definition !!!

#xx <- x[sample.int(nrow(tmp), 10000),c("POINT_X","POINT_Y","eosdbam")]
#plot(xx[,1:2], col=ifelse(x[,3]==0, 1, 2), pch=".")
x <- x[x$eosdbam > 0,]
x <- x[x$COUNTRY != "USA",]
gc()

## HAB: EOSD
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
x$eosdbam[x$eosdbam < 1] <- NA
x$HAB_EOSD2 <- lteosd$Reclass_label2[match(x$eosdbam, lteosd$Value)]
x$HAB_EOSD2 <- relevel(x$HAB_EOSD2, "Conif")
table(x$eosdbam, x$HAB_EOSD2,useNA="a")
x$eosdbam <- NULL
levels(x$HAB_EOSD2) <- c(levels(x$HAB_EOSD2), "0-density")
x$HAB_EOSD2[is.na(x$HAB_EOSD2)] <- "0-density"
table(x$HAB_EOSD2)

## HAB: NALC
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
x$HAB_NALC2 <- ltnalc$Label[match(x$nalcms_05, ltnalc$Value)]
x$HAB_NALC2 <- relevel(x$HAB_NALC2, "Conif")
table(x$nalcms_05, x$HAB_NALC2,useNA="a")
x$nalcms_05 <- NULL
levels(x$HAB_NALC2) <- c(levels(x$HAB_NALC2), "0-density")
x$HAB_NALC2[is.na(x$HAB_NALC2)] <- "0-density"
table(x$HAB_NALC2)

## HAB: LCC
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
x$LCCV1_3Can[x$LCCV1_3Can < 1 | x$LCCV1_3Can > 39] <- NA
x$LCCV1_3Can[x$COUNTRY == "USA"] <- NA
x$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(x$LCCV1_3Can, ltlcc$lcc05v1_2)]
x$HAB_LCC2 <- relevel(x$HAB_LCC2, "Conif")
table(x$LCCV1_3Can, x$HAB_LCC2,useNA="a")
x$LCCV1_3Can <- NULL
levels(x$HAB_LCC2) <- c(levels(x$HAB_LCC2), "0-density")
x$HAB_LCC2[is.na(x$HAB_LCC2)] <- "0-density"
table(x$HAB_LCC2)

## isDM, isNF
## decid + mixed
x$isDM_LCC <- ifelse(x$HAB_LCC2 %in% c("Decid", "Mixed"), 1L, 0L)
x$isDM_EOSD <- ifelse(x$HAB_EOSD2 %in% c("Decid", "Mixed"), 1L, 0L)
x$isDM_NALC <- ifelse(x$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
x$isNF_LCC <- ifelse(x$HAB_LCC2 %in% 
    c("Agr", "Barren", "Burn", "Devel", "Grass", "Wet", "0-density"), 1L, 0L)
x$isNF_EOSD <- ifelse(x$HAB_EOSD2 %in% 
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet", "0-density"), 1L, 0L)
x$isNF_NALC <- ifelse(x$HAB_NALC2 %in% 
    c("Agr", "Barren", "Devel", "Grass", "Shrub", "Wet", "0-density"), 1L, 0L)

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
