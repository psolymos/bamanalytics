library(RColorBrewer)
library(mefa4)
ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/bam/pred-2015"
ROOT3 <- "e:/peter/bam/pred-2016"
source("~/repos/bamanalytics/R/makingsense_functions.R")

PROJECT <- "bam"
Date <- "2016-12-01"

Stage <- 6 # which(names(mods) == "Clim")
BASE_YEAR <- 2012#2002#2012
Percent <- 10
Percent <- as.integer(Percent)

e <- new.env()
load(file.path("e:/peter/bam/Apr2016/out", "data", "pack_2016-12-01.Rdata"), envir=e)
mods <- e$mods
DAT <- e$DAT
DAT$BCR <- DAT$xBCR
DAT$pointid <- NA
DAT$Brandt <- NA
DAT$LEVEL3 <- NA
DAT$POINT_X <- DAT$Xcl
DAT$POINT_Y <- DAT$Ycl
DAT$PKEY <- rownames(DAT)

Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
rm(e)

st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs <- levels(st$LEVEL3)

## includes land cover, height, and fire
hab_col <- c("HABAgr", "HABBarren", "HABDecid", "HABDevel",
    "HABGrass", "HABMixed", "HABShrub", "HABWet", "HABTRAgr", "HABTRBarren",
    "HABTRConifOpen", "HABTRConifSparse", "HABTRDecidDense", "HABTRDecidOpen",
    "HABTRDecidSparse", "HABTRDevel", "HABTRGrass", "HABTRMixedDense",
    "HABTRMixedOpen", "HABTRMixedSparse", "HABTRShrub", "HABTRWetDense",
    "HABTRWetOpen", "HABTRWetSparse",
    "HGT", "HGT2", "HGT05",
    "DTB", "BRN", "LSS", "YSD", "YSF", "YSL",
    "HGT:isDM", "HGT:isWet", "HGT:isDec", "HGT:isMix",
    "HGT2:isDM", "HGT2:isWet", "HGT2:isDec", "HGT2:isMix", "HGT05:isDM",
    "HGT05:isWet", "HGT05:isDec", "HGT05:isMix")

df_sub <- NULL
#regi <- "10.1.1"
for (i in 1:length(regs)) {

regi <- regs[i]
cat(regi, " (", i, "/", length(regs), ")\n", sep="");flush.console()

load(file.path(ROOT3, "chunks", paste0("pgdat-", regi, ".Rdata")))
gc()

#xy <- rbind(df_sub, dat[,c("pointid", "POINT_X", "POINT_Y")])
#save(xy, file="BAM_pred_grid_xy.Rdata")

seq_id <- 1 + (seq_len(nrow(dat)) - 1) %% 100
seq_id <- sample(seq_id)
subset_id <- seq_id <= Percent

dat$HAB <- dat$HAB_NALC2
dat$HABTR <- dat$HAB_NALC1
dat$HGT[dat$HAB %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- 0
dat$HGT2 <- dat$HGT^2
dat$HGT05 <- dat$HGT^0.5

dat$ROAD <- 0L

## YR
dat$YR <- BASE_YEAR - 2001

## disturbance
dat$YearFire[is.na(dat$YearFire)] <- BASE_YEAR - 200
dat$YearLoss[is.na(dat$YearLoss)] <- BASE_YEAR - 200
## backfill means no forest loss, only due to fire

## years since fire
dat$YSF <- BASE_YEAR - dat$YearFire
dat$YSF[dat$YSF < 0] <- 200
## years since loss
dat$YSL <- BASE_YEAR - dat$YearLoss
dat$YSL[dat$YSL < 0] <- 200
## years since most recent burn or loss, with backfill option
dat$YSD <- pmin(dat$YSF, dat$YSL)

## cut at 10 yrs
dat$BRN <- ifelse(dat$YSF <= 10, 1L, 0L)
dat$LSS <- ifelse(dat$YSL <= 10, 1L, 0L)
#dat$LSS[dat$YEAR < 2000] <- NA
dat$DTB <- ifelse(dat$YSD <= 10, 1L, 0L)
#dat$DTB[dat$YEAR < 2000] <- NA

## refining years since variables
AGEMAX <- 50
dat$YSD <- pmax(0, 1 - (dat$YSD / AGEMAX))
dat$YSF <- pmax(0, 1 - (dat$YSF / AGEMAX))
dat$YSL <- pmax(0, 1 - (dat$YSL / AGEMAX))
levels(dat$Brandt) <- c(levels(dat$Brandt), "OUT")
dat$Brandt[is.na(dat$Brandt)] <- "OUT"

dat$subreg <- paste(regi, dat$BCR, dat$JURS, dat$Brandt, sep=" + ")

## these are not considered in the models, thus NAs not tracked down !!!
dat$LIN <- dat$POL <- NULL

dat$SS <- NA
dat$PKEY <- NA

#dat0 <- dat[rowSums(is.na(dat)) == 0,]
#Xn0 <- model.matrix(getTerms(mods[1:Stage], "formula"), dat0)
#colnames(Xn0) <- fixNames(colnames(Xn0))

df_sub <- rbind(df_sub, dat[subset_id,intersect(colnames(DAT), colnames(dat)),drop=FALSE])

}

df_sub$BCR <- as.factor(df_sub$BCR)
DAT$JURS <- factor(as.character(DAT$JURS), levels(df_sub$JURS))
DAT$BCR <- factor(as.character(DAT$BCR), levels(df_sub$BCR))

stacked_data <- rbind(DAT[,colnames(df_sub)], df_sub)
stacked_data$sampled <- c(rep(1, nrow(DAT)), rep(0, nrow(df_sub)))

save(stacked_data, file="BAM_stacked_gap_data.Rdata")

cn <- c("sampled", "HABTR", "HGT", "DTB",
    "BRN", "LSS", "YSD", "YSF", "YSL", "CTI", "SLP",
    "CMIJJA", "DD0", "DD5", "EMT", "MSP", "CMI", "TD")

mod <- glm(sampled ~ ., stacked_data[,cn], family="binomial")
summary(mod)

## you can subset by BCR/JURS to fit regional models

rn <- stacked_data$BCR == "6"
mod <- glm(sampled ~ ., stacked_data[rn,cn], family="binomial")
summary(mod)

library(lattice)
histogram(~ TD | factor(sampled), stacked_data)

ii <- c("HGT", "CTI", "SLP",
    "CMIJJA", "DD0", "DD5", "EMT", "MSP",
    "CMI", "TD", "POINT_X", "POINT_Y")
for (i in ii) {
cat(i, "\n");flush.console()
z <- stacked_data[[i]]
q <- unique(quantile(z, seq(0, 1, by=0.05)))
zc <- cut(z, breaks=q, include.lowest=TRUE)
col <- topo.colors(nlevels(zc))
png(paste0("e:/peter/bam/Apr2016/out/thumbs/", i, ".png"), width=2000, height=1200)
plot(stacked_data$POINT_X, stacked_data$POINT_Y,
    axes=FALSE, ann=FALSE, col=col[zc], pch=".")
title(main=i)
dev.off()
}
