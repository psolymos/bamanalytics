#' ---
#' title: "NWT offsets calculation"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "2019-01-14"
#' output: pdf_document
#' ---
#'
#' ## Preamble
#'
library(mefa4)
library(odbc)
library(DBI)
library(intrval)
library(maptools)
source("~/.ssh/boreal")
source("~/repos/abmianalytics/birds/00-functions.R")
knitr::opts_chunk$set(eval=FALSE)
ROOT <- "d:/bam/BAM_data_v2019"
#'
#' ## Pull data from BU data base
#'
con <- dbConnect(
    odbc::odbc(),
    dsn = "BOREAL",
    database = "EMCLA_SQL_Database",
    driver = "SQL Server",
    uid = ..boreal_db_access$uid,
    pwd = ..boreal_db_access$pwd)
#dbl <- dbListTables(con)
xybu <- dbReadTable(con, "viewCoord_LatLong")
d0 <- dbReadTable(con, "ViewSpecies_LongFormCount")
s0 <- dbReadTable(con, "list Species Code")

dbDisconnect(con)

d0 <- make_char2fact(d0)
s0 <- make_char2fact(s0)
xybu <- make_char2fact(xybu)
#'
#' ## Filtering by PCODE
#'
PCODE_keep <- c(
    "ECYKB",
    "ENWA",
    "NWTI",
    "FR",
    "SH",
    "LV",
    "WR",
    "ARTI")
xybu <- droplevels(xybu[xybu$ProjectID %in% PCODE_keep, ])
d0 <- droplevels(d0[d0$ProjectID %in% PCODE_keep, ])
table(xybu$ProjectID)
table(d0$ProjectID)

xybu$SS <- with(xybu, interaction(
    ProjectID,
    Cluster,
    SITE,
    STATION,
    sep="::", drop=TRUE))

d <- droplevels(d0[d0$Method %in% as.character(c(0, 8, 11:14)) & d0$Replicate == 1,])
d$MAXDUR <- 3
d$MAXDUR[d$Method %in% c("12", "13")] <- 1
d$MAXDUR[d$Method == "0"] <- 10
d$MAXDIS <- Inf
d$SS <- with(d, interaction(
    ProjectID,
    Cluster,
    SITE,
    STATION,
    sep="::", drop=TRUE))
tmp1 <- as.character(d$RECORDING_DATE)
tmp2 <- sapply(strsplit(as.character(d$RECORD_TIME), " "), "[[", 2)
d$YEAR <- as.numeric(substr(tmp1, 1, 4))
#d$DATI <- as.POSIXlt(paste(tmp1, tmp2))
d$DATI <- as.POSIXlt(paste(tmp1, tmp2), tz="America/Denver")
d$PKEY <- as.factor(paste0(
    as.character(d$SS),
    "_",
    tmp1,
    "-",
    tmp2))
d$DATE <- d$RECORDING_DATE
d$SSYR <- paste0(d$SS, "_", d$YEAR)
d$SPECIES <- normalize_sppcode(d$SPECIES)
d$X <- xybu$Longitude[match(d$SS, xybu$SS)]
d$Y <- xybu$Latitude[match(d$SS, xybu$SS)]

pk <- nonDuplicated(d, PKEY, TRUE)
pk$Abundance <- pk$SPECIES <- pk$RCCODE <- pk$Min0 <- pk$Min1 <- pk$Min2 <- pk$maxID <- NULL
#save(pk, file=file.path(ROOT, "nwt", "nwt-pkey-table-2019-01-14.RData"))

lc <- read.csv(file.path(ROOT, "nwt", "NWTpkey.csv"))
pk$landcov <- lc$landcov[match(pk$PKEY, lc$PKEY)]
pk$treecov <- lc$treecov[match(pk$PKEY, lc$PKEY)]
#'
#' ## JDAY/TSSR
#'
dd <- pk
dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
keep <- is.na(dd$DATI) # these will be constant phi
keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
keep[dd$start %[]% c(3, 12)] <- TRUE
dd <- droplevels(dd[keep,])
#' Normalized ordinal day
dd$JDAY <- dd$JULIAN / 365
#' Normalized time since local sunrise
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$TSSR <- (dd$start - dd$srise) / 24
#' Quadratic terms
dd$JDAY2 <- dd$JDAY^2
dd$TSSR2 <- dd$TSSR^2
#'
#' ## LCC/TREE
#'
## Reclass NALCMS
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
dd$NALC <- ltnalc$Label[match(dd$landcov, ltnalc$Value)]
dd$TREE <- dd$treecov / 100
dd$TREE[dd$TREE %)(% c(0, 1)] <- NA

dd$LCC4 <- as.character(dd$NALC)
dd$LCC4[dd$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
dd$LCC4[dd$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
dd$LCC4 <- factor(dd$LCC4,
    c("DecidMixed", "Conif", "Open", "Wet"))
dd$LCC2 <- as.character(dd$LCC4)
dd$LCC2[dd$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
dd$LCC2[dd$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
dd$LCC2 <- factor(dd$LCC2, c("Forest", "OpenWet"))
table(dd$LCC4, dd$LCC2, useNA="a")
#' Maximum counting distance in 100 m units (area in ha)
dd$MAXDIS <- dd$MAXDIS / 100

#'
#' ## Species data
#'
y <- as.matrix(Xtab(Abundance ~ PKEY + SPECIES, d))
dim(y)
#'
#' ## Offsets
#'

Xp <- cbind("(Intercept)"=1, as.matrix(dd[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=dd$TREE,
    LCC2OpenWet=ifelse(dd$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(dd$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(dd$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(dd$LCC4=="Wet", 1, 0))
summary(Xp)
summary(Xq)
#' Use version 3 of BAM QPAD estimates
library(QPAD)
load_BAM_QPAD(version=3)
sppp <- intersect(colnames(y), getBAMspecieslist())
#' Save values in a matrix
off <- matrix(NA, nrow(dd), length(sppp))
rownames(off) <- rownames(dd)
colnames(off) <- sppp
#' Loop for species
for (spp in sppp) {
    ## print out where we are
    cat(spp, "\n");flush.console()
    p <- rep(NA, nrow(dd))
    A <- q <- p
    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC",
        model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    ## design matrices matching the coefs
    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0
    ## calculate p, q, and A based on constant phi and tau for the respective NAs
    p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))
    ## calculate time/lcc varying phi and tau for non-NA cases
    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(dd$MAXDUR[OKp], phi1)
    unlim <- ifelse(dd$MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dd$MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[OKq], tau1))
    ## log(0) is not a good thing, apply constant instead
    ii <- which(p == 0)
    p[ii] <- sra_fun(dd$MAXDUR[ii], cf0[1])
    ## store, next
    off[,spp] <- log(p) + log(A) + log(q)
}
## sanity checks
(Ra <- apply(off, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))

#'
#' ## Save
#'
y <- y[rownames(dd), colnames(off)]
save(y, off, dd, file=file.path(ROOT, "nwt", "nwt-offsets-2019-01-15.RData"))


## WindTrax data
d0 <- read.csv(file.path(ROOT, "nwt", "speciesReport.csv"))
str(d0)
l0 <- read.csv(file.path(ROOT, "nwt", "wildtrax.csv"))
str(l0)

d0$SS <- with(d0, interaction(
    project.name,
    location,
    sep="::", drop=TRUE))
tmp1 <- as.character(d0$recording_date)
tmp2 <- as.character(d0$recording_time)
#tmp2[nchar(tmp2)<8] <- NA
d0$YEAR <- as.numeric(substr(tmp1, 1, 4))
d0$DATI <- as.POSIXlt(paste(tmp1, tmp2), format="%Y-%m-%d %H:%M:%OS", tz="America/Denver")
d0$PKEY <- as.factor(paste0(
    as.character(d0$SS),
    "_",
    tmp1,
    "-",
    tmp2))
d0$DATE <- as.Date(tmp1)
d0$SSYR <- paste0(d0$SS, "_", d0$YEAR)
d0$SPECIES <- normalize_sppcode(d0$species_code)
d0$X <- d0$longitude
d0$Y <- d0$latitude
d0$MAXDUR <- round(d0$recording.length/60)
d0$MAXDIS <- Inf

l0$SS <- with(l0, interaction(
    project.name,
    location,
    sep="::", drop=TRUE))
tmp1 <- as.character(l0$recording_date)
tmp2 <- as.character(l0$recording_time)
l0$PKEY <- as.factor(paste0(
    as.character(l0$SS),
    "_",
    tmp1,
    "-",
    tmp2))

d0$landcov <- l0$landcov[match(d0$PKEY, l0$PKEY)]
d0$treecov <- l0$treecov[match(d0$PKEY, l0$PKEY)]


d0$keep <- TRUE
d0$keep[is.na(d0$DATI)] <- FALSE
d0$keep[d0$abundance != "ONE"] <- FALSE
d0$keep[d0$MAXDUR > 20] <- FALSE

d <- droplevels(d0[d0$keep,])

d[1:25,c("PKEY","SPECIES","species.indiv.number")]
tmp <- paste(d$PKEY, d$SPECIES, d$species.indiv.number)
sum(duplicated(tmp))
d$ABUND <- 1

dd <- nonDuplicated(d, PKEY, TRUE)

dd$JULIAN <- as.POSIXlt(dd$DATE)$yday
dd$start <- as.POSIXlt(dd$DATI)$hour + as.POSIXlt(dd$DATI)$min / 60
keep <- is.na(dd$DATI) # these will be constant phi
keep[dd$JULIAN %[]% c(125, 200)] <- TRUE
keep[dd$start %[]% c(3, 12)] <- TRUE
dd <- droplevels(dd[keep,])
#' Normalized ordinal day
dd$JDAY <- dd$JULIAN / 365
#' Normalized time since local sunrise
Coor <- as.matrix(dd[,c("X", "Y")])
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
dd$srise <- NA
dd$srise[subset] <- sr
dd$TSSR <- (dd$start - dd$srise) / 24
#' Quadratic terms
dd$JDAY2 <- dd$JDAY^2
dd$TSSR2 <- dd$TSSR^2
#'
#' ## LCC/TREE
#'
## Reclass NALCMS
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
dd$NALC <- ltnalc$Label[match(dd$landcov, ltnalc$Value)]
dd$TREE <- dd$treecov / 100
dd$TREE[dd$TREE %)(% c(0, 1)] <- NA

dd$LCC4 <- as.character(dd$NALC)
dd$LCC4[dd$LCC4 %in% c("Decid", "Mixed")] <- "DecidMixed"
dd$LCC4[dd$LCC4 %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- "Open"
dd$LCC4 <- factor(dd$LCC4,
    c("DecidMixed", "Conif", "Open", "Wet"))
dd$LCC2 <- as.character(dd$LCC4)
dd$LCC2[dd$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
dd$LCC2[dd$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
dd$LCC2 <- factor(dd$LCC2, c("Forest", "OpenWet"))
table(dd$LCC4, dd$LCC2, useNA="a")
#' Maximum counting distance in 100 m units (area in ha)
dd$MAXDIS <- dd$MAXDIS / 100


y <- as.matrix(Xtab(ABUND ~ PKEY + SPECIES, d))
dim(y)


Xp <- cbind("(Intercept)"=1, as.matrix(dd[,c("TSSR","JDAY","TSSR2","JDAY2")]))
Xq <- cbind("(Intercept)"=1, TREE=dd$TREE,
    LCC2OpenWet=ifelse(dd$LCC2=="OpenWet", 1, 0),
    LCC4Conif=ifelse(dd$LCC4=="Conif", 1, 0),
    LCC4Open=ifelse(dd$LCC4=="Open", 1, 0),
    LCC4Wet=ifelse(dd$LCC4=="Wet", 1, 0))
summary(Xp)
summary(Xq)
#' Use version 3 of BAM QPAD estimates
library(QPAD)
load_BAM_QPAD(version=3)
sppp <- intersect(colnames(y), getBAMspecieslist())
#' Save values in a matrix
off <- matrix(NA, nrow(dd), length(sppp))
rownames(off) <- rownames(dd)
colnames(off) <- sppp
#' Loop for species
for (spp in sppp) {
    ## print out where we are
    cat(spp, "\n");flush.console()
    p <- rep(NA, nrow(dd))
    A <- q <- p
    ## constant for NA cases
    cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
    ## best model
    mi <- bestmodelBAMspecies(spp, type="BIC",
        model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    ## design matrices matching the coefs
    Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
    OKp <- rowSums(is.na(Xp2)) == 0
    Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
    OKq <- rowSums(is.na(Xq2)) == 0
    ## calculate p, q, and A based on constant phi and tau for the respective NAs
    p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
    unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
    A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
    q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))
    ## calculate time/lcc varying phi and tau for non-NA cases
    phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
    tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
    p[OKp] <- sra_fun(dd$MAXDUR[OKp], phi1)
    unlim <- ifelse(dd$MAXDIS[OKq] == Inf, TRUE, FALSE)
    A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dd$MAXDIS[OKq]^2)
    q[OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[OKq], tau1))
    ## log(0) is not a good thing, apply constant instead
    ii <- which(p == 0)
    p[ii] <- sra_fun(dd$MAXDUR[ii], cf0[1])
    ## store, next
    off[,spp] <- log(p) + log(A) + log(q)
}
## sanity checks
(Ra <- apply(off, 2, range))
summary(t(Ra))
which(!is.finite(Ra[1,]))
which(!is.finite(Ra[2,]))

#'
#' ## Save
#'
y <- y[rownames(dd), colnames(off)]
save(y, off, dd, file=file.path(ROOT, "nwt", "nwt-wildtrax-offsets-2019-01-16.RData"))


