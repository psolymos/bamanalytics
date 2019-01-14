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
    "WR")
xybu <- droplevels(xybu[xybu$ProjectID %in% PCODE_keep, ])
d0 <- droplevels(d0[d0$ProjectID %in% PCODE_keep, ])


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
d$DATI <- as.POSIXlt(paste(tmp1, tmp2))
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
save(pk, file=file.path(ROOT, "nwt", "nwt-pkey-table-2019-01-14.RData"))

## --- I am here



y_bu <- Xtab(Abundance ~ PKEY + SPECIES, d)
d_bu <- make_char2fact(droplevels(nonDuplicated(d, PKEY, TRUE)))
d_bu$VKEY <- d_bu$SS
