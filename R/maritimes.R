##---
##title: "Maritimes data analysis"
##author: "Peter Solymos"
##date: "May 13, 2015"
##output: pdf_document
##---

### Data processing

#### Processing geospatial data for sampling locations

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
library(mefa4)
## Load data files
## use new files
x50 <- read.csv(file.path(ROOT, "maritimes/June19", "avian_50.csv"))
x100 <- read.csv(file.path(ROOT, "maritimes/June19", "avian_100.csv"))
x250 <- read.csv(file.path(ROOT, "maritimes/June19", "avian_250.csv"))
## but CAS data is in old files
x50v <- read.csv(file.path(ROOT, "maritimes", "Avian_50.csv"))
x100v <- read.csv(file.path(ROOT, "maritimes", "Avian_100.csv"))
x250v <- read.csv(file.path(ROOT, "maritimes", "Avian_250.csv"))
x50v <- x50v[match(x50$OBJECTID, x50v$OBJECTID),]
x100v <- x100v[match(x100$OBJECTID, x100v$OBJECTID),]
x250v <- x250v[match(x250$OBJECTID, x250v$OBJECTID),]
stopifnot(all(x50$OBJECTID==x50v$OBJECTID))
stopifnot(all(x100$OBJECTID==x100v$OBJECTID))
stopifnot(all(x250$OBJECTID==x250v$OBJECTID))
x50 <- data.frame(x50, x50v[,setdiff(colnames(x50v), colnames(x50))])
x100 <- data.frame(x100, x100v[,setdiff(colnames(x100v), colnames(x100))])
x250 <- data.frame(x250, x250v[,setdiff(colnames(x250v), colnames(x250))])
rm(x50v,x100v,x250v)
## Checking duplicates
x50[duplicated(x50$SS),1:6]
x100[duplicated(x100$SS),1:6]
x250[duplicated(x250$SS),1:6]
## Dropping duplicates and using SS as row name
x50 <- nonDuplicated(x50, SS, TRUE)
x100 <- nonDuplicated(x100, SS, TRUE)
x250 <- nonDuplicated(x250, SS, TRUE)
## Checking if set of SS values are the same
stopifnot(length(union(rownames(x50), rownames(x100))) ==
    length(intersect(rownames(x50), rownames(x100))))
stopifnot(length(union(rownames(x50), rownames(x250))) ==
    length(intersect(rownames(x50), rownames(x250))))
stopifnot(length(union(rownames(x100), rownames(x250))) ==
    length(intersect(rownames(x100), rownames(x250))))
## Sort the tables by SS
x50 <- x50[rownames(x250),]
x100 <- x100[rownames(x250),]

## CAS based tree species
tree_250 <- c("Abie_bals", "Abie_pice", "Abie_spp", "Acer_rubr",
    "Acer_sacc", "Acer_spp", "Alnu_spp", "Betu_alle", "Betu_papy",
    "Betu_popu", "Betu_spp", "Fagu_gran", "Frax_spp", "Hard_into",
    "Hard_nonc", "Hard_tole", "Hard_unkn", "Lari_deci", "Lari_lari",
    "NOSC_HARD", "NOSC_SOFT", "Pice_abie", "Pice_glau", "Pice_mari",
    "Pice_rube", "Pice_spp", "Pinu_bank", "Pinu_resi", "Pinu_spp",
    "Pinu_stro", "Pinu_sylv", "Popu_balb", "Popu_spp", "Popu_trem",
    "Prun_sero", "Quer_rubr", "Soft_unkn", "Thuj_occi", "Tsug_cana",
    "Ulmu_amer", "Uncl_spp")
## Checking if 50 and 100 m buffer has no tree species
## that that is unique for a <250 m scale (this should not happen, and seems OK)
setdiff(colnames(x50), tree_250)
setdiff(colnames(x100), tree_250)
## CAS related fields from the 250 m buffer
cas_250 <- c(tree_250, "CRCL_AV", "CRCL_STD", "HT_AV", "HT_STD")
## Rows with 0 values for all CAS fields should be excluded
## but only after checking the species data
## (we do not want to throw it away if this is where species live)
x250$keep <- rowSums(x250[,cas_250]) > 0 # drop 1260 rows
## Checking 250 m based percentages for leading tree species
xp <- as.matrix(x250[,tree_250])
range(xp)
## Divide areas (m^2) by buffer size (250^2*pi = 196349.5)
xp <- xp / (250^2*pi)
stopifnot(all(xp <= 1))
## Percentages sometimes add up to >1, but this should not affect
## finding the dominant species (maximum do not change by standardization,
## but we might want to know percentage for comparisons)
hist(rowSums(xp))
summary(rowSums(xp))
sum(rowSums(xp) > 1)
## Standardize with row sums where it is >1
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
summary(rowSums(xp))
## Leading species
tlead <- as.factor(tree_250[apply(xp, 1, which.max)])
data.frame(sort(table(tlead)))
## Percent cover of leading species
plead <- apply(xp, 1, function(z) z[which.max(z)])

## Leading human footprint type in 250 m buffer
hflab <- c("BU","CO","OT","PC","SI","WF")
h <- as.matrix(x250[,colnames(x250) %in% hflab])
hist(rowSums(h))
h <- cbind(Undist=1-rowSums(h),
    groupSums(h, 2, c("Burn","Cut","Other","Cut","Other","Other")))
x250$ldist <- factor(colnames(h)[apply(h, 1, which.max)], colnames(h))
table(x250$ldist)

## Standardization
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "No Data"] <- "Average" # 4 cases
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Below Average"] <- "-1" # 4 cases
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Average"] <- "0" # 4 cases
levels(x250$COMPLEXITY)[levels(x250$COMPLEXITY) == "Above Average"] <- "1" # 4 cases
x250$COMPLEXITY <- as.integer(as.character(x250$COMPLEXITY))

x250$CONNECTEDNESS <- x250$CONNECTEDNESS / 100

x50$CROWNCL_AV <- x50$CROWNCL_AV / 100
x50$CROWNCL_STD <- x50$CROWNCL_STD / 100
x100$CROWNCL_AV <- x100$CROWNCL_AV / 100
x100$CROWNCL_STD <- x100$CROWNCL_STD / 100
x250$CRCL_AV <- x250$CRCL_AV / 100
x250$CRCL_STD <- x250$CRCL_STD / 100

x250$DTW_STD <- x250$DTW_STD / 10

x250$HT_AV <- x250$HT_AV / 10
x250$HT_STD <- x250$HT_STD / 10

x250$HUMAN_FOOTPRINT <- x250$HUMAN_FOOTPRINT / 100

x50$ROAD01 <- ifelse(x50$ROAD <= 400, 1L, 0L)
x100$ROAD01 <- ifelse(x100$ROAD <= 400, 1L, 0L)

x250$WET_LENGTH <- x250$WET_LENGTH / 1000

levels(x50$WET_VEG)[levels(x50$WET_VEG) == "WATER/EXPOSED"] <- "AQUATIC"
levels(x50$WET_VEG)[levels(x50$WET_VEG) %in%
    c("GRAMINOI","GRAMINOID","GRAMNINOID")] <- "GRAMINOID"
levels(x100$WET_VEG)[levels(x100$WET_VEG) == "WATER/EXPOSED"] <- "AQUATIC"
levels(x100$WET_VEG)[levels(x100$WET_VEG) %in%
    c("GRAMINOI","GRAMINOID","GRAMNINOID")] <- "GRAMINOID"

## Reclass leading tree species
lt <- read.csv(file.path(ROOT, "maritimes", "cas-tree-lookup.csv"))
rownames(lt) <- lt$spp
lt <- lt[tree_250,]
xp2 <- groupSums(xp, 2, lt$spp3)
xp2 <- xp2[,colnames(xp2) != "XXX"]
data.frame(tree=sort(table(as.factor(colnames(xp2)[apply(xp2, 1, which.max)]))))
x250$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], c(colnames(xp2), "NOCAS"))
x250$ltree[!x250$keep] <- "NOCAS" # no tree data
## Reclass at 50 m scale
xp <- as.matrix(x50[,colnames(x50) %in% tree_250])
xp <- xp / (50^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x50$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
x50$ltree[!x250$keep] <- "NOCAS" # no tree data
## Reclass at 100 m scale
xp <- as.matrix(x100[,colnames(x100) %in% tree_250])
xp <- xp / (100^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x100$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
x100$ltree[!x250$keep] <- "NOCAS" # no tree data
## Compare scales
table(x50$ltree, x250$ltree)
table(x100$ltree, x250$ltree)
## Get rid of tree cover info
x50 <- x50[,!(colnames(x50) %in% tree_250)]
x100 <- x100[,!(colnames(x100) %in% tree_250)]
x250 <- x250[,!(colnames(x250) %in% tree_250)]
## Define fields that do not change with buffer
asis <- c("OBJECTID", "PCODE", "SS", "POINT_X", "POINT_Y")
## Renaming columns in 50 and 100 m tables to reflect
## 'local' scale (as opposed to 'territory' scale et 250 m)
colnames(x50)[!(colnames(x50) %in% asis)] <- paste("LOC",
    colnames(x50)[!(colnames(x50) %in% asis)], sep="_")
colnames(x100)[!(colnames(x100) %in% asis)] <- paste("LOC",
    colnames(x100)[!(colnames(x100) %in% asis)], sep="_")
xx50 <- data.frame(x50, x250[,!(colnames(x250) %in% asis)])
xx100 <- data.frame(x100, x250[,!(colnames(x250) %in% asis)])


#### Pull in BAM counts and PKEY table

e <- new.env()
load(file.path(ROOT, "out", "data_package_2015-06-19.Rdata"),
    envir=e)
names(as.list(e))
pc <- e$PCTBL
xt <- as.matrix(Xtab(ABUND ~ PKEY + SPECIES, e$PCTBL)[,c("CAWA","RUBL","OSFL")])
pk <- e$PKEY
rm(e)
## Load pre-calculated offsets (QPAD v2)
e <- new.env()
load(file.path(ROOT, "out", "offsets_allspp_BAMBBS_2015-06-19.Rdata"),
    envir=e)
names(as.list(e))
off <- e$OFF[,c("CAWA","RUBL","OSFL")]
rm(e)
## Match SS with buffered Maritimes data sets
dim(pk)
pk <- droplevels(pk[pk$SS %in% rownames(xx50),])
rownames(pk) <- pk$PKEY
dim(pk)
## Inner join for count and PKEY table
pk <- pk[rownames(pk) %in% rownames(xt),]
dim(pk)
## Before and after subsetting the counts
colSums(xt)
apply(xt, 2, table)
xt <- xt[rownames(pk),]
colSums(xt)
apply(xt, 2, table)
## Subsetting offsets
off <- off[rownames(pk),]

#### Packages for each species and save

## use 50 m for CAWA and RUBL
## use 100 m for OSFL
dcawa <- data.frame(
    spp=xt[,"CAWA"],
    off=off[,"CAWA"],
    pk[,c("SITE","STN","ROUND","YEAR","ROAD")],
    xx50[match(pk$SS, xx50$SS),])
drubl <- data.frame(
    spp=xt[,"RUBL"],
    off=off[,"RUBL"],
    pk[,c("SITE","STN","ROUND","YEAR","ROAD")],
    xx50[match(pk$SS, xx50$SS),])
dosfl <- data.frame(
    spp=xt[,"OSFL"],
    off=off[,"OSFL"],
    pk[,c("SITE","STN","ROUND","YEAR","ROAD")],
    xx100[match(pk$SS, xx50$SS),])
## Bootstrap matrix
source("~/repos/detect/R/hbootindex.R")
set.seed(1234)
BB <- hbootindex(groups=dcawa$STN, strata=dcawa$SITE, B=240-1)
## Save the data to be used later
save(dcawa, drubl, dosfl, BB,
    file=file.path(ROOT, "out", "maritimes_3spp.Rdata"))


### Analysis

#### Read in the data set

ROOT <- "c:/bam/May2015"
library(mefa4)
#load(file.path(ROOT, "out", "maritimes_3spp.Rdata"))
load("~/Dropbox/bam/maritimes2015/maritimes_3spp.Rdata")
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/bamanalytics/R/maritimes_mods.R")

## Check if all variables defined in the
## 3 sets of model lists can be found in data

TermsA <- getTerms(modsA, "list")
TermsA <- c("STN", TermsA)
setdiff(TermsA, colnames(dcawa))

TermsB <- getTerms(modsB, "list")
TermsB <- c("STN", TermsB)
setdiff(TermsB, colnames(dcawa))

TermsC <- getTerms(modsC, "list")
TermsC <- c("STN", TermsC)
setdiff(TermsC, colnames(dcawa))

Terms <- unique(c(TermsA, TermsB, TermsC))

## Summary of variables used in modeling
summary(dcawa[,sort(Terms)])

## Model lists
modsA
modsB
modsC

## Settings to test the function below
if (FALSE) {
B <- ncol(BB) - 1
j=1
i="CAWA"
use_wt=TRUE
CAICalpha=1
silent=FALSE
nmax=NULL
}

#### Defining the core functions

## This function is the workhorse.
## Slightly adapted for the Maritimes analysis
do_1spec1run_mar <- function(j, i, mods,
silent=FALSE, use_wt=TRUE,
CAICalpha=0.5, nmax=NULL)
{
    x <- switch(i,
        "CAWA" = dcawa,
        "RUBL" = drubl,
        "OSFL" = dosfl)
    x <- x[BB[,j],]

    y <- x$spp
    off0 <- x$off
    ## spatial weights
    if (use_wt) {
        tmp <- Xtab(~ x$STN + rownames(x), drop.unused.levels=TRUE)
        w0 <- rowSums(tmp)[match(x$STN, rownames(tmp))]
    } else {
        w0 <- rep(1L, length(y))
    }
    if (!is.null(nmax)) {
        if (nmax > length(y))
            stop("nmax > length(y)")
        ss <- sample.int(length(y), nmax, replace=FALSE, prob=1/w0)
        x <- x[ss,]
        x$STN <- droplevels(x$STN)
        y <- y[ss]
        off0 <- off0[ss]
        tmp <- Xtab(~ x$STN + rownames(x), drop.unused.levels=TRUE)
        w0 <- rowSums(tmp)[match(x$STN, rownames(tmp))]
    }
    w0 <- 1/sqrt(w0)
    ## empty objects for storing results
    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
    ## Null
    null <- glm_skeleton(glm(y ~ 1,
        x,
        family=poisson(),
        offset=off0,
        weights=w0,
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    ## Lc tangent is not used
    ip_name <- NULL
    Hi <- NULL
    lam <- NULL
    cv <- NULL
    habmod <- NULL
    ## looping through models list
    for (l1 in 1:nmods) {
        if (nnmods[l1] > 0) {
            mlist <- vector("list", nnmods[l1])
            glist <- vector("list", nnmods[l1])
            for (l2 in 1:nnmods[l1]) {
                mlist[[l2]] <- glm_skeleton(try(update(object=best,
                    formula=mods[[l1]][[l2]]), silent=silent), CAICalpha=CAICalpha)
            }
            mcaic <- sapply(mlist, "[[", "caic")
            attr(mcaic, "StartCAIC") <- best$caic
            for (l2 in 1:length(mlist)) { # check convergence
                if (mlist[[l2]]$class != "try-error" && !mlist[[l2]]$converge)
                    mcaic[l2] <- 2*.Machine$double.xmax^(1/3)
            }
            dcaic <- mcaic - best$caic
            mmid <- which.min(dcaic)
            if (dcaic[mmid] < 0) {
                best <- mlist[[mmid]]
                mid[l1] <- mmid
                gofbest <- glist[[mmid]]
            }
            caicList[[l1]] <- mcaic
        }
        bestList[[l1]] <- best
    }
    ## final assembly
    out <- list(species=i, iteration=j,
        null=null$coef,
        null_caic=null$caic,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
        hi=Hi,
        lc=cv,
        alpha=CAICalpha,
        nmax=nmax,
        use_wt=use_wt,
        habmod=habmod$coef,
        ip_name=ip_name)
    out
}

## This is used to estimate computing times
#system.time(tmp <- do_1spec1run_mar(j=2, i=i, mods=modsA,
#    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL))
#system.time(tmp <- do_1spec1run_mar(j=2, i=i, mods=modsB,
#    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL))
#system.time(tmp <- do_1spec1run_mar(j=2, i=i, mods=modsC,
#    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL))

## This function runs all bootstrap iterations for a species,
## results are returned and optionally saved
wg_fun_mar <- function(i, mods, B,
    output=c("return","rdata","dump"),
    project="", path="", verbose=TRUE, ...)
{
    output <- match.arg(output)
    res <- vector("list", B + 1)
    t0 <- date()
    for (j in seq_len(B+1)) {
        if (verbose) {
            cat(project, i, j, "\n")
            flush.console()
        }
        gc()
        res[[j]] <- try(do_1spec1run_mar(j=j, i=i, mods=mods, ...))
    }
    t1 <- date()
    attr(res, "start_stop") <- c(start=t0, stop=t1)
    if (output != "return") {
        ext <- switch(output,
            "rdata"=".Rdata",
            "dump"=".Rdmp")
        if (project != "")
            project <- paste(project, "_", sep="")
        file <- file.path(path, paste("BirdCoefs_", project,
            i, "_Allrun", ext, sep=""))
        if (output == "rdata")
            save(res, file=file)
        if (output == "dump")
            dump("res", file=file)
        invisible(TRUE)
    }
    res
}

#### Estimation

if (FALSE) {
#B <- 1
B <- ncol(BB) - 1
Output <- "rdata"
Path <- "~/Dropbox/bam/maritimes2015"

res_cawa_A <- wg_fun_mar(i="CAWA", mods=modsA, B=B,
    output=Output, project="MaritimesA", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_cawa_B <- wg_fun_mar(i="CAWA", mods=modsB, B=B,
    output=Output, project="MaritimesB", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_cawa_C <- wg_fun_mar(i="CAWA", mods=modsC, B=B,
    output=Output, project="MaritimesC", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)

res_rubl_A <- wg_fun_mar(i="RUBL", mods=modsA, B=B,
    output=Output, project="MaritimesA", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_rubl_B <- wg_fun_mar(i="RUBL", mods=modsB, B=B,
    output=Output, project="MaritimesB", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_rubl_C <- wg_fun_mar(i="RUBL", mods=modsC, B=B,
    output=Output, project="MaritimesC", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)

res_osfl_A <- wg_fun_mar(i="OSFL", mods=modsA, B=B,
    output=Output, project="MaritimesA", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_osfl_B <- wg_fun_mar(i="OSFL", mods=modsB, B=B,
    output=Output, project="MaritimesB", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
res_osfl_C <- wg_fun_mar(i="OSFL", mods=modsC, B=B,
    output=Output, project="MaritimesC", path=Path,
    silent=FALSE, use_wt=TRUE, CAICalpha=1, nmax=NULL)
}

### Making sense of the results

#### Loading the data and results

load("~/Dropbox/bam/maritimes2015/maritimes_3spp.Rdata")
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/bamanalytics/R/maritimes_mods.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")

load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesA_CAWA_Allrun.Rdata"))
cawa <- list(A=res)
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesB_CAWA_Allrun.Rdata"))
cawa$B <- res
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesC_CAWA_Allrun.Rdata"))
cawa$C <- res

load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesA_RUBL_Allrun.Rdata"))
rubl <- list(A=res)
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesB_RUBL_Allrun.Rdata"))
rubl$B <- res
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesC_RUBL_Allrun.Rdata"))
rubl$C <- res

load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesA_OSFL_Allrun.Rdata"))
osfl <- list(A=res)
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesB_OSFL_Allrun.Rdata"))
osfl$B <- res
load(file.path("~/Dropbox/bam/maritimes2015",
    "BirdCoefs_MaritimesC_OSFL_Allrun.Rdata"))
osfl$C <- res

allres <- list(CAWA=cawa, RUBL=rubl, OSFL=osfl)

#### Comparing A-B-C subsets based on AIC

caic_trajectory <- function(res) {
    f <- function(x) {
        out <- numeric(length(x$mid)+1)
        out[1] <- x$null_caic
        for (i in 1:length(x$mid)) {
            out[i+1] <- if (x$mid[i] > 0)
                x$caic[[i]][x$mid[i]] else out[i]
        }
        out
    }
    t(sapply(res, f))
}


tmp <- lapply(allres$CAWA, caic_trajectory)
tmp <- cbind(A=tmp$A[,9], B=tmp$B[,9], C=tmp$C[,8])
table(colnames(tmp)[apply(tmp, 1, which.min)])

tmp <- lapply(allres$RUBL, caic_trajectory)
tmp <- cbind(A=tmp$A[,9], B=tmp$B[,9], C=tmp$C[,8])
table(colnames(tmp)[apply(tmp, 1, which.min)])

tmp <- lapply(allres$OSFL, caic_trajectory)
tmp <- cbind(A=tmp$A[,9], B=tmp$B[,9], C=tmp$C[,8])
table(colnames(tmp)[apply(tmp, 1, which.min)])


## keep track if Bid, 0+mid vec, aic vec
## adjust level 1&2 for ABC (8 lev for all)
## selectively check sign of terms
## work on design matrix and prediction
## (once pred data comes in -- fast)

f <- function(x) {
        aic <- numeric(length(x$mid)+1)
        aic[1] <- x$null_caic
        for (i in 1:length(x$mid)) {
            aic[i+1] <- if (x$mid[i] > 0)
                x$caic[[i]][x$mid[i]] else aic[i]
        }
        aic
}
fun1 <- function(res) {
    aicA <- t(sapply(res$A, f))[,-2]
    aicB <- t(sapply(res$B, f))[,-2]
    aicC <- t(sapply(res$C, f))
    midA <- t(sapply(res$A, function(z) c(0, z$mid)))
    midA[] <- as.character(midA)
    midA[,3] <- paste(midA[,2], midA[,3], sep="&")
    midA <- midA[,-2]
    midB <- t(sapply(res$B, function(z) c(0, z$mid)))
    midB[] <- as.character(midB)
    midB[,3] <- paste(midB[,2], midB[,3], sep="&")
    midB <- midB[,-2]
    midC <- t(sapply(res$C, function(z) c(0, z$mid)))
    midC[] <- as.character(midC)
    aic <- cbind(A=aicA[,8], B=aicB[,8], C=aicC[,8])
    i <- colnames(aic)[apply(aic, 1, which.min)]
    mid <- t(sapply(1:length(i), function(j)
        switch(i[j],
               "A"=paste("A", midA[j,], sep="-"),
               "B"=paste("B", midB[j,], sep="-"),
               "C"=paste("C", midC[j,], sep="-"))))
    colnames(mid) <- paste0("Stage_", c(0, 1.2, 3:8))
    mid <- mid[,-1]
    colnames(mid)[1] <- "Stage_1&2"
    apply(mid, 2, function(z)
        data.frame(Freq=rev(sort(table(z))),
            Perc=round(100*rev(sort(table(z)))/nrow(mid), 1)))
}

aiccawa <- do.call(cbind, lapply(allres$CAWA, function(z) {
    tmp <- sapply(z, f)
    tmp[nrow(tmp),]
    }))
aicrubl <- do.call(cbind, lapply(allres$RUBL, function(z) {
    tmp <- sapply(z, f)
    tmp[nrow(tmp),]
    }))
aicosfl <- do.call(cbind, lapply(allres$OSFL, function(z) {
    tmp <- sapply(z, f)
    tmp[nrow(tmp),]
    }))
rownames(aiccawa) <- rownames(aicrubl) <- rownames(aicosfl) <- paste0("run", 1:240)

frcawa <- do.call(rbind, fun1(allres$CAWA))
frrubl <- do.call(rbind, fun1(allres$RUBL))
frosfl <- do.call(rbind, fun1(allres$OSFL))

xnA <- model.frame(getTerms(modsA, "formula"), dcawa,
    na.action=na.pass)
XnA <- model.matrix(getTerms(modsA, "formula"), xnA)
colnames(XnA) <- fixNames(colnames(XnA))
Xn <- XnA
cfcawaA <- getEst(allres$CAWA$A)
cfrublA <- getEst(allres$CAWA$A)
cfosflA <- getEst(allres$CAWA$A)

xnB <- model.frame(getTerms(modsB, "formula"), dcawa,
    na.action=na.pass)
XnB <- model.matrix(getTerms(modsB, "formula"), xnB)
colnames(XnB) <- fixNames(colnames(XnB))
Xn <- XnB
cfcawaB <- getEst(allres$CAWA$B)
cfrublB <- getEst(allres$CAWA$B)
cfosflB <- getEst(allres$CAWA$B)

xnC <- model.frame(getTerms(modsC, "formula"), dcawa,
    na.action=na.pass)
XnC <- model.matrix(getTerms(modsC, "formula"), xnC)
colnames(XnC) <- fixNames(colnames(XnC))
Xn <- XnC
cfcawaC <- getEst(allres$CAWA$C)
cfrublC <- getEst(allres$CAWA$C)
cfosflC <- getEst(allres$CAWA$C)

mid <- list(data.frame(spp="CAWA", subset="A", getFancyMidTab(allres$CAWA$A, modsA)),
    data.frame(spp="CAWA", subset="B", getFancyMidTab(allres$CAWA$B, modsB)),
    data.frame(spp="CAWA", subset="C", getFancyMidTab(allres$CAWA$C, modsC)),
    data.frame(spp="RUBL", subset="A", getFancyMidTab(allres$RUBL$A, modsA)),
    data.frame(spp="RUBL", subset="B", getFancyMidTab(allres$RUBL$B, modsB)),
    data.frame(spp="RUBL", subset="C", getFancyMidTab(allres$RUBL$C, modsC)),
    data.frame(spp="OSFL", subset="A", getFancyMidTab(allres$OSFL$A, modsA)),
    data.frame(spp="OSFL", subset="B", getFancyMidTab(allres$OSFL$B, modsB)),
    data.frame(spp="OSFL", subset="C", getFancyMidTab(allres$OSFL$C, modsC)))

mid <- data.frame(do.call(rbind, mid))

write.csv(frcawa, file="~/Dropbox/bam/maritimes2015/res-cawa.csv")
write.csv(frrubl, file="~/Dropbox/bam/maritimes2015/res-rubl.csv")
write.csv(frosfl, file="~/Dropbox/bam/maritimes2015/res-osfl.csv")

write.csv(mid, file="~/Dropbox/bam/maritimes2015/res-all.csv")

write.csv(aiccawa, file="~/Dropbox/bam/maritimes2015/res-cawa-aic.csv")
write.csv(aicrubl, file="~/Dropbox/bam/maritimes2015/res-rubl-aic.csv")
write.csv(aicosfl, file="~/Dropbox/bam/maritimes2015/res-osfl-aic.csv")

write.csv(cfcawaA, file="~/Dropbox/bam/maritimes2015/res-cawa-coef-A.csv")
write.csv(cfcawaB, file="~/Dropbox/bam/maritimes2015/res-cawa-coef-B.csv")
write.csv(cfcawaC, file="~/Dropbox/bam/maritimes2015/res-cawa-coef-C.csv")
write.csv(cfrublA, file="~/Dropbox/bam/maritimes2015/res-rubl-coef-A.csv")
write.csv(cfrublB, file="~/Dropbox/bam/maritimes2015/res-rubl-coef-B.csv")
write.csv(cfrublC, file="~/Dropbox/bam/maritimes2015/res-rubl-coef-C.csv")
write.csv(cfosflA, file="~/Dropbox/bam/maritimes2015/res-osfl-coef-A.csv")
write.csv(cfosflB, file="~/Dropbox/bam/maritimes2015/res-osfl-coef-B.csv")
write.csv(cfosflC, file="~/Dropbox/bam/maritimes2015/res-osfl-coef-C.csv")

### MID figures

pdf("~/Dropbox/bam/maritimes2015/mid-figures.pdf", onefile=TRUE)
plotMid(allres$CAWA$A, modsA, web=TRUE, main="CAWA, A")
plotMid(allres$CAWA$B, modsB, web=TRUE, main="CAWA, B")
plotMid(allres$CAWA$C, modsC, web=TRUE, main="CAWA, C")
plotMid(allres$RUBL$A, modsA, web=TRUE, main="RUBL, A")
plotMid(allres$RUBL$B, modsB, web=TRUE, main="RUBL, B")
plotMid(allres$RUBL$C, modsC, web=TRUE, main="RUBL, C")
plotMid(allres$OSFL$A, modsA, web=TRUE, main="OSFL, A")
plotMid(allres$OSFL$B, modsB, web=TRUE, main="OSFL, B")
plotMid(allres$OSFL$C, modsC, web=TRUE, main="OSFL, C")
dev.off()

pdf("~/Dropbox/bam/maritimes2015/mid-figures-x3.pdf",
    height=5, width=15, onefile=TRUE)
op <- par(mfrow=c(1,3))
plotMid(allres$CAWA$A, modsA, web=TRUE, main="CAWA, A")
plotMid(allres$CAWA$B, modsB, web=TRUE, main="CAWA, B")
plotMid(allres$CAWA$C, modsC, web=TRUE, main="CAWA, C")
plotMid(allres$RUBL$A, modsA, web=TRUE, main="RUBL, A")
plotMid(allres$RUBL$B, modsB, web=TRUE, main="RUBL, B")
plotMid(allres$RUBL$C, modsC, web=TRUE, main="RUBL, C")
plotMid(allres$OSFL$A, modsA, web=TRUE, main="OSFL, A")
plotMid(allres$OSFL$B, modsB, web=TRUE, main="OSFL, B")
plotMid(allres$OSFL$C, modsC, web=TRUE, main="OSFL, C")
par(op)
dev.off()


### Session info

sessionInfo()

## This little utility is used to turn this plain text file into
## a human readable documentation
#mefa4::r2rmd("~/repos/bamanalytics/R/maritimes.R")

