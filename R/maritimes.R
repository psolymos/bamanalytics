##---
##title: "Maritimes data analysis"
##author: "Peter Solymos"
##date: "May 13, 2015"
##output: pdf_document
##---

### Processing geospatial data for sampling locations

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
library(mefa4)
## Load data files
x50 <- read.csv(file.path(ROOT, "maritimes", "Avian_50.csv"))
x100 <- read.csv(file.path(ROOT, "maritimes", "Avian_100.csv"))
x250 <- read.csv(file.path(ROOT, "maritimes", "Avian_250.csv"))
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

## Reclass leading tree species
lt <- read.csv(file.path(ROOT, "maritimes", "cas-tree-lookup.csv"))
rownames(lt) <- lt$spp
lt <- lt[tree_250,]
xp2 <- groupSums(xp, 2, lt$spp3)
xp2 <- xp2[,colnames(xp2) != "XXX"]
data.frame(tree=sort(table(as.factor(colnames(xp2)[apply(xp2, 1, which.max)]))))
x250$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], colnames(xp2))
## Reclass at 50 m scale
xp <- as.matrix(x50[,colnames(x50) %in% tree_250])
xp <- xp / (50^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x50$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
## Reclass at 100 m scale
xp <- as.matrix(x100[,colnames(x100) %in% tree_250])
xp <- xp / (100^2*pi)
xp <- xp / ifelse(rowSums(xp) > 1, rowSums(xp), 1)
xp2 <- groupSums(xp, 2, lt[colnames(xp), "spp3"])
xp2 <- xp2[,colnames(xp2) != "XXX"]
x100$ltree <- factor(colnames(xp2)[apply(xp2, 1, which.max)], levels(x250$ltree))
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

## Pull in BAM counts and PKEY table
e <- new.env()
load(file.path(ROOT, "out", "data_package_2015-05-14.Rdata"),
    envir=e)
names(as.list(e))
pc <- e$PCTBL
pc <- pc[pc$keep,]
xt <- as.matrix(Xtab(ABUND ~ PKEY + SPECIES, pc)[,c("CAWA","RUBL","OSFL")])
pk <- e$PKEY

## data checks
if (FALSE) {
pp <- data.frame(e$PKEY, e$SS[match(e$PKEY$SS, e$SS$SS),])
rownames(pp) <- pp$PKEY
ii <- intersect(rownames(xt),rownames(pp))
xt <- xt[ii,]
pp <- pp[ii,]

table(pp$JURS, xt[,"CAWA"]>0)
with(pp, plot(X, Y, pch="."))
}

rm(e)
## Load pre-calculated offsets (QPAD v2)
e <- new.env()
load(file.path(ROOT, "out", "offsets_allspp_BAMBBS_2015-05-05.Rdata"),
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
    xx50[match(pk$SS, xx50$SS),])

with(dcawa, table(spp, keep))
with(drubl, table(spp, keep))
with(dosfl, table(spp, keep))

save(dcawa, drubl, dosfl, file=file.path(ROOT, "out", "maritimes_3spp.Rdata"))
