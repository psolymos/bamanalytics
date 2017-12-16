## mapping starts here -----------------------------

library(RColorBrewer)
library(mefa4)
ROOT <- "e:/peter/bam/Apr2016"
ROOT2 <- "e:/peter/bam/pred-2015"
ROOT3 <- "e:/peter/bam/pred-2016"

st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs <- levels(st$LEVEL3)

if (FALSE) {
## need to loop over chunks and save a proper XY where rownames match etc
## current load has BCR as NA etc
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs <- levels(st$LEVEL3)

load(file.path(ROOT3, "chunks", paste0("pgdat-", regs[1], ".Rdata")))
XY <- dat[,c("POINT_X","POINT_Y","BCR","JURS","LEVEL3","Brandt")]
for (regi in regs[-1]) {
    load(file.path(ROOT3, "chunks", paste0("pgdat-", regi, ".Rdata")))
    tmp <- dat[,c("POINT_X","POINT_Y","BCR","JURS","LEVEL3","Brandt")]
    XY <- rbind(XY, tmp)
    cat(dim(XY), "\n");flush.console()
}
XY$Brandt <- as.factor(XY$Brandt)
levels(XY$Brandt) <- c(levels(XY$Brandt), "OUT")
XY$Brandt[is.na(XY$Brandt)] <- "OUT"
save(XY, file=file.path(ROOT3, "allXY.Rdata"))
}

## pred grid
## pointid is rownames
load(file.path(ROOT3, "allXY.Rdata"))
XY$subreg <- as.factor(paste(XY$LEVEL3, XY$BCR, XY$JURS, XY$Brandt, sep=" + "))
st <- read.csv(file.path("e:/peter/bam/Apr2016", "BAMCECStudyAreaEcoregionLevel2.csv"))
regs <- levels(st$LEVEL3)
XY$studyarea <- XY$LEVEL3 %in% regs
#summary(XY)
gc()

## Subregions
XY3 <- nonDuplicated(XY, subreg, TRUE)
XY3$POINT_X <- NULL
XY3$POINT_Y <- NULL
rownames(XY3) <- XY3$subreg

source("~/repos/bamanalytics/R/makingsense_functions.R")

## observations
e <- new.env()
#load(file.path(ROOT, "out", "data", "pack_2016-04-18.Rdata"), envir=e)
load(file.path("e:/peter/bam/Apr2016/out", "data", "pack_2016-12-01.Rdata"), envir=e)
mods <- e$mods
Terms <- getTerms(e$mods, "list")
setdiff(Terms, colnames(e$DAT))
xn <- e$DAT[1:500,Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))
yy <- e$YY
xy_p <- e$DAT[,c("Xcl","Ycl")]
load(file.path("e:/peter/bam/Apr2016", "out", "SS-regions-and-xy.Rdata")) # SS01
SS <- SS01[match(e$DAT$SS, SS01$SS),]
jlt <- read.csv("~/repos/bamanalytics/lookup/jurisdictions.csv")
compare_sets(SS$JURSALPHA, jlt$JURS)
SS$JURS <- reclass(SS$JURSALPHA, jlt[,2:1])
SS$LEVEL3 <- SS$CECLEVEL3
SS$Brandt <- SS$BOREALLOC

SS$subreg <- as.factor(paste(SS$LEVEL3, SS$BCR, SS$JURS, SS$Brandt, sep=" + "))
SS$studyarea <- SS$LEVEL3 %in% regs
SS2 <- nonDuplicated(SS, SS, TRUE)[,c("PCODE","SS","X_CLCC","Y_CLCC","X_GEONAD83",
    "Y_GEONAD83","JURS","LEVEL3","Brandt","subreg","studyarea")]
yyss <- groupSums(yy, 1, SS$SS)
yyss[yyss>0] <- 1
mmss <- Mefa(yyss, SS2)
#mmss <- mmss[samp(mmss)$studyarea,]
summary(samp(mmss[samp(mmss)$studyarea,]))
yyl3 <- groupSums(xtab(mmss), 1, samp(mmss)$subreg)
#yyl3 <- yyl3[match(samp(mmss)$LEVEL3, rownames(yyl3)),]
#allSSbyL3 <- table(rep(1, nrow(mmss)), samp(mmss)$LEVEL3)
allSSbySubreg <- table(rep(1, nrow(mmss)), samp(mmss)$subreg)
rm(e)

## Nature Serve range map for CAWA
library(data.table)
getOption("datatable.fread.datatable")
options(datatable.fread.datatable=FALSE)
getOption("datatable.fread.datatable")

PROJECT <- "bam"
Date <- "2016-12-01"
Stage <- 6 # which(names(mods) == "Clim")
BASE_YEAR <- 2012
bfill <- FALSE
spp <- "CAWA"
fo <- paste0(spp, "-", Stage, "-", BASE_YEAR, "-", Date)
cat(fo, "\n");flush.console()

z <- fread("e:/peter/bam/Jan2017/AllPredCAWAallRange.csv")

load(file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "median-pred.Rdata")))
z$D <- x[match(z$PredPointsall_pointid, names(x))]
summary(z$D)
z <- z[!is.na(z$D),]

## range map
png(file.path(ROOT3, "maps", paste0(fo, "-range.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(z[z$PRESENCE==0,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(z$PredPointsall_POINT_X),
    ylim=range(z$PredPointsall_POINT_Y))
points(z[z$PRESENCE>0,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")], col = "tan", pch=".")
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-rangeD.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
Col <- rev(brewer.pal(6, "RdYlBu"))
br <- quantile(z$D, seq(0, 1, by=0.2))
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(z$D)) else as.integer(cut(z$D, breaks=br))
plot(z[,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE, xlim=range(z$PredPointsall_POINT_X),
    ylim=range(z$PredPointsall_POINT_Y))
par(op)
dev.off()

ff <- function(Dt) {
    x <- table(rangemap=z$PRESENCE, density=ifelse(z$D>=Dt, 1, 0))
    sum(diag(x))/sum(x)
}
oo <- optimize(ff, c(0, max(z$D)), maximum=TRUE)

th <- oo$maximum # 0.00724
-log(0.5)/th

sum(z$D >= th) / sum(z$D)

Dv <- exp(seq(-10, log(max(z$D)), len=100))
cmm <- lapply(Dv, function(Dt)
    table(rangemap=z$PRESENCE, density=ifelse(z$D>=Dt, 1, 0)))

a <- sapply(cmm, f)

library(opticut)
l <- lorenz(z$D)
l[which.min(abs(l[,"x"]-th)),]
pp <- t(sapply(Dv, function(i) l[which.min(abs(l[,"x"]-i)),]))

plot(Dv, a, type="l")
plot(-log(0.5)/Dv, a, type="l", xlim=c(0, 1000))
plot(pp[,"L"], a, type="l")
plot(Dv, pp[,"L"], type="l")


Dv <- seq(0, max(z$D), by=max(z$D)/25)
cmm <- lapply(Dv, function(Dt)
    table(rangemap=z$PRESENCE, density=ifelse(z$D>=Dt, 1, 0)))
cmm[[100]] <- cbind(cmm[[100]], "1"=c(0,0))
f <- function(x) sum(diag(x))/sum(x)
sens <- function(x) x[2,2]/sum(diag(x))
spec <- function(x) x[2,1]/sum(diag(t(x)))
a <- sapply(cmm, f)
ase <- sapply(cmm, sens)
asp <- sapply(cmm, spec)
which.max(a)

## todo
## - restrict it to Boreal
## - incorporate uncertainty

plot(Dv, a, type="l")

Dv2 <- seq(0.005087053, 0.010174106, len=5)
cmm2 <- lapply(Dv2, function(Dt)
    table(rangemap=z$PRESENCE, density=ifelse(z$D>=Dt, 1, 0)))
a2 <- sapply(cmm2, f)
which.max(a2)
cmm2[which.max(a2)]
a2[which.max(a2)]
Dv2[which.max(a2)] # 0.00763

plot(Dv2, a2, type="l")


png(file.path(ROOT3, "maps", paste0(fo, "-rangeA.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(z[z$D <= th,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(z$PredPointsall_POINT_X),
    ylim=range(z$PredPointsall_POINT_Y))
points(z[z$D > th,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")], col = "tan", pch=".")
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-rangeAD.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
plot(z[z$D <= th & z$PRESENCE==0,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(z$PredPointsall_POINT_X),
    ylim=range(z$PredPointsall_POINT_Y))
points(z[z$D > th & z$PRESENCE==1,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = Col[4], pch=".")
points(z[z$D <= th & z$PRESENCE==1,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = Col[1], pch=".")
points(z[z$D > th & z$PRESENCE==0,c("PredPointsall_POINT_X", "PredPointsall_POINT_Y")],
    col = Col[5], pch=".")
par(op)
dev.off()








PROJECT <- "bam"
Date <- "2016-12-01"

Stage <- 6 # which(names(mods) == "Clim")
BASE_YEAR <- 2012
bfill <- FALSE

ttt <- list()
brr <- list()

spp <- "CAWA"

gc()
fo <- paste0(spp, "-", Stage, "-", BASE_YEAR, "-", Date)
cat(fo, "\n");flush.console()

fl <- paste0(spp, "-", Stage, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), regs, "-", Date, ".Rdata")

is_null <- integer(length(fl))
names(is_null) <- fl
load(file.path(ROOT3, "species", spp, fl[1]))
if (is.null(lam)) {
    is_null[1] <- 1L
} else {
    plam <- lam
    tlam <- lam_total
}
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT3, "species", spp, fn))
    if (is.null(lam)) {
        is_null[fn] <- 1L
    } else {
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, lam_total)
    }
}
dim(plam)
sum(duplicated(rownames(plam)))

## already done in lamfun()
if (TRUE) {
    q <- quantile(plam[,"Mean"], 0.99)
    plam[plam[,"Mean"] > q,"Mean"] <- q

    q <- quantile(plam[,"Median"], 0.99)
    plam[plam[,"Median"] > q,"Median"] <- q
}

stopifnot(all(rownames(plam) == rownames(XY)))
gc()
#compare_sets(rownames(plam), rownames(XY2))

XY2all <- as.matrix(XY[,c("POINT_X","POINT_Y")])
## ideally, this is empty
#XY2miss <- XY[rownames(XY) %notin% rownames(plam),]
#XY2miss <- XY2miss[XY2miss$studyarea,]
#dim(XY2miss)

x <- plam[,"Median"]
save(x, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "median-pred.Rdata")))
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
br <- Lc_quantile(x, probs=probs, type="L")
if (!is.finite(br[length(br)]))
    br[length(br)] <- 1.01* max(x, na.rm=TRUE)
brr[[fo]] <- br
ttt[[fo]] <- tlam

## total million males
tlam0 <- tlam
tlam[!is.finite(tlam)] <- 0
for (i in 1:nrow(tlam)) {
    q <- quantile(tlam[i,], 0.99)
    tlam[i, tlam[i,] > q] <- q
}
tlam[!is.finite(tlam)] <- max(tlam[is.finite(tlam)])
summary(colSums(tlam))
fstat <- function(x, level=0.95) {
    c(Mean=mean(x), Median=median(x), quantile(x, c((1-level)/2, 1 - (1-level)/2)))
}
fstat(colSums(tlam/10^6), 0.9)
fstat(colSums(tlam0/10^6), 0.9)
## quick numbers
100*sum(plam[,"Mean"])/10^6
100*sum(plam[,"Median"])/10^6


XY3s <- droplevels(XY3[rownames(tlam),])
colnames(tlam) <- paste0(spp, "_run", 1:ncol(tlam))
XY3s$nSSinSubreg <- 0
for (i in colnames(allSSbySubreg))
    if (i %in% levels(XY3s$subreg))
        XY3s$nSSinSubreg[XY3s$subreg == i] <- allSSbySubreg[1,i]
XY3s$nDETinSubreg <- 0
for (i in colnames(allSSbySubreg))
    if (i %in% rownames(yyl3))
        XY3s$nDETinSubreg[XY3s$subreg == i] <- yyl3[i,spp]
ddd <- data.frame(XY3s, tlam)
rownames(ddd) <- rownames(tlam)

write.csv(ddd, row.names=FALSE, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "totals.csv")))
save(XY3s, tlam, file=file.path(ROOT3, "maps",
    paste0(fo, "-", BASE_YEAR, ifelse(bfill, "-bf-", "-"), "totals.Rdata")))


png(file.path(ROOT3, "maps", paste0(fo, "-mean-", BASE_YEAR, ifelse(bfill, "bf", ""), ".png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
Col <- rev(brewer.pal(6, "RdYlBu"))
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(x)) else as.integer(cut(x, breaks=br))
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
points(xy_p[yy[,spp] > 0,c("Xcl","Ycl")], pch=19, cex=0.1, col=1)
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "mean abundance"))
    #title=paste(spp, "median abundance"))
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-cov.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
TEXT <- paste0(br[-length(br)], "-", br[-1])
TEXT[length(TEXT)] <- paste0(">", br[length(br)-1])
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "darkgrey", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "SD / mean"))
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-sd.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
CoV <- plam[,"SD"] / mean(plam[,"Mean"])
zval <- cut(CoV, breaks=br)
br <- round(br*mean(plam[,"Mean"]), 4)
TEXT <- paste0(br[-length(br)], "-", br[-1], "%")
TEXT[length(TEXT)] <- paste0(">", br[length(br)-1], "%")
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "darkgrey", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "SD"))
par(op)
dev.off()




CAN <- c("ALBERTA", "BRITISH COLUMBIA", "MANITOBA",
    "NEW BRUNSWICK", "NEWFOUNDLAND",
    "NORTHWEST TERRITORIES", "NOVA SCOTIA", "NUNAVUT",
    "ONTARIO", "PRINCE EDWARD ISLAND", "QUEBEC", "SASKATCHEWAN", "YUKON")

## pop size in full study area
fstat(colSums(tlam)/10^6, 0.9)
## pop size in Canada
ss <- XY3s$JURS %in% CAN
fstat(colSums(tlam[ss,])/10^6, 0.9)
## pop size in Brandt boreal
ss <- XY3s$Brandt != "OUT"
fstat(colSums(tlam[ss,])/10^6, 0.9)
## pop size in Canada/Boreal
ss <- XY3s$JURS %in% CAN & XY3s$Brandt != "OUT"
fstat(colSums(tlam[ss,])/10^6, 0.9)

## by state/prov/terr
by_jurs <- data.frame(t(apply(groupSums(tlam/10^6, 1, XY3s$JURS), 1, fstat)))
by_jurs$perc <- by_jurs[,2] * 100 / sum(by_jurs[,2])
by_jurs <- by_jurs[order(by_jurs$perc),]
round(by_jurs, 4)
## by bcr
by_bcr <- data.frame(t(apply(groupSums(tlam/10^6, 1, XY3s$BCR), 1, fstat)))
by_bcr$perc <- by_bcr[,2] * 100 / sum(by_bcr[,2])
by_bcr <- by_bcr[order(by_bcr$perc),]
round(by_bcr, 4)

chfun <- function(Na, Nb, ta, tb) {
    100 * ((Nb/Na)^(1/(tb-ta)) - 1)
}
chfun(5.729, 5.704, 2002, 2012)

## difference maps

e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-bf-median-pred.Rdata", envir=e)
x0 <- e$x
e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2002-2016-12-01-2002-median-pred.Rdata", envir=e)
x1 <- e$x
e <- new.env()
load("e:/peter/bam/pred-2016/maps/CAWA-6-2012-2016-12-01-2012-median-pred.Rdata", envir=e)
x2 <- e$x
x1 <- x1[names(x2)]
x0 <- x0[names(x2)]

## thresholds for TF maps
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")

br0 <- Lc_quantile(x0, probs=probs, type="L")
if (!is.finite(br0[length(br0)]))
    br0[length(br0)] <- 1.01* max(x0, na.rm=TRUE)
br1 <- Lc_quantile(x1, probs=probs, type="L")
if (!is.finite(br1[length(br1)]))
    br1[length(br1)] <- 1.01* max(x1, na.rm=TRUE)
br2 <- Lc_quantile(x2, probs=probs, type="L")
if (!is.finite(br2[length(br2)]))
    br2[length(br2)] <- 1.01* max(x2, na.rm=TRUE)

dd <- data.frame(UpperCutpoint=probs, bf2012=br0, nobf2002=br1, nobf2012=br2)
dd[1,] <- 0
write.csv(dd, row.names=FALSE, file="CAWA-density-breaks.csv")

#pcawa <- data.frame(pointid=names(x2), median_2012=x2, median_2002=x1, median_backfilled=x0)
#write.csv(pcawa, row.names=FALSE, file="w:/bam-cawa/cawa-pred-med-2016-12-13.csv")

br <- c(0.05, 0.25, 0.5, 0.95, 1/rev(c(0.05, 0.25, 0.5, 0.95)))
d21 <- x2/x1
d20 <- x2/x0
c21 <- cut(d21, c(-1, br, Inf))
c20 <- cut(d20, c(-1, br, Inf))

if (FALSE) {

all(names(d21)==names(d20))
tow <- data.frame(pointid=names(d21),
    diff_2012div2002=d21,
    diff_2012div2012bf=d20,
    index_2012div2002=c21,
    index_2012div2012bf=c20)
tow$index_2012div2002 <- as.integer(tow$index_2012div2002)
tow$index_2012div2012bf <- as.integer(tow$index_2012div2012bf)
write.csv(tow, row.names=FALSE, file="w:/bam-cawa/difference-maps-20170106.csv")
}

png(file.path(ROOT3, "maps", paste0(fo, "-diff-2002-2012.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
Col <- brewer.pal(9, "RdYlGn")
zval <- as.integer(c21)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "tan", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
TEXT <- paste0(round(br[-length(br)],2), "-", round(br[-1],2), "x")
TEXT <- c(paste0("0-", br[1], "x"), TEXT, paste0(">", br[length(br)], "x"))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff 2002-2012"))
par(op)
dev.off()

png(file.path(ROOT3, "maps", paste0(fo, "-diff-bfill-2012.png")),
    width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
Col <- brewer.pal(9, "RdYlGn")
zval <- as.integer(c20)
plot(XY[!XY$studyarea,1:2], col = "lightgrey", pch=".",
    ann=FALSE, axes=FALSE, xlim=range(XY$POINT_X), ylim=range(XY$POINT_Y))
#points(XY2miss[,c("POINT_X","POINT_Y")], col = "tan", pch=".")
points(XY2all[,c("POINT_X","POINT_Y")], col = Col[zval], pch=".")
TEXT <- paste0(round(br[-length(br)],2), "-", round(br[-1],2), "x")
TEXT <- c(paste0("0-", br[1], "x"), TEXT, paste0(">", br[length(br)], "x"))
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    title=paste(spp, "diff bfill-2012"))
par(op)
dev.off()


## --

if (FALSE) {
br <- c(0, 0.4, 0.8, 1.2, 1.6, Inf)
Col <- rev(brewer.pal(5, "RdYlGn"))
TEXT <- paste0(100*br[-length(br)], "-", 100*br[-1], "%")
TEXT[length(TEXT)] <- paste0(">", 100*br[length(br)-1], "%")
#CoV <- plam[,"SD"] / plam[,"Mean"]
CoV <- plam[,"IQR"] / plam[,"Median"]
zval <- cut(CoV, breaks=br)
if (ids$SEXT[fid] == "nam") {
    plot(XYfull[rownames(plam),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
} else {
    plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
}
points(xy1, pch=19, cex=2)
legend("topright", bty = "n", legend=rev(TEXT),
    fill=rev(Col), border=1, cex=3,
    #title=paste(spp, "SD / mean"))
    title=paste(spp, "IQR / median"))
}

tlam <- data.frame(t(tlam))
write.csv(tlam, row.names=FALSE,
    file=file.path(ROOT2, "species", "cawa-nmbca-tabs", paste0("byregion-", fo, ".csv")))
plam <- data.frame(id=rownames(plam), median=plam[,"Median"], cov=CoV)
write.csv(plam, row.names=FALSE,
    file=file.path(ROOT2, "species", "cawa-nmbca-tabs", paste0("bypoint-", fo, ".csv")))


rm(plam)

## brr: list of Lc based breaks
## ttt: list of BCR/prov x B matrices
#save(brr, ttt, file=file.path(ROOT2, "species", "tlam-CAWA.Rdata"))


## visualizing climate

xx <- "int"
x <- DAT$CMIJJA * DAT$DD0
br <- unique(quantile(x, seq(0, 1, 0.2), na.rm=TRUE))
br[1] <- -Inf
br[length(br)] <- Inf
zval <- as.integer(cut(x, br))
COL <- rev(brewer.pal(length(br)-1, "RdYlBu"))

op <- par(mar=c(1,1,1,1)+0.1)
with(DAT, plot(Xcl, Ycl, pch=".", col=COL[zval], main=xx))
par(op)




load(file.path(ROOT2, "species", "tlam-CAWA.Rdata"))

f <- function(x) {
    q <- unname(quantile(x, c(0.5, 0.05, 0.95)))
    c(Mean=mean(x), SD=sd(x), Median=q[1], LCL90=q[2], UCL90=q[3], IQR=q[3]-q[2])
}
for (i in 1:length(ttt)) {

tmp <- t(apply(ttt[[i]], 1, f))
write.csv(tmp, row.names=TRUE,
    file=file.path(ROOT2, "species", "cawa-nmbca-tabs",
    paste0("summary-by-region-", names(ttt)[i], ".csv")))
}



## dealing with outliers
ttt2 <- ttt
for (i in 1:16) {
    tmp <- ttt[[i]]
    if (any(tmp > 10^10)) {
        for (j in 1:nrow(tmp)) {
            vv <- tmp[j,]
            vv <- vv[is.finite(vv)]
            q <- quantile(vv, 0.99)
            cat(j, q, "\n")
            tmp[j, tmp[j,] > q] <- q
        }
    }
    ttt2[[i]] <- tmp
}

t(sapply(ttt, function(z)
        c(mean=mean(colSums(z)/10^6), median=median(colSums(z)/10^6))))
t(sapply(ttt2, function(z)
        c(mean=mean(colSums(z)/10^6), median=median(colSums(z)/10^6))))


df <- data.frame(ids[rep(1:6, each=2),1:3],
    Year=c(2002, 2012),
    t(sapply(ttt[1:12], function(z)
        c(mean=mean(colSums(z)/10^6), median=median(colSums(z)/10^6)))))
df$change <- NA
for (i in 1:6) {
    N0 <- df$median[i*2-1]
    N10 <- df$median[i*2]
    df$change[i*2] <- 100 * ((N10/N0)^(1/10) - 1)
}

write.csv(df, file=file.path(ROOT2, "species", "popsize.csv"))

pe <- data.frame(ids[1:6,1:3],
    Year=BASE_YEAR,
    t(sapply(ttt, function(z)
        c(mean=mean(colSums(z)/10^6), median=median(colSums(z)/10^6)))))
save(ttt, file=file.path(ROOT, "out", "figs", "nmbca2",
    paste0(paste0("popsize-", spp, "-", Stage, "-", BASE_YEAR, "-", Date), ".Rdata")))

## 2003

pe03 <- structure(list(TEXT = structure(c(1L, 2L, 1L, 2L, 1L, 2L), .Label = c("gfw",
    "fre"), class = "factor"), SEXT = structure(c(1L, 1L, 1L, 1L,
    1L, 1L), .Label = c("can", "nam"), class = "factor"), LCTU = structure(c(1L,
    1L, 2L, 2L, 3L, 3L), .Label = c("nlc", "lcc", "eos"), class = "factor"),
        Year = c(2003, 2003, 2003, 2003, 2003, 2003), mean = c(7.04675889420143,
        8.1948604617475, 7.30885952161788, 8.12954353917993, 7.32127463952691,
        8.04625186575813), median = c(7.03908748814014, 7.90960556092234,
        7.28424257344771, 7.97434394721777, 7.29854788856675, 7.70893150271355
        )), .Names = c("TEXT", "SEXT", "LCTU", "Year", "mean", "median"
    ), row.names = c(NA, 6L), class = "data.frame")

pe13 <- structure(list(TEXT = structure(c(1L, 2L, 1L, 2L, 1L, 2L), .Label = c("gfw",
    "fre"), class = "factor"), SEXT = structure(c(1L, 1L, 1L, 1L,
    1L, 1L), .Label = c("can", "nam"), class = "factor"), LCTU = structure(c(1L,
    1L, 2L, 2L, 3L, 3L), .Label = c("nlc", "lcc", "eos"), class = "factor"),
        Year = c(2013, 2013, 2013, 2013, 2013, 2013), mean = c(7.03823934905092,
        8.2201686110189, 7.29573582539286, 8.15710339698905, 7.2890423740622,
        8.0574300369338), median = c(7.02512525548913, 7.93361330635701,
        7.27357416838527, 7.98576510700306, 7.2667175762312, 7.70826606796275
        )), .Names = c("TEXT", "SEXT", "LCTU", "Year", "mean", "median"
    ), row.names = c(NA, 6L), class = "data.frame")






## -- old

spp <- "CAWA"
wid <- 1
Stage <- 6
BASE_YEAR <- 2015
regs <- gsub(".Rdata", "",
    gsub("pgdat-", "", list.files(file.path(ROOT2, "chunks"))))

xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]

fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
load(file.path(ROOT2, "species", paste0(spp, "-ver3"), fl[1]))
plam <- lam
tlam <- attr(lam, "total")
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT2, "species", paste0(spp, "-ver3"), fn))
    plam <- rbind(plam, lam)
    tlam <- rbind(tlam, attr(lam, "total"))
}
dim(plam)
sum(duplicated(rownames(plam)))

hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)
mean(colSums(tlam)/10^6)


plotfun <- function(plam, what="Median") {
    z <- plam[,what]
    #z <- plam[,"IQR"]/plam[,"Median"]
    Col <- rev(brewer.pal(10, "RdBu"))
    cz <- cut(z, breaks = c(min(z)-1, quantile(z, seq(0.1, 1, 0.1))))
    plot(XYeosd[rownames(plam),], col = Col[cz], pch=".")
    invisible(NULL)
}

plotfun(z)

z <- dat1$HAB
Col <- brewer.pal(9, "Set3")
plot(dat[,2:3], col = Col[z], pch=".")

hs <- rowSums(pg4x4$nalc[rownames(dat0), c("Decid", "Mixed")])
plot(hs, aa)
table(unlist(lapply(sres$nalc, "[[", "hi")))







png(file.path(ROOT, "out", "figs", "CAWA-ver3-2015.png"), width = 2000, height = 2000)
op <- par(mfrow=c(2,1), mar=c(1,1,1,1)+0.1)

library(RColorBrewer)
x <- plam[,"Mean"]
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))
br <- Lc_quantile(x, probs=probs, type="L")
zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(x)) else as.integer(cut(x, breaks=br))
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)


br <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col <- rev(brewer.pal(10, "RdYlGn"))
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br)
plot(XYeosd[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=2)

par(op)
dev.off()

yrs <- c(2001, 2005, 2009, 2013)
gfw <- list()
for (BASE_YEAR in yrs) {
    fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
    load(file.path(ROOT2, "species", spp, fl[1]))
    plam <- lam
    tlam <- attr(lam, "total")
    for (fn in fl[-1]) {
        cat("loading", fn, "\n");flush.console()
        load(file.path(ROOT2, "species", spp, fn))
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, attr(lam, "total"))
    }
    gfw[[as.character(BASE_YEAR)]] <- tlam
}

hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)

ps <- sapply(gfw, function(z) {
    tmp <- colSums(z)/10^6
    c(Mean=mean(tmp), quantile(tmp, c(0.025, 0.975)))
    })

ps <- rbind(ps, Percent=(100 * ps[1,] / ps[1,1]) - 100)
ps <- rbind(ps, Decadal=10*ps[4,]/c(yrs-yrs[1]))


#spp <- "CAWA"

getDistPred <- function(spp) {
    xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]
    fd <- function(i) {
        xy <- XYeosd[i, ]
        min(sqrt((xy[1] - xy1[,1])^2 + (xy[2] - xy1[,2])^2)) / 1000
    }
    pbsapply(seq_len(nrow(XYeosd)), fd)
}

dp_cawa <- getDistPred("CAWA")
save(dp_cawa, file=file.path(ROOT2, "dp-cawa.Rdata"))


png(file.path(ROOT, "out", "figs", "CAWA-d-2015.png"), width = 2000, height = 1000)
op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)

col <- rev(brewer.pal(9, "Reds"))
z <- cut(dp_cawa, c(-1, 0, 5, 10, 50, 100, 200, 500, 1000, Inf))
plot(XYeosd, pch=".", col=col[z],
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)

par(op)
dev.off()


spp <- "CAWA"
Stage <- 7
BASE_YEAR <- 2015
est <- list()
map <- list()
for (wid in c(0,2,3)) {
    fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
    load(file.path(ROOT2, "species", spp, fl[1]))
    plam <- lam
    tlam <- attr(lam, "total")
    for (fn in fl[-1]) {
        cat("loading", fn, "\n");flush.console()
        load(file.path(ROOT2, "species", spp, fn))
        plam <- rbind(plam, lam)
        tlam <- rbind(tlam, attr(lam, "total"))
    }
    est[[as.character(wid)]] <- tlam
    map[[as.character(wid)]] <- plam
}

ps <- sapply(est, function(z) {
    tmp <- colSums(z)/10^6
    c(Mean=mean(tmp), quantile(tmp, c(0.025, 0.975)))
    })
colnames(ps) <- c("NALC","EOSD","LCC")
ps

png(file.path(ROOT, "out", "figs", "CAWA-all-7-2015.png"), width = 2000, height = 3000)
op <- par(mfrow=c(3,1), mar=c(1,1,1,1)+0.1)

for (i in 1:3) {
    x <- map[[i]][,"Mean"]
    probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
    TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
    Col <- rev(brewer.pal(5, "RdYlBu"))
    br <- Lc_quantile(x, probs=probs, type="L")
    zval <- if (length(unique(round(br,10))) < 5)
        rep(1, length(x)) else as.integer(cut(x, breaks=br))
    plot(XYeosd[rownames(map[[i]]),], col = Col[zval], pch=".",
        ann=FALSE, axes=FALSE)
    points(xy1, pch=19, cex=1.2)
}

par(op)
dev.off()


library(RColorBrewer)
ROOT <- "c:/bam/May2015"
load("e:/peter/bam/pred-2015/pg-clim.Rdata")
load("e:/peter/bam/pred-2015/pg-loss.Rdata")
clim$YearFire <- loss$YearFire[match(clim$pointid, loss$pointid)]
clim$YearLoss <- loss$YearLoss[match(clim$pointid, loss$pointid)]

for (fn in c("CMIJJA", "CMI", "TD", "DD0",
    "DD5", "EMT", "MSP", "CTI", "SLP")) {

    png(file.path(ROOT, "out", "figs", paste0("x-", fn,".png")), width = 2000, height = 1000)
    op <- par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)

    Col <- brewer.pal(5, "OrRd")
    z <- cut(clim[,fn], breaks=quantile(clim[,fn], seq(0,1,0.2), na.rm=TRUE))
    plot(clim[,c("POINT_X","POINT_Y")], col = Col[z], pch=".",
        ann=FALSE, axes=FALSE)
    title(main=fn)
    legend("bottomleft", bty = "n", legend=rev(levels(z)), fill=rev(Col))


    par(op)
    dev.off()

}


## mapping, fid = 4

library(RColorBrewer)
library(mefa4)
library(pbapply)
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/pred-2015"
source("~/repos/bamanalytics/R/makingsense_functions.R")
source("~/repos/bamanalytics/R/analysis_mods.R")
load(file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))
load(file.path(ROOT2, "XYfull.Rdata"))

spp <- "CAWA"
wid <- 1
Stage <- 6
BASE_YEAR <- 2015
regs <- gsub(".Rdata", "",
    gsub("pgdat-full-", "", list.files(file.path(ROOT2, "chunks2"))))

xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]

fl <- paste0(paste(spp, wid, Stage, BASE_YEAR, regs, sep="-"), ".Rdata")
load(file.path(ROOT2, "species", paste0(spp, "-4"), fl[1]))
plam <- lam
tlam <- attr(lam, "total")
for (fn in fl[-1]) {
    cat("loading", fn, "\n");flush.console()
    load(file.path(ROOT2, "species", paste0(spp, "-4"), fn))
    plam <- rbind(plam, lam)
    tlam <- rbind(tlam, attr(lam, "total"))
}
dim(plam)
sum(duplicated(rownames(plam)))


hist(colSums(tlam)/10^6, col=3)
sum(tlam[,1])/10^6
median(colSums(tlam)/10^6)
mean(colSums(tlam)/10^6)

z <- ifelse(rownames(XYfull) %in% names(Brandt), 2, 1)
plot(XYfull, col = z, pch=".", ann=FALSE, axes=FALSE)

z <- plam[,"Median"]
z[z >= 100] <- max(z[z < 100])
probs <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
TEXT <- paste0(100*probs[-length(probs)], "-", 100*probs[-1], "%")
Col <- rev(brewer.pal(5, "RdYlBu"))
br <- Lc_quantile(z, probs=probs, type="L")

png(file.path(ROOT, "out", "figs", "CAWAfid4-2015.png"), width = 2000, height = 2000)
op <- par(mfrow=c(2,1), mar=c(1,1,1,1)+0.1)

zval <- if (length(unique(round(br,10))) < 5)
    rep(1, length(z)) else as.integer(cut(z, breaks=br))
plot(XYfull[rownames(plam),], col = Col[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)


br2 <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, Inf)
Col2 <- rev(brewer.pal(10, "RdYlGn"))
CoV <- plam[,"SD"] / plam[,"Mean"]
zval <- cut(CoV, breaks=br2)
plot(XYfull[rownames(plam),], col = Col2[zval], pch=".",
    ann=FALSE, axes=FALSE)
points(xy1, pch=19, cex=1.2)

par(op)
dev.off()

## marginal plots

ROOT <- "e:/peter/bam/pred-2015"
library(mefa4)

load(file.path(ROOT, "pg-main.Rdata"))
x <- x[x$EOSD_COVER == 1,]
rownames(x) <- x$pointid

x <- x[rownames(plam),]
load(file.path(ROOT, "pg-loss.Rdata"))
ii <- loss$YearFire >= 9000 & !is.na(loss$YearFire)
loss$YearFire[ii] <- loss$YearFire[ii] - 8000
x$YearFire <- loss$YearFire[match(x$pointid, loss$pointid)]
x$YearLoss <- loss$YearLoss[match(x$pointid, loss$pointid)]

i <- sample.int(nrow(x), 5000)
tc <- c(rgb(1,0,0,alpha=0.2), rgb(0,0,1,alpha=0.2), rgb(0,1,0,alpha=0.2), rgb(0,0,0,alpha=0.2))
j <- as.integer(x$HAB_NALC2[i])
j[] <- 4
j[x$HAB_NALC2[i] == "Decid"] <- 1
j[x$HAB_NALC2[i] == "Mixed"] <- 2
j[x$HAB_NALC2[i] == "Conif"] <- 3

boxplot(plam[i,"Mean"] ~ x$HAB_NALC2[i], col="gold", range=0, main="Habitat", ylab="D")
boxplot(plam[i,"Mean"] ~ x$TR3[i], col="gold", range=0, main="Tree", ylab="D")
plot(plam[i,"Mean"] ~ jitter(x$HGT[i]), pch=19, cex=1, col=tc[j], main="Height", ylab="D")
legend("topleft", pch=19, col=tc, legend=c("Dec","Mix","Con","Else"))

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(x$LIN[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(x$POL[i]), pch=19, cex=1, col=tc[j], ylab="D")


load(file.path(ROOT, "pg-clim.Rdata"))
rownames(clim) <- clim$pointid
clim <- clim[match(x$pointid, clim$pointid),4:14]
clim <- clim[rownames(plam),]

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(clim$CTI[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(clim$SLP[i]), pch=19, cex=1, col=tc[j], ylab="D")

par(mfrow=c(2,1))
plot(plam[i,"Mean"] ~ jitter(2015-x$YearFire[i]), pch=19, cex=1, col=tc[j], ylab="D")
plot(plam[i,"Mean"] ~ jitter(2015-x$YearLoss[i]), pch=19, cex=1, col=tc[j], ylab="D")
