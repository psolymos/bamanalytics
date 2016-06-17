library(mefa4)
library(detect)
library(QPAD)
load_BAM_QPAD(version=3)

d <- read.csv(paste0("c:/Users/Peter/Dropbox/bam/ARU/samuel/",
    "qryBBSRoadData_pairedDistance_PC_Solymos_VegInfo_2016-06-02.csv"))
sptab <- read.csv(paste0("c:/Users/Peter/Dropbox/bam/ARU/samuel/",
    "qryBBS_DirectedStudy_SpeciesList.csv"))

## classify habitats
levels(d$Additional)[levels(d$Additional) == "Coniter"] <- "Conifer"
table(d$Veg.Type..W.side.of.road., d$Additional)
## there are unknown veg cases: shelf the habitat effect

## BAM based EDR estimates (no habitat effect)
SPP_QPAD <- getBAMspecieslist()
compare_sets(sptab$SPECIES, SPP_QPAD)
EDRm1 <- 100*exp(sapply(sort(intersect(sptab$SPECIES, SPP_QPAD)),
    function(i) unname(coefBAMspecies(i)$edr)))

## estimate EDR for species that have data and part of spp list
compare_sets(sptab$SPECIES, d$SpeciesCode)
SPP <- sort(intersect(sptab$SPECIES, d$SpeciesCode))

d$PKEY <- with(d, interaction(Project, Cluster, Site, Station, sep=":", drop=TRUE))
table(d$DistanceDescrip, d$Distance)

table(d$BehaviourDescrip)
dd <- d[d$BehaviourDescrip %in% c("Calling, Sex unknown", "Counter-singing Male",
    "Singing male"),]

xt <- Xtab(Abundance ~ PKEY + DistanceDescrip + SpeciesCode, dd)
xt <- xt[SPP]
for (i in 1:length(xt))
    xt[[i]] <- xt[[i]][,c("0-50m", "50-100m", ">100m")]

D <- matrix(c(0.5, 1, Inf), nrow(xt[[1]]), 3, byrow=TRUE)

mods0 <- list()
for (spp in SPP) {
    Y <- as.matrix(xt[[spp]])
    mods0[[spp]] <- try(cmulti(Y|D ~ 1, type="dis"))
}

mods <- mods0[!sapply(mods0, inherits, "try-error")]
EDRm2 <- 100*exp(sapply(mods, coef))
names(EDRm2) <- names(mods)
n <- sapply(mods, nobs)

sptab$EDRm_BAM <- EDRm1[match(sptab$SPECIES, names(EDRm1))]
sptab$EDRm_BBS <- EDRm2[match(sptab$SPECIES, names(EDRm2))]
sptab$n <- n[match(sptab$SPECIES, names(n))]

LIM <- range(sptab$EDRm_BAM, sptab$EDRm_BBS, na.rm=TRUE)
with(sptab, plot(EDRm_BAM, EDRm_BBS, ylim=LIM, xlim=LIM,
    cex=log(n+1)))
abline(0,1, lty=2)
abline(lm(EDRm_BBS ~ EDRm_BAM, sptab, weights=sqrt(sptab$n)), col=2)

ss <- !is.na(sptab$n) & sptab$n >= 20
with(sptab[ss, ], plot(EDRm_BAM, EDRm_BBS, ylim=LIM, xlim=LIM,
    cex=log(n+1)))
abline(0,1, lty=2)
abline(lm(EDRm_BBS ~ EDRm_BAM, sptab[ss,], weights=sqrt(sptab$n[ss])), col=2)

write.csv(sptab, row.names=FALSE,
    file=paste0("c:/Users/Peter/Dropbox/bam/ARU/samuel/",
    "BBS_DirectedStudy_EDRestimates.csv"))

