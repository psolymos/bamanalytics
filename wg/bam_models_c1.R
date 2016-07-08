#### preliminaries ####

if (!interactive()) {
.Last <- function() {
    if (getOption("CLUSTER_ACTIVE")) {
        stopCluster(cl)
        cat("active cluster stopped by .Last\n")
    } else {
        cat("no active cluster found\n")
    }
}
options("CLUSTER_ACTIVE" = FALSE)
}
library(snow)
if (!interactive())
    library(Rmpi)
library(MASS)
library(ResourceSelection)
library(mefa4)
source("~/repos/bragging/R/glm_skeleton.R")
source("~/repos/bamanalytics/R/analysis_functions.R")
if (!interactive())
    (args <- commandArgs(trailingOnly = TRUE))

TEST <- interactive()

ROOT <- "~/bam"
CA <- 1 # 1=AIC, 0=BIC

#### setup ####

## arg1: nodes, arg2: species, arg3: text, arg4: sext, arg5: lctu
nodes <- if (interactive())
    5 else as.numeric(args[1])
BBB <- if (TEST) 2 else 240 # = 4*5*12
ncl <- if (TEST) 2 else nodes*12

#### load all object on the master ####

if (interactive())
    #setwd("e:/peter/bam/Apr2016/out")
    setwd(ROOT)

Date <- "2016-04-18"
fn <- paste0("pack_", Date, ".Rdata")
load(file.path("data", fn))
#mods <- mods[c("Clim", "Hab", "Road", "Hgt", "Dist", "Wet", "Year")]

#DAT$CMI2 <- DAT$CMI^2
#DAT$CMIJJA2 <- DAT$CMIJJA^2
mods$Clim <-list(
    . ~ . + CMIJJA + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP,
    . ~ . + CMI + DD0 + DD5 + EMT + MSP + DD02 + DD52 + CMI2 +
        CMI:DD0 + CMI:DD5 + EMT:MSP,
    . ~ . + CMI + CMIJJA + DD0 + MSP + TD + DD02 + CMI2 + CMIJJA2 +
        CMI:DD0 + CMIJJA:DD0 + MSP:TD,
    . ~ . + CMI + CMIJJA + DD5 + MSP + TD + DD52 + CMI2 + CMIJJA2 +
        CMI:DD5 + CMIJJA:DD5 + MSP:TD,
    . ~ . + CMIJJA + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + MSP:TD + MSP:EMT,
    . ~ . + CMI + DD0 + DD5 + EMT + TD + MSP + DD02 + DD52 + CMI2 +
        CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT,
    . ~ . + CMI + CMIJJA + DD0 + MSP + TD + EMT + DD02 + CMI2 + CMIJJA2 +
        CMI:DD0 + CMIJJA:DD0 + MSP:TD + MSP:EMT,
    . ~ . + CMI + CMIJJA + DD5 + MSP + TD + EMT + DD52 + CMI2 + CMIJJA2 +
        CMI:DD5 + CMIJJA:DD5 + MSP:TD + MSP:EMT)

if (TEST)
    mods <- mods[1:3]

#load(file.path("data", "analysis_package_distances.Rdata"))

#### spawning the slaves ####

cl <- if (interactive())
    makeCluster(ncl) else makeMPIcluster(ncl)
if (!interactive())
    options("CLUSTER_ACTIVE" = TRUE)

#### loading packages on slaves ####

tmpcl <- clusterEvalQ(cl, library(ResourceSelection))
tmpcl <- clusterEvalQ(cl, library(MASS))
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, source("~/repos/bragging/R/glm_skeleton.R"))

#### load all the objects on the slaves ####

tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("e:/peter/bam/Apr2016/out"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

#### project identifier ####

PROJECT <- if (TEST)
    "bam-test" else "bam"

spp <- if (interactive()) # CAWA OSFL RUBL WEWP OVEN MOWA
    "CAWA" else as.character(args[2])

#system.time(aaa <- do_1spec1run_noW(1, i=spp, mods=mods, hsh_name=NA, CAICalpha=CA))

res <- parLapply(cl, 1:BBB, do_1spec1run_noW, i=spp, mods=mods,
    hsh_name=NA, CAICalpha=CA)

fout <- paste0(PROJECT, "_", spp, "_", Date, "_c1.Rdata")
save(res, file=file.path("results", fout))


#### shutting down ####

stopCluster(cl)
if (!interactive()) {
    options("CLUSTER_ACTIVE" = FALSE)
    mpi.quit("no")
} else {
    quit("no")
}
