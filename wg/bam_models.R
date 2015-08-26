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

#### setup ####

nodes <- if (interactive())
    5 else as.numeric(args[1])
BBB <- if (TEST) 2 else 240 # = 4*5*12
ncl <- if (TEST) 2 else nodes*12

#### load all object on the master ####

if (interactive())
    setwd("c:/bam/May2015/out")
fid <- if (interactive())
    1 else as.numeric(args[2])
fl <- paste0(c("analysis_package_gfwfire-nalc-",
    "analysis_package_gfwfire-eosd-",
    "analysis_package_gfwfire-lcc-",
    "analysis_package_fire-nalc-"), "2015-08-26.Rdata")
fn <- fl[fid]
load(file.path("data", fn))
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
    tmpcl <- clusterEvalQ(cl, setwd("c:/bam/May2015/out"))
tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))

#tmpcl <- clusterEvalQ(cl, load(file.path("data", "analysis_package_distances.Rdata")))

#### project identifier ####

PROJECT <- if (TEST)
    paste0("bam-", fid, "-test") else paste0("bam-", fid) 


#### checkpoint ####
if (FALSE) { # this is for multi species runs -------------- !!!
DONE <- substr(sapply(strsplit(list.files("results"), "_"), "[[", 3), 1, 4)
SPP <- setdiff(SPP, DONE)
if (TEST)
    SPP <- SPP[1:2]

for (SPP1 in SPP) {
    cat(SPP1, date(), "\n")
    res <- parLapply(cl, 1:BBB, wg_fun2, i=SPP1, mods=mods, 
        output="return", path="results", project=PROJECT, 
        use_wt=TRUE, ip_name=ip_name, CAICalpha=CAICalpha, nmax=nmax)
    #res <- wg_fun2(1, i=SPP[1], mods=mods, 
    #    output="return", path="results", project=PROJECT, 
    #    use_wt=TRUE, ip_name=ip_name, CAICalpha=CAICalpha, nmax=nmax)
    save(res, file=paste("results/birds_", PROJECT, "_", SPP1, ".Rdata", sep=""))
}
}

#hsh_name <- "HAB"
hsh_name <- NA # no landscape level effects
CAICalpha <- 1
spp <- if (interactive()) # CAWA OSFL RUBL WEWP
    "CAWA" else as.character(args[3])

#DAT$ND2 <- -(d_all[match(DAT$SS, rownames(d_all)),spp] / 1000)^2
#nd2 <- DAT$ND2
#tmpcl <- clusterExport(cl, "nd2")
##clusterEvalQ(cl, summary(DAT$ND2))
#tmpcl <- clusterEvalQ(cl, DAT$ND2 <- nd2)
##clusterEvalQ(cl, summary(DAT$ND2))

#system.time(aaa <- do_1spec1run_noW(1, i=spp, mods=mods, hsh_name=hsh_name, CAICalpha=CAICalpha))

res <- parLapply(cl, 1:BBB, do_1spec1run_noW, i=spp, mods=mods, 
    hsh_name=hsh_name, CAICalpha=CAICalpha)
attr(res, "fid") <- fid
attr(res, "spp") <- spp
attr(res, "hsh_name") <- hsh_name
attr(res, "CAICalpha") <- CAICalpha

fout <- paste0(PROJECT, "_", spp, ".Rdata", sep="")
save(res, file=file.path("results", fout))


#### shutting down ####

stopCluster(cl)
if (!interactive()) {
    options("CLUSTER_ACTIVE" = FALSE)
    mpi.quit("no")
} else {
    quit("no")
}
