library(mefa4)

load("e:/peter/bam/bcr4/bcr4-data.RData")
source("~/repos/bamanalytics/projects/bcr4/bcr4-models.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-functions.R")
source("~/repos/mep/R/diagnostics-functions.R")

#setwd("your/path/here")
#load("bcr4-data.RData")
#source("bcr4-functions.R")
#source("diagnostics-functions.R")
#source("bcr4-models.R")

#m <- do_1spec1run(1, "OSFL", mods[c(4,1,2,3)], CAICalpha=0, return_best=TRUE)
#m <- do_1spec1run(1, "OSFL", mods[c(1,2,3,4)], CAICalpha=0, return_best=TRUE)
#summary(m$best)
#mep(m$best)

## here is how you can loop over species
## best model is saved only for the 1st run

B <- 25 # number of bootstrap runs you want
alpha <- 0

SPP2 <- colnames(YY[,colSums(YY[DAT$xBCR == 4,]>0) > 999])
SPP1 <- colnames(YY[,colSums(YY[DAT$xBCR == 4,]>0) > 99])

SPP <- setdiff(SPP1, SPP2)
#SPP <- c("OVEN", "OSFL") # comment this out for all species

for (i in SPP) {
    out <- list()
    for (j in 1:B) {
        cat("Species:", i, "- Run:", j, "\n")
        flush.console()
        out[[j]] <- try(do_1spec1run(j, i, mods, CAICalpha = alpha,
            return_best = j==1))
    }
    save(out, file=paste0("e:/peter/bam/bcr4/results/", i, ".RData"))
}
