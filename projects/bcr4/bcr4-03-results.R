library(mefa4)

source("~/repos/mep/R/diagnostics-functions.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-models.R")
#source("~/repos/bamanalytics/projects/bcr4/bcr4-models2.R")
source("~/repos/bamanalytics/projects/bcr4/bcr4-functions.R")
source("~/repos/bamanalytics/R/makingsense_functions.R")
load("e:/peter/bam/bcr4/bcr4-data.RData")

Terms <- getTerms(mods, "list")
setdiff(Terms, colnames(DAT))
xn <- DAT[BB[,1],Terms]
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))

spp <- "RCKI"
load(paste0("e:/peter/bam/bcr4/results/", spp, ".RData"))
100 * sum(getOK(out)) / length(out)
est1 <- getEst(out, stage = 2, X=Xn)
est2 <- getEst(out, stage = 4, X=Xn)
m <- out[[1]]$best
visualize_ysd2()

## explore the 1st run
summary(m)
#map(m)
boxplot(fitted(m) ~ hab, xn, range=0, ylab="Density")

## checking results
#getMid(out, mods)
getFancyMidTab(out, mods)

par(mfrow=c(1,3))
plotMid(out, mods)
visualize_road()
visualize_ysd()


SPP <- gsub("\\.RData", "", list.files("e:/peter/bam/bcr4/results/"))

pdf("e:/peter/bam/bcr4/results.pdf", onefile=TRUE, height=5, width=14)
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("e:/peter/bam/bcr4/results/", spp, ".RData"))
    est1 <- getEst(out, stage = 2, X=Xn)
    est2 <- getEst(out, stage = 4, X=Xn)
    par(mfrow=c(1,3))
    plotMid(out, mods)
    visualize_road()
    visualize_ysd2()
}
dev.off()


visualize_ysd2 <- function(use_climate=TRUE, ysd=0:100) {
    k <- nlevels(DAT$hab)
    Xn0 <- Xn[1:length(ysd),]
    rownames(Xn0) <- ysd
    Xn0[] <- 0
    Xn0[,1] <- 1
    if ("ysd" %in% colnames(Xn)) {
        Xn0[,"ysd"] <- pmin(100, ysd)
        Xn0[,"ysd05"] <- sqrt(Xn0[,"ysd"])
        Xn0[,"ysd2"] <- Xn0[,"ysd"]^2
        Xn0[,"ysd3"] <- Xn0[,"ysd3"]^2
    }
    qq <- list()
    tmp <- table(DAT$hab, DAT$isForest)
    est <- if(use_climate) est2 else est1
    for (i in 1:k) {
        Xn1 <- Xn0
        Xn1[,i] <- 1
        if (tmp[i,"1"] > 0) {
            Xn1[,"isForest:ysd"] <- ysd
            Xn1[,"isForest:ysd05"] <- sqrt(Xn1[,"isForest:ysd"])
            Xn1[,"isForest:ysd2"] <- Xn1[,"isForest:ysd"]^2
            Xn1[,"isForest:ysd3"] <- Xn1[,"isForest:ysd"]^3
        }
        pr <- apply(est, 1, function(z) Xn1 %*% z)
        qq[[i]] <- t(apply(exp(pr), 1, quantile, c(0.5, 0.05, 0.95)))
    }
    names(qq) <- levels(DAT$hab)

    mat <- sapply(qq, function(z) z[,1])
    matplot(mat, lty=1, type="l", lwd=2, ylim=c(0, 1.2*max(mat)),
        ylab="Relative abundance", xlab="Years since last disturbance",
        main=paste0(spp, " (B=", nrow(est1), ")"))
    legend("topright", lty=1, col=1:k, lwd=2, legend=names(qq))

    invisible(qq)
}
