ROOT <- "c:/bam/May2015"
library(mefa4)
library(pbapply)
load(file.path(ROOT, "out", "analysis_package_YYSS.Rdata"))

XYSS <- XYSS[order(XYSS[,"Xcl"], XYSS[,"Ycl"]),]
YYSS <- YYSS[rownames(XYSS),]

#spp <- "CAWA"

getDist <- function(spp) {
    xy1 <- XYSS[YYSS[,spp] > 0, c("Xcl","Ycl")]
    fd <- function(i) {
        xy <- XYSS[i, c("Xcl","Ycl")]
        min(sqrt((xy[1] - xy1[,1])^2 + (xy[2] - xy1[,2])^2)) / 1000
    }
    nd <- numeric(nrow(YYSS))
    i0 <- which(YYSS[,spp] == 0)
    d0 <- pbsapply(i0, fd)
    nd[i0] <- d0
    nd
}
#d_cawa <- getDist("CAWA")

#SPP <- colnames(YYSS)
SPP <- c("CAWA","OSFL","WEWP","RUBL")
d_all <- matrix(0, nrow(YYSS), length(SPP))
rownames(d_all) <- rownames(YYSS)
colnames(d_all) <- SPP

for (spp in SPP) {
    cat(spp, "\n");flush.console()
    d_all[,spp] <- getDist(spp)
}
save(d_all, file=file.path(ROOT, "out", "analysis_package_distances.Rdata"))

library(RColorBrewer)
col <- rev(brewer.pal(9, "Reds"))
z <- cut(d_cawa, c(-1, 0, 5, 10, 50, 100, 200, 500, 1000, Inf))
plot(XYSS[,c("Xcl","Ycl")], pch=".", col=col[z])
