## sparse-to-dense approximation
#library(ineq)

## currently this only includes longest ts at each location
## but there might be interruptions (e.g. 2 ts at same loc)
sparsity <- function(x, overlap=TRUE) {
#    x <- ifelse(x==0, 0L, 1L)
    z <- x
    z[] <- 0
    #x[,-1] * x[,-ncol(x)]
#    for (j in 1:nrow(x)) {
        for (t in 2:ncol(x)) {
#            z[j,t] <- sign(x[j,t-1]) * sign(x[j,t]) * (z[j,t-1]+1)
            z[,t] <- sign(x[,t-1]) * sign(x[,t]) * (z[,t-1]+1)
            #z[j,t] <- (z[j,t-1]+1)
        }
#    }
    v <- apply(z, 1, max)
    v_max <- max(v)

    vs <- sort(v[v>0])
    G <- sum(vs * 1:length(vs))
    G <- 2 * G/(length(vs) * sum(vs))
    G <- G - 1 - (1/length(vs))
    e <- 1 - G
    #e <- 1-Gini(Max[Max>0])
    if (is.na(e))
        e <- 0

    ts_len <- v_max / (ncol(x) - 1)

    is_ts <- sum(v > 0) / nrow(x)

    ol <- 0
    if (any(v>0) && overlap) {
        t <- z[v>0,]
        t[t>0] <- 1L
        v0 <- v[v>0]
        t0 <- t
        t0[] <- 0L
        for (j in 1:nrow(t))
            t0[j,1:v0[j]] <- 1L
        co <- t %*% t(t)
        sco <- sum(co[lower.tri(co)])
        co0 <- t0 %*% t(t0)
        sco0 <- sum(co0[lower.tri(co0)])
        ol <- sco / sco0
    }

    list(ts_mat=z, v=v, v_max=v_max, 
        evenness=e, proportion=is_ts, length=ts_len, overlap=ol,
        T_density=ts_len*e*is_ts, T_density2=ts_len*e*is_ts*ol)
}

plot.mat <- function(x, text="") {
    sp <- sparsity(x)
    xx <- 1 - t(x)
    image(xx, ylim=c(1+0.5/nrow(x),0-0.5/nrow(x)), 
        ann=FALSE, axes=FALSE, col=c("grey","white"))
    box()
    axis(side=3, at=seq(0,1,len=ncol(x)), labels=colnames(x))
    axis(side=2, at=seq(1,0,len=nrow(x)), labels=rev(rownames(x)), las=1)
    #title(sub=paste0("T-density = ", round(sp$T_density, 2)), xlab=text)
    txt <- paste0(text,"\n\nT-density = ", round(sp$T_density, 2),
        "\n\nproportion = ", round(sp$proportion, 2),
        "\nlength = ", round(sp$length, 2),
        "\nevenness = ", round(sp$evenness, 2),
        "\noverlap = ", round(sp$overlap, 2))
    text(0.5,0.5,txt)
    invisible(x)
}

## quantifying sparsity

set.seed(123)

m <- 10
T <- 5

## T-dense, S-dense
x1 <- matrix(1, 10, 5)
dimnames(x1) <- list(paste0("j=", 1:m), paste0("t=", 1:T))

## T-dense, S-sparse
x2 <- x1
x2[c(2,4,9),] <- 0

## T-sparse, S-dense
x3 <- x1
x3[,c(2,4)] <- 0

## T-sparse, S-sparse
x4 <- matrix(rbinom(50,1,0.5), 10, 5)
dimnames(x4) <- dimnames(x1)

par(mfrow=c(2,2))
plot.mat(x1, "T-dense, S-dense")
plot.mat(x2, "T-dense, S-sparse")
plot.mat(x3, "T-sparse, S-dense")
plot.mat(x4, "T-sparse, S-sparse")

x5 <- x4
x5[6:10,] <- 1

x6 <- x1-1
x6[1:5,1] <- 1
x6[1:5,3] <- 1
x6[1:5,5] <- 1
x6[6:10,2] <- 1
x6[6:10,4] <- 1

par(mfrow=c(1,2))
plot.mat(x5, "Anchor locations")
plot.mat(x6, "Rotating panel")

## BAM data 2014

load("c:/bam/Feb2014/data_ec_2014.Rdata")
library(mefa4)
d <- DAT[,c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND","YEAR",
    "BCR","PROV")]
rm(DAT, xy)
d$BBS <- substr(as.character(d$PCODE),1,3) == "BBS"
d1 <- droplevels(d[d$BBS,])
d2 <- droplevels(d[!d$BBS,])
d2$bcrprov <- interaction(d2$BCR, d2$PROV, drop=TRUE)

x1 <- Xtab(~SS+YEAR,d1)
x1[x1>0] <- 1
x2 <- Xtab(~SS+YEAR,d2)
x2[x2>0] <- 1

x3 <- Xtab(~bcrprov+YEAR,d2)
x3[x3>0] <- 1

sp1 <- sparsity(x1, overlap=FALSE)
sp2 <- sparsity(x2, overlap=FALSE)
sp3 <- sparsity(x3, overlap=FALSE)
str(sp1)
str(sp2)
str(sp3)


## BAM data 2015

load("c:/bam/May2015/out/data_package_2015-08-26.Rdata")
library(mefa4)
#d <- PKEY[,c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND","YEAR")]
isBBS <- substr(as.character(PKEY$PCODE),1,3) == "BBS"

d1 <- PKEY[PKEY$YEAR %in% c(2002:2012),]
d1$YEAR <- as.factor(d1$YEAR)
d2 <- d1[isBBS,]
d2$SS <- droplevels(d2$SS)
d3 <- d1[!isBBS,]
d3$SS <- droplevels(d3$SS)

x1 <- Xtab(~SS+YEAR,d1)
x1[x1>0] <- 1
sp1 <- sparsity(x1, overlap=FALSE)

x2 <- Xtab(~SS+YEAR,d2)
x2[x2>0] <- 1
sp2 <- sparsity(x2, overlap=FALSE)

x3 <- Xtab(~SS+YEAR,d3)
x3[x3>0] <- 1
sp3 <- sparsity(x3, overlap=FALSE)

str(sp1)
str(sp2)
str(sp3)

dim(x1)
dim(x2)
dim(x3)
levels(d1$YEAR)
