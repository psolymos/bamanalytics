library(mefa4)
ROOT <- "c:/bam/May2015"
source("~/repos/bragging/R/glm_skeleton.R")

fid <- 1
fl <- c("analysis_package_gfwfire-nalc-2015-07-24.Rdata",
    "analysis_package_gfwfire-eosd-2015-07-24.Rdata",
    "analysis_package_gfwfire-lcc-2015-07-24.Rdata",
    "analysis_package_fire-nalc-2015-07-24.Rdata")
load(file.path(ROOT, "out", fl[fid]))



if (FALSE) {
B <- ncol(BB) - 1
j=1
i="CAWA"
w_id="gridcode"
CAICalpha=1
hsh_name="HAB"
silent=FALSE
nmax=25000
}


do_1spec1run <- function(j, i, mods, 
silent=FALSE, w_id=NA, 
hsh_name=NA, CAICalpha=0.5, nmax=NULL) 
{
    select_hsh <- !is.na(hsh_name)
    use_wt <- !is.na(w_id)
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    if (select_hsh)
        hsh <- HSH[BB[,j],]
    off <- OFF[BB[,j], i]
    ## spatial weights
    if (use_wt) {
        tmp <- Xtab(~ x[[w_id]] + rownames(x), drop.unused.levels=TRUE)
        w <- rowSums(tmp)[match(x[[w_id]], rownames(tmp))]
    } else {
        w <- rep(1L, length(y))
    }
    if (!is.null(nmax)) {
        if (nmax > length(y))
            stop("nmax > length(y)")
        ss <- sample.int(length(y), nmax, replace=FALSE, prob=1/w)
        x <- x[ss,]
        x[[w_id]] <- droplevels(x[[w_id]])
        y <- y[ss]
        off <- off[ss]
        if (select_hsh)
            hsh <- hsh[ss,]
        tmp <- Xtab(~ x[[w_id]] + rownames(x), drop.unused.levels=TRUE)
        w <- rowSums(tmp)[match(x[[w_id]], rownames(tmp))]
    }
    w <- 1/sqrt(w)
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
        offset=off, 
        weights=w,
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    ## Lorenz-tangent approach for core habitat delineation
    if (select_hsh) {
        HABV <- x[,hsh_name]
        habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD,
            x,
            family=poisson(), 
            offset=off, 
            weights=w,
            x=FALSE, y=FALSE, model=FALSE), silent=silent), CAICalpha=CAICalpha)
        ## need to correct for linear effects
        ## so that we estimate potential pop in habitats (and not realized)
        XHSH <- model.matrix(~ HABV + ROAD, x)
        XHSH[,"ROAD"] <- 0 # not predicting edge effects
        ## some levels might be dropped (e.g. Marsh)
        XHSH <- XHSH[,names(habmod$coef)]
        lam <- exp(drop(XHSH %*% habmod$coef))
        cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
        Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
        Prob <- Freq[,"1"] / rowSums(Freq)
        ## missing/dropped levels are NaN=0/0
        Prob[is.na(Prob)] <- 0
        Hi <- names(Prob)[Prob > 0.5]
        #tb <- ifelse(Freq > 0, 1, 0)
        #Hi <- rownames(tb)[tb[,"1"] > 0]
        #Lo <- rownames(tb)[tb[,"1"] == 0]
        x$HSH <- unname(rowSums(hsh[, Hi, drop=FALSE]))
        x$HSH2 <- x$HSH^2
    } else {
        Hi <- NULL
        lam <- NULL
        cv <- NULL
        habmod <- NULL
    }
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
        w_id=w_id,
        habmod=habmod$coef,
        hsh_name=hsh_name)
    out
}


system.time(res <- do_1spec1run(1, "CAWA", mods, 
    w_id="gridcode", hsh_name="HAB", CAICalpha=1, nmax=NULL))
data.frame(id=structure(res$mid,names=names(mods)))







## ---- packages: save data packages for separate runs

## hab=HAB31 = NALCMS (all NAm)
## hab=HAB12 = LCC05v2 (Can)
## hab=HAB22 = EOSD (Can Boreal)
## hab=HAB41 = CASFRI (Can Boreal+?, CAS area)

library(mefa4)
load("c:/bam/Feb2014/data_ec_2014.Rdata")
if (FALSE) {
    DAT <- DAT[sample(1:nrow(DAT), 10^4),]
    boxplot(HEIGHT ~ HAB31, DAT)
ABdat <- droplevels(DAT[!is.na(DAT$PROV) & DAT$PROV=="AB",
    c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND","YEAR")])

}
DAT0 <- DAT
source("c:/Dropbox/bam/EC_contract_2014/BAM-EC_02_analysis_functions.R")
source("c:/Dropbox/bam/EC_contract_2014/BAM-EC_03_making_sense_functions.R")

hab_list <- c("HAB12","HAB12","HAB22","HAB22","HAB31","HAB31","HAB41","HAB41","HAB31")
ClimBCR_list <- c("BCR", "Clim", "BCR", "Clim", "BCR", "Clim", "BCR", "Clim", "Clim")
DoCell_list <- c(rep("Pt", 8), "Cl")
#hab <- "HAB31" # HAB12, HAB22, HAB31, HAB41
#ClimBCR <- Clim or BCR

for (zz in 1:9) {

    hab <- hab_list[zz]
    ClimBCR <- ClimBCR_list[zz]
    DoCell <- DoCell_list[zz]

    ## reset data because of BCR filtering
    DAT <- DAT0

    cat(hab, ClimBCR, DoCell, "\n");flush.console()
    source("c:/Dropbox/bam/EC_contract_2014/BAM-EC_02_analysis_mods.R")

    ## SPATIAL switch: TRUE=BCR interaction + exlude USeast (& Arctic for CAS)
    if (ClimBCR == "BCR") {
        DAT$xBCR[DAT$xBCR == "USeast"] <- NA
        if (hab == "HAB41")
            DAT$xBCR[DAT$xBCR == "Arctic"] <- NA
    }
    print(dim(DAT));flush.console()

    NA_ACTION <- na.exclude # na.pass
    Terms <- getTerms(c(mods, m0=mod_0, m1=mod_1), "list")
    xn <- model.frame(getTerms(c(mods, m0=mod_0, m1=mod_1), "formula"), DAT[,Terms], 
        na.action=NA_ACTION)
    print(dim(xn));flush.console()
    cat("\n\n")

    xn$HAB <- xn[,hab]
    xn$isDecid <- ifelse(xn$HAB == "Decid", 1L, 0L)
    xn$isMixed <- ifelse(xn$HAB == "Mixed", 1L, 0L)
    #xn$isForest <- ifelse(xn$HAB %in% c("Conif","Decid","Mixed"), 1L, 0L)
    xn$isNonforest <- ifelse(xn$HAB %in% c("Conif","Decid","Mixed"), 0L, 1L)
    xn[,hab] <- NULL

    #xn <- xn[1:200,]
    #save(xn, mods, Terms, 
    #    file=paste0("c:/Dropbox/bam/EC_contract_2014/res2/Xn_", hab, "_", 
    #    ClimBCR, "_", DoCell, "_EC_2014.Rdata"))

    save(xn, mods, hab, ClimBCR,
        do_1spec1run, glm_skeleton, do_sample_0, do_sample, Lc_cut,
        file=paste0("c:/bam/Feb2014/data_", hab, "_", ClimBCR,
        "_", DoCell, "_EC_2014.Rdata"))

    write.csv(getFancyModsTab(mods), 
        file=paste0("c:/bam/Feb2014/ModelID_", hab, "_", ClimBCR, "_", DoCell, ".csv"))

}

## ---- analysis

library(mefa4)
hab_list <- c("HAB12","HAB22","HAB31","HAB41")
SPATIAL_list <- c(TRUE, FALSE)
n <- 10
B <- 2

hab <- "HAB31" # HAB12, HAB22, HAB31, HAB41
SPATIAL <- TRUE

load(paste0("c:/bam/Feb2014/data_", hab, "_", 
    ifelse(SPATIAL, "Clim", "BCR"), "_EC_2014.Rdata"))

res_cawa <- lapply(1:B, do_1spec1run, 
    i="CAWA", mods=mods, xn=xn, hab=hab, n=n, use_wt=TRUE, silent=TRUE)
save(res_cawa, file=paste0("res_CAWA_", hab, "_", 
    ifelse(SPATIAL, "Clim", "BCR"), "_EC_2014.Rdata"))
res_osfl <- lapply(1:B, do_1spec1run, 
    i="OSFL", mods=mods, xn=xn, hab=hab, n=n, use_wt=TRUE, silent=TRUE)
save(res_osfl, file=paste0("res_OSFL_", hab, "_", 
    ifelse(SPATIAL, "Clim", "BCR"), "_EC_2014.Rdata"))
res_coni <- lapply(1:B, do_1spec1run, 
    i="CONI", mods=mods, xn=xn, hab=hab, n=n, use_wt=TRUE, silent=TRUE)
save(res_coni, file=paste0("res_CONI_", hab, "_", 
    ifelse(SPATIAL, "Clim", "BCR"), "_EC_2014.Rdata"))












if (FALSE) {
## BCR interaction
x <- glm(CAWA ~ HAB31*xBCR+offset(QPAD_CAWA), DAT[1:70000,], family=poisson)
## becareful with BCR interaction: NAs are not to be treated as 0
## exclude BCRs up fron to avoid this:
with(DAT, table(xBCR, HAB31, useNA="a"))
apply(with(DAT, table(xBCR, HAB31)),1,min)
## HAB31: exclude Arctic and USeast
## HAB41: exclude 13, 
}
hab <- "HAB31" # HAB12, HAB22, HAB31, HAB41
library(mefa4)
load(paste0("c:/bam/Feb2014/data", hab, "_ec_2014.Rdata"))

system.time(res <- do_1spec1run(1, "CAWA", mods, xn=xn, hab=hab, n=10))
data.frame(id=structure(res$mid,names=names(mods)))


## -- prediction: becareful with BCR interaction: NAs are not to be treated as 0

source("c:/Dropbox/bam/EC_contract_2014/BAM-EC_02_analysis_functions.R")
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))

cf <- res$coef[[1]]
setdiff(names(cf), colnames(Xn))
Xn2 <- Xn[,names(cf)]
pr <- exp(drop(Xn2 %*% cf))


hab
dim(xn)
plot(DAT$POINT_X, DAT$POINT_Y, pch=19, cex=0.1)
points(xn$POINT_X, xn$POINT_Y, pch=19, cex=0.1, col=2)

#set.seed(1234)
#Bm <- pblapply(1:200, do_sample, g=xn$bootg, x=xn, n=5)



if (FALSE) {
iii <- sample(nrow(DAT), 5000)
dat <- DAT[iii,]
pdf("c:/bam/Feb2014/na-map.pdf", onefile=TRUE)
for (i in 1:ncol(DAT)) {
tmp <- is.na(dat[,i])
with(dat,plot(POINT_X,POINT_Y,pch=19,cex=0.2,col=1,
    main=colnames(DAT)[i]))
with(dat[tmp,],points(POINT_X,POINT_Y,pch=19,cex=0.3,col=2))
}
dev.off()
}


system.time(res <- do_1spec1run(13, "CAWA", mods[c(1:10)], xn=xn, hab=hab, n=5))
data.frame(id=structure(res$mid,names=names(mods)))

res_cawa <- do_1spec1run(1, "CAWA", mods[1:3], xn=xn, hab=hab, n=1)
res_osfl <- do_1spec1run(1, "OSFL", mods)
res_coni <- do_1spec1run(1, "CONI", mods)
res_cawa$coef[[10]]
res_osfl$coef[[10]]
res_coni$coef[[10]]
res_cawa$mid
res_osfl$mid
res_coni$mid
res_cawa$hi
res_osfl$hi
res_coni$hi

save(res_cawa, res_osfl, res_coni, file="c:/bam/Feb2014/prelim_ec_2014.Rdata")

m <- glm(CAWA ~ 1, data=DAT, family=poisson, offset=DAT$qpad_cawa)
m1 <- glm(CAWA ~ HAB, data=DAT, family=poisson, offset=DAT$qpad_cawa)
m2 <- glm(CAWA ~ HAB + HEIGHT, data=DAT, family=poisson, offset=DAT$qpad_cawa)
m3 <- glm(CAWA ~ HAB * HEIGHT, data=DAT, family=poisson, offset=DAT$qpad_cawa)
BIC(m, m1, m2, m3)

HAB/AGE
ROAD
CTI
DISTURB
CLIMATE


Terms <- getTerms(mods, "list")
if (hab=="HAB31") {
    DAT$HAB <- DAT[,hab]
    DAT$isDecid <- ifelse(DAT$HAB == "Decid", 1L, 0L)
    DAT$isMixed <- ifelse(DAT$HAB == "Mixed", 1L, 0L)
    DAT$isForest <- ifelse(DAT$HAB %in% c("Conif","Decid","Mixed"), 1L, 0L)
    DAT$isNonforest <- ifelse(!(DAT$HAB %in% c("Conif","Decid","Mixed")), 1L, 0L)
}
xn <- model.frame(getTerms(mods, "formula"), DAT[,Terms], 
    na.action=na.pass)
Xn <- model.matrix(getTerms(mods, "formula"), xn)
colnames(Xn) <- fixNames(colnames(Xn))

Xnj <- Xn[jj,]

with(DAT, plot(POINT_X, POINT_Y, pch=19, cex=0.1, col=DAT$RANGE_CAWA+1))


#rn <- sample(rownames(xn), nrow(xn))
rn <- sample(rownames(xn), 5000)
xnqs <- xn[rn,]
Xnqs <- Xn[rn,]





plotfun <- function(i) {
    x <- if (is.factor(xy[,i]))
        table(xy[,i]) else density(xy[,i], na.rm=TRUE)
    plot(x, main=colnames(xy)[i])
    invisible(NULL)
}

pdf("vars.pdf", onefile=TRUE)
for (i in 1:ncol(xy))
    plotfun(i)
dev.off()






## reclass LCC etc

## process CAS

## figure out ranges/distributions -- transform/rescale

## figure out meaning (metadata)


ff <- ~ (NORM_6190_CMDpm  +  NORM_6190_CMIJJApm +
 NORM_6190_CMIpm  +  NORM_6190_DD0    +  NORM_6190_DD5  +   
 NORM_6190_EMT    +  NORM_6190_FFP    +  NORM_6190_MAP  +   
 NORM_6190_MAT    +  NORM_6190_MCMT   +  NORM_6190_MSP  +   
 NORM_6190_MWMT   +  NORM_6190_NFFD   +  NORM_6190_PAS  +   
 NORM_6190_PET    +  NORM_6190_PPT_sm +  NORM_6190_PPT_wt + 
 NORM_6190_TD )

x <- DAT[Bm[[1]],]

m0 <- glm(CAWA ~ 1, data=x, offset=x$qpad_cawa, family=poisson)
m1 <- step(m0, scope=ff, direction="f")

m01 <- glm(OSFL ~ 1, data=x, offset=x$qpad_osfl, family=poisson)
m11 <- step(m01, scope=ff, direction="f")

m02 <- glm(CONI ~ 1, data=x, offset=x$qpad_coni, family=poisson)
m12 <- step(m02, scope=ff, direction="f")

anova(m1)
anova(m11)
anova(m12)


cm <- model.frame(ff,x)
colnames(cm) <- sub("NORM_6190_","",colnames(cm))
rownames(cc) <- NULL

pc <- princomp(cm)
cor(cbind(pc$scores[,1], cm))[,1]
