offset_fun <- function(j, i, x) {
    if (i %in% BAMspp) { ## spp with offsets (>=75 obs)
        type <- "full"
    } else {
        if (i %in% BAMCOEFS25$spp) { ## SPP based on >=25 obs
            type <- "const"
        } else {
            stop("no offsets available")
        }
    }
    if (type == "full") {
        xna <- as.data.frame(lapply(x[,c("LCC_combo","TREE","JDAY","TSSR")], 
            function(z) as.integer(is.na(z))))
        xna <- with(xna, interaction(LCC_combo,TREE,JDAY,TSSR, drop=TRUE, sep=""))
        xx <- lapply(levels(xna), function(z) x[xna==z,])
        id <- lapply(levels(xna), function(z) which(xna==z))
        names(xx) <- levels(xna)
        off <- lapply(xx, function(z) offset_fun0(j, i, x=z, type=type))
        OFF <- numeric(nrow(x))
        for (k in seq_len(length(off)))
            OFF[id[[k]]] <- off[[k]]
        gc <- globalBAMcorrections(i, 
            r=x$MAXDIS[!is.finite(OFF)],
            t=x$MAXDUR[!is.finite(OFF)])
        gc <- rowSums(log(as.matrix(gc)))
        OFF[!is.finite(OFF)] <- gc
    } else {
        OFF <- offset_fun0(j, i, x=x, type=type)
    }
    OFF
}

offset_fun0 <- function(j, i, x, type=c("full","const")) {
    BOOT <- j != 1
    if (type == "full") { ## spp with offsets (>=75 obs)
        ml1 <- getBAMmodellist()$sra
        ml2 <- getBAMmodellist()$edr
        if (any(is.na(x$JDAY)))
            ml1 <- ml1[!grepl("JDAY", ml1)]
        if (any(is.na(x$TSSR)))
            ml1 <- ml1[!grepl("TSSR", ml1)]
        if (any(is.na(x$TREE)))
            ml2 <- ml2[!grepl("TREE", ml2)]
        if (any(is.na(x$LCC_combo)))
            ml2 <- ml2[!grepl("LCC", ml2)]
        best <- bestmodelBAMspecies(i, 
            model.sra=names(ml1), model.edr=names(ml2),
            type=ifelse(BOOT, "multi", "BIC"))
        out <- with(x, localBAMcorrections(i,
            r=MAXDIS,
            t=MAXDUR,
            jday=JDAY, 
            tssr=TSSR, 
            tree=TREE, 
            lcc=LCC_combo,
            model.sra=best$sra, 
            model.edr=best$edr,
            boot=BOOT, ## MVN approach of j != 1
            ver=1))
    }
    if (type == "const") { ## SPP based on >=25 obs
        if (BOOT) {
            PHI <- exp(rnorm(1, 
                mean = BAMCOEFS25$sra_estimates[[i]][["0"]]$coef[1],
                sd = BAMCOEFS25$sra_estimates[[i]][["0"]]$vcov[1,1]))
            TAU <- exp(rnorm(1, 
                mean = BAMCOEFS25$edr_estimates[[i]][["0"]]$coef[1],
                sd = BAMCOEFS25$edr_estimates[[i]][["0"]]$vcov[1,1]))
        } else {
            PHI <- exp(BAMCOEFS25$sra_estimates[[i]][["0"]]$coef[1])
            TAU <- exp(BAMCOEFS25$edr_estimates[[i]][["0"]]$coef[1])
        }
        out <- customBAMcorrections(r=x$MAXDIS,
            t=x$MAXDUR,
            phi=PHI,
            tau=TAU)
    }
    #out
    rowSums(log(as.matrix(out)))
}

## Function to compare sets (factors are left untouched)
compare.sets <- function(x, y) {
    x <- as.factor(x)
    y <- as.factor(y)
    xl <- levels(x)
    yl <- levels(y)
    xa <- levels(droplevels(x))
    ya <- levels(droplevels(y))
    lab <- c(xlength=length(xl), ylength=length(yl),
        intersect=length(intersect(xl, yl)),
        union=length(union(xl, yl)),
        xbutnoty=length(setdiff(xl, yl)),
        ybutnotx=length(setdiff(yl, xl)))
    act <- c(xlength=length(xa), ylength=length(ya),
        intersect=length(intersect(xa, ya)),
        union=length(union(xa, ya)),
        xbutnoty=length(setdiff(xa, ya)),
        ybutnotx=length(setdiff(ya, xa)))
    rbind(labels=lab, unique=act)
}

arrange.intervals <-
function(x, sep="-")
{
    x <- ifelse(x > 0, 1L, 0L)
    ## start/stop
    ss <- strsplit(colnames(x), sep)
    names(ss) <- colnames(x)
    ss <- lapply(ss, as.numeric)
    ss <- do.call(rbind, ss)
    nr <- nrow(x)
    nc <- ncol(x)
    End <- Id <- array(NA, dim(x))
    rownames(End) <- rownames(Id) <- rownames(x)
    for (i in seq_len(nr)) {
        id <- which(x[i,] > 0)
        #cn <- colnames(x)[id]
        endv <- ss[id,2]
        id <- id[order(endv)]
        endv <- endv[order(endv)]
        End[i, seq_len(length(endv))] <- endv
        Id[i, seq_len(length(endv))] <- id
    }
    lc <- which(rev(cumsum(rev(nr - colSums(is.na(Id))))) < 1)[1L] - 1L
    list(x=x, end=End[,seq_len(lc)], id=Id[,seq_len(lc)])
}
