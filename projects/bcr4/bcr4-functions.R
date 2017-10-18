'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1,
        nobs = 1, class = "logLik")
}

glm_skeleton <- function(object, ..., CAICalpha=0.5) {
    if (inherits(object, "try-error"))
        return(structure(as.character(object), class="try-error"))
    out <- structure(list(
        call=object$call,
        formula=formula(object),
        coef=coef(object),
        converge=object$converge,
        logLik=as.numeric(logLik(object)),
        df=attr(logLik(object), "df"),
        nobs=nobs(object)), class="glm_skeleton")
    out$class0 <- class(object)[1L]
    out$aic <- -2*out$logLik + 2*out$df
    out$bic <- -2*out$logLik + log(out$nobs)*out$df
    out$caic <- CAICalpha * out$aic + (1-CAICalpha) * out$bic
    out
}

## j is the bootstrap run (survey row IDs taken from jth column of BB matrix)
## i is the species (ith column in YY matrix, and offset as well)
## mods is the list of model formulas to be used in updates
## CAICalpha=1 means AIC is used, 0 means BIC is used
do_1spec1run <- function(j, i, mods, silent=FALSE, CAICalpha=1, return_best=FALSE)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- OFF[BB[,j], i]
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
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    ## looping through models list
    for (l1 in 1:nmods) {
        if (nnmods[l1] > 0) {
            mlist <- vector("list", nnmods[l1])
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
        alpha=CAICalpha)
    if (return_best)
        out$best <- eval(best$call)
    out
}

visualize_road <- function(use_climate=TRUE) {
    k <- nlevels(DAT$hab)
    Xn0 <- Xn[1:k,]
    rownames(Xn0) <- levels(DAT$hab)
    Xn0[] <- 0
    diag(Xn0) <- 1
    Xn0[,1] <- 1
    Xn1 <- Xn0
    Xn1[,"ROAD"] <- 1
    tmp <- table(DAT$hab, DAT$isForest)
    Xn1[tmp[,"1"] > 0,"isForest:ROAD"] <- 1
    xx <- rbind(Xn0, Xn1)
    rownames(xx) <- paste(rep(c("OffR", "OnR"), each=k), rownames(Xn0))
    est <- if(use_climate) est2 else est1
    pr <- apply(est, 1, function(z) xx %*% z)
    rownames(pr) <- rownames(xx)
    q <- t(apply(exp(pr), 1, quantile, c(0.5, 0.05, 0.95)))
    mat <- t(matrix(q[,1], k, 2))
    pos <- barplot(mat, beside=TRUE, names.arg=rownames(Xn0),
        legend.text=c("Off road", "On road"), ylim=c(0, max(q)*1.2),
        ylab="Density in undisturbed habitat (ysd=100)", main=spp)
    xpos <- c(pos[1,], pos[2,])
    for (i in 1:length(xpos))
        segments(x0=xpos[i], y0=q[i,2], y1=q[i,3], lwd=3, col=1)
    invisible(q)
}

visualize_ysd <- function(use_climate=TRUE) {
    k <- nlevels(DAT$hab)
    ysd <- 0:120
    Xn0 <- Xn[1:length(ysd),]
    rownames(Xn0) <- ysd
    Xn0[] <- 0
    Xn0[,1] <- 1
    Xn0[,"ysd"] <- pmin(100, ysd)
    Xn0[,"ysd2"] <- Xn0[,"ysd"]^2
    Xn0[,"ysd05"] <- sqrt(Xn0[,"ysd"])
    Xn0[ysd <= 10,"ysd10"] <- 1
    Xn0[,"ysd60"] <- pmax(0, 1 - ysd / 60)
    Xn0[,"ysd100"] <- pmax(0, 1 - ysd / 100)
    Xn0[ysd > 10 & ysd < 100,"ysdmid"] <- 1
    qq <- list()
    tmp <- table(DAT$hab, DAT$isForest)
    est <- if(use_climate) est2 else est1
    for (i in 1:k) {
        Xn1 <- Xn0
        Xn1[,i] <- 1
        if (tmp[i,"1"] > 0) {
            Xn1[,"isForest:ysd"] <- Xn1[,"ysd"]
            Xn1[,"isForest:ysd2"] <- Xn1[,"ysd2"]
            Xn1[,"isForest:ysd05"] <- Xn1[,"ysd05"]
            Xn1[,"isForest:ysd10"] <- Xn1[,"ysd10"]
            Xn1[,"isForest:ysd60"] <- Xn1[,"ysd60"]
            Xn1[,"isForest:ysd100"] <- Xn1[,"ysd100"]
            Xn1[,"isForest:ysdmid"] <- Xn1[,"ysdmid"]
        }
        if (FALSE) {
        if (levels(DAT$hab)[i] == "Conif") {
            Xn1[,"isConif:ysd"] <- Xn1[,"ysd"]
            Xn1[,"isConif:ysd2"] <- Xn1[,"ysd2"]
            Xn1[,"isConif:ysd05"] <- Xn1[,"ysd05"]
            Xn1[,"isConif:ysd10"] <- Xn1[,"ysd10"]
            Xn1[,"isConif:ysd60"] <- Xn1[,"ysd60"]
            Xn1[,"isConif:ysd100"] <- Xn1[,"ysd100"]
            Xn1[,"isConif:ysdmid"] <- Xn1[,"ysdmid"]
        }
        if (levels(DAT$hab)[i] == "Decid") {
            Xn1[,"isDecid:ysd"] <- Xn1[,"ysd"]
            Xn1[,"isDecid:ysd2"] <- Xn1[,"ysd2"]
            Xn1[,"isDecid:ysd05"] <- Xn1[,"ysd05"]
            Xn1[,"isDecid:ysd10"] <- Xn1[,"ysd10"]
            Xn1[,"isDecid:ysd60"] <- Xn1[,"ysd60"]
            Xn1[,"isDecid:ysd100"] <- Xn1[,"ysd100"]
            Xn1[,"isDecid:ysdmid"] <- Xn1[,"ysdmid"]
        }
        }
        pr <- apply(est, 1, function(z) Xn1 %*% z)
        qq[[i]] <- t(apply(exp(pr), 1, quantile, c(0.5, 0.05, 0.95)))
    }
    names(qq) <- levels(DAT$hab)

    mat <- sapply(qq, function(z) z[,1])
    matplot(mat, lty=1, type="l", lwd=2, ylim=c(0, 1.2*max(mat)),
        ylab="Relative abundance", xlab="Years since last disturbance", main=spp)
    legend("topright", lty=1, col=1:k, lwd=2, legend=names(qq))

    invisible(qq)
}
