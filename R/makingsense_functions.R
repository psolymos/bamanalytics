getOK <- function(res) {
    sapply(res, class) != "try-error"
}

fix_underscore <- function(x) {
        x <- lapply(x, strsplit, "_")
        x <- sapply(x, function(z) paste(z[[1]], collapse="\\_", sep=""))
        x
}

getTerms <- function(mods, type=c("formula", "list"), intercept=TRUE) {
    type <- match.arg(type)
    x <- unlist(lapply(unlist(mods), function(z) as.character(z)[3]))
#    x <- unname(substr(x, 5, nchar(x)))
    x <- gsub(". + ", "", x, fixed=TRUE)
    x <- unlist(strsplit(x, "+", fixed=TRUE))
    x <- unlist(strsplit(x, "*", fixed=TRUE))
    if (type == "list")
        x <- unlist(strsplit(x, ":", fixed=TRUE))
    x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
    x <- unique(x)
    if (type == "formula") {
        x <- paste("~", paste(x, collapse=" + ", sep=""))
        if (!intercept)
            x <- paste(x, "- 1")
        x <- as.formula(x)
    }
    x
}

fixNames <- function(x, sep=":") {
    unlist(lapply(x, function(z) {
        paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
    }))
}

Grepl <- function(pattern, x) {
    rowSums(sapply(pattern, function(z) grepl(z, x))) > 0
}

getEst <- function(res, stage=NULL, na.out=TRUE, X) {
    if (!missing(X))
        Xn <- X
    OK <- !sapply(res, inherits, "try-error")
    if (any(!OK))
        warning(paste("try-error found:", sum(!OK)))
    ii <- sapply(res[OK], "[[", "iteration")
    est <- Xn[1:length(ii),,drop=FALSE]
    rownames(est) <- ii
    est[] <- 0
    if (is.null(stage))
        stage <- length(res[[ii[1]]]$coef)
    if (stage > 0) {
        for (i in 1:length(ii)) {
            tmp <- res[[ii[i]]]$coef[[stage]]
            names(tmp) <- fixNames(names(tmp))
            est[i,match(names(tmp), colnames(est))] <- tmp
        }
    } else {
        for (i in 1:length(ii)) {
            est[i,1] <- res[[ii[i]]]$null
        }
    }
    if (any(!OK) && na.out) {
        nas <- matrix(NA, sum(!OK), ncol(est))
        rownames(nas) <- which(!OK)
        est <- rbind(est, nas)
    }
    est
}

getCaic <- function(res, stage=NULL, na.out=TRUE) {
    OK <- !sapply(res, inherits, "try-error")
    if (is.null(stage))
        stage <- length(res[[which(OK)[1]]]$coef)
    caic <- numeric(length(OK))
    caic[!OK] <- NA
    for (run in which(OK)) {
        if (stage == 0) {
            cc <- attr(res[[run]]$caic[[1]], "StartCAIC")
        } else {
            cc <- res[[run]]$caic[[stage]]
            cc <- cc[which.min(cc)]
        }
        caic[run] <- cc
    }
    if (!na.out)
        caic <- caic[OK]
    caic
}

getSummary <- function(res, stage=NULL) {
    est <- getEst(res, stage=stage)
    est <- est[,colSums(abs(est), na.rm=TRUE) > 0]
    fr <- colMeans(abs(est) > 0, na.rm=TRUE)
    cf <- colMeans(est, na.rm=TRUE)
    se <- apply(est, 2, sd, na.rm=TRUE)
    z <- cf/se
    p <- 2 * pnorm(-abs(z))
    cmat <- cbind(cf, se, fr, z, p)
    colnames(cmat) <- c("Estimate", "Std. Error", "Freq.", "z value", "Pr(>|z|)")
    cmat
}
#printCoefmat(getSummary(res))

getVcov <- function(res) {
    est <- getEst(res)
    est <- est[,colSums(abs(est), na.rm=TRUE) > 0]
    cov(est)
}
getConfint <- function(res, level=0.95, type=c("tboot","quantile")) {
    type <- match.arg(type)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    s <- getSummary(res)
    parm <- rownames(s)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%", sep="")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    if (type == "tboot") {
        fac <- qnorm(a)
        ci[] <- s[,1] + s[,2] %o% fac
    } else {
        est <- getEst(res)
        est <- est[,colSums(abs(est)) > 0]
        cii <- t(apply(est, 2, quantile, probs=a))
        rownames(cii) <- parm
        ci[] <- cii[parm,,drop=FALSE]
    }
    return(ci)
}
#getConfint(res, type="tboot")

getMidPure <- function(res, mods) {
    OK <- !sapply(res, inherits, "try-error")
    mid <- data.frame(t(sapply(res[OK], function(z) unlist(z$mid))))
    rownames(mid) <- sapply(res[OK], "[[", "iteration")
    colnames(mid) <- names(mods)
    mid
}
getMid <- function(res, mods, use_rmax=FALSE) {
    OK <- !sapply(res, inherits, "try-error")
    mid <- data.frame(t(sapply(res[OK], function(z) unlist(z$mid))))
    rownames(mid) <- sapply(res[OK], "[[", "iteration")
    colnames(mid) <- names(mods)
    mid0 <- do.call(cbind,
        lapply(1:ncol(mid), function(i) 
            apply(mid[,1:i,drop=FALSE], 1, paste, collapse=" ")))
    if (use_rmax)
        require(ade4)
    Rao <- sapply(1:ncol(mid), function(i) {
        est <- suppressWarnings(getEst(res, stage=i, na.out=FALSE))
        est0 <- ifelse(est>0, 1, 0)
    #    d <- dist(est0, "manhattan") / ncol(est)
        d <- dist(est0)
        rmax <- if (use_rmax)
            suppressWarnings(divcmax(d))$value else 1
        PI <- data.matrix(rep(1/nrow(est), nrow(est)))
        drop(t(PI) %*% as.matrix(d) %*% PI) / rmax
    })
    attr(Rao, "use_rmax") <- use_rmax
    attr(mid, "Rao") <- Rao
    Gini <- 1 - apply(mid0, 2, function(z) sum((table(z)/length(z))^2))
    attr(mid, "Gini") <- Gini
    mid
}

getFancyMid <- function(res, mods) {
    mid <- getMidPure(res, mods)
    out <- lapply(mid, function(z) data.frame(table(Prop=z)))
    for (i in 1:length(out)) {
        tmp <- out[[i]]
        rownames(tmp) <- tmp$Prop
        Terms <- c("NULL", sapply(mods[[i]], function(z) as.character(z)[3]))
        names(Terms) <- 0:(length(Terms) - 1)
        tmp$Terms <- Terms[match(tmp$Prop, names(Terms))]
        tmp$Prop <- round(tmp$Freq/sum(tmp$Freq), 4)
        out[[i]] <- tmp
    }
    out
}
#getFancyMid(res, mods)
#getFancyMid(allres, mods, allspp=TRUE)
getFancyMidTab <- function(res, mods, ...) {
    x <- getFancyMid(res, mods, ...)
    for (i in 1:length(x)) {
        x[[i]] <- data.frame(Stage=rep(names(x)[i], nrow(x[[i]])),
            Model_ID=rownames(x[[i]]), 
            Full_ID=paste(i, rownames(x[[i]]), sep="."),
            x[[i]])
    }
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
}
getFancyModsTab <- function(mods, chmax=60) {
    out <- mods
    for (i in 1:length(out)) {
        Terms <- c("NULL", sapply(mods[[i]], function(z) as.character(z)[3]))
#        if (any(nchar(as.character(Terms)) > chmax))
#            Terms <- compress(Terms)
        tmp <- data.frame(Stage=rep(names(mods)[i], length(Terms)),
            Model_ID=0:(length(Terms) - 1), 
            Full_ID=paste(i, 0:(length(Terms) - 1), sep="."),
            Terms=Terms)
        out[[i]] <- tmp
    }
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out
}

## bootstrap based model support

midfig <- function(mid, m=apply(mid, 2, max), ...) {
    k <- ncol(mid)
    n <- nrow(mid)
    ylim <- c(0, k+1)
    xlim <- c(-(2+max(m)/2), (2+max(m)/2))
    pt <- lapply(1:k, function(i) {
        tt <- table(mid[,i])
        tt <- tt[match(0:m[i], names(tt))]
        tt[is.na(tt)] <- 0
        names(tt) <- 0:m[i]
        tt/n
    })
    yt <- lapply(1:k, function(i) rep(i, m[i]+1))
    xt <- lapply(1:k, function(i) 0:m[i] - m[i]/2)
    wl <- lapply(1:(k-1), function(i) {
        tt <- table(interaction(mid[,i], mid[,i+1]))
        tt[tt==0] <- NA
        tt/n
    })
    yl <- lapply(1:(k-1), function(i) {
        matrix(c(i, i+1), length(wl[[i]]), 2, byrow=TRUE)
    })
    xl <- lapply(1:(k-1), function(i) {
        z <- strsplit(names(wl[[i]]), "\\.")
        z <- matrix(as.integer(unlist(z)), length(z), 2, byrow=TRUE)
        z[,1] <- z[,1] - m[i]/2
        z[,2] <- z[,2] - m[i+1]/2
        z
    })
    wl <- do.call(c, wl)
    yl <- do.call(rbind, yl)
    xl <- do.call(rbind, xl)
    rr <- order(wl)
    wl <- wl[rr]
    yl <- yl[rr,]
    xl <- xl[rr,]
    CEX <- ifelse(unlist(pt)>0,1,0)*1+2.5*unlist(pt)/max(unlist(pt))
#    CEX <- 2+2.5*unlist(pt)/max(unlist(pt))
    LWD <- 3+5*wl/max(wl,na.rm=TRUE)
    COL <- grey(ifelse(is.na(wl), 0, wl/max(wl,na.rm=TRUE)))
    plot(xlim, ylim, type="n", xlab="", ylab="", axes=FALSE, ...)
    segments(-m/2, 1:k, m/2, 1:k)
    segments(xl[,1], yl[,1],xl[,2], yl[,2], 
        col=ifelse(is.na(wl), NA, 1), lwd=LWD)
    segments(xl[,1], yl[,1],xl[,2], yl[,2], 
        col=ifelse(is.na(wl), NA, COL), lwd=LWD-3)
    points(unlist(xt), unlist(yt), 
        cex=CEX,
        pch=19, col=grey(unlist(pt)/max(unlist(pt))))
    points(unlist(xt), unlist(yt), 
        cex=CEX,
        pch=21)
    text(unlist(xt), unlist(yt), 
        unlist(lapply(m, function(z) 0:z)))
    invisible(NULL)
}

getModelVariation <- function(res, mods, use_rmax=FALSE) {
    mid <- getMid(res, mods, use_rmax)
    out <- cbind(Gini=attr(mid, "Gini"),
        Rao=attr(mid, "Rao"))
    rownames(out) <- names(mods)
    out
}

plotMid <- function(res, mods, web=TRUE, ...) {
    if (is.character(res))
        res <- allres[[res]]
    mid <- getMidPure(res, mods)
    if (web) {
        opar <- par(mai=0.1*c(1,1,2,1))
        midfig(mid, sapply(mods, length), ...)
            #main=paste(taxa(mm)[res[[1]]$species,"CommonName"], " (", res[[1]]$species, ")", sep="")
            #main=as.character(taxa(mm)[res[[1]]$species,"English_Name"])
            #, ...)
        text(rep(-0.5-0.5*max(sapply(mods, length)), ncol(mid)), 
            seq_len(ncol(mid))+0.25, names(mods))
        par(opar)
    } else {
        rc <- c(2,2)
        if (ncol(mid) != 4)
            rc[1] <- ceiling(ncol(mid)/2)
        opar <- par(mfrow=rc, mai=0.5*c(1,5,1,1), las=1)
        for (i in 1:ncol(mid)) {
            tmp <- table(mid[,i])
            aa <- rep(0, length(mods[[i]]) + 1)
            names(aa) <- 0:length(mods[[i]])
            aa[names(tmp)] <- tmp
            names(aa) <- c(".", mods[[i]])
            col <- grey(aa/sum(aa))
            barplot(rev(aa), main=names(mods)[i], col=rev(col), horiz=TRUE)
        }
        par(opar)
    }
    invisible(NULL)
}

## --- prediction

## add in option for offsets ???
getDataPred <- function(res, ...) {
    est <- getEst(res, ...)
    tX <- t(Xn[,colnames(est)])
    mu <- apply(est, 1, function(z) crossprod(tX, z))
    mu
}
#mu <- getDataPred(res)
#summary(apply(exp(mu), 1, median, na.rm=TRUE))

## ---------------- predict: natural veg and age [AUPRF]

predStat <- function(X, est, level=0.95, n=0, ci=TRUE, raw=FALSE) {
    tX <- t(X)
    mu <- apply(est, 1, function(z) crossprod(tX, z))
    if (raw)
        return(mu)
    pr <- exp(mu)
    out <- matrix(NA, nrow(pr), ifelse(ci, 6, 2))
    rownames(out) <- rownames(X)
    out[,1] <- rowMeans(pr, na.rm=TRUE)
    out[,2] <- apply(pr, 1, median, na.rm=TRUE)
    if (!ci) {
        colnames(out) <- c("Mean", "Median")
    } else {
        colnames(out) <- c("Mean", "Median", "CL1n", "CL2n", "CL1q", "CL2q")
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        if (n > 0) {
            tmp <- exp(replicate(n, rnorm(nrow(mu), rowMeans(mu), apply(mu, 1, sd))))
            out[,3:4] <- t(apply(tmp, 1, quantile, probs=a, na.rm=TRUE))
        }
        out[,5:6] <- t(apply(pr, 1, quantile, probs=a, na.rm=TRUE))
    }
    out
}
