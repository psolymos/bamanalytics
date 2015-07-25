#### step 2 -------------------------------------

glm_skeleton <- function(object, ...) {
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
    out$caic <- 0.5 * (out$aic + out$bic)
    out
}

getTerms <- function(mods, type=c("formula", "list")) {
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
    if (type == "formula")
        x <- as.formula(paste("~", paste(x, collapse=" + ", sep="")))
    x
}

fixNames <- function(x, sep=":") {
    unlist(lapply(x, function(z) {
        paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
    }))
}

## x=DAT, n is max no of pts in grid
do_sample_0 <- function(x, n=10) {
    require(mefa4)
    ## random PKEY at each SS -- THIS IS THE KEY!!!
    xx <- nonDuplicated(x[sample.int(nrow(x)),,drop=FALSE], SS)
    xt <- Xtab(~SS + GRID_ID, xx)
    ## smaller no of SS in grid is kept
    xt1 <- xt[,colSums(xt) <= n]
    ## >n no of SS is shuffled and 1st n is kept
    xt2 <- xt[,colSums(xt) > n]
    xt2 <- xt2[sample.int(nrow(xt2)),]
    for (cc in seq_len(ncol(xt2))) {
        rr <- which(unname(xt2[,cc]) > 0)
        #xt[sample(rr, length(rr) - n),cc] <- 0
        xt2[rr[-seq_len(n)],cc] <- 0
    }
    m1 <- Melt(xt1)
    m2 <- Melt(xt2)
    m <- rbind(m1, m2)
    xxx <- xx[xx$SS %in% m$rows,]
    out <- x$PKEY %in% xxx$PKEY
    out
}
## j is placeholder for Boot ID, g is grouping vector
do_sample <- function(j, g, x, n=10) {
    g <- as.integer(droplevels(as.factor(g)))
    out <- logical(length(g))
    for (k in seq_len(max(g))) {
        kk <- g==k
        dat <- droplevels(x[kk,,drop=FALSE])
        out[kk] <- do_sample_0(dat, n=n)
    }
    which(out)
}
#system.time(jj <- do_sample(1, DAT$bootg, DAT, n=5))
#Bm <- list(jj)

Lc_cut <-
function (lam, transform=FALSE) 
{
    if (transform)
        lam <- 1-exp(-lam)
    o <- order(lam)
    x <- lam[o]
    p <- seq_len(length(x))/sum(length(x))
    L <- cumsum(x)/sum(x)
    p <- c(0, p)
    L <- c(0, L)
    J <- p - L

    G <- sum(x * 1:length(x))
    G <- 2 * G/(length(x) * sum(x))
    G <- G - 1 - (1/length(x))

    m1 <- which.max(J)
    list(lam=unname(ifelse(transform, -log(1-x[m1]), x[m1])), 
        L=unname(L[m1+1]), 
        p=unname(p[m1+1]), S=unname(L[m1+1]+p[m1+1]), 
        G=G, J=max(p - L))
}


do_1spec1run <- function(j, i, mods, xn, hab, n=10, use_wt=TRUE, silent=TRUE) {
    if (length(j)==1) {
        set.seed(1000+j)
        jj <- do_sample(j, g=xn$bootg, x=xn, n=n)
    } else {
        jj <- j
    }

    x <- xn[jj,]
    y <- as.integer(x[,i])
    off0 <- as.numeric(x[,paste0("QPAD_", i)])
    if (use_wt) {
        tmp <- Xtab(~ x$GRID_ID + rownames(x), drop.unused.levels=TRUE)
        w0 <- rowSums(tmp)[match(x$GRID_ID, rownames(tmp))]
        w0 <- 1/sqrt(w0)
        rm(tmp)
    } else {
        w0 <- rep(1L, length(y))
    }


    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
#    gofList <- vector("list", nmods)

    ## ++++++++++++++ Null

    null <- glm_skeleton(glm(y ~ 1, 
        x, 
        family=poisson(), 
        offset=off0, 
        weights=w0,
        x=FALSE, y=FALSE, model=FALSE))
    best <- null
#    gof0 <- gof_fun(j, i, object=null, x0=x0, off00=off00, ss_ext=ss_ext)
#    gofbest <- gof0

    select_ip <- "IP" %in% names(mods)
    if (select_ip) {
        habmod <- glm_skeleton(try(glm(y ~ HAB,
            x,
            family=poisson(), 
            offset=off0, 
            weights=w0,
            x=FALSE, y=FALSE, model=FALSE), silent=silent))
        lam <- exp(drop(model.matrix(~ HAB, x) %*% habmod$coef))
        cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
        tb <- ifelse(table(hab=x$HAB[!is.na(x$HAB)], 
            lc=ifelse(lam >= cv$lam, 1, 0))>0, 1, 0)
        Hi <- rownames(tb)[tb[,"1"] > 0]
#        Lo <- rownames(tb)[tb[,"1"] == 0]
        x$IP <- rowSums(x[,paste0("Cell_", Hi),drop=FALSE])
    } else {
        Hi <- NULL
        cv <- NULL
    }

    for (l1 in 1:nmods) {
#        gc()
        mlist <- vector("list", nnmods[l1])
        glist <- vector("list", nnmods[l1])
        for (l2 in 1:nnmods[l1]) {
            mlist[[l2]] <- glm_skeleton(try(update(object=best, 
                formula=mods[[l1]][[l2]]), silent=silent))
#            glist[[l2]] <- gof_fun(j, i, object=mlist[[l2]], x0=x0, off00=off00, ss_ext=ss_ext)
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
#            gofbest <- glist[[mmid]]
        }
        bestList[[l1]] <- best
        caicList[[l1]] <- mcaic
#        gofList[[l1]] <- gofbest
    }
    tmp <- list(species=i, 
        iteration=if (length(j)==1) j else NULL,
        #bootid=jj,
        null=null$coef,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
#        gof0=gof0,
#        gof=gofList,
#        n_ext=length(ss_ext),
        hi=Hi,
        lc=cv,
        hab=hab,
        habmod=if (select_ip) habmod$coef else NULL)
    tmp
}
