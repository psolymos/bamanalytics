do_1spec1run_noW <- function(j, i, mods, 
silent=FALSE, hsh_name=NA, CAICalpha=1) 
{
    select_hsh <- !is.na(hsh_name)
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    if (select_hsh)
        hsh <- HSH[BB[,j],]
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
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE), CAICalpha=CAICalpha)
    best <- null
    ## Lorenz-tangent approach for core habitat delineation
    if (select_hsh) {
        HABV <- x[,hsh_name]
        habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD,
            x,
            family=poisson(), 
            offset=off, 
            #weights=w,
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
        #nmax=nmax,
        #w_id=w_id,
        habmod=habmod$coef,
        hsh_name=hsh_name)
    out
}
