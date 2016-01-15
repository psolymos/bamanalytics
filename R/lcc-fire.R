library(mefa4)
load("c:/bam/May2015/out/data/pack_LCC05-fire-Canada_toNicole_2016-01-14.Rdata")

## do something with LCC05 here (treat it as factor, reclass, etc)
## fire year and size is also kept FYI
## bootg is the blocking factor for bootstrap (close to BCR/PROV)

DAT$LCC <- as.factor(DAT$LCC05)
table(DAT$LCC)
DAT <- droplevels(DAT)

## one more thing to note: be careful when subsetting!
## for DAT, YY, OFF it works fine (TAX is just for species names etc)
## you might have to redo BB after subsetting the rest

## this ensures that errors will not break the process
## failed model gets ridiculously high AIC
'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1, 
        nobs = 1, class = "logLik")
}

do_1spec1run_noW_simplified <- function(j, i, ...) 
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- OFF[BB[,j], i]
    ## Null model
    null <- glm(y ~ 1, 
        data=x, 
        family=poisson(), 
        offset=off, 
        x=FALSE, y=FALSE, model=FALSE, ...)
    best <- null
    ## Alternative models
    mlist <- list()
    cat("\n", i, j, "model: 0");flush.console()
    mlist[[1]] <- glm(y ~ 1, 
        data=x, 
        family=poisson(), 
        offset=off, 
        x=FALSE, y=FALSE, model=FALSE)
    gc()
    cat(" 1");flush.console()
    mlist[[2]] <- try(update(object=mlist[[1]], 
        formula=.~.+LCC))
    gc()
    cat(" 2");flush.console()
    mlist[[3]] <- try(update(object=mlist[[1]], 
        formula=.~.+bootg))
    gc()
    cat(" 3");flush.console()
    mlist[[4]] <- try(update(object=mlist[[1]], 
        formula=.~.+LCC+bootg))
    gc()
    cat(" 4");flush.console()
    mlist[[5]] <- try(update(object=mlist[[1]], 
        formula=.~.+LCC*bootg))
    gc()
    aic <- sapply(mlist, AIC)
    mlist[[which.min(aic)]]
}

## I suggest run j=1 for multiple species 1st
## so that you get an idea of how thing will look eventually
## and go for the bootstrap if you have to
m <- do_1spec1run_noW_simplified(j=1, i="OVEN")

## make a data frame to supply as newdata for prediction
ndf <- DAT
## do prediction by hand (predict method expects an offset as part of newdata)
f <- formula(m)
f[[2]] <- NULL
X <- model.matrix(f, ndf)
pr <- exp(drop(X[,names(coef(m))] %*% coef(m)))


