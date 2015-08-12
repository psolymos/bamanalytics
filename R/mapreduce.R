library(mefa4)
set.seed(1234)
y <- Matrix(rpois(50, 0.5), 20, 10)
dimnames(y) <- list(letters[1:20], LETTERS[1:10])
x <- Melt(y)
x <- x[sample.int(nrow(x)),]
x <- data.frame(id=1:nrow(x), x)

file <- "trydata.csv"
write.csv(x, file, row.names=FALSE)
FUN <- function(x) return(x)
REDUCE <- rbind
nrows <- 20

nlines <- function(file) {
    ## http://r.789695.n4.nabble.com/Fast-way-to-determine-number-of-lines-in-a-file-td1472962.html
    ## needs Rtools on Windows
    if (.Platform$OS.type == "windows") { 
        nr <- as.integer(strsplit(system(paste("/RTools/bin/wc -l", 
            file), intern=TRUE), " ")[[1]][1])
    } else {
        nr <- as.integer(strsplit(system(paste("wc -l", 
            file), intern=TRUE), " ")[[1]][1])
    }
    nr
}

MapReduce_function <- function(file, nrows, FUN, REDUCE, ...) {
    ## Map
    nr <- nlines(file)
    m <- floor((nr-1) / nrows)
    mm <- (nr-1) %% nrows
    if (mm > 0)
        m <- m+1
    ## Reduce
    tmp0 <- read.csv(file, nrows=2, skip=0, header=TRUE, ...)
    cn <- colnames(tmp0)
    res <- list()
    for (i in 1:m) {
        tmp <- read.csv(file, nrows=nrows, skip=(i-1)*nrows+1, 
            header=FALSE, ...)
        colnames(tmp) <- cn
        res[[i]] <- FUN(tmp)
    }
    out <- do.call(REDUCE, res)
}

out <- MapReduce_function(file, nrows, FUN, REDUCE)
fff <- Xtab(value ~ rows + cols, out)
fff
y[rownames(fff),]

