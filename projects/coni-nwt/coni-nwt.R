#' ---
#' title: "Developing ARU based QPAD offsets for CONI in AB+NWT"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "`r as.Date(Sys.time())`"
#' output: pdf_document
#' ---
#+ echo=FALSE
#knitr::opts_chunk$set(eval=FALSE)
ROOT <- "~/GoogleWork/collaborations/coni-nwt"
par(las=1)
#'
#' # Intro
#'
#' The goal is to develop ARU/recognizer based model to estimate
#' - CONI availability as a function of date and time,
#' - and CONI EDR based on sound pressure level and measured distances.
#'
#' # Preamble
library(mefa4)
library(maptools)
library(survival)
#'
#' # Recognizer data
#'
#' Read NWT data of distances
x_dis <- read.csv(file.path(ROOT, "data", "EDRCalculationData_65Score.csv"))
str(x_dis)
x_dis <- x_dis[rowSums(is.na(x_dis)) == 0,]
#' Read AB data of 1st detections
x_ab <- read.csv(file.path(ROOT, "data", "BorealDensity_Alberta_Clean.csv"))
str(x_ab)
#' Read NWT data of 1st detections
x_nwt <- read.csv(file.path(ROOT, "data", "BorealDensity_NWT_Clean.csv"))
str(x_nwt)
#' Compare fields nd adjust
setdiff(colnames(x_ab), colnames(x_nwt))
setdiff(colnames(x_nwt), colnames(x_ab)) # all 1s
x_ab$count <- 1 # add 1s to AB data
#'
#' # EDR estimation
#'
#' Binomial model: $(Y \mid p) \sim Bernoulli(p(d))$,
#' probability following half-normal distance curve $p(d) = exp(-d^2 / \tau^2)$.
#' We fit the loglinear model
#' $log(p(d)) = -d^2 / \tau^2 = x \frac{1}{\tau^2} = 0 + x \beta$
#' with the log (can be unstable) or cloglog (more robust) link
md1 <- glm(hit ~ I(-distance^2) - 1, x_dis, family=binomial("log"))
md2 <- glm(hit ~ I(-distance^2) - 1, x_dis, family=binomial("cloglog"))
#' Both were successful, let's see which fit better
AIC(md1, md2)
#' Here are EDRs from the 2 models
data.frame(EDR_in_m=c(log_link=unname(sqrt(1/coef(md1))),
    cloglog_link=unname(sqrt(1/coef(md2)))))
#' We pick `md1` over `md2` based on AIC
EDR <- unname(sqrt(1/coef(md1)))
#' Get affine transformation of the distance function
#' that fits the level-distance scatterplot best
f <- function(x) {
    p <- x[1] + x[2] * exp(-x_dis$distance^2/EDR^2)
    sum((x_dis$level - p)^2)
}
opt <- optim(c(100, 1), f, method="Nelder-Mead")
#' Plot level vs distance with the transformed distance function
with(x_dis, plot(distance, level,
    col=c("lightblue", "blue")[hit+1],
    ylim=c(0, max(x_dis$level))))
rug(x_dis$distance[x_dis$hit == 1], side=3, col="blue")
rug(x_dis$distance[x_dis$hit == 0], side=1, col="lightblue")
curve(opt$par[1] + opt$par[2] * exp(-x^2/EDR^2), add=TRUE)
abline(v=EDR, h=c(opt$par[1], sum(opt$par)), lty=2)
#'
#' # Availability estimation
#'
#' Join the 2 tables
stopifnot(all(colnames(x_nwt) == colnames(x_ab))) # consistency check once more
x <- rbind(x_ab, x_nwt)
with(x, plot(long, lat, pch="."))
#'
#' Determining date/time and related variables
x$datetime <- as.POSIXlt(
    with(x, paste0(year, "-", month, "-", day, " ", hr, ":", min, ":00")))
#' Ordinal day
x$JULIAN <- x$datetime$yday
x$JDAY <- x$JULIAN / 365
#' Time since local sunrise
x$start <- x$datetime$hour + x$datetime$min / 60
x$srise <- sunriset(
    as.matrix(x[,c("long", "lat")]),
    as.POSIXct(x$datetime, tz="America/Edmonton"),
    direction="sunrise", POSIXct.out=FALSE) * 24
x$TSSR <- (x$start - x$srise) / 24
#' Transformed variables
x$TOD <- x$start / 24
x$JDAY2 <- x$JDAY^2
x$sin <- sin(x$TOD * 2 * pi)
x$cos <- cos(x$TOD * 2 * pi)
x$sin2 <- sin(x$TOD * 2 * pi)^2
x$cos2 <- cos(x$TOD * 2 * pi)^2
#' check if all these triginometric functions make sense given
#' the distribution of the data
op <- par(mfrow=c(2,2))
plot(0:240/10, sin(pi*2*(0:240)/240), ylab="sin(TOD)", xlab="TOD", type="l")
points(x$TOD*24, x$sin)
plot(0:240/10, sin(pi*2*(0:240)/240)^2, ylab="sin2(TOD)", xlab="TOD", type="l")
points(x$TOD*24, x$sin2)
plot(0:240/10, cos(pi*2*(0:240)/240), ylab="cos(TOD)", xlab="TOD", type="l")
points(x$TOD*24, x$cos)
plot(0:240/10, cos(pi*2*(0:240)/240)^2, ylab="cos2(TOD)", xlab="TOD", type="l")
points(x$TOD*24, x$cos2)
par(op)
#' Well, not quite. Let's just use `sin` and `cos` as these seems to capture
#' the expected patterns. Also: we cannot really assess interactions
#' between `TOD` and `JDAY` because of the large gap in the daytime.
#'
#' Subset based on `station`:
#' we drop stations where the species was never observed.
#' If it was observed at least once, we treat 0s as false negatives
st1 <- unique(as.character(x$station)[x$detection > 0])
#' This is the proportion of stations with >0 detections
length(st1)/nlevels(x$station)
#' `x1` is the subset
x1 <- x[x$station %in% st1,]
#' Proportion of detection in the assumed present locations
mean(x1$detection)
#'
#' Making the data survival model compatible:
#' we can go over rows of `x1` because files are not duplicated
#' (i.e. no multiple detections from same file).
stopifnot(!any(duplicated(x1$file)))
#' Events are coded into a survival object:
#' we treat nondetections as censored events at 600 sec,
x1$time <- ifelse(is.na(x1$timeofdetection), 600, x1$timeofdetection)
#' time cannot be 0, so we use 1 sec instead
x1$time[x1$time == 0] <- 1
#' We give time in minutes, so we get rate as events/min.
x1$sv <- Surv(x1$time/60, x1$detection)
#' Fit a series of survival models
mods <- list(m0 = survreg(sv ~ 1, x1, dist="exponential"),
    m1 = survreg(sv ~ JDAY, x1, dist="exponential"),
    m2 = survreg(sv ~ JDAY + JDAY2, x1, dist="exponential"),
    m3 = survreg(sv ~ cos, x1, dist="exponential"),
    m4 = survreg(sv ~ cos + sin, x1, dist="exponential"),
    m5 = survreg(sv ~ JDAY + cos, x1, dist="exponential"),
    m6 = survreg(sv ~ JDAY + cos + sin, x1, dist="exponential"),
    m7 = survreg(sv ~ JDAY + JDAY2 + cos, x1, dist="exponential"),
    m8 = survreg(sv ~ JDAY + JDAY2 + cos + sin, x1, dist="exponential"))

aic <- data.frame(df=sapply(mods, function(z) length(coef(z))), AIC=sapply(mods, AIC))
aic$dAIC <- aic$AIC - min(aic$AIC)
aic
#' `survreg` fits accelerated failure models, not proportional
#' hazards models, so the coefficients are logarithms of ratios of
#' survival times, and a positive coefficient means longer survival.
mb <- mods[[which.min(aic$AIC)]]
summary(mb)
#' Survival times
summary(predict(mb))
#' Event rate per unit (1 min) time
summary(1/predict(mb))
#' Probability of at least 1 event per 10 time units (mins)
summary(1-exp(-(1/predict(mb))*10))
#'
#' Visualize the date-time relationships:
#' set up a date-time grid
vjd <- seq(min(x$JDAY), max(x$JDAY), len=51*10)
vtd <- seq(0, 23/24, len=24*10)
pr <- expand.grid(JDAY=vjd, TOD=vtd)
pr$JDAY2 <- pr$JDAY^2
pr$cos <- cos(pr$TOD * 2 * pi)
pr$sin <- sin(pr$TOD * 2 * pi)
#' make predictions for the grid
fit <- 1/predict(mb, newdata=pr)
summary(fit)
summary(1-exp(-fit*10))
#' we plot contours for P(detection in 10 min)
z <- matrix(1-exp(-fit*10), length(vjd), length(vtd))
plot(TOD*24 ~ jitter(JULIAN), x, pch=19, cex=0.6, col="#80808020",
    main="Probability of availability (10 min)",
    xlab="Ordinal day", ylab="Hour")
points(TOD*24 ~ JULIAN, x[x$detection > 0,], pch=19, cex=0.6,
    col="#FF0040")
#contour(vjd*365, vtd*24, z, add=TRUE, col=2)
l <- c(0.001, 0.01, 0.1, 0.2, 0.4)
for (i in seq_along(l)) {
    col <- colorRampPalette(c("blue", "red"))( length(l) )[i]
    contour(vjd*365, vtd*24, z, add=TRUE, col=col, levels=l[i],
        labcex=0.8, lwd=1)
}
#'
#' # Calculating offsets
#'
#' - A (effective area): EDR is given in 100 m units si that area refers to ha
#' - q = 1 because we used EDR for unlimited distances
#' - p is based on survival model predictions and 10 min duration
OUT <- data.frame(
    A = (EDR/100)^2*pi,
    q = 1,
    p = 1-exp(-(1/predict(mb, newdata = x))*10))
OUT$Offset <- rowSums(log(OUT))
summary(OUT)
if (FALSE) {
fn <-  "CONI-AB-NWT-withOffsets.csv"
od <- setwd(ROOT)
write.csv(data.frame(x, OUT), row.names=FALSE, file=fn)
zip(gsub(".csv", "", fn), fn)
unlink(fn)
setwd(od)
}
#'
#' The End / Fin
