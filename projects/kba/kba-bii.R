#' ---
#' title: BII calculation for KBA
#' output: pdf_document
#' ---
#'
#' # Objective
#'
#' We are using species specific predicted density values and footprint
#' maps to come up with intactness index.
#' In the absence of 'reference state', we are having trouble
#' defining the index for a single species.
#' (The multi-species version should ideally be an aggregation of
#' single species intactness indices/layers).
#' In this write-up, I am trying to conceptualize a way of defining
#' the index, also making the algorithm explicit.
#'
#' # Notation
#'
#' * Density: $D_{i}$ density of a species in the $i$th pixel/watershed (non-negative number)
#' * Footprint: $H_{i}$ human footprint value in pixel/watershed $i$ (0--1 value, binary or proportion)
#'
#' # Preamble
library(ggplot2)
library(mvtnorm)
set.seed(1)
Theme <- theme_minimal() + theme(
    panel.grid=element_blank(),
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank())
#' # Dummy data
#'
#' Let's define a grid
## function to make a plot
Extent <- 100
df <- expand.grid(x=1:Extent, y=1:Extent)
#' Reference density surface, i.e. $D$ given $H=0$
df$D_ref <- exp(sin(df$x/100) + cos(df$y/100) -
    sin(df$x/50)^2 + cos(df$y/50)^2 -
    cos(df$x*df$y/2000)) / 100
Sigma1 <- matrix(c(1, 0.5, 0.5, 1)*500, 2, 2)
Sigma2 <- matrix(c(1, -0.5, -0.5, 1)*250, 2, 2)
Sigma3 <- diag(1,2) * 100
df$D_ref <- 1000*(dmvnorm(df[,c("x", "y")], c(25, 50), Sigma1) +
    dmvnorm(df[,c("x", "y")], c(75, 25), Sigma2) +
    dmvnorm(df[,c("x", "y")], c(80, 60), Sigma1))
ggplot(df, aes(x, y, z=D_ref, fill = D_ref)) +
    geom_raster() +
    geom_contour(colour = "white") +
    Theme
#' Human footprint
df$H <- plogis(1 - 10 * df$y/Extent + 6*df$y*df$x/Extent^2)
ggplot(df, aes(x, y, z=H, fill = H)) +
    geom_raster() +
    geom_contour(colour = "white") +
    Theme
#' Now a relationship between $H$ and $D$ to create the
#' density that is $D'=D \mid H$, i.e. conditional on footprint.
#' This is the density we are modelling spatially,
#' and trying to infer about $D$ (which is $D$ given $H=0$).
#' `Effect` is the effect of footprint.
#' We add error to ut (contour shows the expected value, the
#' salt-and-pepper patterns is due to variation)
Effect <- -2
df$Expected <- df$D_ref * exp(Effect * df$H) # expected value
df$D <- df$Expected * exp(rnorm(Extent^2, 0, 0.25)) # error
ggplot(df, aes(x, y, z=Expected, fill = D)) +
    geom_raster() +
    geom_contour(colour = "white") +
    Theme
#' Let's see how the conditional density relates to footprint
ggplot(df, aes(H, D)) +
    geom_bin2d() +
    theme_minimal()
#' and how it relates to the reference density
ggplot(df, aes(D_ref, D)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    theme_minimal()
#' # Model
#'
#' The model really is:
#' $D = D_{ref} e^{\beta H}$, or on the log scale:
#' $log(D) = log(D_{ref}) + \beta H + error$.
#'
#' We build a model between the predicted density (in the simulation, this
#' is the conditional density `D`) and footprint (`H`)
m <- lm(log(D) ~ H, df)
b <- coef(m)
#' We expect the intercept to be close to the mean of the log reference
#' density (this is what we are trying to estimate after all)
b[1] # intercept
mean(log(df$D_ref)) # mean log reference density
#' The slope should be close to `Effect` ($\beta$) defined above
Effect
b[2]
#' # Intactness
#'
#' Now we want to approximate `D_ref` based on the coefficients, `D`, and `H`.
#' First, let's predict the expected value based on the log-linear model
df$Pred <- exp(b[1] + b[2] * df$H)
ggplot(df, aes(x, y, z=Pred, fill = Pred)) +
    geom_raster() +
    geom_contour(colour = "white") +
    Theme
#' this surface as we see is just a weighted version of `H` because
#' we ignored all other spatial variation.
#' To bring that in, we use the conditional density that we
#' call adjusted density ($D_{adj} = D / e^{\beta H}$)
df$Dadj <- df$D / exp(b[2] * df$H)
ggplot(df, aes(x, y, z=Dadj, fill = Dadj)) +
    geom_raster() +
    Theme
#' this looks very similar to $D_{ref}$
ggplot(df, aes(D_ref, Dadj)) +
    geom_bin2d() +
    geom_abline(slope=1, col=2) +
    theme_minimal()
#'
#' # Intactness
#'
#' There are 2 ways of expressing intactness:
#'
#' 1. spatial representation: the ratio of conditional vs. reference (adjusted density is our estimator for reference density)
#' 2. overall index: this takes the ratio of the sums of the 2 layers
#'
#' Here is the spatial version first: calculate min/max
#' (this forces the index to be < 1 even when $D > D_{ref}$, i.e.
#' the species responds positively to footprint)
df$BI <- pmin(df$D, df$Dadj) / pmax(df$D, df$Dadj)
ggplot(df, aes(x, y, z=BI, fill = BI)) +
    geom_raster() +
    geom_contour(colour = "white") +
    Theme
#' Hmm. This seems a bit curcular. Hope it will work better when $D$ is unknown...
#'
#' Anyways, the regional or overall index is calculated by taking the ratio
#' of the sums (or means):
Nc <- mean(df$D)
Nc
Nr <- mean(df$Dadj)
Nr
min(Nc, Nr) / max(Nc, Nr)
#' This can be calculated as a scaled footprint
mean(exp(b[2] * df$H))
#' A way around the fact that the spatial (pixel level) representation
#' will be a scaled version of the footprint could be to
#' calculate the scaled footprint , $mean(e^{\beta H})$, in a watershed,
#' and display that on a map. This should take into account of how
#' effectively footprint is influencing density, which is how
#' we define the opposite of intact.


