% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdmModelling.R
\name{sdmModelling}
\alias{sdmModelling}
\title{A function to carry out the sdm using random forest.}
\usage{
sdmModelling(samples, envStack, modFormula, ntrees, plot = FALSE)
}
\arguments{
\item{samples}{data frame of the samples for model building (e.g. samples
taken from species generated through \code{\link{create_sp}}). The data
frame needs to have columns 'x' and 'y' for the coordinates, and 'presence'
with the presence (1) and absence (0) of the species}

\item{envStack}{raster stack of environmental variables used for the model
(e.g output from \code{\link{create_env_nsets}})}

\item{modFormula}{model formula to be applied}

\item{ntrees}{number of trees for generating the random forest}

\item{plot}{TRUE or FALSE whether to plot the predicted distribution. For
the plot the probability of occurrence map is converted to presence-absence
by thresholding at 0.5 which is a sensible guess for random forests (i.e.
50% of trees predicted presence and 50% absence). Sampled coordinates are
also plotted (absence = black, presences = red)}
}
\value{
A raster of predicted probability of occurrences
}
\description{
This function uses random forest, but you could replace it with
 your own if you want to use a different algorithm/model/ensemble method.
 Anything can be used as long as the output is a single raster object of the
 predicted probability of occurrences.
}
\examples{

### generate a set of environmental variables
set.seed(9999)
envSet <- create_env_nsets(cellDims = c(100, 100),
                           sets     = c(4, 4, 3, 1),
                           model    = "Sph",
                           psill    = 1.5,
                           dep1     = 1,
                           rangeFun = function() exp(runif(1, 1, 6)),
                           propSamp = 0.25)

### generate a virtual species from the variables
sp <- create_sp(envStack = envSet,
                spFun    = "x[1] * x[5] * x[9]",
                spModel  = "Sph",
                spPsill  = 1,
                spRange  = 500,
                propSamp = 0.5,
                prev     = 0.1)

### an initial 'sample' of the species (assuming perfect detection)
sampPts <- data.frame(sampleRandom(sp$presence, 50, na.rm = TRUE, xy = TRUE))

### a formula to fit to random forest (additive for all vars + quadratics)
(form <- paste0("presence ~ ", paste(names(envSet), collapse = " + "), " + I(",
                paste(names(envSet), collapse = " ^ 2) + I("), " ^ 2)"))

### run the model
spMod <- sdmModelling(samples = sampPts,
                      envStack = envSet,
                      modFormula = form,
                      ntrees = 500,
                      plot = TRUE)

### compare the model against the true species distribution
par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3, 3))
plot(sp$prob, axes = FALSE, main = "True prob. of occ.")
plot(sp$presence, axes = FALSE, main = "True distribution")
plot(spMod, axes = FALSE, main = "SDM prediction")

}
\author{
Charlie Marsh (charlie.marsh@mailbox.org) & Yoni Gavish
}
