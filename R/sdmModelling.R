##############################################################################
### sdmModelling.R
### 01/2022

#' A function to carry out the sdm using random forest.
#'
#' @description This function uses random forest, but you could replace it with
#'  your own if you want to use a different algorithm/model/ensemble method.
#'  Anything can be used as long as the output is a single raster object of the
#'  predicted probability of occurrences.
#'
#' @author Charlie Marsh (charlie.marsh@@mailbox.org) & Yoni Gavish
#'
#' @param samples data frame of the samples for model building (e.g. samples
#'  taken from species generated through \code{\link{create_sp}}). The data
#'  frame needs to have columns 'x' and 'y' for the coordinates, and 'presence'
#'  with the presence (1) and absence (0) of the species
#' @param envStack raster stack of environmental variables used for the model
#'  (e.g output from \code{\link{create_env_nsets}})
#' @param modFormula model formula to be applied
#' @param ntrees number of trees for generating the random forest
#' @param plot TRUE or FALSE whether to plot the predicted distribution. For
#'  the plot the probability of occurrence map is converted to presence-absence
#'  by thresholding at 0.5 which is a sensible guess for random forests (i.e.
#'  50% of trees predicted presence and 50% absence). Sampled coordinates are
#'  also plotted (absence = black, presences = red)
#'
#' @return A raster of predicted probability of occurrences
#'
#' @examples
#'
#' ### generate a set of environmental variables
#' set.seed(9999)
#' envSet <- create_env_nsets(cellDims = c(100, 100),
#'                            sets     = c(4, 4, 3, 1),
#'                            model    = "Sph",
#'                            psill    = 1.5,
#'                            dep1     = 1,
#'                            rangeFun = function() exp(runif(1, 1, 6)),
#'                            propSamp = 0.25)
#'
#' ### generate a virtual species from the variables
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = "x[1] * x[5] * x[9]",
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 500,
#'                 propSamp = 0.5,
#'                 prev     = 0.1)
#'
#' ### an initial 'sample' of the species (assuming perfect detection)
#' sampPts <- data.frame(sampleRandom(sp$presence, 50, na.rm = TRUE, xy = TRUE))
#'
#' ### a formula to fit to random forest (additive for all vars + quadratics)
#' (form <- paste0("presence ~ ", paste(names(envSet), collapse = " + "), " + I(",
#'                 paste(names(envSet), collapse = " ^ 2) + I("), " ^ 2)"))
#'
#' ### run the model
#' spMod <- sdmModelling(samples = sampPts,
#'                       envStack = envSet,
#'                       modFormula = form,
#'                       ntrees = 500,
#'                       plot = TRUE)
#'
#' ### compare the model against the true species distribution
#' par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3, 3))
#' plot(sp$prob, axes = FALSE, main = "True prob. of occ.")
#' plot(sp$presence, axes = FALSE, main = "True distribution")
#' plot(spMod, axes = FALSE, main = "SDM prediction")
#'
################################################################################

#' @export
#' @importFrom raster extract raster res extent cellFromXY plot as.data.frame
#' @importFrom randomForest randomForest
#' @importFrom sp SpatialPoints
#' @importFrom graphics points
#' @importFrom stats as.formula
sdmModelling <- function(samples,
                         envStack,
                         modFormula,
                         ntrees,
                         plot = FALSE) {
  requireNamespace("randomForest")

  ### prepare sampled coordinates
  coordsSample <- samples
  presence <- coordsSample$presence
  presence[presence > 1] <- 1
  presence <- as.factor(presence)

  ### extract environmental values for points
  coords <- SpatialPoints(coordsSample[, c("x", "y")])
  envValues <- as.data.frame(extract(envStack, coords))
  envLayers <- as.data.frame(envStack)

  ### modelling and predict
  mod <- randomForest(as.formula(modFormula), data = envValues, ntree = ntrees)

  ### generate prediction raster
  pred <- predict(mod, newdata = envLayers, type = "prob")[, "1"]
  pred <- raster(ext  = extent(envStack[[1]]),
                 res  = res(envStack[[1]]),
                 vals = as.vector(pred))
  cellsSample <- cellFromXY(pred, coordsSample[, c("x", "y")])
  pred[cellsSample] <- predict(mod, type = "prob")[, "1"]

  ### optional plotting function
  if(plot == TRUE){
    predMap <- pred
    predMap[predMap  < 0.5] <- 0
    predMap[predMap >= 0.5] <- 1
    plot(predMap,
         colNA = "dark grey",
         main = paste("Samples =", length(presence)))
    points(coords[presence == 0], pch = 1, cex=1)
    points(coords[presence == 1], pch = 16, col = "red", cex=1)
  }

  ### return prediction raster
  return(pred)
}
