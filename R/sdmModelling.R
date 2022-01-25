####################################################################################################
### sdmModelling.R
### 01/2022
### Charlie Marsh (charlie.marsh@mailbox.org)
###
### A function to carry out the sdm. This function uses random forest, but you could replace it with
### your own if you want to use a different algorithm/model/ensemble method. Anything can be used as
### long as the output is a single raster object of the predicted probability of occurrences.
###
### Args:
###   coords      = coordinates of samples,
###   envStack    = stack of environmental variables used for the model,
###   mod_formula = regform,
###   ntrees      = number of trees for randome forest
###
### Output:
###   pred  = a raster of predicted probability of occurrences
###
####################################################################################################

sdmModelling <- function(coords,
                         envStack,
                         mod_formula,
                         plot = FALSE,
                         ntrees) {
  require(raster)
  require(randomForest)

  ### prepare sampled coordinates
  coordsSample <- coords
  presence <- coordsSample[, 1]
  presence[presence > 1] <- 1
  presence <- as.factor(presence)

  ### extract environmental values for points
  coords <- SpatialPoints(coords.sample[, c("x", "y")])
  envValues <- as.data.frame(extract(envStack, coords))
  envLayers <- as.data.frame(envStack)

  ### modelling and predict
  mod <- randomForest(as.formula(mod_formula), data = envValues, ntree = ntrees)

  ### generate prediction raster
  pred <- predict(mod, newdata = envLayers, type = "prob")[, "1"]
  pred <- raster(ext  = extent(envStack[[1]]),
                 res  = res(envStack[[1]]),
                 vals = as.vector(pred))
  pred[cellFromXY(pred, coords.sample[, c("x", "y")])] <-  predict(mod, type = "prob")[, "1"]

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
