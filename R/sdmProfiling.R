################################################################################
# sdmProfiling
# 01/2022
#
# Args:
#   unsampled_coords  = data frame of unsampled cells with column names "x", "y" and "ID"
#   sampled_coords    = data frame of sampled cells with column names "presence", "x" and "y"
#   original_SDM      = raster file of SDM probability of occurrences using sampled cells
#   envStack          = stack of raster files for environmental variables for building SDM
#   regform           = model formula if using built-in random forest model
#   ntrees            = number of trees for generating forest if using built-in random forest model
#   parallel          = TRUE or FALSE to whether to run profiling on multiple cores.
#   mc.cores          = if parallel = TRUE then the number of cores to share jobs to
#   ...               = further arguments to be passed to the sdmModelling function
#
# Output: list of two matrices, diff0 and diff1 where each row is an unsampled
#   cell, the first three columns are "x", "y" and "ID", the rest of the columns
#   are the predicted probabilities for all non-NA cells in the SDM raster when
#   that unsampled cell was made an absence (diff0) and a presence (diff1)
#
################################################################################

sdmProfiling <- function(unsampled_coords,
                         sampled_coords,
                         original_SDM,
                         envStack,
                         regform = NULL,
                         ntrees = 500,
                         parallel = FALSE,
                         mc.cores = mc.cores,
                         ...) {
  ### read in script to carry out SDMs using random forest
  source("sdmModelling.R")

  ### tidy up names
  names(sampled_coords)   <- tolower(names(sampled_coords))
  names(unsampled_coords) <- tolower(names(unsampled_coords))

  ## data storage
  diff0 <- rep(NA, nrow(unsampled_coords))
  diff1 <- rep(NA, nrow(unsampled_coords))

  if(parallel == FALSE) {
    pb <- txtProgressBar(min = 0, max = nrow(unsampled_coords), style = 3)

    for(cell in 1:nrow(unsampled_coords)) {
      ## change point to absence and run new SDM
      coords0 <- rbind(sampled_coords,
                       data.frame(presence = 0,
                                  x        = unsampled_coords[cell, "x"],
                                  y        = unsampled_coords[cell, "y"],
                                  id       = unsampled_coords[cell, "id"]))

      pred0 <- sdmModelling(coords      = coords0,
                            envStack   = envStack,
                            mod_formula = regform,
                            ntrees      = ntrees)

      ## change point to presence and run new SDM
      coords1 <- rbind(sampled_coords,
                       data.frame(presence = 1,
                                  x        = unsampled_coords[cell, "x"],
                                  y        = unsampled_coords[cell, "y"],
                                  id       = unsampled_coords[cell, "id"]))

      pred1 <- sdmModelling(coords      = coords1,
                            envStack   = envStack,
                            mod_formula = regform,
                            ntrees      = ntrees)

      dev0 <- abs(pred0@data@values - original_SDM@data@values)
      dev1 <- abs(pred1@data@values - original_SDM@data@values)
      dev0 <- dev0[!is.na(dev0)]
      dev1 <- dev1[!is.na(dev1)]
      diff0[cell] <- sum(as.vector(dev0), na.rm = TRUE)
      diff1[cell] <- sum(as.vector(dev1), na.rm = TRUE)

      setTxtProgressBar(pb, cell)
    }
  }

  if(parallel == TRUE) {
    source("parallel_profiling.R")
    require(parallel)
    require(parallelsugar)

    all.cells <- 1:length(unsampled_coords[, 1])
    both <- mclapply(all.cells,
                     function(x) {parallel_profiling(cell             = x,
                                                     unsampled_coords = unsampled_coords,
                                                     sampled_coords   = sampled_coords,
                                                     original_SDM     = original_SDM,
                                                     envStack        = envStack,
                                                     regform          = regform,
                                                     ntrees           = ntrees)},
                     mc.cores = mc.cores)
    both <- unlist(both)
    diff0 <- both[names(both) == "diff0"]
    diff1 <- both[names(both) == "diff1"]
  }

  ### standardise differences
  diffmax <- max(c(diff0, diff1), na.rm = TRUE)
  diffmin <- min(c(diff0, diff1), na.rm = TRUE)

  diff0 <- (diff0 - diffmin) / (diffmax- diffmin)
  diff1 <- (diff1 - diffmin) / (diffmax- diffmin)

  ### compile and return data frame
  return(data.frame(x     = unsampled_coords$x,
                    y     = unsampled_coords$y,
                    id    = unsampled_coords$id,
                    diff0 = diff0,
                    diff1 = diff1))
}
