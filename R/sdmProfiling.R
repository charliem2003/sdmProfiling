################################################################################
### sdmProfiling
### 01/2022
#'
#' Run SDM profiling for unsampled sites
#'
#' @author Charlie Marsh (charlie.marsh@@mailbox.org) & Yoni Gavish
#'
#' @param unsampledCoords data frame of unsampled cells for which we wish to
#'  run SDM profiling on, with column names "x", "y" and "ID"
#' @param sampledCoords data frame of sampled cells used to build original SDM,
#'  with column names "presence" (0 = absence, 1 = presence), "x" and "y"
#' @param origSDM raster file of the probability of occurrences of the original
#'  SDM built using the sampled cells (i.e. those in \code{'sampledCoords'})
#' @param envStack raster stack of environmental variables for building SDM
#' @param sdmFun name of modelling function to use. The default,
#'  \code{'sdmModelling'}, is the built-in random forest model. The user can
#'  use their own defined function as long as the output is a single raster
#'  object of the predicted probability of occurrences
#' @param sdmFunArgs a list with the arguments for the sdmFun (which is called
#'  through \code{do.call}). The first argument should be the coordinates with
#'  to build the model, and will be automatically generated and replaced the
#'  first list item. See the examples for how it works with the default random
#'  forest SDM function.
#' @param modFormula model formula if using built-in random forest model
#' @param ntrees number of trees for generating forest if using built-in
#'  random forest model
#' @param parallel TRUE or FALSE to whether to run profiling on multiple cores
# #' @param mc.cores if parallel = TRUE then the number of cores to share jobs to
#'
#' @return A data frame containing results of the SDM profile where each row is
#'  an unsampled cell and eight columns:
#'  the 'x' and 'y' coordinates of the unsampled cells we profiled;
#'  the sum standardised channge in probability of occurrencces when the
#'  unsampled point was changed to an absence ('diff0') or presence ('diff1');
#'  the leverage values ('redundancy', 'presence_leverage', 'absence_leverage'
#'  and  'dual_leverage') measured as the distance to the four corners of the
#'  profile plot.
#'
#' @examples
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
#' form <- paste0("presence ~ ", paste(names(envSet), collapse = " + "), "+ I(",
#'                 paste(names(envSet), collapse = " ^ 2) + I("), " ^ 2)")
#'
#' ### run the initial model
#' spMod <- sdmModelling(samples = sampPts,
#'                       envStack = envSet,
#'                       modFormula = form,
#'                       ntrees = 500,
#'                       plot = FALSE)
#'
#' ### a random set of 500 points to profile
#' unsampPts <- data.frame(x = runif(500, 1, 100), y = runif(500, 1, 100))
#' unsampPts <- unsampPts[!paste(unsampPts$x, unsampPts$y) %in%
#'                          paste(sampPts$x, sampPts$y), ]
#'
#' profile <- sdmProfiling(unsampledCoords = unsampPts,
#'                         sampledCoords   = sampPts,
#'                         origSDM         = spMod,
#'                         envStack        = envSet,
#'                         sdmFun          = "sdmModelling",
#'                         sdmFunArgs      = list(samples    = NULL,
#'                                                envStack   = envSet,
#'                                                modFormula = form,
#'                                                ntrees     = 50))
#'
################################################################################

#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
sdmProfiling <- function(unsampledCoords,
                         sampledCoords,
                         origSDM,
                         envStack,
                         sdmFun = "sdmModelling",
                         sdmFunArgs = list(samples    = NULL,
                                           envStack   = NULL,
                                           modFormula = NULL,
                                           ntrees     = 500),
                         parallel = FALSE#,
                         # mc.cores = mc.cores,
                         ) {

  ##############################################################################
  ### tidy up names
  names(sampledCoords) <- tolower(names(sampledCoords))
  sampledCoords <- sampledCoords[, c("presence", "x", "y")]

  names(unsampledCoords) <- tolower(names(unsampledCoords))
  unsampledCoords <- unsampledCoords[, c("x", "y")]

  ##############################################################################
  ## data storage - vector of sum standardised difference from original SDM
  diff0 <- rep(NA, nrow(unsampledCoords))
  diff1 <- rep(NA, nrow(unsampledCoords))

  ##############################################################################
  ### run profiling on unsampled points
  if(parallel == FALSE) {
    pb <- txtProgressBar(min = 0, max = nrow(unsampledCoords), style = 3)

    ### loop through unsampled coords and run SDM each as presence and absence
    for(cell in 1:nrow(unsampledCoords)) {

      ##########################################################################
      ### change point to absence and run new SDM

      ### add unsampled cell as absence to the data frame
      coords0 <- rbind(sampledCoords,
                       data.frame(presence = 0,
                                  x = unsampledCoords[cell, "x"],
                                  y = unsampledCoords[cell, "y"]))

      ### change the first list element to the new coordinates
      sdmFunArgs0 <- sdmFunArgs
      sdmFunArgs0[[1]] <- coords0

      ### call SDM function
      pred0 <- do.call(sdmFun, sdmFunArgs0)

      ##########################################################################
      ## change point to presence and run new SDM

      ### add unsampled cell as presence to the data frame
      coords1 <- rbind(sampledCoords,
                       data.frame(presence = 1,
                                  x = unsampledCoords[cell, "x"],
                                  y = unsampledCoords[cell, "y"]))

      ### change the first list element to the new coordinates
      sdmFunArgs1 <- sdmFunArgs
      sdmFunArgs1[[1]] <- coords1

      ### call SDM function
      pred1 <- do.call(sdmFun, sdmFunArgs1)

      ##########################################################################
      ### calculate deviance from original SDM and save to storage

      dev0 <- abs(pred0@data@values - origSDM@data@values)
      dev1 <- abs(pred1@data@values - origSDM@data@values)
      dev0 <- dev0[!is.na(dev0)]
      dev1 <- dev1[!is.na(dev1)]

      ##########################################################################
      ### Save to storage

      diff0[cell] <- sum(as.vector(dev0), na.rm = TRUE)
      diff1[cell] <- sum(as.vector(dev1), na.rm = TRUE)

      setTxtProgressBar(pb, cell)
    }
  }
  #
  # if(parallel == TRUE) {
  #   source("parallel_profiling.R")
  #   require(parallel)
  #   require(parallelsugar)
  #
  #   all.cells <- 1:length(unsampledCoords[, 1])
  #   both <- mclapply(all.cells,
  #                    function(x) {parallel_profiling(cell             = x,
  #                                                    unsampledCoords = unsampledCoords,
  #                                                    sampledCoords   = sampledCoords,
  #                                                    original_SDM     = original_SDM,
  #                                                    envStack        = envStack,
  #                                                    regform          = regform,
  #                                                    ntrees           = ntrees)},
  #                    mc.cores = mc.cores)
  #   both <- unlist(both)
  #   diff0 <- both[names(both) == "diff0"]
  #   diff1 <- both[names(both) == "diff1"]
  # }

  ##############################################################################
  ### standardise differences

  diffmax <- max(c(diff0, diff1), na.rm = TRUE)
  diffmin <- min(c(diff0, diff1), na.rm = TRUE)

  diff0 <- (diff0 - diffmin) / (diffmax- diffmin)
  diff1 <- (diff1 - diffmin) / (diffmax- diffmin)

  ##############################################################################
  ### calculate leverage values

  ### redundancy (bottom left corner)
  redund <- apply(cbind(diff0, diff1), 1,
                  function(x) sqrt(2) - dist(rbind(c(0, 0), c(x[1], x[2]))))

  ### presence-leverage (top left corner)
  pres <-   apply(cbind(diff0, diff1), 1,
                  function(x) sqrt(2) - dist(rbind(c(0, 1), c(x[1], x[2]))))

  ### absence-leverage (bottom right cornre)
  abs <-    apply(cbind(diff0, diff1), 1,
                  function(x) sqrt(2) - dist(rbind(c(1, 0), c(x[1], x[2]))))

  ### dual leverage (top right corner)
  dual <-   apply(cbind(diff0, diff1), 1,
                  function(x) sqrt(2) - dist(rbind(c(1, 1), c(x[1], x[2]))))

  ##############################################################################
  ### compile and return data frame
  return(data.frame(x                 = unsampledCoords$x,
                    y                 = unsampledCoords$y,
                    diff0             = diff0,
                    diff1             = diff1,
                    redundancy        = redund,
                    presence_leverage = pres,
                    absence_leverage  = abs,
                    dual_leverage     = dual))
}
