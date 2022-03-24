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
#' @param noCores if parallel = TRUE then the number of cores to share jobs to.
#'  If noCores is larger than values retrieved through parallel::detectCores(),
#'  then the values from parallel::detectCores() (the total number of cores on
#'  the machine) will be used. Be wary of memory requirements when thinking
#'  about how many cores to use/
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
#'                                                ntrees     = 50),
#'                         parallel = FALSE)
#'
################################################################################

#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach `%dopar%`
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
sdmProfiling <- function(unsampledCoords,
                         sampledCoords,
                         origSDM,
                         envStack,
                         sdmFun = "sdmModelling",
                         sdmFunArgs = list(samples    = NULL,
                                           envStack   = NULL,
                                           modFormula = NULL,
                                           ntrees     = 500),
                         parallel = FALSE,
                         noCores = 2
) {
  if(parallel == TRUE) {
    requireNamespace("parallel")
    requireNamespace("foreach")
    requireNamespace("doSNOW")
    # requireNamespace("doParallel")
  }

  ##############################################################################
  ### tidy up names
  names(sampledCoords) <- tolower(names(sampledCoords))
  sampledCoords <- sampledCoords[, c("presence", "x", "y")]

  names(unsampledCoords) <- tolower(names(unsampledCoords))
  unsampledCoords <- unsampledCoords[, c("x", "y")]

  ##############################################################################
  ### run profiling on unsampled points
  if(parallel == FALSE) {
    pb <- txtProgressBar(min = 0, max = nrow(unsampledCoords), style = 3)

    ##############################################################################
    ## data storage - vector of sum standardised difference from original SDM
    diff0 <- rep(NA, nrow(unsampledCoords))
    diff1 <- rep(NA, nrow(unsampledCoords))

    ### loop through unsampled coords and run SDM each as presence and absence
    for(cell in 1:nrow(unsampledCoords)) {
      ##########################################################################
      ### profile cell with point changed point to absenceand presence

      diffs <- profileCell(cell            = cell, #rownumber of unsampledCoords
                           sampledCoords   = sampledCoords,
                           unsampledCoords = unsampledCoords,
                           origSDM         = origSDM,
                           sdmFunArgs      = sdmFunArgs,
                           sdmFun          = sdmFun)

      ##########################################################################
      ### change point to presence, run new SDM and calculate absolute deviance

      diff0[cell] <- diffs$cellDiff0
      diff1[cell] <- diffs$cellDiff1

      setTxtProgressBar(pb, cell)
    }
  }

  if(parallel == TRUE) {
    ## calculate number of cores and set as appropriate
    noCoresToUse <- parallel::detectCores()
    noCoresToUse <- ifelse(noCoresToUse < noCores, noCoresToUse, noCores)

    ### initiate cluster
    cl <- parallel::makeCluster(noCoresToUse)
    doSNOW::registerDoSNOW(cl)

    ### set up progress bar using doSNOW (from stackoverflow.com
    ### how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r)
    pb <- txtProgressBar(1, nrow(unsampledCoords), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    diffs <- foreach::foreach(
      cell = 1:nrow(unsampledCoords),
      .combine = rbind,
      .options.snow = opts,
      .packages = c("sdmProfiling")) %dopar% {
        profileCell(cell            = cell, # rownumber of unsampledCoords
                    sampledCoords   = sampledCoords,
                    unsampledCoords = unsampledCoords,
                    origSDM         = origSDM,
                    sdmFunArgs      = sdmFunArgs,
                    sdmFun          = sdmFun)
      }

    ### stop cluster
    close(pb)
    parallel::stopCluster(cl)


    # ### initiate cluster
    # cl <- parallel::makeCluster(noCoresToUse)
    # doParallel::registerDoParallel(cl)
    #
    # diffs <- foreach::foreach(
    #   cell = 1:nrow(unsampledCoords),
    #   .combine = rbind,
    #   .packages = c("sdmProfiling")) %dopar% {
    #     profileCell(cell            = cell, # rownumber of unsampledCoords
    #                 sampledCoords   = sampledCoords,
    #                 unsampledCoords = unsampledCoords,
    #                 origSDM         = origSDM,
    #                 sdmFunArgs      = sdmFunArgs,
    #                 sdmFun          = sdmFun)
    #   }
    #
    # ### stop cluster
    # parallel::stopCluster(cl)


    ### extract results
    diff0 <- as.vector(unlist(diffs[, "cellDiff0"]))
    diff1 <- as.vector(unlist(diffs[, "cellDiff1"]))

    # all.cells <- 1:length(unsampled_coords[, 1])
    # both <- mclapply(all.cells,
    #                  function(x) {parallel_profiling(cell             = x,
    #                                                  unsampled_coords = unsampledCoords,
    #                                                  sampled_coords   = sampledCoords,
    #                                                  original_SDM     = origSDM,
    #                                                  env.stack        = envStack,
    #                                                  regform          = modFormula,
    #                                                  ntrees           = ntrees)},
    #                  mc.cores = mc.cores)
    # both <- unlist(both)
    # diff0 <- both[names(both) == "diff0"]
    # diff1 <- both[names(both) == "diff1"]
  }

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
