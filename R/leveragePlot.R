################################################################################
### leveragePlot
### 01/2022
#'
#' Plot a map with leverage values as RGB colours. The leverage values
#'  outputted from SDM profilinf for unsampled cells for three leverage
#'  categories are mapped on to the RGB values
#'
#' @author Charlie Marsh (charlie.marsh@@mailbox.org) & Yoni Gavish
#'
#' @param profile the output from \code{\link{sdmProfiling}}
#' @param origSDM a raster file of the probability of occurrences of the original
#'  SDM built using the sampled cells (i.e. those in \code{'sampledCoords'})
#' @param leverage_type which leverage type to plot. One of \code{"redundancy},
#'  \code{presence_leverage}, \code{absence_leverage} or \code{dual_leverage"}
#'  (default)
#' @param col the colour ramp to use (default = "red"). Low values are white
#' @param plot_type either \code{"points"} (default) or \code{"raster"}. Use
#'  \code{"points"} where you have profiled a subsample of unsampled cells. Each
#'  profiled cell will be represented by a point where point colour and size is
#'  the leverage value, and the background will be the original SDM map. If all
#'  unsampled cells have been profiled use \code{"raster"} will plot the
#'  profiled cells as a raster map, with cell colour representing the
#'  leverage value
#' @param save_plots TRUE or FALSE (default) whether to save the final plot to
#'  a ggplot object
#'
#' @returns Automatically generates a map of leverage values. If  Use
#'  \code{plot_type = "points"} each profiled cell will be represented by a
#'  point with colour and size relative to the leverage value, and the
#'  background will be the original SDM map. If all unsampled cells have been
#'  profiled use \code{"raster"} will plot the profiled cells as a raster map,
#'  with cell colour representing the leverage value.
#'
#'  Optionally, if save_plot = TRUE will output the ggplot object.
#'
#' @examples
#' set.seed(111)
#' envSet <- create_env_nsets(cellDims = c(25, 25),
#'                            sets     = c(2, 2, 2, 1),
#'                            model    = "Sph",
#'                            psill    = 1.5,
#'                            dep1     = 1,
#'                            rangeFun = function() exp(runif(1, 1, 6)),
#'                            propSamp = 0.1)
#'
#' ### generate a virtual species from the variables
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = "x[1] * x[3] * x[5]",
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 500,
#'                 propSamp = 0.5,
#'                 prev     = 0.25)
#'
#' ### an initial 'sample' of the species (assuming perfect detection)
#' sampPts <- data.frame(sampleRandom(sp$presence, 15, na.rm = TRUE, xy = TRUE))
#'
#' ### a formula to fit to random forest (additive for all vars + quadratics)
#' form <- paste0("presence ~ ", paste(names(envSet), collapse = " + "), "+ I(",
#'                paste(names(envSet), collapse = " ^ 2) + I("), " ^ 2)")
#'
#' ### run the initial model
#' spMod <- sdmModelling(samples = sampPts,
#'                       envStack = envSet,
#'                       modFormula = form,
#'                       ntrees = 500,
#'                       plot = FALSE)
#'
#' ### the full set of unsampled points to profile
#' unsampPts <- expand.grid(x = 1:25, y = 1:25)
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
#'                                                ntrees     = 100))
#'
#' ###
#' leveragePlot(profile = profile, origSDM = spMod, plot_type = "raster",
#'              leverage_type = "dual_leverage", col = "forest green")
#'
#' leveragePlot(profile = profile, origSDM = spMod, plot_type = "raster",
#'              leverage_type = "presence_leverage", col = "red")
#'
#'
#' #' ### the full set of unsampled points to profile
#' unsampPts <- data.frame(x = sample(1:25, 25), y = sample(1:25, 25))
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
#'                                                ntrees     = 100))
#'
#' leveragePlot(profile = profile, origSDM = spMod,
#'             leverage_type = "presence_leverage",
#'             col = "red", plot_type = "points")


################################################################################

#' @export
#' @importFrom raster ncell cellFromXY
#' @importFrom viridis viridis
#' @import ggplot2

leveragePlot <- function(profile,
                         origSDM,
                         leverage_type = "dual-leverage",
                         col = "red",
                         plot_type = "points",
                         save_plots = FALSE
) {
  requireNamespace("ggplot2")

  ### leverage values and get cell IDs
  levs <- profile[, c("x", "y", "redundancy", "presence_leverage",
                      "absence_leverage", "dual_leverage")]
  levs$id <- cellFromXY(origSDM, levs[, c("x", "y")])

  ### plot as points
  if(plot_type == "points") {
    sdm <- data.frame(xyFromCell(origSDM, 1:ncell(origSDM)),
                      prob = origSDM@data@values)

    pLev <- ggplot(profile, aes(x, y)) +
      theme(panel.background = element_rect(fill = "grey30", color = "black"),
            panel.grid = element_blank()) +
      coord_fixed(expand = FALSE) +
      labs(fill = "Prob. of occ.", colour = leverage_type, size = leverage_type) +
      scale_fill_gradientn(colours = viridis(100), limits = c(0, 1)) +
      geom_raster(data = sdm, aes(fill = prob))
    pLev <- pLev +
      scale_colour_gradientn(colours = c("white", col), limits = c(0, sqrt(2))) +
      geom_point(aes(size = get(leverage_type), colour = get(leverage_type)),
                 alpha = 1, shape = 19)#, stroke = 1)
    pLev
  }

  ### plot as raster
  if(plot_type == "raster") {
    pLev <- ggplot(profile, aes(x, y)) +
      theme(panel.background = element_rect(fill = "grey30", color = "black"),
            panel.grid = element_blank()) +
      coord_fixed(expand = FALSE) +
      labs(fill = leverage_type) +
      scale_fill_gradientn(colours = c("white", col), limits = c(0, sqrt(2))) +
      geom_raster(aes(fill = get(leverage_type)))
  }

  ### generate final plot
  print(pLev)

  ### save plot
  if(save_plots == TRUE) { return(pLev) }
}
