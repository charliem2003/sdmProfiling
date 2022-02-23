################################################################################
### profilePlot
### 01/2022
#'
#' Make a 'profile plot' - the Diff0 values plotted against Diff1 values
#'  outputted from SDM profiling for unsampled cells
#'
#' @author Charlie Marsh (charlie.marsh@@mailbox.org) & Yoni Gavish
#'
#' @param profile the output from \code{\link{sdmProfiling}}
#' @param points TRUE (default) or FALSE whether to plot the individual points
#' @param contours TRUE or FALSE (default) whether to plot contours of point
#'  density
#' @param density TRUE or FALSE (default) whether to plot a raster of point
#'  density (blue to red)
#' @param save_plot TRUE or FALSE (default) whether to save the ggplot object
#'
#' @returns Automatically generates a profile plot.
#'
#'  Optionally, if save_plot = TRUE will output the ggplot object.
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
#'                                                ntrees     = 100))
#'
#' ### simple plot of points
#' profilePlot(profile, points = TRUE, contours = FALSE, density = FALSE)
#'
#' ### also include a background of point density and contours
#' profilePlot(profile, points = TRUE, contours = TRUE, density = TRUE)
#'
################################################################################

#' @export
#' @importFrom grDevices colorRampPalette
#' @import ggplot2

profilePlot <- function(profile,
                        points    = TRUE,
                        contours  = FALSE,
                        density   = FALSE,
                        save_plot = FALSE
) {
  requireNamespace("ggplot2")

  basePlot <- ggplot(data = profile, aes(diff0, diff1)) +
    theme(panel.background = element_rect(fill = NA, color = "black", size = 1),
          axis.ticks = element_line(color = "black", size = 1),
          panel.grid = element_blank(),
          legend.position = "none") +
    expand_limits(x = c(0, 1), y = c(0, 1)) +
    coord_fixed(expand = FALSE) +
    labs(x = bquote(italic(Diff)[0]), y = bquote(italic(Diff)[1]))

  if(density == TRUE) {
    pal <- colorRampPalette(c("white", "white", "white", "dark blue",
                              "green", "yellow", "red", "dark red"))(100)
    basePlot <- basePlot +
      scale_fill_gradientn(colours = pal) +
      stat_density_2d(aes(fill = ..density.. ^0.2),
                      geom = "tile", contour = FALSE, n = 250)
  }

  if(contours == TRUE) {
    basePlot <- basePlot +
      geom_density_2d(col = "red")
  }

  if(points == TRUE) {
    basePlot <- basePlot +
      geom_point(alpha = 0.15)
  }

  basePlot <- basePlot +
    geom_abline(slope = 1, intercept = 0, lty = 2)

  ### generate final plot
  print(basePlot)

  ### save plot
  if(save_plot == TRUE) { return(basePlot) }
}


