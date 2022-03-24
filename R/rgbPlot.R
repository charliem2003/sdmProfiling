################################################################################
### rgbPlot
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
#' @param dual_col The colour (either \code{"green"} (default) or
#'  \code{"purple"} for the dual-leverage (top-right) corner
#' @param absence_col  The colour (either \code{"red"} (default) or
#'  \code{"blue"} for the absence-leverage (bottom-right) corner
#' @param save_plots TRUE or FALSE whether to save the final plots to a
#'  TableGrob object
#'
#' @returns a tableGrob of two raster plots:
#'  1) a legend plot that includes the profiling values as well as a background
#'     of the RGB colours
#'  2) a raster of the RGB colours for the profilied cells
#'
#'  Optionally, if save_plots = TRUE will output the tableGrob.
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
#'                                                ntrees     = 100),
#'                         parallel = FALSE)
#'
#' ### simple plot of points
#' rgbPlot(profile = profile, origSDM = spMod)
#'
#' ### you can specify how to map the red-green-blue colours
#' ### the top right (dual_col) can be either green (default) or purple
#' rgbPlot(profile = profile, origSDM = spMod,
#'         dual_col = "purple", absence_col = "red")
#'
#' ### the bottom right (absence_col) can be either red (default) or blue
#' rgbPlot(profile = profile, origSDM = spMod,
#'         dual_col = "purple", absence_col = "blue")

################################################################################

#' @export
#' @importFrom raster brick raster cellFromXY
#' @importFrom stats dist
#' @importFrom scales rescale
#' @importFrom RStoolbox ggRGB
#' @importFrom gridExtra grid.arrange
#' @import ggplot2

rgbPlot <- function(profile,
                    origSDM,
                    dual_col = "green",
                    absence_col = "red",
                    save_plots = FALSE
) {
  requireNamespace("ggplot2")

  ### leverage values and get cell IDs
  levs <- profile[, c("x", "y", "redundancy", "presence_leverage",
                      "absence_leverage", "dual_leverage")]
  levs$id <- cellFromXY(origSDM, levs[, c("x", "y")])

  ### leverage values are measured as closeness to the corner
  ### reverse this for rgb plot to make it distance from the corner
  levs$redundancy <- sqrt(2) - levs$redundancy
  levs$absence_leverage <- sqrt(2) - levs$absence_leverage
  levs$presence_leverage <- sqrt(2) - levs$presence_leverage
  levs$dual_leverage <- sqrt(2) - levs$dual_leverage

  ### rescale to between 0 and 1
  levs[, 3:6] <- apply(levs[, 3:6], 2,
                       function(x) scales::rescale(x, to = c(0, 1),
                                                   from = c(0, sqrt(2))))

  ### make raster for each leverage type with values as rgb intensity
  redund <- pres <- abs <- dual <- raster(ext = extent(origSDM),
                                          resolution = res(origSDM),
                                          vals = NA)
  pres[levs$id] <- levs$presence_leverage * 255
  abs[levs$id]  <- levs$absence_leverage * 255
  redund[levs$id] <- levs$redundancy * 255
  dual[levs$id] <- levs$dual_leverage * 255

  ### background plot
  pRGB <- ggplot() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          panel.grid = element_blank()) +
    labs(title = "RGB map") +
    lims(x = c(extent(origSDM)[1], extent(origSDM)[2]),
         y = c(extent(origSDM)[3], extent(origSDM)[4])) +
    coord_fixed(expand = FALSE)

  ### rgb raster - depending on colour choices
  if(dual_col == "green" & absence_col == "red") {
    pRGB <- pRGB + ggRGB(brick(pres, redund, abs), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "green" & absence_col == "blue") {
    pRGB <- pRGB + ggRGB(brick(abs, redund, pres), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "purple" & absence_col == "red") {
    pRGB <- pRGB + ggRGB(brick(pres, dual, abs), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "purple" & absence_col == "blue") {
    pRGB <- pRGB + ggRGB(brick(abs, dual, pres), 1, 2, 3, ggLayer = TRUE)
  }

  ##############################################################################
  ### create a legend plot

  ### generate 50 x 50 grid of each leverage measure
  ### note: measured as distance from corner like above
  width <- 50
  dd <- data.frame(x = rep(seq(1, 50, length.out = width), width),
                   y = rep(rev(seq(0, 50, length.out = width)), each = width),
                   diff0 = rep(seq(0, 1, length.out = width), width),
                   diff1 = rep(rev(seq(0, 1, length.out = width)), each = width))
  dd$red <- apply(cbind(dd$diff0, dd$diff1), 1,
                  function(x) dist(rbind(c(0, 0), c(x[1], x[2]))))
  dd$pres <- apply(cbind(dd$diff0, dd$diff1), 1,
                   function(x) dist(rbind(c(0, 1), c(x[1], x[2]))))
  dd$abs <- apply(cbind(dd$diff0, dd$diff1), 1,
                  function(x) dist(rbind(c(1, 0), c(x[1], x[2]))))
  dd$dual <- apply(cbind(dd$diff0, dd$diff1), 1,
                   function(x) dist(rbind(c(1, 1), c(x[1], x[2]))))

  ### widthcale to between 0 and 1 and convert to rgb intensity
  dd[, 5:8] <- apply(dd[, 5:8], 2,
                     function(x) scales::rescale(x, to = c(0, 1),
                                                 from = c(0, sqrt(2)))) * 255

  ### generate rasters for each leverage value
  redLg  <- raster(nrow = width, ncol = width,
                   xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = dd$red)
  presLg <- raster(nrow = width, ncol = width,
                   xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = dd$pres)
  absLg  <- raster(nrow = width, ncol = width,
                   xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = dd$abs)
  dualLg <- raster(nrow = width, ncol = width,
                   xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = dd$dual)

  ### base plot
  pLeg <- ggplot() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, color = "black"),
          panel.grid = element_blank()) +
    labs(x = bquote(italic(Diff)[0]), y = bquote(italic(Diff)[1]),
         title =  "Profile plot") +
    lims(x = c(0, 1), y = c(0, 1)) +
    coord_fixed(expand = FALSE)

  ### add in rgb raster - depending on colour choices
  if(dual_col == "green" & absence_col == "red") {
    pLeg <- pLeg + ggRGB(brick(presLg, redLg, absLg), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "green" & absence_col == "blue") {
    pLeg <- pLeg + ggRGB(brick(absLg, redLg, presLg), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "purple" & absence_col == "red") {
    pLeg <- pLeg + ggRGB(brick(presLg, dualLg, absLg), 1, 2, 3, ggLayer = TRUE)
  }
  if(dual_col == "purple" & absence_col == "blue") {
    pLeg <- pLeg + ggRGB(brick(absLg, dualLg, presLg), 1, 2, 3, ggLayer = TRUE)
  }

  ### add in points from profile plot
  pLeg <- pLeg +
    geom_point(data = profile, aes(x = diff0, y = diff1),
               cex = 1.25, alpha = 0.5)

  ### generate final plot
  pBoth <- grid.arrange(pLeg, pRGB, nrow = 1)
  pBoth

  ### save plot
  if(save_plots == TRUE) { return(pBoth) }
}
