################################################################################
### create_sp
### 01/2022

#' Generate a virtual species based on a set of environmental variables
#'
#' @description Takes a raster stack of environmental variables and fits a
#'  user-defined function to a (subset) of these. This is then subsampled to
#'  generate a new conditional Gaussian simulation to generate a probability of
#'  occurrence surface, and converted to presence-absence for a defined
#'  prevalence. The idea is that the virtual species is related only to a
#'  subset of modelling variables, and that some modelling variables are
#'  somewhat related to those actually important to the species.
#'
#' @details The function first generates a base probability surface based upon
#'  a user-specified relationship with the variables. It then randomly
#'  subsamples this surface and generates a new conditional Gaussian simulation
#'  based on this subsample. Finally, the cells with the highest values are
#'  defined as presences, based on the desired overall prevalence.
#'
#' If the variables are generated using \code{\link{create_env_nsets}}, the
#'  idea is that the species is related to a subset of the those variables
#'  (although not precisely as we subsample the base probability surface).
#'  However, when we run the sdm we will be using a number of variables
#'  unrelated to the generation of the species probability surface, although
#'  some of those will be somewhat related to the relevant variables
#'
#'  So we are attempting to include several pitfalls of an SDM on a real
#'  species, that generally are excluded when generating virtual species which
#'  could lead to unwarranted confidence on model performance in the real world:
#'
#'  \enumerate{
#'    \item There are unknown processes unrelated to environmental variables
#'     not included in the model (hence the extra conditional Gaussian
#'     simulation step) and so we can never perfectly recover the distribution
#'
#'    \item We will likely include environmental variables unrelated to the
#'     species distribution in the model
#'
#'    \item Some of these unnecessary variables will be somewhat correlated
#'     with the variable that is important (and thus a model may select the
#'     'wrong' one)
#'
#'    \item There may be other spatial surfaces that may mislead the model,
#'     such as spatial variation in sampling effort
#'  }
#'
#' @author Charlie Marsh (charlie.marsh@@mailbox.org) & Yoni Gavish, based on
#'  the original code from
#'  http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
#'
#' @param envStack stack of environmental variables to generate base
#'  probability surface (i.e. the output from \code{\link{create_env_nsets}} or
#'  \code{\link{create_env_set}}, although any raster stack will work)
#' @param spFun a function used in \code{raster::calc} (in quotes) for
#'  generating the base probability surface from the raster stack of
#'  environmental variables, which is then subsampled for the variogram that
#'  generates the final surface. The function should take the form of the layer
#'  number in the raster stack, and how they should be treated. E.g.
#'  \code{"x[1] * x[6] * x[8]"} would mean multiplying layers 1, 6 and 8
#' @param spModel model for variogram; default = "Sph" (see
#'  \code{\link[gstat]{vgm}} for options)
#' @param spPsill (partial) sill of the variogram model; default = 1.5 (see
#'  \code{\link[gstat]{vgm}} for options)
#' @param spRange range parameter of the variogram; default = 100 (see
#'  \code{\link[gstat]{vgm}} for details)
#' @param propSamp proportion of cells to sample in each submodel (default = 0.5)
#' @param prev proportion of cells to be defined as a presence
#'
#' @return A raster stack of the virtual species. The first layer 'prob' is the
#'  probability of occurrence (range = 0-1), and 'pa' is the presence-absence
#'  distribution based on the desired prevalence.
#'
#' @return A raster stack of environmental variables. For each variable all
#'  values are standardised between 0 and 1.
#'
#' @references Variations on this method have been used to generate virtual
#'  species in:
#'
#' Gavish, Y., Marsh, C.J., Kuemmerlen, M., Stoll, S., Haase, P., Kunin, W.E.,
#'  2017. Accounting for biotic interactions through alpha-diversity
#'  constraints in stacked species distribution models. Methods in Ecology and
#'  Evolution 8, 1092–1102. https://doi.org/10.1111/2041-210X.12731
#'
#' Marsh, C.J., Gavish, Y., Kunin, W.E., Brummitt, N.A., 2019. Mind the gap:
#'  Can downscaling Area of Occupancy overcome sampling gaps when assessing
#'  IUCN Red List status? Diversity and Distributions 025, 1832–1845.
#'  https://doi.org/10.1111/ddi.12983
#'
#' @examples
#'
#' ### first generate sets of related environmental variables
#'
#' set.seed(9999)
#' envSet <- create_env_nsets(cellDims = c(100, 100),
#'                            sets     = c(4, 4, 3, 1),
#'                            model    = "Sph",
#'                            psill    = 1.5,
#'                            dep1     = 1,
#'                            rangeFun = function() exp(runif(1, 1, 6)),
#'                            propSamp = 0.25)
#'
#' plot(envSet)
#'
#' # the species relationship to the environmental variables is defined through
#' # spFun (although that surface is then subsampled for a new variogram).
#'
#' # In spFun the bracketed numbers refers to the layer number of the variable
#' # In this example we multiply the 1st variable of each environmental set, but
#' # the last set is not used (perhaps it is a nuisance sampling effort surface)
#'
#' # if you want you can define the function as it's own variable, or you can
#' # specify it within the function argument (as in the later examples)
#' envFun <- "x[1] * x[5] + x[9]"
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = envFun,
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 50,
#'                 propSamp = 0.5,
#'                 prev     = 0.1)
#'
#' # the output is a raster stack which can be plotted
#' plot(sp)
#'
#' # spatial autocorrelation can be increased with a higher spRange value
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = "x[1] * x[5] * x[9]",
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 500,
#'                 propSamp = 0.5,
#'                 prev     = 0.1)
#' plot(sp)
#'
#' # the higher the propSamp the less the species will deviate from the surface
#' # defined in spFun. E.g.
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = "x[1] * x[5] * x[9]",
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 50,
#'                 propSamp = 0.95,
#'                 prev     = 0.1)
#' # add the spFun defined surface
#' sp$env <- calc(envSet, function(x) x[1] * x[5] * x[9])
#' plot(sp)
#'
#' # prevalence in the final presence-absence map can be controlled using prev
#' sp <- create_sp(envStack = envSet,
#'                 spFun    = "x[1] * x[5] * x[9]",
#'                 spModel  = "Sph",
#'                 spPsill  = 1,
#'                 spRange  = 500,
#'                 propSamp = 0.5,
#'                 prev     = 0.5)
#' plot(sp)
#'
################################################################################

#' @export
#' @importFrom raster calc ncell stack xyFromCell
#' @importFrom sp "gridded<-"
#' @importFrom stats predict
create_sp <- function(envStack,
                      spFun    = "x[1] * x[6] * x[8]",
                      spModel  = "Sph",
                      spPsill  = 1,
                      spRange  = 100,
                      propSamp = 0.5,
                      prev     = 0.1) {

  ### calculate base prob layer using spFun
  spFun <- as.function(list(str2lang(spFun)))
  formals(spFun) <- alist(x = )
  baseProb <- calc(envStack, fun = spFun)

  ### create random subset to base kriging on
  cells <- sample.int(n = ncell(baseProb),
                      size = round(ncell(baseProb) * propSamp, 0))

  ### convert to data frame
  baseProb <- data.frame(xyFromCell(baseProb, cells),
                         sim1 = baseProb@data@values[cells])
  xy <- expand.grid(1:dim(envStack)[1], 1:dim(envStack)[2])
  names(xy) <- c("x","y")

  ##############################################################################
  ### create probability of occurrenec surface
  sp <- gstat(formula = sim1 ~ 1,
              locations = ~x + y,
              dummy = FALSE,
              beta = 1,
              model = vgm(psill = spPsill,
                          model = spModel,
                          range = spRange),
              nmax = 20,
              data = baseProb)

  ### predict for the full grid
  sp <- predict(sp, newdata = xy, nsim = 1)

  ### standardize between 0 and 1
  sp$sim1 <- (sp$sim1 - min(sp$sim1, na.rm = TRUE)) /
    (max(sp$sim1, na.rm = TRUE) - min(sp$sim1, na.rm = TRUE))

  ##############################################################################
  ### convert to presence-absence

  sp$cellID <- 1:nrow(sp)

  ### number of cells to occupy
  ncells <- ncell(envStack)
  ncells <- ncells * prev

  ### top x cells
  paCells <- sp[order(sp$sim1, decreasing = TRUE), ]
  paCells <- paCells[1:ncells, ]

  ### add pa info to sp
  sp$pa <- 0
  sp$pa[sp$cellID %in% paCells$cellID] <- 1

  ##############################################################################
  ### convert to raster

  ### grid
  gridded(sp) = ~x + y

  ### convert to raster
  spRas <- stack(raster(sp["sim1"]),
                 raster(sp["pa"]))

  ### name and save
  names(spRas) <- c("prob", "presence")
  return(spRas)
}
