####################################################################################################
### create_sp
### Generate a virtual species based on a set of environmental variables. If these variables are
###
###
### The function first
### generates a base probability surface based upon a user-specified relationship with the
### variables. It then randomly subsamples this surface and generates a new variogram based on this
### subsample. Finally, the cells with the highest values are defined as presences, based on the
### desired overall prevalence.
###
### If the variables are generated using create_env_nsets, the idea is that the species is related
### to a subset of the those variables (although not precisely as we subsample the base probability
### surface). However, when we run the sdm we will be using a number of variables unrelated to the
### generation of the surface, although some of those will be somewhat related to the relevant
### variables
###
### 01/2022
### Charlie Marsh (charlie.marsh@mailbox.org)
###
### envStack  = stack of environmental variables to generate base probability surface (i.e. the
###             output from create_env_nsets or create_env_set, although any raster stack will work)
### spFun     = a function used in raster::calc (in quotes) for generating the base probability
###             surface from the raster stack of environmental variables, which is then subsampled
###             for the variogram that generates the final surface. The function should take the
###             form of the layer number in the raster stack, and how they should be treated.
###             E.g. "x[1] * x[6] * x[8]" would mean multiplying layers 1, 6 and 8 of the stack
### spModel   = model for variogram; default = "Sph" (see ?vgm for options)
### spPsill   = (partial) sill of the variogram model; default = 1.5 (see ?vgm for options)
### spRange   = range parameter of the variogram; default = 100 (see ?vgm for options)
### propSamp  = propotion of cells to sample in each submodel (default = 0.5)
### prev      = proportion of cells to be defined as a presence
###
### output:
###   a raster stack of the virtual species. The first layer 'prob' is the probability of occurrence
###   (range = 0-1), and 'pa' is the presence-absence distribution based on the desired prevalence
###
####################################################################################################

create_sp <- function(envStack = envStack,
                      spFun    = "x[1] * x[6] * x[8]",
                      spModel  = "Sph",
                      spPsill  = 1,
                      spRange  = 100,
                      propSamp = 0.5,
                      prev     = 0.1) {
  require(gstat)
  require(sp)
  require(raster)

  ### calculate base prob layer using spFun
  spFun <- as.function(list(str2lang(spFun)))
  formals(spFun) <- alist(x = )
  baseProb <- calc(envStack, fun = spFun)

  ### create random subset to base kriging on
  cells <- sample.int(n = ncell(baseProb), size = round(ncell(baseProb) * propSamp, 0))

  ### convert to data frame
  baseProb <- data.frame(xyFromCell(baseProb, cells), sim1 = baseProb@data@values[cells])
  xy <- expand.grid(1:dim(envStack)[1], 1:dim(envStack)[2])
  names(xy) <- c("x","y")

  ##################################################################################################
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

  ##################################################################################################
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

  ##################################################################################################
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
