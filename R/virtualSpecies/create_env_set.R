####################################################################################################
### create_env_set
### Generate a set of environmental variables (can be just 1 variable). The first variable is
### generated using create_env and then subsequent variables are generated based on sampling a
### subset of the 1st variable. Each new variable is generated using a random range parameter
### generated from the rangeFun. Therefore the output is a set of variables that are related but
### with some random variation between them.
###
### 01/2022
### Charlie Marsh (charlie.marsh@mailbox.org)
###
### cellDims  = vector containing variable dimensions (n cell width, n cell height)
### nPerSet   = the total number of environmental variables to generate
### model     = model for variogram; default = "Sph" (see ?vgm for options)
### psill     = (partial) sill of the variogram model; default = 1.5 (see ?vgm for options)
### rangeFun  = name of function (in quotes!) to generate random number for the range parameter of
###             the variogram model. E.g. function() exp(runif(1, 1, 6))
### propSamp  = propotion of cells to sample in each submodel (default = 0.5)
### dep1      = multiplied by the psill in submodels (default = 1)
###
### output:
###   a raster stack of environmental variable set
###
####################################################################################################

create_env_set <- function(cellDims = c(100, 100),
                           nPerSet  = 5,
                           model    = "Sph",
                           psill    = 1.5,
                           rangeFun = "vrangeFun",
                           propSamp = 0.25,
                           dep1     = 1) {
  require(gstat)
  require(sp)
  require(raster)

  ### generate base variable
  envBase <- create_env(cellDims = cellDims,
                        model    = model,
                        psill    = psill,
                        rangeFun = rangeFun)

  ### generate variables based on subsamples of base variable and save as list
  if(nPerSet > 1) {
    envSet <- list(envBase)
    for(i in 2:nPerSet) {
      envSet[[i]] <- create_env_samp(envBase,
                                     model    = model,
                                     psill    = psill,
                                     dep1     = dep1,
                                     rangeFun = rangeFun,
                                     propSamp = propSamp)
    }

    ### stanrdise and convert to raster
    for(i in 1:nPerSet) {
      ### stanardise variables between 0 and 1
      envSet[[i]]$sim1 <- (envSet[[i]]$sim1 - min(envSet[[i]]$sim1, na.rm = TRUE)) /
        (max(envSet[[i]]$sim1, na.rm = TRUE) - min(envSet[[i]]$sim1, na.rm = TRUE))

      ### grid
      gridded(envSet[[i]]) = ~x + y

      ### convert to raster
      envSet[[i]] <- raster(envSet[[i]])
    }
  }

  ### if only one variable then return the base variable
  if(nPerSet == 1) {
    envSet <- envBase

    ### standardise variables between 0 and 1
    envSet$sim1 <- (envSet$sim1 - min(envSet$sim1, na.rm = TRUE)) /
      (max(envSet$sim1, na.rm = TRUE) - min(envSet$sim1, na.rm = TRUE))

    ### grid
    gridded(envSet) = ~x + y

    ### convert to raster
    envSet <- raster(envSet)
  }

  ### stack and rename
  envSet <- stack(envSet)
  names(envSet) <- paste0("var_", 1:nPerSet)
  return(envSet)
}
