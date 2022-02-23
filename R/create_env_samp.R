####################################################################################################
### create_env_samp
### Generate an environmental variable based on a subsample of another variable
###
### 01/2022
### Charlie Marsh (charlie.marsh@mailbox.org)
###
### envBase   = the data frame for the base variable that will be subsampled (output of create_env)
### model     = model for variogram; default = "Sph" (see ?vgm for options)
### psill     = (partial) sill of the variogram model; default = 1.5 (see ?vgm for options)
### rangeFun  = name of function (in quotes!) to generate random number for the range parameter of
###             the variogram model. E.g. function() exp(runif(1, 1, 6))
### propSamp  = propotion of cells to sample in each submodel (default = 0.5)
### dep1      = multiplied by the psill in submodels (default = 1)
###
### output:
###   data frame of 'x' and 'y' coordinates and (non-standardised) environmental values ('sim1')
###
####################################################################################################

#' @importFrom gstat gstat vgm
#' @importFrom stats predict
create_env_samp <- function(envBase,
                            model    = "Sph",
                            psill    = 1.5,
                            rangeFun = "vrangeFun",
                            propSamp = 0.5,
                            dep1     = 1) {

  ### sumEnv is for error catching when NAs are sometimes predicted
  sumEnv <- NA
  while(is.na(sumEnv)) {
    ### create random subset to base kriging on
    subSamp <- envBase[sample.int(n = nrow(envBase),
                                  size = round(nrow(envBase) * propSamp, 0)), ]

    ### create a few other random variable based on the subset of the base variable
    gEnvSamp <- gstat(formula = sim1 ~ 1,   # note that sim1 is the dependent
                      locations = ~x + y,
                      dummy = FALSE,
                      beta = 1,
                      model = vgm(psill = psill * dep1,
                                  model = model,
                                  range = round(do.call(rangeFun, list()))),
                      nmax = 20,
                      data = subSamp) # the data is the subset

    ### predict for the full grid
    suppressWarnings(env <- predict(gEnvSamp, newdata = envBase,
                                    nsim = 1, debug.level = 0))
    sumEnv <- sum(env$sim1)
  }

  return(env)
}
