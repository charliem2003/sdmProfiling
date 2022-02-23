####################################################################################################
### create_env
### Generate single environmental variable
###
### 01/2022
### Charlie Marsh (charlie.marsh@mailbox.org)
###
### cellDims  = vector containing variable dimensions (n cell width, n cell height)
### model     = model for variogram; default = "Sph" (see ?vgm for options)
### psill     = (partial) sill of the variogram model; default = 1.5 (see ?vgm for options)
### vrange    = name of function (in quotes!) to generate random number for the range parameter of
###             the variogram model. E.g. function() exp(runif(1, 1, 6))
###
### output:
###   data frame of 'x' and 'y' coordinates and (non-standardised) environmental values ('sim1')
###
####################################################################################################

#' @importFrom gstat gstat vgm
#' @importFrom stats predict
create_env <- function(cellDims = c(100, 100),
                       model    = "Sph",
                       psill    = 1.5,
                       rangeFun = "vrangeFun") {
  ### create structure
  xy <- expand.grid(1:cellDims[1],
                    1:cellDims[2])
  names(xy) <- c("x","y")

  ### sumEnv is for error catching when NAs are sometimes predicted
  sumEnv <- NA
  while(is.na(sumEnv)) {
    ### the gstat model and prediction for the entire grid
    gEnv <- gstat(formula = z ~ 1,
                  locations = ~x + y,
                  dummy = TRUE,
                  beta = 1,
                  model = vgm(psill = psill,
                              model = model,
                              range = round(do.call(rangeFun, list()))),
                  nmax = 20)

    ### predict across whole grid
    suppressWarnings(env <- predict(gEnv, newdata = xy,
                                    nsim = 1, debug.level = 0))
    sumEnv <- sum(env$sim1)
  }

  return(env)
}
