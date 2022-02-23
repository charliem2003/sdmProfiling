library(testthat)
library(sdmProfiling)

### generate vars
set.seed(999)
sets <- c(5, 4, 3, 1)
envSet <- create_env_nsets(cellDims = c(100, 100),
                           sets     = sets,
                           model    = "Sph",
                           psill    = 1.5,
                           dep1     = 1,
                           rangeFun = function() exp(runif(1, 1, 6)),
                           propSamp = 0.25)

####################################################################################################
### raster characteristics

### is it a raster stack
test_that("Produces a raster stack", {
  expect_s4_class(envSet, "RasterStack")
})

### are the names correct
varNames <- unlist(sapply(1:length(sets), FUN = function(x) paste0("var_", x, ".", 1:sets[x])))
test_that("Variable names are correct", {
  expect_equal(names(envSet), varNames)
})

####################################################################################################
### checking layer values

### function to loop through layers and apply required function
layerLoop <- function(envSet, fun) {
  unlist(lapply(1:nlayers(envSet),
                function(x) do.call(fun, list(envSet@layers[[x]]@data@values))))
}

### are there NAs
test_that("No layers contain NAs", {
  layerNAs <- layerLoop(envSet, is.na)
  expect_false(any(layerNAs),
               info = paste0())
})

### are all values between 0 and 1
test_that("Each layer has min = 0 and min = 1", {
  expect_equal(layerLoop(envSet, min), rep(0, nlayers(envSet)))
  expect_equal(layerLoop(envSet, max), rep(1, nlayers(envSet)))
})

### sum of values across all rasters
test_that("Sum across all layer values = expected", {
  layerSums <- layerLoop(envSet, sum)
  expect_equal(round(sum(layerSums), 2), 63396.23)
})

