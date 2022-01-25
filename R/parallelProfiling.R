

parallel_profiling <- function(cell,
                               unsampled_coords = unsampled_coords,
                               sampled_coords   = sampled_coords,
                               original_SDM     = original_SDM,
                               env.stack        = env.stack,
                               regform          = regform,
                               ntrees           = ntrees) {

  ## change point to absence
  coords.0 <- rbind(sampled_coords,
                    data.frame(presence = 0,
                               x = unsampled_coords[cell, "x"],
                               y = unsampled_coords[cell, "y"],
                               id = unsampled_coords[cell, "id"]))
  pred.0 <- sdmModelling(coords = coords.0,
                         env.stack = env.stack,
                         mod_formula = regform,
                         ntrees = ntrees)

  ## change point to presence
  coords.1 <- rbind(sampled_coords,
                    data.frame(presence = 1,
                               x = unsampled_coords[cell, "x"],
                               y = unsampled_coords[cell, "y"],
                               id = unsampled_coords[cell, "id"]))
  pred.1 <- sdmModelling(coords = coords.1,
                         env.stack = env.stack,
                         mod_formula = regform,
                         ntrees = ntrees)

  dev0 <- abs(pred.0@data@values - original_SDM@data@values)
  dev1 <- abs(pred.1@data@values - original_SDM@data@values)
  dev0 <- dev0[!is.na(dev0)]
  dev1 <- dev1[!is.na(dev1)]

  return(data.frame(diff0 = sum(as.vector(dev0), na.rm = TRUE),
                    diff1 = sum(as.vector(dev1), na.rm = TRUE)))
}
