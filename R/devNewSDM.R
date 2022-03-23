################################################################################
## devNewSDM
## Charlie Marsh (charlie.marsh@mailbox.org)
## 18/03/2022
##
## Runs a new SDM with the unsampled cell added as either a presence or absence.
## Returns absolute (non-standardised) deviance from the original SDM for the
## new SDM
##
################################################################################

devNewSDM <- function(cellState, # either 0 (absence) or 1 (presence)
                      cell,      # rownumber of unsampledCoords to profile
                      sampledCoords   = sampledCoords,
                      unsampledCoords = unsampledCoords,
                      origSDM         = origSDM,
                      sdmFunArgs      = sdmFunArgs,
                      sdmFun          = sdmFun)
{

  ##########################################################################
  ### change point to absence and run new SDM

  ### add unsampled cell as absence to the data frame
  coordsNew <- rbind(sampledCoords,
                     data.frame(presence = cellState,
                                x = unsampledCoords[cell, "x"],
                                y = unsampledCoords[cell, "y"]))

  ### change the first list element to the new coordinates
  sdmFunArgsNew <- sdmFunArgs
  sdmFunArgsNew[[1]] <- coordsNew

  ### call SDM function
  predNew <- do.call(sdmFun, sdmFunArgsNew)

  ##########################################################################
  ### calculate deviance from original SDM and save to storage

  dev <- abs(predNew@data@values - origSDM@data@values)
  dev <- dev[!is.na(dev)]

  ##########################################################################
  ### Save to storage

  diff <- sum(as.vector(dev), na.rm = TRUE)
  return(diff)
}
