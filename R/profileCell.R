################################################################################
## profileCell
## Charlie Marsh (charlie.marsh@mailbox.org)
## 18/03/2022
##
## Profiles an unsampled cell with the cell as presence and absence. Runs a new
## SDM with the unsampled cell added as either a presence or absence Returns a
## 2 item list of the absolute (non-standardised) deviance from the  original
## SDM for the new SDM for an absence (cellDiff0) and presence (cellDiff1)
##
################################################################################

profileCell <- function(cell,      # rownumber of unsampledCoords to profile
                        sampledCoords   = sampledCoords,
                        unsampledCoords = unsampledCoords,
                        origSDM         = origSDM,
                        sdmFunArgs      = sdmFunArgs,
                        sdmFun          = sdmFun)
{

  ##########################################################################
  ### change point to absence, run new SDM and calculate absolute deviance

  cellDiff0 <- devNewSDM(cellState       = 0,
                         cell            = cell,
                         sampledCoords   = sampledCoords,
                         unsampledCoords = unsampledCoords,
                         origSDM         = origSDM,
                         sdmFunArgs      = sdmFunArgs,
                         sdmFun          = sdmFun)

  ##########################################################################
  ### change point to presence, run new SDM and calculate absolute deviance

  cellDiff1 <- devNewSDM(cellState       = 1,
                         cell            = cell,
                         sampledCoords   = sampledCoords,
                         unsampledCoords = unsampledCoords,
                         origSDM         = origSDM,
                         sdmFunArgs      = sdmFunArgs,
                         sdmFun          = sdmFun)

  return(list(cellDiff0 = cellDiff0,
              cellDiff1 = cellDiff1))
}
