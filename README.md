# sdmProfiling

An R package for carrying out SDM profiling,a spatially-explicit sensitivity analysis for species distribution models which seeks to assess the leverage that sampled or unsampled cells have on the overall model.

Details can be found in the manuscript (in review):

**SDM profiling: a spatially-explicit species distribution model evaluation tool and its potential applications**
*Charles J. Marsh, Yoni Gavish, Mathias Kuemmerlen, Stefan Stoll, Peter Haase, William E. Kunin*


## Installation

`## If you don't have it installed then install devtools`

`install.packages(devtools)`

`devtools::install_github("charliem2003/sdmProfiling")`

## Details

Detailed tutorial to come but see details and examples in help pages.

SDM profiling for a given species and presence-absence sdm is done through the `sdmProfiling` function. A simple function for carrying out a random forest model is availabe with `sdmModelling`, although the user may provide their own modelling function inside `sdmProfiling` (see help page).

Results from the sdm profiling can be visualised using `profilePlot` (scatterplot where each cell is plotted as the value for standardised difference when the cell was added as a presence (diff1), and when the cell was added as an absence (diff0)), `leveragePlot` (map where colour ramp is based on values of one of the leverage measures) and `rgbPlot` (map where the values of 3 leverage measures have been visualsied in RGB space).

Use `create_sp` for generating virtual species using sets of environmental layers generated through `create_env_nsets` and `create_env_set`.

