# comixr
This is an R package for fitting gaussian mixture models with shared components across 1-dimensional data sets.

## Installation

```R
# install.packages("devtools")
devtools::install_github("daniel-wells/comixr")
```

For examples of how to use the function please see the vignette.

Functions to fit the models and plot resulting components are in the [R/](R/) directory.

[vignettes/](vignettes/) contains examples of how to use the GMM-EM functions ([vignette.Rmd](vignettes/vignette.Rmd)) as well as a full mathematical description of the model and the update equations ([GMM-EM.Rmd](vignettes/GMM-EM.Rmd)).