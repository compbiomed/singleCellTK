
## Code taken from:
## https://rstudio.github.io/reticulate/articles/package.html
## Meant to delay loading of python environmnet so user can set the python environment

# python modules to use 
scrublet <- NULL
scipy <- NULL

.onLoad <- function(libname, pkgname) {
  # delay load foo module (will only be loaded when accessed via $)
  scrublet <<- import("scrublet", delay_load = TRUE)
  scipy <<- import("scipy", delay_load = TRUE)
}