
## Code taken from:
## https://rstudio.github.io/reticulate/articles/package.html
## Meant to delay loading of python environmnet so user can set the python environment

# python modules to use 
scrublet <- NULL
scipy <- NULL
sparse <- NULL
numpy <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scrublet <<- reticulate::import("scrublet", delay_load = TRUE)
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
  sparse <<- reticulate::import("scipy.sparse", delay_load = TRUE)
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
}


#' @name sctk_py_install
#' @title Installs Python packages
#' @description Install all Python packages used in the \code{\link{singleCellTK}} package
#' using \code{\link[reticulate]{py_install}} from package \code{\link{reticulate}}.
#' @param packages List of packages to install
#' @param ... Other parameters to pass to \code{\link[reticulate]{py_install}}. Can be useful for
#' selecting which environment to install these packages into. If no parameters are supplied,
#' then packages will be installed into default Python environment chosen by \code{\link{reticulate}}.
#' @export
sctk_py_install = function(packages = c("scipy", "scrublet"), ...) {
  reticulate::py_install(packages, ...)
}