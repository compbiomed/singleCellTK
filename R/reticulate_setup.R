
## Code taken from:
## https://rstudio.github.io/reticulate/articles/package.html
## Meant to delay loading of python environmnet so user can set the python environment

# python modules to use
scrublet <- NULL
scipy <- NULL
sparse <- NULL
numpy <- NULL
scnrm <- NULL
sc <- NULL
bbknn <- NULL
pkg_resources <- NULL
ad <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scrublet <<- reticulate::import("scrublet", delay_load = TRUE)
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
  sparse <<- reticulate::import("scipy.sparse", delay_load = TRUE)
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
  scnrm <<- reticulate::import("scanorama", delay_load = TRUE)
  sc <<- reticulate::import("scanpy", delay_load = TRUE)
  bbknn <<- reticulate::import("bbknn", delay_load = TRUE)
  pkg_resources <<- reticulate::import('pkg_resources',delay_load = TRUE)
  ad <<- reticulate::import('anndata',delay_load = TRUE,convert = FALSE)
}

#' @name sctkPythonInstallConda
#' @title Installs Python packages into a Conda environment
#' @description Install all Python packages used in the \code{\link{singleCellTK}} package
#' using \code{\link[reticulate]{conda_install}} from package \code{\link{reticulate}}. This
#' will create a new Conda environment with the name \code{envname} if not already present.
#' Note that Anaconda or Miniconda already need to be installed on the local system.
#' @param envname Character. Name of the conda environment to create.
#' @param conda Character. Path to conda executable. Usue "auto" to find conda using the PATH and other conventional install locations. Default 'auto'.
#' @param packages Character Vector. List of packages to install from Conda.
#' @param pipPackages Character Vector. List of packages to install into the Conda environment using 'pip'.
#' @param selectConda Boolean. Run \code{\link[singleCellTK]{selectSCTKConda}} after installing all packages to select the Conda environment. Default TRUE.
#' @param forge Boolean. Include the Conda Forge repository.
#' @param pipIgnoreInstalled Boolean. Ignore installed versions when using pip. This is TRUE by default so that specific package versions can be installed even if they are downgrades.
#'        The FALSE option is useful for situations where you don't want a pip install to attempt an overwrite of a conda binary package (e.g. SciPy on Windows which is very difficult
#'        to install via pip due to compilation requirements).
#' @param pythonVersion Passed to \code{python_version} variable in \code{\link[reticulate]{conda_install}}. Default NULL.
#' @param ... Other parameters to pass to \code{\link[reticulate]{conda_install}}.
#' @return None. Installation of Conda environment.
#' @examples
#' \dontrun{
#' sctkPythonInstallConda(envname = "sctk-reticulate")
#' }
#' @seealso See \code{\link[reticulate]{conda_create}} for more information on creating a Conda environment.
#' See \code{\link[reticulate]{conda_install}} for more description of the installation parameters.
#' See \url{https://rstudio.github.io/reticulate/} for more information on package \code{\link{reticulate}}.
#' See \code{\link[singleCellTK]{selectSCTKConda}} for reloading the Conda environment if R is restarted without
#' going through the whole installation process again.
#' See \url{https://docs.conda.io/en/latest/} for more information on Conda environments.
#' @export
sctkPythonInstallConda <- function(envname = "sctk-reticulate",
                                   conda = "auto",
                                   packages = c("scipy", "numpy", "astroid", "six"),
                                   pipPackages = c("scrublet", "scanpy", "bbknn", "scanorama", "anndata"),
                                   selectConda = TRUE,
                                   forge = FALSE,
                                   pipIgnoreInstalled = TRUE,
                                   pythonVersion = NULL,
                                   ...) {

  path <- reticulate::conda_create(envname = envname, packages = "python", conda = conda)

  for(i in packages) {
    reticulate::conda_install(envname = envname, packages = i, conda = conda,
                  pip = FALSE, pip_ignore_installed = pipIgnoreInstalled,
                  python_version = pythonVersion, ...)
  }

  reticulate::conda_install(envname = envname, packages = pipPackages,
                pip = TRUE, pip_ignore_installed = pipIgnoreInstalled,
                python_version = pythonVersion, ...)

  if(isTRUE(selectConda)) selectSCTKConda(envname = envname)

  invisible(path)
}



#' @name sctkPythonInstallVirtualEnv
#' @title Installs Python packages into a virtual environment
#' @description Install all Python packages used in the \code{\link{singleCellTK}} package
#' using \code{\link[reticulate]{virtualenv_install}} from package \code{\link{reticulate}}. This
#' will create a new virtual environment with the name \code{envname} if not already present.
#' @param envname Character. Name of the virtual environment to create.
#' @param packages Character Vector. List of packages to install.
#' @param selectEnvironment Boolean. Run \code{\link[singleCellTK]{selectSCTKVirtualEnvironment}} after installing all packages to select the virtual environment. Default TRUE.
#' @param python The path to a Python interpreter, to be used with the created virtual environment. When NULL, the Python interpreter associated with the current session will be used. Default NULL.
#' @return None. Installation of virtual environment.
#' @examples
#' \dontrun{
#' sctkPythonInstallVirtualEnv(envname = "sctk-reticulate")
#' }
#' @seealso See \code{\link[reticulate]{virtualenv_create}} for more information on creating a Conda environment.
#' See \code{\link[reticulate]{virtualenv_install}} for more description of the installation parameters.
#' See \url{https://rstudio.github.io/reticulate/} for more information on package \code{\link{reticulate}}.
#' See \code{\link[singleCellTK]{selectSCTKVirtualEnvironment}} for reloading the virtual environment if R is restarted without
#' going through the whole installation process again.
#' @export
sctkPythonInstallVirtualEnv <- function(envname = "sctk-reticulate",
                                        packages = c("scipy", "numpy", "astroid", "six", "scrublet", "scanpy", "scanorama", "bbknn", "anndata"),
                                        selectEnvironment = TRUE,
                                        python = NULL) {

  path <- reticulate::virtualenv_create(envname = envname, python = python)

  for(i in packages) {
    reticulate::virtualenv_install(envname = envname, packages = i, ignore_installed = TRUE)
  }

  if(isTRUE(selectEnvironment)) selectSCTKVirtualEnvironment(envname = envname)

  invisible(path)
}


#' @name selectSCTKConda
#' @title Selects a Conda environment
#' @description Selects a Conda environment with Python packages used in \code{\link{singleCellTK}}.
#' @param envname Character. Name of the conda environment to activate.
#' @return None. Selects Conda environment.
#' @examples
#' \dontrun{
#' sctkPythonInstallConda(envname = "sctk-reticulate", selectConda = FALSE)
#' selectSCTKConda(envname = "sctk-reticulate")
#' }
#' @seealso \code{\link[reticulate]{conda-tools}} for more information on using Conda environments with package \code{\link{reticulate}}.
#' See \url{https://rstudio.github.io/reticulate/} for more information on package \code{\link{reticulate}}.
#' @export
#' @seealso See \code{\link[singleCellTK]{sctkPythonInstallConda}} for installation of Python modules into a Conda environment.
#' See\code{\link[reticulate]{conda-tools}} for more information on using Conda environments with package \code{\link{reticulate}}.
#' See \url{https://rstudio.github.io/reticulate/} for more information on package \code{\link{reticulate}}.
#' See \url{https://docs.conda.io/en/latest/} for more information on Conda environments.
selectSCTKConda <- function(envname = "sctk-reticulate") {
  condaList <- reticulate::conda_list()
  ix <- condaList$name == envname

  if(!any(ix)) {
    stop(paste0("Environment '", envname, "', not found. Run sctkPythonInstallConda(envname = '", envname, "') to install Python packages into a conda environmanet with this name."))
  }
  if(sum(ix) > 1) {
    envs <- paste(condaList[ix,"python"], collapse="\n")
    warning(paste0("More than one Conda environment detected with the name '", envname, "'. Selecting the first one in the list:\n", envs))
  }

  reticulate::use_condaenv(condaenv = envname, required = TRUE)
}



#' @name selectSCTKVirtualEnvironment
#' @title Selects a virtual environment
#' @description Selects a virtual environment with Python packages used in \code{\link{singleCellTK}}
#' @param envname Character. Name of the virtual environment to activate.
#' @return None. Selects virtual environment.
#' @examples
#' \dontrun{
#' sctkPythonInstallVirtualEnv(envname = "sctk-reticulate", selectEnvironment = FALSE)
#' selectSCTKVirtualEnvironment(envname = "sctk-reticulate")
#' }
#' @seealso See \code{\link[singleCellTK]{sctkPythonInstallVirtualEnv}} for installation of Python modules into a virtual environment.
#' See\code{\link[reticulate]{virtualenv-tools}} for more information on using virtual environments with package \code{\link{reticulate}}.
#' See \url{https://rstudio.github.io/reticulate/} for more information on package \code{\link{reticulate}}.
#' @export
selectSCTKVirtualEnvironment <- function(envname = "sctk-reticulate") {
  res <- reticulate::virtualenv_list()
  ix <- res == envname

  if(!any(ix)) {
    stop(paste0("Environmnet '", envname, "', not found. Run selectSCTKVirtualEnvironment(envname = '", envname, "') to install Python packages into a virtual environmanet with this name."))
  }
  if(sum(ix) > 1) {
    warning(paste0("More than one virtual environment detected with the name '", envname, "'. Selecting the first one in the list."))
  }

  reticulate::use_virtualenv(res[which(ix)[1]], required = TRUE)
}
