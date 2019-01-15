# Single Cell TK

[![Travis build status](https://travis-ci.org/compbiomed/singleCellTK.svg?branch=master)](https://travis-ci.org/compbiomed/singleCellTK)
[![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)
[![BioC status](https://www.bioconductor.org/shields/build/release/bioc/singleCellTK.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/singleCellTK)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

## Installation

### Release Version

You can download the release version of the Single Cell Toolkit in
[Bioconductor v3.7](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html):

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("singleCellTK")
```

### Devel Version

You can download the development version of the Single Cell Toolkit in
[Bioconductor v3.8](https://bioconductor.org/packages/devel/bioc/html/singleCellTK.html)
or from this repository:

```r
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK")
```

### R 3.4 Version

If you are still running an earlier version of R instead of 3.5, you can install
the following version from this repository:

```r
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK", ref="r_3_4")
```

#### Troubleshooting Installation

For the majority of users, the commands above will install the latest version
of the singleCellTK without any errors. Rarely, you may encounter an error due
to previously installed versions of some packages that are required for the
singleCellTK. If you encounter an error during installation, use the commands
below to check the version of Bioconductor that is installed:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::version()
```

If the version number is not 3.6 or higher, you must upgrade Bioconductor to
install the toolkit:

```r
BiocManager::install()
```

After you install Bioconductor 3.6 or higher, you should be able to install the
toolkit using `devtools::install_github("compbiomed/singleCellTK")`. If you
still encounter an error, ensure your Bioconductor packages are up to date by
running the following command.

```r
BiocManager::valid()
```

If the command above does not return `TRUE`, run the following command to
update your R packages:

```r
BiocManager::install()
```

Then, try to install the toolkit again:

```r
devtools::install_github("compbiomed/singleCellTK")
```

If you still encounter an error, please [contact us](mailto:dfj@bu.edu) and
we'd be happy to help.

## Develop singleCellTK

To contribute to singleCellTK, follow these steps:

__Note__: singleCellTK is developed to eventually be added to bioconductor. The current
version of bioconductor works best with R >= 3.3.1 

1. Fork the repo using the "Fork" button above.
2. Download a local copy of your forked repository "```git clone https://github.com/{username}/singleCellTK.git```"
3. Open Rstudio
4. Go to "File" -> "New Project" -> "Existing Directory" and select your git repository directory

You can then make your changes and test your code using the Rstudio build tools.
There is a lot of information about building packages available here: http://r-pkgs.had.co.nz/.

Information about building shiny packages is available here: http://shiny.rstudio.com/tutorial/.

When you are ready to upload your changes, commit them locally, push them to your
forked repo, and make a pull request to the compbiomed repository.

Report bugs and request features on our [GitHub issue tracker](https://github.com/compbiomed/singleCellTK/issues).

Join us on [slack](https://compbiomed.slack.com/)!
