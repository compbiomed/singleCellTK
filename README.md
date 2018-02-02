[![Build Status](https://travis-ci.org/compbiomed/singleCellTK.svg?branch=master)](https://travis-ci.org/compbiomed/singleCellTK)

# Single Cell TK

## Installation

singleCellTK is under development. You can install the development version from github:

```r
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK")
```
> Note: Some package dependencies require Bioconductor v3.6,
> https://bioconductor.org/install/

#### Troubleshooting Installation

For the majority of users, the commands above will install the latest version
of the singleCellTK without any errors. Rarely, you may encounter an error due
to previously installed versions of some packages that are required for the
singleCellTK. If you encounter an error during installation, use the commands
below to check the version of Bioconductor that is installed:

```r
source("https://bioconductor.org/biocLite.R")
biocVersion()
```

If the version number is not 3.6 or higher, you must upgrade Bioconductor to
install the toolkit:

```r
biocLite("BiocUpgrade")
```

After you install Bioconductor 3.6 or higher, you should be able to install the
toolkit using `devtools::install_github("compbiomed/singleCellTK")`. If you
still encounter an error, ensure your Bioconductor packages are up to date by
running the following command.

```r
biocValid()
```

If the command above does not return `TRUE`, run the following command to
update your R packages:

```r
biocLite()
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
