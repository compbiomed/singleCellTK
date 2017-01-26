# Single Cell TK (name will change)

## Installation

singleCellTK is under development. You can install the development version from github:

```
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK")
```

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
