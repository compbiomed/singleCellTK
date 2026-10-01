# Canonical commands for singleCellTK: one answer to "how do I ...?" for
# humans, agents, and CI. Agents may run these targets but not edit this file.
#
# _R_CHECK_FORCE_SUGGESTS_=false mirrors .BBSoptions: the Suggests list is
# large, and some entries are optional.

PKG      := singleCellTK
# Build artifacts go to a fresh, uniquely named directory outside the repo.
# Nothing is ever cleaned or deleted.
BUILD_DIR = $(shell mktemp -d "$${TMPDIR:-/tmp}/sctk-build.XXXXXX")

.PHONY: help test check bioccheck build docs lint site app test-app

help:        ## list targets
	@grep -E '^[a-z-]+:.*## ' $(MAKEFILE_LIST) | sed 's/:.*## /\t/'

test:        ## fast loop: run after every change
	Rscript -e 'devtools::test()'

check:       ## full R CMD check: run before opening a PR
	_R_CHECK_FORCE_SUGGESTS_=false Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual"), build_args = "--no-manual", error_on = "warning", check_dir = tempfile("sctk-check-"))'

build:       ## build the source tarball (printed path is outside the repo)
	@d="$(BUILD_DIR)"; cd "$$d" && R CMD build "$(CURDIR)" && echo "Tarball: $$d/$$(ls *.tar.gz)"

bioccheck:   ## BiocCheck on the built tarball, plus BiocCheckGitClone on the checkout
	@d="$(BUILD_DIR)"; cd "$$d" && R CMD build "$(CURDIR)" && \
	  Rscript -e 'BiocCheck::BiocCheck(Sys.glob("*.tar.gz"), `quit-with-status` = FALSE); BiocCheck::BiocCheckGitClone("$(CURDIR)")'; \
	  echo "BiocCheck output: $$d"

docs:        ## regenerate man/ and NAMESPACE from roxygen
	Rscript -e 'devtools::document()'

lint:        ## lint R/ and inst/ (including inst/shiny) using .lintr
	Rscript -e 'lintr::lint_package()'

site:        ## full local pkgdown build into docs/ (site owner only)
	Rscript -e 'pkgdown::build_site()'

app:         ## launch the Shiny app from the working tree
	Rscript -e 'devtools::load_all(); singleCellTK()'

test-app:    ## shinytest2 smoke suite (slow; excluded from make test)
	@if [ -d tests/app ]; then \
	  Rscript -e 'testthat::test_dir("tests/app")'; \
	else \
	  echo "No shinytest2 suite yet (tests/app/ does not exist). See AGENTS.md."; \
	fi
