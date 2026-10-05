# Makefile for singleCellTK, an r-bioc-dev-standards package. Run `make help` to list
# targets. Recipe lines must start with a tab.
#
# The standard targets (test, check, bioccheck, ...) live in the shared
# standards.mk in r-bioc-dev-standards, not here. It's downloaded to a cache
# on first use and refreshed at every Claude session start.

# Package settings (see standards.mk).
# The Suggests list is large and some entries are optional; this mirrors
# _R_CHECK_FORCE_SUGGESTS_ in .BBSoptions.
FORCE_SUGGESTS = FALSE

# Extra targets that only people may run (list them in AGENTS.md too), e.g.
# PEOPLE_ONLY := clean site-deploy
PEOPLE_ONLY := site

# Guards. Settings allow `make test-one` with any arguments, so this file
# checks them itself, before anything else runs:
# - FILTER may contain only letters, digits, '.', '_' and '-' (checked on
#   its raw value, before make expands it);
# - no other variable may be set on the command line, since one could
#   change which makefile, shell, or source is used (environment variables
#   still work, e.g. R_BIOC_STANDARDS_REF=<branch> make help);
# - test-one and article must be the only target;
# - PEOPLE_ONLY targets refuse to run from Claude Code, which sets
#   CLAUDECODE in its shell (settings deny rules match only exact commands).
_ok_chars := a b c d e f g h i j k l m n o p q r s t u v w x y z \
  A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 0 1 2 3 4 5 6 7 8 9 . _ -
_strip_ok = $(if $(2),$(call _strip_ok,$(subst $(firstword $(2)),,$(1)),$(wordlist 2,$(words $(2)),$(2))),$(1))
_filter_error := FILTER may contain only letters, digits, '.', '_' and '-'
ifneq ($(findstring $$,$(value FILTER)),)
  $(error $(_filter_error))
endif
ifneq ($(words $(value FILTER)),$(if $(value FILTER),1,0))
  $(error $(_filter_error))
endif
ifneq ($(call _strip_ok,$(value FILTER),$(_ok_chars)),)
  $(error $(_filter_error))
endif
_cmdline_vars := $(filter-out FILTER,$(foreach v,$(.VARIABLES),$(if $(filter command line,$(origin $(v))),$(v))))
ifneq ($(_cmdline_vars),)
  $(error Only FILTER=<pattern> may be set on the make command line (found: $(_cmdline_vars)); set other variables in the environment or above)
endif
ifneq ($(filter test-one article,$(MAKECMDGOALS)),)
  ifneq ($(words $(MAKECMDGOALS)),1)
    $(error make $(firstword $(filter test-one article,$(MAKECMDGOALS))) must be run on its own)
  endif
endif
ifneq ($(CLAUDECODE),)
  ifneq ($(filter $(PEOPLE_ONLY),$(MAKECMDGOALS)),)
    $(error $(filter $(PEOPLE_ONLY),$(MAKECMDGOALS)) is for people only; ask the developer to run it)
  endif
endif

R_BIOC_STANDARDS_REF ?= v1
R_BIOC_STANDARDS_BASE ?= https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)
STANDARDS_MK := $(or $(XDG_CACHE_HOME),$(HOME)/.cache)/r-bioc-dev-standards/$(R_BIOC_STANDARDS_REF)/standards.mk

include $(STANDARDS_MK)

$(STANDARDS_MK):
	@mkdir -p "$(@D)"
	curl -fsSL --max-time 30 "$(R_BIOC_STANDARDS_BASE)/shared/standards.mk" -o "$@.tmp"
	@mv -f "$@.tmp" "$@"

# Extra targets for this package go below, each listed in AGENTS.md. To
# replace a standard target, define it here; make warns that it overrides
# the shared recipe.

PKG := singleCellTK
# Build artifacts go to a fresh, uniquely named directory outside the repo.
# Nothing is ever cleaned or deleted.
BUILD_DIR = $(shell mktemp -d "$${TMPDIR:-/tmp}/sctk-build.XXXXXX")

.PHONY: build site app test-app

build:  ## Build the source tarball (printed path is outside the repo)
	@d="$(BUILD_DIR)"; cd "$$d" && R CMD build "$(CURDIR)" && echo "Tarball: $$d/$$(ls *.tar.gz)"

site:  ## Full local pkgdown build into docs/ (site owner only)
	Rscript -e 'pkgdown::build_site()'

app:  ## Launch the Shiny app from the working tree
	Rscript -e 'devtools::load_all(); singleCellTK()'

test-app:  ## shinytest2 smoke suite (slow; excluded from make test)
	@if [ -d tests/app ]; then \
	  Rscript -e 'testthat::test_dir("tests/app")'; \
	else \
	  echo "No shinytest2 suite yet (tests/app/ does not exist). See AGENTS.md."; \
	fi
