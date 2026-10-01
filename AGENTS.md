# AGENTS.md

## Campbell Lab Playbook (common across lab packages — v2.0, do not edit per-repo)

### Common commands
make test / make check / make bioccheck / make docs / make lint / make site
(See Makefile for definitions. These are the ONLY sanctioned entry points.
Run `make test` after every change; `make check` before opening a PR.)

### Git and PR workflow
- Branch from devel; all work lands via PR. Never push to devel or master.
- Use plan mode for any non-trivial change.
- Run /code-review before requesting human review.
- Every user-facing change gets a NEWS.md entry.

### Coding conventions
- Style enforced by lintr/styler (config in repo); <= 80-char lines (BiocCheck).
- roxygen2 owns man/ and NAMESPACE — NEVER hand-edit them.
- Use accessor functions, not @ slot access, outside class definition files.

### Documentation (pkgdown)
- The website is GENERATED. Improve docs by editing roxygen comments, vignettes,
  and _pkgdown.yml — never files under docs/ or the gh-pages branch.
- New exported functions MUST be added to the _pkgdown.yml reference index;
  verify with pkgdown::check_pkgdown().
- To preview one changed page: pkgdown::build_article("<name>") or
  build_reference_index(). NEVER run a full build_site() as verification —
  full site builds/deploys are a local maintainer action (make site-deploy).
- Files under vignettes/articles/ are pkgdown-only and NOT checked by
  R CMD check — knit locally when you edit them.

### Shiny app rules (packages with inst/shiny only)
- The app contains NO analysis logic. Server code only wires inputs to
  exported package functions and renders results. New app features are
  implemented as tested, exported functions first.
- Reactive logic is tested with shiny::testServer(); the golden path is
  covered by a small shinytest2 smoke suite (make test-app).
- UI changes are verified with a screenshot of the RUNNING app
  (make app + browser), not just passing tests.
- inst/ code is invisible to R CMD check — tests and lintr are the only
  guards; inst/shiny is included in the lint paths.

### Versioning and releases
- Bioconductor even/odd x.y.z scheme; releases ~April and ~October.
- Follow dev/RELEASE.md for the release checklist.

### Safety rules
- No structural refactors (file splits, DESCRIPTION dependency changes,
  class redesign) without an approved ADR — propose via a GitHub issue.
- Never commit secrets, tokens, or absolute local paths.
- Architectural decisions are recorded in dev/adr/ (see template and index
  there). Never store anything in docs/ — that is pkgdown build output.
- Maintainer docs (release, roadmap, audits) live in dev/, not the root.


## This package: singleCellTK

### Project overview
singleCellTK (SCTK) is a Bioconductor package for single-cell RNA-seq
analysis. It covers import, QC, doublet detection, ambient RNA removal,
normalization, batch correction, dimensionality reduction, clustering,
markers, differential expression, cell type labeling, and pathway analysis.
It offers three interfaces to the same functions: the R console, an
interactive Shiny GUI, and a command-line QC pipeline. It also produces
HTML reports via R Markdown.

### Repository map
- `R/`: all analysis logic. By prefix: `import*.R` (data import),
  `run*.R` (analysis wrappers), `plot*.R` (plotting),
  `*_doubletDetection.R`, `seuratFunctions.R`, `scanpyFunctions.R`
  (Python via reticulate). `R/singleCellTK.R` launches the app.
- `inst/shiny/`: the Shiny GUI (`server.R`, `ui.R`, one `ui_NN_*.R` file per
  tab, a few modules, `www/`).
- `inst/rmarkdown/`: HTML report templates. `inst/extdata/`: small example
  datasets for examples and tests.
- `exec/SCTK_runQC.R`: command-line QC pipeline.
- `vignettes/singleCellTK.Rmd`: the Bioconductor vignette.
  `vignettes/articles/`: pkgdown-only tutorials.
- `tests/testthat/`: one test file per feature area.
- `dev/`: maintainer docs and ADRs. `.github/CONTRIBUTING.md` has the full
  layout.

### Object model
- Every analysis function takes and returns a `SingleCellExperiment`
  (argument `inSCE`, assay chosen with `useAssay`).
- Results are stored on the object: new assays (`outAssayName`),
  `reducedDims`, `altExps`, `colData` columns, or `metadata`. `run*`
  wrappers don't return bare matrices.
- Use Bioconductor accessors (`assay()`, `reducedDim()`, `colData()`,
  `metadata()`), not `@` slot access.
- Naming is camelCase (`runNormalization`, `useAssay`); functions are verbs
  (`run*`, `plot*`, `import*`, `get*`).

### Environment setup
- R and Bioconductor must match the target cycle in
  <https://bioconductor.org/config.yaml>. Never assume the pairing.
- `BiocManager::install("singleCellTK", dependencies = TRUE)`, or install
  from `DESCRIPTION` for development. The dependency list is large; see
  `DESCRIPTION`.
- Python-backed features (scrublet, scanpy, AnnData export) need a Python
  environment via reticulate. Tests for them skip when it's unavailable.
- macOS needs `fftw` for some dependencies.

### Package-specific overrides of the common section
These differ from the common playbook for this repo, by maintainer decision:
- **Branches:** commits go on `bioc_release_2026_09`; PRs target `devel`.
  Nothing is committed until a person has reviewed the changes.
- **Website:** `docs/` stays committed on the main branches (no `gh-pages`
  migration). Never edit it by hand. Only the site owner rebuilds it
  (`make site`). There is no `make site-deploy` target.
- **Lint hook:** report-only lint of the edited file (`dev/hooks/`). Don't
  auto-restyle files. The existing lint backlog is fixed deliberately, never
  mixed into functional changes.
- **Agents may not edit** `Makefile`, `.claude/settings.json`, or
  `dev/hooks/`. People change them via a reviewed PR.
- **Never** run `rm -rf`, `rm -r`, `find -delete`, or `git clean`. Ask a
  person to delete things.

### Package-specific notes
- Shiny app in scope for agent-driven changes: **TODO (maintainer
  decision).** Until decided, propose app changes rather than making them.
- No `shinytest2` suite exists yet. The old `inst/shiny/tests/shinytest.R`
  uses the deprecated `shinytest` package. `make test-app` reports this.
- No compiled code (no `src/`).
- CI: which jobs are required, the coverage threshold, and the BiocCheck
  container/cron setup are all **TODO**. The existing workflows in
  `.github/workflows/` are unchanged.
- Test fixtures: `data/` (for example `scExample`, `sceBatches`) and
  `inst/extdata/`. The full test suite is slow.
- Known debt (backlog, don't fix in unrelated changes): `make lint` reports
  about 33,700 existing lints, roughly half in `R/` and half in
  `inst/shiny/`. Most are indentation (the code uses 2 spaces; `.lintr`
  enforces Bioconductor's 4), line length, and trailing whitespace. Only
  lines you change need to be lint-clean.
