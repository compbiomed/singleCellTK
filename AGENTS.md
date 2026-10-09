# singleCellTK: notes for coding agents

The shared development standards (r-bioc-dev-standards) load automatically
at session start in Claude Code. This file adds only what is specific to
this package. Other agents (for example through GEMINI.md) don't get them
automatically: read
https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/v1/standards.md
(or the cached copy in `~/.cache/r-bioc-dev-standards/v1/`) before
starting, and follow it; those agents also aren't bound by
`.claude/settings.json`.

## About

singleCellTK (SCTK) is a Bioconductor package for single-cell RNA-seq
analysis: import, QC, doublet detection, ambient RNA removal,
normalization, batch correction, dimensionality reduction, clustering,
markers, differential expression, cell type labeling, and pathway analysis.
The same functions are reached three ways: the R console, an interactive
Shiny GUI, and a command-line QC pipeline. It also produces HTML reports
via R Markdown.

## Layout

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
- `dev/`: maintainer docs, ADRs, plans, and hooks. `.github/CONTRIBUTING.md`
  has the full layout.
- No compiled code (no `src/`).

## Object model

No classes of its own. Everything is built on `SingleCellExperiment`:

- Every analysis function takes and returns a `SingleCellExperiment`
  (argument `inSCE`, assay chosen with `useAssay`).
- Results are stored on the object: new assays (`outAssayName`),
  `reducedDims`, `altExps`, `colData` columns, or `metadata`. `run*`
  wrappers don't return bare matrices.
- Use `assay()`, `reducedDim()`, `altExp()`, `colData()`, `metadata()`,
  never `@`.
- Naming is camelCase (`runNormalization`, `useAssay`); functions are verbs
  (`run*`, `plot*`, `import*`, `get*`).

## Tests

- Fixtures: `data/` (for example `scExample`, `sceBatches`) and
  `inst/extdata/`.
- The full suite is slow. While developing, use `make test-one` on the
  file for the feature area (for example `FILTER=qc`).
- Python-backed tests (scrublet, scanpy, AnnData export) skip when no
  Python environment is available through reticulate.
- There is no `shinytest2` suite yet. The old `inst/shiny/tests/shinytest.R`
  uses the deprecated `shinytest` package; `make test-app` reports this.

## Extra make targets

- `build`: builds the source tarball in a new temporary folder outside the
  repo and prints its path. Safe.
- `app`: launches the Shiny app from the working tree. Safe; run it in the
  background and stop it when done.
- `test-app`: runs the `shinytest2` suite in `tests/app/` (not yet present).
  Safe.
- `site`: full pkgdown build into the committed `docs/`. People only (the
  site owner).

## Setup in a new worktree

- R and Bioconductor must match the target cycle in
  <https://bioconductor.org/config.yaml>. Never assume the pairing.
- Install dependencies with
  `BiocManager::install("singleCellTK", dependencies = TRUE)`, or from
  `DESCRIPTION`. The dependency list is large.
- Python-backed features need a Python environment via reticulate.
- macOS needs `fftw` for some dependencies.
- `renv.lock` and `renv_R4.0.lock` are not used for development; don't run
  `renv::restore()`.

## Related packages

- `decontX` (campbio/decontX): singleCellTK imports it and calls
  `decontX::decontX()` and `decontXcounts()` in `runDecontX()`
  (`R/celda_decontX.R`), which `runCellQC()` runs by default. Changes to
  these functions in decontX, or to how singleCellTK calls them, should be
  checked in both packages.
- `celda` (campbio/celda): singleCellTK imports it for
  `celda::distinctColors()` (plotting helpers). Since celda 1.23.0,
  `celda::decontX()` only forwards to the decontX package.

## Overrides

- **Shiny app:** whether agents may change `inst/shiny/` is undecided
  (see `dev/ROADMAP.md`). Until the maintainer decides, propose app changes
  in the plan rather than making them. Reason: the app has almost no test
  coverage, so unreviewed changes there are hard to catch.

## Package notes

- CI: `R-CMD-check.yaml` (macOS, Windows, Ubuntu) and `BioC-check.yaml`
  (Ubuntu) install Bioconductor devel packages via `R_BIOC_VERSION`, to
  match the Bioconductor builders (ADR 0007). Which jobs are required and
  the coverage threshold are still TODO.
- Documentation is generated with roxygen2 8.1, which writes NAMESPACE in
  a multi-line `importFrom()` format. Use roxygen2 8.1 or later for
  `make docs`; an older roxygen2 rewrites the whole NAMESPACE.
- The website is served from `docs/`, which is committed. Never edit it;
  only the site owner rebuilds it with `make site`.
- Known debt (backlog, don't fix in unrelated changes): `make lint` reports
  about 20,300 existing lints (October 2026), roughly half in `R/` and half
  in `inst/`. The most common are indentation inside calls, line length,
  trailing whitespace, braces, and infix spacing. The code uses 2-space
  indentation, which `.lintr` matches. Only lines you change need to be
  lint-clean.
