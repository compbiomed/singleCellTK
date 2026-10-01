# 0002. Support both harmony interfaces in runHarmony() by version dispatch

- **Status:** Proposed
- **Date:** 2026-09-20
- **Deciders:** TODO (maintainer approval required)

## Context

`runHarmony()` called `harmony::HarmonyMatrix()`. harmony removed that
function in version 1.0.0 (current release: 2.0.5), so the function fails for
anyone with a current harmony installed. The Bioconductor 3.24 build report
and the GitHub Actions checks both reported it:

- `checking dependencies in R code ... WARNING: Missing or unexported object:
  'harmony::HarmonyMatrix'`
- `checking Rd cross-references ... WARNING: Missing link(s) in Rd file
  'runHarmony.Rd': '[harmony]{HarmonyMatrix}'`

Because `R-CMD-check.yaml` uses `error_on = "warning"`, these two warnings
were failing CI on every platform. The example is wrapped in `\dontrun{}`, so
the breakage never surfaced as a test failure.

The replacement, `harmony::RunHarmony()`, differs in three ways that matter:
`max.iter.harmony` is now `max_iter`; it no longer performs PCA, so
`do_pca`/`npcs` are gone; and tuning arguments such as `epsilon.cluster` moved
into `.options = harmony_options()`.

The lab requires that users on older harmony installations keep working.

## Decision

`runHarmony()` checks `utils::packageVersion("harmony")` at run time and
dispatches to one of two internal helpers:

- `.runHarmonyLegacy()` (harmony < 1.0.0) calls `HarmonyMatrix()` with exactly
  the arguments used before this change, so behaviour on old installations is
  unchanged.
- `.runHarmonyCurrent()` (harmony >= 1.0.0) calls `RunHarmony()` with
  `max_iter`, and computes the PCA with `scater::calculatePCA()` when the user
  passes a full-size assay, since harmony no longer does it. That call passes
  `ntop = Inf, scale = TRUE` so the PCA uses every feature of the selected
  assay, scaled, as the legacy `HarmonyMatrix(do_pca = TRUE)` did. scater's
  defaults (`ntop = 500`, `scale = FALSE`) would silently drop features and
  produce an embedding matching neither the legacy path nor the package's own
  `scaterPCA()`.

The legacy helper looks the function up with
`utils::getFromNamespace("HarmonyMatrix", "harmony")` rather than writing
`harmony::HarmonyMatrix`. This is required: `R CMD check` inspects code
statically, so a literal namespaced call to a function that no longer exists
still raises the WARNING, and CI would stay red.

The same run-time lookup is deliberately **not** applied to
`harmony::RunHarmony()`, which is written as a normal namespaced call. On a
machine with harmony < 1.0.0 installed, `R CMD check` therefore reports the
mirror-image `Missing or unexported object: 'harmony::RunHarmony'`. That is
accepted: checks are run by developers and CI against current dependencies,
nobody is expected to run `R CMD check` against a 2020-era harmony, and
keeping the normal `::` form for the supported path preserves the static
checking that would otherwise catch a misspelled name. Run-time behaviour on
old harmony is unaffected, since that branch is never evaluated there.

The user-facing arguments of `runHarmony()` are unchanged. The documentation
states which extra `...` arguments apply to which harmony version, and that
corrected embeddings are not expected to match between versions.

## Alternatives considered

- **Replace `HarmonyMatrix()` outright (no legacy path).** One code path and
  simpler to maintain, and it matches Bioconductor's expectation that a
  package works with current dependencies. Rejected: it breaks existing lab
  environments pinned to harmony 0.1.x.
- **Detect old harmony and stop with an upgrade message.** Also one code path,
  with no silent failure. Rejected for the same reason.
- **`suppressWarnings()` or removing the Rd link only.** Hides the problem;
  `runHarmony()` would stay broken on current harmony.

## Consequences

- Results change for users who upgrade harmony: the algorithm differs between
  0.1.x and 2.x, so corrected embeddings will not reproduce older runs. This
  is recorded in `NEWS.md`.
- The legacy path cannot be exercised in CI, which installs current harmony,
  so it can rot unnoticed. It was verified once during this change by
  installing a stub package presenting the 0.1.1 interface, confirming the
  arguments passed are identical to the previous implementation. Repeat that
  check before relying on the legacy path again.
- With harmony >= 1.0.0, the PCA for the assay input is now computed by
  singleCellTK (`scater::calculatePCA`, all features, scaled) rather than by
  harmony's internal scaling plus `irlba`. The inputs match what the legacy
  path used, but the implementations differ, so embeddings from the assay path
  are not expected to be identical across versions.
- Using every feature is slower than scater's 500-feature default on large
  datasets. That cost is accepted to keep the assay path behaving as before.
- `ntop = Inf` also makes the `nComponents` guard correct. With scater's
  500-feature default the guard compared `nComponents` against
  `min(dim(mat))` while the PCA ran on at most 500 features, so a request for
  more components than the retained features produced a silently smaller
  embedding: 1000 features x 800 cells with `nComponents = 600` returned 499
  columns, with no warning from the wrapper. Using every feature makes
  `min(dim(mat))` the real bound, and the same call now returns the requested
  600 columns. The unit test asserts the column count for this reason.
- No dependency changes: `scater` and `utils` were already used, and `harmony`
  stays in Suggests.
- `docs/reference/runHarmony.html` was regenerated for this change with
  `pkgdown::build_reference(topics = "runHarmony")` rather than a full
  `make site`, so that page is built by pkgdown 2.2.1 while the rest of the
  committed site is 2.1.1: its version badge and `lang` attribute differ from
  the other pages. Accepted so the published help page matches the code; the
  inconsistency clears the next time the site owner rebuilds the whole site.
- Both CI warnings are cleared on machines with harmony >= 1.0.0, which is
  what CI and the Bioconductor builders use: `R CMD check` drops from
  3 WARNINGs to 1 (the remaining one is the unrelated scuttle/scran
  deprecation warning). Checking the package on an old-harmony machine still
  reports one WARNING, as described above.
