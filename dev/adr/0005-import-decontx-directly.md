# 0005. Import decontX directly instead of through celda

- **Status:** Proposed
- **Date:** 2026-10-06
- **Deciders:** TODO (maintainer approval required)

## Context

`runDecontX()` called `celda::decontX()` and `celda::decontXcounts()`. The
DecontX algorithm now lives in its own Bioconductor package, `decontX`
(1.10.0 in Bioc 3.23, 1.11.1 in devel). Starting with celda 1.23.0 (Bioc
3.24 devel), `celda::decontX()` is a thin wrapper that requires `decontX`,
but celda lists `decontX` only in Suggests.

On the Bioconductor devel builders, `decontX` is therefore not installed
for singleCellTK, and `runDecontX()` stops with "The 'decontX' package is
required for this function". This causes a check ERROR (the
`detectCellOutlier` example) and two test failures (`test-decontX.R`, and
`test-qc.R` through `runCellQC()`). It blocks the Bioc 3.24 release.

celda and decontX are both lab packages (see AGENTS.md, Related packages).

## Decision

We will add `decontX` to Imports and call `decontX::decontX()` and
`decontX::decontXcounts()` directly. `celda` stays in Imports, because
singleCellTK still uses `celda::distinctColors()`.

`runDecontX()` records `packageDescription("decontX")$Version` in its run
metadata, so the record names the package that actually ran.

## Consequences

- The Bioc devel ERROR is fixed, without depending on celda's Suggests.
- One more direct dependency. decontX brings in rstan and scrapper, but
  celda already needs decontX on devel, so the installed set is effectively
  unchanged.
- Results on devel are unchanged by this decision, because celda's wrapper
  forwards to the same `decontX::decontX()`. Separately, decontX 1.11
  initializes its clusters with scrapper (normalization, HVGs, PCA, UMAP)
  rather than scater and scran, so contamination estimates in Bioc 3.24 can
  differ from 3.23 whichever package singleCellTK calls. That deserves a NEWS
  note.
- On Bioc 3.23 (release), decontX 1.10.0 still calls the deprecated
  `scater::logNormCounts()`, so `runDecontX()` emits deprecation warnings
  there until decontX 1.11 is released.
- The metadata field `packageVersion` now holds the decontX version
  (1.x) instead of the celda version.

## Alternatives considered

- **Keep calling celda and ask celda to move decontX to Imports.** Rejected:
  it fixes singleCellTK only once celda changes, adds a needless
  indirection, and the celda wrapper exists only for backward compatibility.
- **Put decontX in Suggests and check for it at run time.** Rejected:
  `runCellQC()` runs decontX by default, so it is core functionality.
