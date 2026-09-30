# 0003. Move scMerge from Imports to Suggests

- **Status:** Proposed
- **Date:** 2026-09-22
- **Deciders:** TODO (maintainer approval required)

## Context

The Bioconductor 3.24 build report for singleCellTK shows an ERROR on the
`taishan` builder (openEuler Linux, aarch64):

```
ERROR: dependency 'scMerge' is not available for package 'singleCellTK'
```

The package therefore fails to install on that platform and its check is
skipped entirely: no QC, clustering, import or Shiny functionality is
verified there.

The cause is three levels down and outside this package:

```
singleCellTK -> scMerge -> proxyC -> Intel TBB
```

scMerge's own install log on that builder ends with

```
unable to load shared object '.../proxyC/libs/proxyC.so':
undefined symbol: _ZTIN3tbb4taskE
```

which is a toolchain problem with TBB on that machine, not a defect in
scMerge, proxyC or singleCellTK. It is not a general aarch64 problem: proxyC
and scMerge install and load correctly on aarch64 macOS. Fixing it properly
belongs to the Bioconductor build administrators or the proxyC maintainer,
on their timeline.

scMerge is used by exactly one function, `runSCMerge()`
(`R/runBatchCorrection.R`), one of several batch-correction methods. Being
listed in `Imports` makes it a hard install-time requirement for the whole
package.

## Decision

Move `scMerge (>= 1.2.0)` from `Imports` to `Suggests`, and have
`runSCMerge()` check for it at run time with `requireNamespace()`, stopping
with an informative install message when it is missing. This matches the
existing pattern used for `harmony` in the same file.

The documentation for `runSCMerge()` states that scMerge must be installed
separately and why. The example was already wrapped in `\dontrun{}`, so it
does not run where scMerge is absent. No test uses scMerge.

## Alternatives considered

- **Report the TBB problem upstream and wait.** Correct, and still worth
  doing, but it leaves singleCellTK uninstallable on that platform for an
  unknown period. Not mutually exclusive with this decision.
- **Vendor or replace the scMerge functionality.** Far more work than the
  feature warrants, and would duplicate a maintained package.
- **Pin an older proxyC.** Not possible: the Bioconductor builders control
  their own environment.

## Consequences

- singleCellTK installs on platforms where scMerge cannot be installed;
  everything except `runSCMerge()` works there.
- Users of `runSCMerge()` must install scMerge themselves. They get a clear
  message rather than a failed package install.
- The Shiny app still offers scMerge as a batch-correction method. Selecting
  it without scMerge installed now surfaces that error message; the app does
  not check availability before listing the option.
- `R CMD check` after the change: still 1 WARNING (the unrelated scuttle and
  scran deprecation warnings), with no new notes about conditional use of a
  suggested package.
- If the upstream TBB problem is fixed, this decision does not need to be
  reversed: Suggests remains appropriate for a dependency used by one
  optional method.
