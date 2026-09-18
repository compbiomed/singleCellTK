# singleCellTK Dependency and Deprecation Audit

A periodic audit (suggested: once per Bioconductor cycle, before the release
checklist in `dev/RELEASE.md`) to catch deprecated functions and dependency
problems before they break the build.

**This audit does not change code.** Findings become GitHub issues, and ADRs
where the fix is structural.

## Prompt

Give an agent this prompt from the repo root:

> Run `BiocCheck::BiocCheck()` and `BiocManager::valid()` on this package.
> Use r-lib lifecycle practices to find deprecated functions or S4 methods
> and report replacements compatible with the current Bioconductor release.
> Write findings to `dev/agent-log.md`. Do not change code — findings become
> GitHub issues (and ADRs where structural).

Useful skills: `build-check-bioccheck`, `bioc-pkg-dev`, `r-lib:lifecycle`,
`security-audit-r-package`.

## Log

| Date | Bioconductor version | Run by | Findings (issues/ADRs) |
|---|---|---|---|
| TODO | | | |
