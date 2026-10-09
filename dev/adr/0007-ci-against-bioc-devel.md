# 0007. Run GitHub CI against Bioconductor devel

- **Status:** Proposed
- **Date:** 2026-10-08
- **Deciders:** TODO (maintainer approval required)

## Context

The `devel` branch is the Bioconductor devel version of singleCellTK, but
GitHub CI installed the Bioconductor *release* packages that match the
runner's R version (Bioc 3.23 for R 4.6). CI therefore tested devel code
against older dependencies than the Bioconductor builders use. It failed
on warnings the builders don't see, such as decontX 1.10's deprecated
scuttle calls, which decontX 1.11 in devel no longer makes, and it could
miss problems they do see.

The macOS jobs never reached the check: `alabaster.base` (needed by
celldex and scRNAseq) has no macOS binary for the current release, and its
source build fails to link OpenSSL. The BiocCheck job ran on macOS, so it
failed the same way.

## Decision

We will:

- Set `R_BIOC_VERSION` to the current Bioconductor devel version ("3.24")
  in `R-CMD-check.yaml` and `BioC-check.yaml`. pak, used by
  `setup-r-dependencies`, then installs Bioconductor devel packages.
- Before installing dependencies on macOS, install Homebrew's
  `openssl@3` and add its `lib` and `include` directories to `LDFLAGS` and
  `CPPFLAGS` in `~/.R/Makevars`.
- Run BiocCheck on `ubuntu-latest` and install it with the other
  dependencies (`bioc::BiocCheck`), so it also comes from devel.
- Run BiocCheck on the built tarball rather than on the Git checkout.
  BiocCheck 1.49 (Bioc 3.24) errors on Git-tracked `.claude/` files, and the
  lab standards commit `.claude/settings.json`. The tarball excludes it
  through `.Rbuildignore`, and the Bioconductor builders also check
  tarballs. Whether the standards should stop tracking it is raised in
  campbio/r-bioc-dev-standards.
- Wrap singleCellTK's calls into batchelor (`fastMNN`, `reducedMNN`,
  `mnnCorrect`) in `.muffleUpstreamDeprecations()`, which muffles only
  base R deprecation warnings raised inside batchelor.

## Consequences

- CI tests what the Bioconductor devel builders build, so a green CI
  should predict a clean build report.
- `R_BIOC_VERSION` must be updated at each Bioconductor release (April and
  October). `dev/RELEASE.md` should gain that step. If it is forgotten, CI
  silently tests against an old devel.
- The first CI runs build many Bioconductor packages from source on Linux
  until the cache is warm.
- The macOS OpenSSL step and the batchelor muffling are workarounds. Remove
  them once binaries exist or batchelor updates.

## Alternatives considered

- **Keep testing against release.** Rejected: it tests the wrong
  dependencies and fails on warnings that are already fixed in devel.
- **Use the `bioconductor/bioconductor_docker:devel` container.** This
  matches the builders' system libraries more closely, but only runs on
  Linux and would replace the existing r-lib/actions setup. It is worth
  revisiting when CI is modernized (ROADMAP).
- **Drop macOS from CI.** Rejected: one of the Bioconductor builders is
  macOS, and the OpenSSL fix is small.
