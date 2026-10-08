# singleCellTK Release Checklist

Bioconductor releases twice a year, around **April** and **October**. This
checklist runs once per cycle. Sections marked `TODO` need a maintainer
decision. It adds singleCellTK's specifics to the shared standards'
[Checks and releases](https://github.com/campbio/r-bioc-dev-standards#checks-and-releases)
and
[Release day](https://github.com/campbio/r-bioc-dev-standards#branches-releases-and-tags)
steps.

Authoritative sources:
- Release schedule and freeze dates: <https://bioconductor.org/developers/release-schedule/>
- Current release/devel versions and matching R: <https://bioconductor.org/config.yaml>
- Build reports: <https://bioconductor.org/checkResults/>

## 0. This cycle

| Item | Value |
|---|---|
| Target Bioconductor release | 3.24 (October 2026) |
| Package freeze date | TODO (from the release schedule) |
| Release date | TODO |
| Release owner | TODO |
| Site deploy owner | TODO |

## 1. Prepare (about 6 weeks before freeze)

- [ ] Confirm the local R and Bioconductor versions match the target in
      `config.yaml` (never assume the pairing).
- [ ] Sync with Bioconductor: `git fetch bioc` (the remote whose URL is
      `git.bioconductor.org`), merge `bioc/devel` into `devel`, and push
      `devel` to `compbiomed`. Bioconductor's own commits (for example the
      version bumps at each release) must be merged before any push to
      `bioc`.

## 2. Check and fix

- [ ] `make check-full`: `R CMD check` with vignettes and `\donttest`
      examples, no ERROR or WARNING.
- [ ] `make bioccheck`: BiocCheck on the tarball plus `BiocCheckGitClone()`.
      Tarball under 10 MB (printed), no single file over 5 MB.
- [ ] Triage every finding into a fix plan (use the `build-check-bioccheck`
      skill): real defects vs. environment gaps vs. known false positives.
- [ ] Fix in focused PRs, one category per PR.
- [ ] Re-run `make check-full` and `make bioccheck` until clean.
- [ ] `make test` passes; `make coverage` hasn't dropped since the last
      release; `make lint` shows no new lints.
- [ ] Deprecations advanced one stage (`.Deprecated()` → `.Defunct()` →
      removed).
- [ ] `/security-review` run on the release diff.

## 3. Documentation

- [ ] `NEWS.md` updated for this version (use the `update-r-news` skill).
- [ ] `make site-check` passes (every export is in `_pkgdown.yml`).
- [ ] Articles under `vignettes/articles/` that changed this cycle render
      with `make article FILTER=<name>` (`R CMD check` does not run them).
- [ ] Site rebuilt with `make site` and the updated `docs/` committed by the
      site deploy owner. The site is served from `docs/` on the main branch.
      TODO: confirm the owner and when this happens in the cycle.

## 4. Version bump

- [ ] Follow Bioconductor's `x.y.z` scheme: `y` is odd in devel and even in
      release, and `z` is bumped on every commit pushed to Bioconductor. At
      release, Bioconductor itself bumps devel to the next odd `y`.
- [ ] Record any structural decisions made this cycle in `dev/adr/`.

## 5. After release

- [ ] Check the Bioconductor build report for singleCellTK on all platforms.
- [ ] Update `R_BIOC_VERSION` in `.github/workflows/R-CMD-check.yaml` and
      `BioC-check.yaml` to the new devel version (ADR 0007).
- [ ] Fix any platform-specific failures on `devel` first, then port them
      to `RELEASE_x_y` on their own branch by cherry-pick with a release z
      bump. Never merge `devel` into a release branch.
- [ ] Confirm the "Sync stable branch" Action updated `master` and tagged
      the new release version.
- [ ] Update `dev/ROADMAP.md` for the next cycle.
- [ ] Fill in section 0 for the next cycle.

## TODO (maintainer decisions)

- TODO: freeze and release dates for each cycle.
- TODO: who owns the release, and who rebuilds and commits the site.
- TODO: whether the twice-yearly check → fix → PR loop is automated (scheduled
  agent) or run by hand.
