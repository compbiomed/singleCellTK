# singleCellTK Release Checklist

Bioconductor releases twice a year, around **April** and **October**. This
checklist runs once per cycle. Sections marked `TODO` need a maintainer
decision.

Authoritative sources:
- Release schedule and freeze dates: <https://bioconductor.org/developers/release-schedule/>
- Current release/devel versions and matching R: <https://bioconductor.org/config.yaml>
- Build reports: <https://bioconductor.org/checkResults/>

## 0. This cycle

| Item | Value |
|---|---|
| Target Bioconductor release | TODO |
| Package freeze date | TODO (from the release schedule) |
| Release date | TODO |
| Release owner | TODO |
| Site deploy owner | TODO |

## 1. Prepare (about 6 weeks before freeze)

- [ ] Confirm the local R and Bioconductor versions match the target in
      `config.yaml` (never assume the pairing).
- [ ] Sync with Bioconductor: `git fetch upstream` from `git.bioconductor.org`
      (only `devel` and `RELEASE_x_y` branches accept pushes there).
- [ ] Merge any Bioconductor-side changes into the working branch.

## 2. Check and fix

- [ ] `make check`: `R CMD check` with no ERROR or WARNING.
- [ ] `make bioccheck`: BiocCheck on the tarball plus `BiocCheckGitClone()`.
- [ ] Triage every finding into a fix plan (use the `build-check-bioccheck`
      skill): real defects vs. environment gaps vs. known false positives.
- [ ] Fix in focused PRs, one category per PR.
- [ ] Re-run `make check` and `make bioccheck` until clean.
- [ ] `make test` and `make lint` pass.

## 3. Documentation

- [ ] `NEWS.md` updated for this version (use the `update-r-news` skill).
- [ ] `pkgdown::check_pkgdown()` passes (every export is in `_pkgdown.yml`).
- [ ] Articles under `vignettes/articles/` that changed this cycle were knit
      locally (`R CMD check` does not run them).
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
- [ ] Fix any platform-specific failures on the `RELEASE_x_y` branch and
      merge the fixes to `devel`.
- [ ] Update `dev/ROADMAP.md` for the next cycle.
- [ ] Fill in section 0 for the next cycle.

## TODO (maintainer decisions)

- TODO: freeze and release dates for each cycle.
- TODO: who owns the release, and who rebuilds and commits the site.
- TODO: whether the twice-yearly check → fix → PR loop is automated (scheduled
  agent) or run by hand.
