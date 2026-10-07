# 0004. Adopt the shared r-bioc-dev-standards

- **Status:** Proposed
- **Date:** 2026-10-05
- **Deciders:** TODO (maintainer approval required)

## Context

singleCellTK set up its AI-agent tooling (September 2026) before the lab
had a shared home for it. The package therefore carried its own full copy
of the lab's rules (the "Campbell Lab Playbook" section of `AGENTS.md`), its
own Makefile recipes, lint hook, Claude Code permissions, and PR template.
The other lab packages carry similar copies, and each copy drifts as it is
fixed in one place but not the others.

The lab now keeps these in one repository,
<https://github.com/campbio/r-bioc-dev-standards>, published through the
`v1` tag. It provides condensed standards loaded at every Claude Code
session start, the standard `make` targets, the lint hook, and shared
GitHub workflows that keep the stable branch matching the current
Bioconductor release.

Separately, the GitHub `master` branch had fallen behind the Bioconductor
release (2.18.0 while release 3.23 ships 2.22.0), because nothing updated
it automatically.

## Decision

We will adopt the shared standards as described in their `ADOPTING.md`:

- `dev/hooks/load-standards.sh` (SessionStart) downloads `standards.md`,
  the lint hook, and `standards.mk` from tag `v1`, and loads the standards
  into each session. `dev/hooks/lint-changed.sh` becomes a stub that runs
  the downloaded lint hook.
- The Makefile includes the shared `standards.mk` for the standard targets
  (`test`, `test-one`, `check`, `check-full`, `bioccheck`, `docs`, `lint`,
  `coverage`, `site-check`, `article`). It keeps only singleCellTK's
  settings (`FORCE_SUGGESTS = FALSE`, mirroring `.BBSoptions`) and extra
  targets (`build`, `app`, `test-app`, and the people-only `site`).
- `.claude/settings.json` follows the shared template, so `git push` and
  `gh pr create` always ask, and guardrail files can't be edited by agents.
- `AGENTS.md` keeps only package-specific notes and overrides; the only
  override is that Shiny app changes are proposed, not made, until the
  maintainer decides the app's scope.
- `.lintr` indentation is set to 2 spaces to match the existing code.
  Converting the code to Bioconductor's recommended 4 spaces is a separate
  decision.
- `sync-stable.yaml` and `pr-base-devel.yaml` call the shared workflows,
  with `master` as the stable branch.
- The dedicated `bioc_release_2026_09` working branch is retired; work uses
  `fix/<topic>` and `feature/<topic>` branches from `devel`.

## Consequences

- Fixes to the rules, make targets, lint hook, and workflows reach
  singleCellTK when the lab moves the `v1` tag, without editing this repo.
- Sessions need network access to GitHub to refresh the shared files; if it
  is unreachable, the last cached copies are used and Claude says so.
- `make check` now skips vignettes; the full check that rebuilds them is
  `make check-full`, which is required before a PR.
- `master` becomes CI-maintained. Nobody commits to it by hand, and a
  ruleset blocks force-pushes and deletion.
- The lint backlog reported by `make lint` shrinks, since indentation no
  longer conflicts with the code.
- Follow-up: maintainer pushes `RELEASE_3_23` to GitHub and protects
  `master` (ADOPTING.md step 7); decide Shiny agent scope and how the
  pkgdown site is deployed (`dev/ROADMAP.md`).

## Alternatives considered

- **Keep the local copies.** Rejected: every fix would have to be repeated
  in each lab package, and the copies had already started to differ.
- **Git submodule or vendored copy of the standards.** Rejected: it still
  needs a commit in every package to pick up a change, and submodules
  complicate cloning and Bioconductor's git server.
