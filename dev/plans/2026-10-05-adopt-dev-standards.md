# Finish moving singleCellTK to the shared r-bioc-dev-standards

## Context

singleCellTK began its AI-agent setup before the lab's shared repo existed
(https://github.com/campbio/r-bioc-dev-standards, tag `v1`). As a result, it
holds its own full copies of rules and tooling that now live centrally:

- `AGENTS.md` repeats a "Campbell Lab Playbook v2.0" that the central
  `standards.md` replaces.
- `Makefile`, `dev/hooks/lint-changed.sh`, `.claude/settings.json` and the PR
  template are local copies that will drift from the shared versions.
- There is no `load-standards.sh` SessionStart hook, so the shared standards
  are never loaded.
- There are no `sync-stable` or `pr-base-devel` workflows.

The goal is to follow `ADOPTING.md`: shared rules arrive through the hook and
the `v1` tag, and the package keeps only what's specific to it.

While checking the repo I also found it isn't in sync with Bioconductor:

- Local `devel` is at 2.21.1, but `bioc/devel` is at 2.23.0. The two
  RELEASE_3_23 version-bump commits, `1867491f` and `3ac722fc`, were never
  merged.
- GitHub `master` is still at 2.18.0.
- `compbiomed` has no `RELEASE_3_23` branch.
- Bioconductor 3.24 is due this month. `config.yaml` lists release 3.23 and
  devel 3.24, and local R is 4.6.1.

Decisions made:

- I stage the guardrail files and you copy them into place yourself.
- `.lintr` uses `indent = 2L`.
- The "Overrides" section keeps only "Shiny: propose, don't change".
  - The `bioc_release_2026_09` branch rule is dropped.
  - The "no commits before review" rule is dropped in favour of the standard
    local commits plus a hand-off before anything is pushed.
  - There is no `docs/` override. The standards already forbid hand edits and
    agent site builds, and `site` becomes a people-only target.
- The release-sync steps are included, as steps for you to do.

## Phase 0: Sync with Bioconductor (maintainer)

1. `git fetch bioc && git fetch compbiomed`, then check out `devel` and run
   `git merge bioc/devel`.
2. Resolve the conflict in the DESCRIPTION `Version` field to `2.23.0`. Bioc
   owns x.y.
3. Retitle the NEWS.md 2.21.1 section so it falls under 2.23.x. Note also that
   the 2.21.1 date is wrong: it says `2026-01-1`.
4. Push `devel` to `compbiomed`. Push it to `bioc` only when a z-bumped change
   is ready.

## Phase 1: Agent-editable files

Work on branch `feature/adopt-dev-standards`, created from the freshly fetched
`compbiomed/devel`.

1. **`AGENTS.md`**: rewrite on `templates/AGENTS.md`.
   - Delete the whole "Campbell Lab Playbook" section.
   - Keep the template's header paragraph pointing non-Claude agents
     (`GEMINI.md`) to the `v1` `standards.md` URL and its cache copy.
   - Fill the template sections from the current content:
     - **About**: the current "Project overview", condensed.
     - **Layout**: the current "Repository map".
     - **Object model**: the current text (`inSCE`, `useAssay`, results stored
       on the SCE, accessors, camelCase verbs).
     - **Tests**: fixtures in `data/` and `inst/extdata/`; the full suite is
       slow; Python tests skip without reticulate; no shinytest2 yet; the old
       `shinytest` file is deprecated.
     - **Extra make targets**: `build` (safe, tarball goes outside the repo),
       `app` (safe, launches the app from the working tree), `test-app` (safe),
       `site` (people only, rebuilds the committed `docs/`).
     - **Setup in a new worktree**: the Bioc/R pairing comes from
       `config.yaml`; install with BiocManager; Python via reticulate; macOS
       needs `fftw`. `renv.lock` is not used for development (confirm).
     - **Related packages**: `celda` (campbio; SCTK wraps decontX and celda),
       described the same way in celda's AGENTS.md.
     - **Overrides**: Shiny is propose-only until the maintainer decides.
   - Keep the package notes: no `src/`; CI required jobs, coverage threshold
     and the BiocCheck container are TODO; the lint backlog, with the count
     re-measured after `indent = 2L`.
2. **`CLAUDE.md` and `GEMINI.md`**: no change (`@AGENTS.md`).
3. **`.github/CONTRIBUTING.md`**:
   - Remove `bioc_release_2026_09`. Use `fix/<topic>` or `feature/<topic>`
     branches from `devel`, with PRs to `devel`.
   - Link the central README for the six-step process.
   - Update the layout table: shared standards loaded by
     `dev/hooks/load-standards.sh`, `dev/plans/`, and the two new workflows.
   - Update the command table: `test-one`, `check-full`, `coverage`,
     `site-check` and `article`, plus the extras.
   - Change the style text from "4-space" to "2-space, matching the existing
     code".
4. **`.github/PULL_REQUEST_TEMPLATE.md`**: replace the content with
   `templates/github/pull_request_template.md`. Keep the existing file name
   to avoid a case-only rename on macOS.
5. **`.github/workflows/`**: add `sync-stable.yaml` with
   `stable-branch: master`, and add `pr-base-devel.yaml`, both from the
   templates. Leave `R-CMD-check.yaml` and `BioC-check.yaml` alone; CI
   modernization is a separate ROADMAP item.
6. **`dev/adr/`**:
   - Add `0004-adopt-shared-dev-standards.md`. Changing build, test, release
     and agent tooling needs an ADR. Base it on the existing `template.md`.
   - Update the `README.md` index.
   - Align the `README.md` process with the central one: issue → Proposed →
     approved → implement.
   - Drop the reference to the nonexistent `adr-author` skill.
7. **`dev/plans/`**: create it and save this plan there as
   `2026-10-05-adopt-dev-standards.md`. Plans are committed.
8. **`dev/RELEASE.md`**:
   - Change `git fetch upstream` to `bioc`.
   - Use `make check-full`, `make bioccheck` (which prints tarball size),
     `make site-check`, and `make article FILTER=` in place of hand-knitting.
   - Add the standard release items: deprecations advanced, `/security-review`
     of the release diff, tarball under 10 MB.
   - Link the README "Release day" steps.
   - Fill the target as Bioc 3.24 and leave the dates and owners TODO.
9. **`dev/AUDIT.md`**:
   - Use `make bioccheck` in place of raw BiocCheck.
   - Replace the nonexistent skill names (`bioc-pkg-dev`,
     `security-audit-r-package`) with `build-check-bioccheck` and
     `r-lib:lifecycle`.
   - Write the audit output outside the repo or to `dev/`, as now.
10. **`dev/ROADMAP.md`**: add these items:
    - Decide on Shiny agent scope.
    - Decide on pkgdown deploy (CI vs committed `docs/`).
    - Decide on 2→4-space conversion (its own PR).
11. **`.Rbuildignore`**: add `^\.worktrees$`. The other required entries are
    already there.
12. **`.gitignore`**: add `.worktrees/` and `.claude/settings.local.json`.
13. **`_pkgdown.yml`**: add `url: https://www.camplab.net/sctk/`, which
    `site-check` needs and which matches DESCRIPTION `URL`. Confirm that this
    is the pkgdown site's URL.
14. **NEWS.md**: no entry. This change isn't user-facing and `dev/` and
    `.github/` aren't in the build. No version bump.

## Phase 2: Guardrail files (I stage them, you copy them)

I write each file to the session scratchpad and give you one `cp` block. You
review and copy. I then `git add` and commit them. Committing is allowed,
since the deny rules cover only Edit and Write.

1. **`dev/hooks/load-standards.sh`**: new, a verbatim copy of central
   `hooks/load-standards.sh`.
2. **`dev/hooks/lint-changed.sh`**: replace it with the central stub. The
   shared hook lints any edited `.R` file, so `inst/shiny` is still covered.
3. **`Makefile`**: start from `templates/Makefile` and set:
   - `FORCE_SUGGESTS = FALSE`, which mirrors `.BBSoptions`.
   - `PEOPLE_ONLY := site`.
   - Below the include, move over only `build`, `site`, `app` and `test-app`,
     with the `BUILD_DIR` helper.
   - Drop the local `help`, `test`, `check`, `bioccheck`, `docs` and `lint`
     targets.
4. **`.claude/settings.json`**: `templates/settings.json` merged with these
   package additions:
   - allow: `make help`, `make build`, `make app`, `make test-app`,
     `git fetch compbiomed`, `git blame:*`.
   - deny: `make site`, `git push --force-with-lease:*`, `rm -fr:*`.
   - Remove the overly broad `git branch:*` allow. The template's safe forms
     replace it.
   - Hooks: SessionStart (`load-standards.sh`) plus PostToolUse
     (`lint-changed.sh`), using the template's relative paths.
5. **`.lintr`**: change `indentation_linter(indent = 2L)` and add `"man"` to
   the exclusions. Keep the other linter settings as they are.

## Phase 3: Verification

- `make help` lists the standard targets and the 4 extras. The first run
  downloads `standards.mk` into `~/.cache/r-bioc-dev-standards/v1/`.
- `CLAUDECODE=1 make site` refuses to run (people-only guard).
  `make test-one FILTER=qc` runs; `make test-one FILTER='a;b'` is rejected.
- `make lint`: record the new total, expected to be well under 33,700, and
  update the AGENTS.md debt note.
- `make site-check` passes, as does `BiocCheckGitClone`, via `make bioccheck`.
- `make check` passes, or any failures already exist on `devel`. Compare the
  results against `devel` before treating anything as new.
- Start a new `claude` session from the repo root, approve the changed hooks,
  and ask: "Which remote do the standards say never to push to, and when may a
  PR be opened?" Claude should answer from the standards without reading
  files.
- `grep -rn "Playbook\|bioc_release_2026_09"` returns nothing outside `docs/`.

## Phase 4: Review and hand-off

- Run requesting-code-review against this plan, then `/code-review`.
- Show `git log devel..HEAD` and the full diff. Push and open a PR to `devel`
  only after you approve.

## Phase 5: Stable branch (maintainer, after the PR merges)

Follow `ADOPTING.md` step 7 with `<shared>` = `compbiomed`:

1. Run `git push compbiomed bioc/RELEASE_3_23:refs/heads/RELEASE_3_23`.
2. Cherry-pick the workflow and PR-template commit onto `RELEASE_3_23` and
   push it. That runs the first sync, so `master` → 2.22.0 and gets the
   `v2.22.0` tag.
3. `master` is already the default branch. Add a ruleset on it with
   **Restrict deletions** and **Block force pushes**.
4. Confirm the "Sync stable branch" Action is green.
5. On Bioc 3.24 release day, follow the README "Release day" steps.
