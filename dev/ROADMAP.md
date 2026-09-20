# singleCellTK Roadmap

> **Status: SKELETON.** Sections marked `TODO` need input from the maintainer
> (Joshua D. Campbell) and the lab. Agents: treat `TODO` sections as *unknown
> direction*. Ask before proposing work there, and don't infer priorities from
> them.

This file tells contributors (human and agent) which direction proposals
should lean. It is not a changelog (see `NEWS.md`) or a decision record (see
`dev/adr/`).

## Near term (current Bioconductor cycle)

Known items from the AI-agent setup (2026-09-18):

- [ ] Clear `R CMD check` and BiocCheck ERRORs/WARNINGs for the current
      Bioconductor release (triage first, then fix in focused PRs).
- [ ] Establish a lintr baseline for `R/` and `inst/shiny/`. New code must be
      clean, and existing debt is reduced deliberately, not in bulk.
- [ ] Add test coverage for the Shiny app (`shiny::testServer()` for reactive
      logic, and a small `shinytest2` smoke suite). The old `shinytest` test is
      deprecated.
- [ ] Modernize CI workflows (BiocCheck in the Bioconductor container, and
      `check_pkgdown()`).
- [ ] TODO: other near-term goals.

## Medium term (next 1–2 releases)

- [ ] TODO: new analysis methods or integrations planned.
- [ ] TODO: dependency reduction (the package has ~80 Imports). Which can move
      to Suggests?
- [ ] TODO: Shiny app direction (modules / restructuring). Refactor territory,
      so an ADR is required first.

## Long term

- [ ] TODO: long-term vision for SCTK.

## Out of scope / not planned

Proposals in these areas should not be pursued without maintainer approval.

- TODO: list directions the lab has decided against.

## Open questions for the maintainer

- TODO: Is the Shiny app currently in scope for agent-driven changes?
- TODO: Which deprecated methods or wrappers should be retired?
