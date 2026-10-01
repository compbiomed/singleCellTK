# 0001. Record architecture decisions

- **Status:** Proposed
- **Date:** 2026-09-18
- **Deciders:** TODO (maintainer approval required)

## Context

singleCellTK is maintained by a changing group of lab members, increasingly
helped by AI coding agents. The reasons behind past choices (dependencies,
interfaces, workarounds for Bioconductor checks) live in people's heads or in
scattered commit messages, and are lost when people leave. Agents also need
a written record so they don't propose reversing decisions already made.

## Decision

- Significant decisions are recorded as numbered ADRs in `dev/adr/`, using
  `dev/adr/template.md`. They are never stored in `docs/`, which is pkgdown
  build output.
- An ADR is required for decisions that are hard to reverse or that future
  contributors will question (see `dev/adr/README.md`).
- ADRs are proposed in the PR that makes the change, and accepted by the
  maintainer when that PR is approved.
- Numbering is sequential and zero padded. ADRs are append-only: decisions
  are never edited, and a reversal is a new ADR that marks the old one
  `Superseded by NNNN`.

## Consequences

- Contributors and agents must add an ADR in the same PR as any qualifying
  change. `AGENTS.md` and the PR template point to this.
- `dev/` is excluded from the package build via `.Rbuildignore`.
- A small, ongoing cost per structural PR, in exchange for a durable record.
