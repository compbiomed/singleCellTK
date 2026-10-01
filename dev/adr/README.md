# Architecture Decision Records

This folder records **why** significant decisions about singleCellTK were
made, so future contributors (people and agents) don't have to guess or
re-argue them. `NEWS.md` records *what* changed for users; git history
records *which lines* changed.

## When an ADR is required

Write one for decisions that are **hard to reverse** or that **future
contributors will question**:
- adding, removing, or replacing a dependency (DESCRIPTION changes)
- module or file structure (for example, splitting files, restructuring
  the Shiny app)
- class design (new S4 classes, changing where results are stored)
- build, CI, release, or deploy machinery
- deliberately ignoring a BiocCheck or style recommendation

**Not** for routine choices: bug fixes, version bumps, docs typos, or anything
self-evident from the diff.

## How

1. Copy `template.md` to `NNNN-short-title.md` (next unused number, zero
   padded). The `adr-author` skill can draft it.
2. Open it as `Proposed` in the same PR as the change, and add it to the index
   below.
3. The maintainer accepts it (status `Accepted`) when the PR is approved.
4. **Decisions are never edited.** To reverse or replace one, write a new ADR
   and mark the old one `Superseded by NNNN`.

## Index

| # | Title | Status | Date |
|---|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Proposed | 2026-09-18 |
| [0002](0002-harmony-version-dispatch.md) | Support both harmony interfaces in runHarmony() by version dispatch | Proposed | 2026-09-20 |
| [0003](0003-scmerge-to-suggests.md) | Move scMerge from Imports to Suggests | Proposed | 2026-09-22 |
