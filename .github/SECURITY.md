# Security Policy

## Reporting a vulnerability

Please report security issues **privately** to the package maintainer:

**Joshua David Campbell** — <camp@bu.edu>

Do not open a public GitHub issue for security problems. Include a
description, steps to reproduce, and the affected version. You should receive
a response within TODO business days.

## Scope

singleCellTK is an R/Bioconductor package with a Shiny app and a command-line
QC pipeline. Relevant issues include code execution through crafted input
files, unsafe handling of downloaded data, and exposure of local files
through the Shiny app.

## Rules for automated tools and AI agents

Agents working in this repository must:
- **Never read, print, or commit credentials**: tokens, API keys, `.Renviron`,
  `.Rprofile` secrets, SSH keys, or cloud credentials.
- **Never exfiltrate data**: do not upload, paste, or send data files,
  including anything under `data/`, `inst/extdata/`, or user data loaded in a
  session, to external services.
- **Never commit absolute local paths** or machine-specific configuration.
- Report suspected vulnerabilities to the maintainer instead of publishing
  fixes that disclose them.

See also `AGENTS.md` (safety rules) and `.claude/settings.json` (enforced
permissions for Claude Code).
