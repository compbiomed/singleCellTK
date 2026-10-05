#!/usr/bin/env bash
# PostToolUse hook stub for Claude Code. It runs the shared lint hook, which
# dev/hooks/load-standards.sh downloads at each session start. Does nothing
# if no copy has been downloaded yet.
#
# Copy this file to dev/hooks/lint-changed.sh in each package and register
# it in .claude/settings.json (see ADOPTING.md).
# Source: https://github.com/campbio/r-bioc-dev-standards

hook="${XDG_CACHE_HOME:-$HOME/.cache}/r-bioc-dev-standards/${R_BIOC_STANDARDS_REF:-v1}/lint-changed.sh"
[ -f "$hook" ] && exec bash "$hook"
exit 0
