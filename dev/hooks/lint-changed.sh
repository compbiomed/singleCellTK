#!/usr/bin/env bash
# PostToolUse hook (see .claude/settings.json): report-only lint of the file
# an agent just edited. It never modifies files, only reports lints back to
# the agent, and always exits 0 so it never blocks work.
#
# Scope: .R files under R/ and inst/shiny/ (inst/ is invisible to R CMD check).

set -u
input="$(cat)"
file="$(printf '%s' "$input" | python3 -c 'import sys, json
try:
    d = json.load(sys.stdin)
    print(d.get("tool_input", {}).get("file_path", ""))
except Exception:
    print("")')"

[ -n "$file" ] || exit 0
root="${CLAUDE_PROJECT_DIR:-$(pwd)}"
rel="${file#"$root"/}"

case "$rel" in
  R/*.R|R/*.r|inst/shiny/*.R|inst/shiny/*.r) ;;
  *) exit 0 ;;
esac
[ -f "$root/$rel" ] || exit 0

report="$(cd "$root" && Rscript -e 'l <- lintr::lint(commandArgs(TRUE)[1]); if (length(l)) print(l)' "$rel" 2>/dev/null | head -n 60)"
[ -n "$report" ] || exit 0

# Hand the lint report back to the agent as context (report-only).
printf '%s' "$report" | python3 -c 'import sys, json
msg = "lintr (report-only) for the file you just edited. Fix lints in lines you changed; do not restyle untouched code:\n" + sys.stdin.read()
print(json.dumps({"hookSpecificOutput": {"hookEventName": "PostToolUse", "additionalContext": msg}}))'
exit 0
