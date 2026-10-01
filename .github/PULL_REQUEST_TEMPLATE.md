## Summary

<!-- What does this PR change, and why? Link related issues. -->

## Checklist

- [ ] `make test` passes
- [ ] `make check` passes (no ERRORs or WARNINGs)
- [ ] `make lint`: no new lints in changed lines
- [ ] `NEWS.md` updated (user-facing changes)
- [ ] roxygen edited and `make docs` run; no hand edits to `man/`, `NAMESPACE`, or `docs/`
- [ ] New exports added to `_pkgdown.yml` (`pkgdown::check_pkgdown()` passes)
- [ ] `/code-review` run on the diff
- [ ] ADR added in `dev/adr/` and linked below, if this is a structural or dependency decision
- [ ] UI changes: screenshot of the running app attached (`make app`)

ADR: <!-- dev/adr/NNNN-...md or "n/a" -->

## Scientific correctness (human judgment)

<!-- Required for changes that affect analysis results. Who checked that the
results are scientifically correct, and how (for example, compared outputs on
a reference dataset)? Tests passing is not enough. -->

- [ ] A person has checked that analysis results are scientifically correct, or this PR does not affect results
