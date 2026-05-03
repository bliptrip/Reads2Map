---
name: MAPpoly 0.4.0 API quirks relevant to scripts/mappoly_report.R
description: Non-obvious MAPpoly behaviors that affect refactoring of the report script
type: reference
---

MAPpoly 0.4.0 quirks observed while refactoring `scripts/mappoly_report.R`:

- `read_vcf(filter.non.conforming = TRUE)` is the only path that populates
  `dat$chisq.pval`. The script intentionally uses `filter.non.conforming = FALSE`
  to preserve `dat$geno` (needed by `est_full_hmm_with_prior_prob`), then
  re-derives `chisq.pval` manually using `mappoly:::mrk_chisq_test` (note the
  `:::` — `mrk_chisq_test` is internal, not exported).

- `mappoly:::mrk_chisq_test` will crash on offspring dosage calls that violate
  Mendelian expectations given the parents (e.g. dosage=2 progeny from a 0x1
  cross). The fix in this script is to scrub bad dosages to `pl+1` (missing)
  for the chi-square computation only, leaving `dat$geno.dose`/`dat$geno`
  intact — pattern lives in `mappoly_preprocess()`.

- `est_pairwise_rf` and PSOCK clusters: each PSOCK worker gets a serialized
  copy of `seq.init`, `dat`, and the helper closures. `clusterExport()` reads
  from a specified `envir`, so when refactoring this code into functions,
  bind helpers/data to the local function env first then pass
  `envir = environment()`.

- `summary_maps()` returns a data.frame; the script tags it post-hoc with
  `parent`, `prob_source`, `global_error`, `data` columns.

**How to apply:** When extending this report or building a new MAPpoly task,
prefer the manual chisq.pval pattern over `filter.non.conforming = TRUE` if
prior-probability mapping is on the table; and always pass closures via a
`helpers = list(...)` arg for testability.
