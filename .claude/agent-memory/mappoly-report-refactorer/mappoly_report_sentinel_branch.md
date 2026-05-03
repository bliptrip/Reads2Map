---
name: mappoly_report.R skip-branch sentinel writes 1 to output RDS
description: Behavioral quirk in scripts/mappoly_report.R that any refactor must preserve
type: project
---

When `length(seq.init$seq.mrk.names) <= sample_size`, the bootstrap module
overwrites `summaries`, `info`, `dat`, `mat` to the integer `1` (sentinel)
instead of running the bootstrap. The downstream RDS export then writes those
sentinel `1` values to `_summaries.rds`, `_info.rds`, `_dat.rds`, `_mat2.rds`.

The export block then unconditionally recomputes `tpt` and `mat` from
`seq.init` and overwrites the sentinel `mat <- 1`. So the actual on-disk
contents in the skip branch are: summaries=1, info=1, dat=1,
mat=<real rf matrix>, maps=<empty boot.refactored>.

**Why:** Downstream tasks in `tasks/JointReports.wdl` and the per-population
aggregator distinguish "bootstrap ran" vs. "skipped" by checking for the
sentinel `1` in the RDS bundle.

**How to apply:** Any refactor of this script MUST keep both the sentinel
overwrite AND the unconditional `mat` recomputation. The refactored
`mappoly_bootstrap()` returns `mat = NULL` when bootstrap runs (because
`export_mappoly_results()` computes `mat` itself either way) and `mat = 1`
in the skip branch — but `export_mappoly_results()` then overwrites with the
real matrix before saveRDS, matching the original byte-for-byte.
