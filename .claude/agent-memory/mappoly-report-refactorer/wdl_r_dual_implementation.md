---
name: MAPpoly task has dual implementation (WDL heredoc + standalone script)
description: tasks/mappoly.wdl inlines R via heredoc; scripts/mappoly_report.R is a parallel standalone version
type: project
---

`tasks/mappoly.wdl` task `MappolyReport` inlines its R code as a `R --vanilla
--no-save <<RSCRIPT ... RSCRIPT` heredoc with WDL-interpolated `~{var}`
placeholders. `scripts/mappoly_report.R` is a parallel standalone version of
the same logic, designed to be run via `Rscript scripts/mappoly_report.R
'name=value' ...` inside the `cristaniguti/reads2map:0.0.8` container.

**Why:** The WDL heredoc form is what currently runs in production. The
standalone script appears to be a portable / locally-debuggable mirror — the
header comment-block calls it "Standalone implementation of the MappolyReport
WDL task". As of 2026-05-01, no WDL task references `scripts/mappoly_report.R`
directly (grep confirms zero hits in tasks/ and subworkflows/).

**How to apply:** Edits to one should be mirrored to the other to keep them
in sync. The refactor of `scripts/mappoly_report.R` into four functions does
NOT auto-propagate to `tasks/mappoly.wdl`. If/when that mirroring matters,
either inline the new functions into the WDL heredoc, or change the WDL task
command to `Rscript scripts/mappoly_report.R 'name=value' ...` and ship the
script via `localization_optional` or a Docker layer.

The WDL heredoc uses `[[]]` accessors (not `$`) because the unquoted bash
heredoc would otherwise expand `$var` as shell vars — a constraint the
standalone script does not have but adopts anyway for consistency.
