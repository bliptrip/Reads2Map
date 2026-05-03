#!/usr/bin/env Rscript
#
# mappoly_report.R
#
# Standalone implementation of the MappolyReport WDL task (tasks/mappoly.wdl).
# Builds bootstrap and (optionally) full-marker linkage maps with MAPpoly for a
# single VCF, writes per-step diagnostic data.frames and map objects to .rds,
# and tars them into <SNPCall>_<GenoCall>_<Counts>[<chr_suffix>]_poly_results.tar.gz.
#
# Arguments are passed as `name=value` R expressions, parsed via
# eval(parse(text=...)). Defaults match the WDL task defaults.
#
# Usage:
#   Rscript scripts/mappoly_report.R \
#     'vcf_file="/path/to/in.vcf.gz"' \
#     'SNPCall_program="gatk"' \
#     'GenotypeCall_program="updog"' \
#     'CountsFrom="vcf"' \
#     'parent1="P1"' \
#     'parent2="P2"' \
#     'ploidy=4' \
#     'max_cores=4' \
#     'global_errors="0.05,0.1"' \
#     'chromosome="chr1"' \
#     'build_full_map=TRUE' \
#     'output_dir="./out"' \
#     'verbose="FALSE"'
#
# Required: vcf_file, SNPCall_program, GenotypeCall_program, CountsFrom,
#           parent1, parent2, max_cores, ploidy
#
# Optional (default): prob_thres (0.8), repetitions (100), sample_size (30),
#                     filt_segr ("TRUE"), global_errors ("0.05"),
#                     chromosome (""), build_full_map (TRUE), output_dir (".")
#
# global_errors accepts either a comma-separated string ("0.05,0.1") or a
# character vector (c("0.05","0.1")).

suppressPackageStartupMessages({
  library(mappoly)
  library(parallel)
  library(vcfR)
})


# ── Helper functions ──────────────────────────────────────────────────────────
gen_map_sequential <- function(seq.init, max_cores = 4, verbose = FALSE) {
  tpt <- est_pairwise_rf(seq.init, ncpus = max_cores)
  map.seq <- est_rf_hmm_sequential(input.seq = seq.init,
                                    start.set = 5,
                                    thres.twopt = 10,
                                    thres.hmm = 10,
                                    extend.tail = 30,
                                    info.tail = TRUE,
                                    twopt = tpt,
                                    phase.number.limit = 10,
                                    reestimate.single.ph.configuration = TRUE,
                                    tol = 10e-2,
                                    tol.final = 10e-4,
                                    verbose = verbose,
                                    detailed.verbose = verbose #Set to the same value as verbose
                                  )
  map <- filter_map_at_hmm_thres(map.seq, thres.hmm = 0.0001)
  return(list(map = map, tpt = tpt))
}

gen_map_global <- function(boot_seq, error, verbose = FALSE) {
  # Global error
  map <- boot_seq[["map"]]
  tpt <- boot_seq[["tpt"]]
  map2 <- est_full_hmm_with_global_error(map, error = error, tol = 10e-3, verbose = verbose)
  # split_and_rephase can fail with "subscript out of bounds" when the combined-parent
  # ("all") sequence in a diploid cross produces joint-parent phase-config strings at
  # block boundaries that aren't covered by the transition matrix A (which is built for
  # single-parent marker types). Fall back to map2 directly so the final HMM pass still runs.
  map3 <- tryCatch({
    split_and_rephase(map2, gap.threshold = 5, size.rem.cluster = 3, twopt = tpt)
    est_full_hmm_with_global_error(map3, error = error, tol = 10e-4, verbose = verbose)},
    error = function(e) {
      message("split_and_rephase (global) failed: ", conditionMessage(e), " — skipping rephase")
      map2
    }
  )
  return(map3)
}

gen_map_prob <- function(boot_seq, dat, verbose = FALSE) {
  # Prob error
  map <- boot_seq[["map"]]
  tpt <- boot_seq[["tpt"]]
  map.prob <- NULL
  if (!is.null(dat[["geno"]])) {
    map2 <- est_full_hmm_with_prior_prob(map, tol = 10e-3, verbose = verbose)
    map3 <- tryCatch({
      split_and_rephase(map2, gap.threshold = 5, size.rem.cluster = 3, twopt = tpt)
      est_full_hmm_with_prior_prob(map3, tol = 10e-4, verbose = verbose)},
      error = function(e) {
        message("split_and_rephase (prob) failed: ", conditionMessage(e), " — skipping rephase")
        map2
      }
    )
  }
  return(map3)
}

run_bootstraps <- function(seed, s_all, sample_size, dat, max_cores, verbose) {
  maps.global <- maps_init_global
  maps.prob <- maps_init

  anchor_err <- global_errors[[1]]
  for (p in map_parentals) {
    set.seed(seed[[1]])
    # Anchor completed-rep count on the first error rate; all errors advance together.
    while (length(maps.global[[p]][[anchor_err]]) < length(seed)) {
      tryCatch({
        mrk.id <- sample(s_all[["seq.mrk.names"]], sample_size)
        s <- suppressWarnings(get_genomic_order(make_seq_mappoly(s_all, mrk.id)))
        s.p <- make_seq_mappoly(s, info.parent = p)
        maps.seq <- gen_map_sequential(s.p, max_cores = max_cores, verbose = verbose)
        for (err in global_errors) {
          maps.global[[p]][[err]][[length(maps.global[[p]][[err]]) + 1]] <-
            gen_map_global(maps.seq, error = as.numeric(err), verbose = verbose)
        }
      }, error = function(e) {})
    }
    if (!is.null(dat[["geno"]])) { # Only generate map on prior genotype probabilities if they exist
      set.seed(seed[[1]])
      while (length(maps.prob[[p]]) < length(seed)) {
        tryCatch({
          mrk.id <- sample(s_all[["seq.mrk.names"]], sample_size)
          s <- suppressWarnings(get_genomic_order(make_seq_mappoly(s_all, mrk.id)))
          s.p <- make_seq_mappoly(s, info.parent = p)
          maps.seq <- gen_map_sequential(s.p, max_cores = max_cores, verbose = verbose)
          maps.prob[[p]][[length(maps.prob[[p]]) + 1]] <- gen_map_prob(maps.seq, dat, verbose = verbose)
        }, error = function(e) {})
      }
    }
  }

  maps.ret <- maps_init
  for (p in map_parentals) {
    maps.ret[[p]][["map.global"]] <- maps.global[[p]]   # named list keyed by error
    maps.ret[[p]][["map.prob"]]   <- maps.prob[[p]]
  }
  return(maps.ret)
}

# ── Module 1: Preprocessing ───────────────────────────────────────────────────
#' Read VCF, build mappoly.data, filter markers, accumulate diagnostic info
#'
#' @description
#' Encapsulates the data-loading + filtering steps that precede mapping:
#'   1. read VCF via vcfR, scrub `NA` chromosome rows, write `bugfix.vcf`
#'   2. read into a `mappoly.data` object (without `filter.non.conforming`)
#'   3. compute `chisq.pval` locally (working around MAPpoly 0.4.0 behaviour)
#'   4. filter by missingness, then by segregation, then optionally by chromosome
#'   5. accumulate a per-step diagnostic data.frame `info`
#'
#' This function intentionally writes the temporary `bugfix.vcf` to the current
#' working directory to preserve the original behaviour.
#'
#' @param vcf_file              Path to input VCF.
#' @param SNPCall_program       SNP caller name (used in `data_label`).
#' @param GenotypeCall_program  Genotype caller name (drives prob.thres relaxation
#'                              for `supermassa`).
#' @param CountsFrom            Counts-source label (used in `data_label`).
#' @param parent1,parent2       Parent IDs as they appear in the VCF.
#' @param prob_thres            Probability threshold for genotype calls.
#' @param ploidy                Population ploidy (integer).
#' @param filt_segr             Logical (or "TRUE"/"FALSE" string) — apply
#'                              chi-square segregation filter.
#' @param chromosome            Optional chromosome to subset to ("" = all).
#' @param global_errors         Character vector of error rates (used to
#'                              pre-seed the `maps_init_global` containers).
#'
#' @return A named list with keys: `dat`, `seq.init`, `info`, `data_label`,
#'         `chr_suffix`, `prob.thres`, `map_parentals`, `maps_init`,
#'         `maps_init_global`.
mappoly_preprocess <- function(vcf_file,
                               SNPCall_program,
                               GenotypeCall_program,
                               CountsFrom,
                               parent1,
                               parent2,
                               prob_thres,
                               ploidy,
                               filt_segr,
                               chromosome,
                               global_errors) {
  # ── Per-parent control structures ─────────────────────────────────────────
  map_parentals <- c("p1", "p2")
  if (ploidy <= 2) {
    map_parentals <- c(map_parentals, "all")
  }
  maps_init <- setNames(vector("list", length(map_parentals)), map_parentals)
  maps_init_global <- setNames(lapply(map_parentals, function(p) {
    setNames(lapply(global_errors, function(e) list()), global_errors)
  }), map_parentals)
  stopifnot(is.character(vcf_file), length(vcf_file) == 1, nzchar(vcf_file))
  stopifnot(file.exists(vcf_file))
  stopifnot(is.character(parent1), is.character(parent2))
  stopifnot(is.numeric(prob_thres), length(prob_thres) == 1)
  stopifnot(is.numeric(ploidy) || is.integer(ploidy))

  # ── Read VCF and prep dat ──────────────────────────────────────────────────
  if (GenotypeCall_program == "supermassa") {
    prob.thres <- prob_thres - prob_thres * 0.3
  } else {
    prob.thres <- prob_thres
  }

  bugfix <- read.vcfR(vcf_file)
  idx <- which(bugfix@fix[, 1] == "NA")

  if (length(idx) > 0) {
    bugfix@fix <- bugfix@fix[-idx, ]
    bugfix@gt  <- bugfix@gt[-idx, ]
  }

  write.vcf(bugfix, file = "bugfix.vcf")

  dat <- read_vcf(file = "bugfix.vcf",
                  parent.1 = parent1,
                  parent.2 = parent2,
                  verbose = FALSE,
                  read.geno.prob = TRUE,
                  prob.thres = prob.thres,
                  ploidy = ploidy,
                  filter.non.conforming = FALSE)
  # MAPpoly functions call get(data.name, pos=1) to look up the data object in
  # .GlobalEnv. Sync dat there so those lookups succeed when we're inside a function.
  assign("dat", dat, envir = .GlobalEnv)

  # MAPpoly 0.4.0 only computes chisq.pval inside filter.non.conforming=TRUE;
  # we skip that filter intentionally (to preserve dat$geno probabilities for
  # est_full_hmm_with_prior_prob) but still need chisq.pval for filter_segregation.
  # Some offspring dosage calls violate Mendelian expectations given the parents
  # (e.g. dosage=2 progeny from a 0x1 cross — genotyping errors). mrk_chisq_test
  # crashes on these because seg.exp drops zero-expectation classes then name-lookup
  # returns NA. Locally scrub those calls to missing (pl+1) only for the chi-square
  # computation, leaving dat$geno.dose and dat$geno intact.
  if (is.null(dat[["chisq.pval"]])) {
    pl <- dat[["ploidy"]]
    Ds <- array(NA, dim = c(pl + 1, pl + 1, pl + 1))
    for (i in 0:pl) for (j in 0:pl)
      Ds[i + 1, j + 1, ] <- segreg_poly(ploidy = pl, dP = i, dQ = j)
    Dpop <- cbind(dat[["dosage.p1"]], dat[["dosage.p2"]])
    M <- t(apply(Dpop, 1, function(x) Ds[x[1] + 1, x[2] + 1, ]))
    dimnames(M) <- list(dat[["mrk.names"]], c(0:pl))
    allowed <- M != 0
    geno_scrub <- as.matrix(dat[["geno.dose"]])
    for (i in seq_len(nrow(geno_scrub))) {
      bad <- !(geno_scrub[i, ] %in% c(which(allowed[i, ]) - 1, pl + 1))
      if (any(bad)) geno_scrub[i, bad] <- pl + 1
    }
    M <- cbind(M, geno_scrub)
    dat[["chisq.pval"]] <- apply(M, 1, mappoly:::mrk_chisq_test, ploidy = pl)
    assign("dat", dat, envir = .GlobalEnv)
  }

  # ── Diagnostic info accumulator ────────────────────────────────────────────
  data_label <- paste0(SNPCall_program, "_", GenotypeCall_program, "_", CountsFrom)

  info <- data.frame(dat = data_label,
                     step = "raw",
                     n.markers = dat[["n.mrk"]],
                     mis.perc = round((sum(dat[["geno.dose"]] == dat[["ploidy"]] + 1) /
                                         (nrow(dat[["geno.dose"]]) * ncol(dat[["geno.dose"]]))) * 100, 2),
                     n.redundant = round(100 * (nrow(dat[["elim.correspondence"]]) /
                                                  (length(dat[["kept"]]) + nrow(dat[["elim.correspondence"]]))), 2),
                     n.ind = dat[["n.ind"]])

  dat <- filter_missing(input.data = dat, type = "marker",
                        filter.thres = 0.25, inter = FALSE)
  assign("dat", dat, envir = .GlobalEnv)

  info_temp <- data.frame(dat = data_label,
                          step = "miss filtered",
                          n.markers = dat[["n.mrk"]],
                          mis.perc = round((sum(dat[["geno.dose"]] == dat[["ploidy"]] + 1) /
                                              (nrow(dat[["geno.dose"]]) * ncol(dat[["geno.dose"]]))) * 100, 2),
                          n.redundant = round(100 * (nrow(dat[["elim.correspondence"]]) /
                                                       (length(dat[["kept"]]) + nrow(dat[["elim.correspondence"]]))), 2),
                          n.ind = dat[["n.ind"]])
  info <- rbind(info, info_temp)

  if (as.logical(filt_segr)) {
    pval.bonf <- 0.05 / dat[["n.mrk"]]
    mrks.chi.filt <- filter_segregation(dat,
                                        chisq.pval.thres = pval.bonf,
                                        inter = FALSE)
    seq.init <- make_seq_mappoly(mrks.chi.filt)
  } else {
    seq.init <- make_seq_mappoly(dat, "all")
  }

  info_temp <- data.frame(dat = data_label,
                          step = "segr filtered",
                          n.markers = length(seq.init[["seq.mrk.names"]]),
                          mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% seq.init[["seq.mrk.names"]]), ] == dat[["ploidy"]] + 1) /
                                              (ncol(dat[["geno.dose"]]) * length(seq.init[["seq.mrk.names"]]))) * 100, 2),
                          n.redundant = round(100 * (nrow(dat[["elim.correspondence"]]) /
                                                       (length(dat[["kept"]]) + nrow(dat[["elim.correspondence"]]))), 2),
                          n.ind = dat[["n.ind"]])
  info <- rbind(info, info_temp)

  # Filter to a specific chromosome when requested
  chr_suffix <- ""
  if (nzchar(chromosome)) {
    chr_suffix <- paste0("_", chromosome)
    chr_idx <- which(seq.init[["chrom"]] == chromosome)
    if (length(chr_idx) > 0) {
      chr_names <- seq.init[["seq.mrk.names"]][chr_idx]
      seq.init <- make_seq_mappoly(seq.init, chr_names)
    } else {
      # No markers on this chromosome — downstream size check will skip mapping
      seq.init[["seq.mrk.names"]] <- character(0)
    }
    info_temp <- data.frame(dat = data_label,
                            step = "chromosome filtered",
                            n.markers = length(seq.init[["seq.mrk.names"]]),
                            mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% seq.init[["seq.mrk.names"]]), ] == dat[["ploidy"]] + 1) /
                                                (ncol(dat[["geno.dose"]]) * length(seq.init[["seq.mrk.names"]]))) * 100, 2),
                            n.redundant = round(100 * (nrow(dat[["elim.correspondence"]]) /
                                                         (length(dat[["kept"]]) + nrow(dat[["elim.correspondence"]]))), 2),
                            n.ind = dat[["n.ind"]])
    info <- rbind(info, info_temp)
  }

  # Calculate some additional informational stats for the parentals
  for (p in map_parentals) {
    s.p <- make_seq_mappoly(seq.init, info.parent = p)
    info_temp <- data.frame(dat = data_label,
                            step = p,
                            n.markers = length(s.p[["seq.mrk.names"]]),
                            mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% s.p[["seq.mrk.names"]]), ] == dat[["ploidy"]] + 1) /
                                                (ncol(dat[["geno.dose"]]) * length(s.p[["seq.mrk.names"]]))) * 100, 2),
                            n.redundant = round(100 * (nrow(dat[["elim.correspondence"]]) /
                                                         (length(dat[["kept"]]) + nrow(dat[["elim.correspondence"]]))), 2),
                            n.ind = dat[["n.ind"]])
    info <- rbind(info, info_temp)
  }

  list(dat = dat,
       seq.init = seq.init,
       info = info,
       data_label = data_label,
       chr_suffix = chr_suffix,
       prob.thres = prob.thres,
       map_parentals = map_parentals,
       maps_init = maps_init,
       maps_init_global = maps_init_global)
}

# ── Module 2: Bootstrap ───────────────────────────────────────────────────────
#' Run parallel bootstrap mapping with global-error and prior-prob arms
#'
#' @description
#' Runs `repetitions` bootstrap iterations split across `max_cores` PSOCK
#' workers. Each worker subsamples `sample_size` markers, runs sequential
#' phasing, then estimates per-error global-error maps and (optionally) a
#' prior-probability map. Results are flattened per (parent, error) and
#' summarized into a tidy `summaries` data.frame.
#'
#' If `length(seq.init$seq.mrk.names) <= sample_size`, the bootstrap is skipped
#' and the `summaries` / `info` / `dat` / `mat` slots in the return value are
#' set to the sentinel value `1` — preserving the original script's behavior
#' so the downstream RDS exports remain byte-identical.
#'
#' @param seq.init            mappoly.sequence object from preprocessing.
#' @param dat                 mappoly.data object from preprocessing.
#' @param info                Diagnostic data.frame from preprocessing.
#' @param data_label,chr_suffix  Naming components carried from preprocessing.
#' @param sample_size,repetitions,max_cores   Bootstrap params.
#' @param map_parentals       Per-parent labels (from `mappoly_preprocess` return).
#' @param maps_init,maps_init_global  Pre-shaped result containers (from
#'                            `mappoly_preprocess` return).
#' @param global_errors       Character vector of error rates.
#'
#' @return A named list with keys: `boot.refactored`, `summaries`, `info`,
#'         `dat`, `mat`, `ran_bootstrap`. The first five mirror the script's
#'         original variables (and inherit the sentinel-`1` overwrites in
#'         the skip branch). `ran_bootstrap` is a logical flag for the caller.
#'
#' @note RNG: `set.seed()` is invoked inside `run_bootstraps` (driven by
#'       per-worker seed slices `seeds[[i]][[1]]`), so reproducibility is
#'       identical to the original script as long as `repetitions` and
#'       `max_cores` are unchanged.
#' @note Helper functions (`gen_map_sequential`, `gen_map_global`,
#'       `gen_map_prob`, `run_bootstraps`) are resolved from `.GlobalEnv`
#'       at call time. Source the full `mappoly_report.R` before calling
#'       this function externally.
mappoly_bootstrap <- function(seq.init,
                              dat,
                              info,
                              data_label,
                              chr_suffix,
                              sample_size,
                              repetitions,
                              max_cores,
                              map_parentals,
                              maps_init,
                              maps_init_global,
                              global_errors,
                              verbose = FALSE) {
  stopifnot(is.list(seq.init), !is.null(seq.init[["seq.mrk.names"]]))

  boot.refactored <- maps_init
  summaries <- data.frame()
  mat <- NULL
  ran_bootstrap <- FALSE

  if (length(seq.init[["seq.mrk.names"]]) > sample_size) {
    ran_bootstrap <- TRUE
    seeds <- split(1:repetitions, rep(1:max_cores, each = repetitions / max_cores))

    cat("<MapPolyReport() Starting cluster and run_bootstraps()>\n")
    clust <- makeCluster(max_cores)

    clusterExport(clust, c("seq.init", "sample_size", "dat", "run_bootstraps",
                           "maps_init", "maps_init_global", "map_parentals",
                           "global_errors", "max_cores",
                           "gen_map_sequential",
                           "gen_map_global",
                           "gen_map_prob"))
    results <- parLapply(clust, seeds, function(x) {
      suppressPackageStartupMessages(library(mappoly))
      run_bootstraps(seed = x, s_all = seq.init, sample_size = sample_size, dat = dat, max_cores = max_cores, verbose = verbose)
    })
    stopCluster(clust)
    cat("</MapPolyReport() Stopping cluster and run_bootstraps()>\n")

    # Flatten per-worker bootstrap maps into a single list per (parent, error), drop NULLs,
    # summarize once per (parent, prob_source, global_error) tuple, and accumulate into a
    # tidy data.frame.
    for (p in map_parentals) {
      boot.refactored[[p]][["map.global"]] <- setNames(vector("list", length(global_errors)), global_errors)
      for (err in global_errors) {
        global_maps <- unlist(lapply(results, function(r) r[[p]][["map.global"]][[err]]), recursive = FALSE)
        global_maps <- Filter(Negate(is.null), global_maps)
        boot.refactored[[p]][["map.global"]][[err]] <- global_maps
        if (length(global_maps) > 0) {
          s <- summary_maps(global_maps)
          s[["parent"]] <- p
          s[["prob_source"]] <- "global"
          s[["global_error"]] <- as.numeric(err)
          summaries <- rbind(summaries, s)
        }
      }

      prob_maps <- unlist(lapply(results, function(r) r[[p]][["map.prob"]]), recursive = FALSE)
      prob_maps <- Filter(Negate(is.null), prob_maps)
      boot.refactored[[p]][["map.prob"]] <- prob_maps
      if (length(prob_maps) > 0) {
        s <- summary_maps(prob_maps)
        s[["parent"]] <- p
        s[["prob_source"]] <- "genoprobs"
        s[["global_error"]] <- NA_real_
        summaries <- rbind(summaries, s)
      }
    }
    if (nrow(summaries) > 0) {
      summaries[["data"]] <- paste0(data_label, chr_suffix)
    }
  } else {
    # Sentinel overwrites preserve byte-identical RDS contents in the
    # "too few markers" branch of the original script.
    summaries <- 1
    info <- 1
    dat <- 1
  }

  list(boot.refactored = boot.refactored,
       summaries = summaries,
       info = info,
       dat = dat,
       mat = mat,
       ran_bootstrap = ran_bootstrap)
}

# ── Module 3: Full linkage map ────────────────────────────────────────────────
#' Build the optional full-marker production linkage map per parent
#'
#' @description
#' Production full-marker linkage map per parent. Independent of the
#' bootstrap; wrapped in tryCatch so HMM failures (non-convergence on very
#' large marker sets, especially the SNPCaller-direct arm) do not destroy
#' the bootstrap output.
#'
#' Returns `list(full_map = <list or NULL>)` without writing any files;
#' file export is handled by `export_mappoly_results()`.
#'
#' @param build_full_map   Logical gate. When FALSE, the function returns
#'                         `list(full_map = NULL)` immediately.
#' @param seq.init         mappoly.sequence object from preprocessing.
#' @param dat              mappoly.data object from preprocessing.
#' @param data_label,chr_suffix  Naming components carried from preprocessing.
#' @param map_parentals    Per-parent labels (from `mappoly_preprocess` return).
#' @param maps_init,maps_init_global  Pre-shaped result containers (from
#'                         `mappoly_preprocess` return).
#' @param global_errors    Character vector of error rates.
#' @param verbose          Logical indicating verbose setting for mappoly map generating functions.
#'                         Defaults to FALSE.
#'
#' @note Helper functions (`gen_map_sequential`, `gen_map_global`,
#'       `gen_map_prob`) are resolved from `.GlobalEnv` at call time. Source
#'       the full `mappoly_report.R` before calling this function externally.
#'
#' @return `list(full_map = <list or NULL>)`.
mappoly_full_map <- function(build_full_map,
                             seq.init,
                             dat,
                             data_label,
                             chr_suffix,
                             map_parentals,
                             maps_init,
                             maps_init_global,
                             global_errors,
                             verbose = FALSE) {
  stopifnot(is.list(seq.init), !is.null(seq.init[["seq.mrk.names"]]))

  map.full.global <- maps_init_global
  map.full.prob   <- maps_init
  map.full.n_mrk  <- maps_init
  s <- suppressWarnings(get_genomic_order(seq.init))
  for (p in map_parentals) {
    s.p <- make_seq_mappoly(s, info.parent = p)
    map.full.n_mrk[[p]] <- n_mrk <- length(s.p[["seq.mrk.names"]])
    cat("<MapPolyReport() Building full map for ", p, ", n = ", n_mrk, ">\n")
    if (n_mrk >= 3) {
      cat("<MapPolyReport() gen_map_sequential() for ", p, ">\n")
      tryCatch({
        map.seq <- gen_map_sequential(s.p, verbose = verbose)
        for (err in global_errors) {
          cat("<MapPolyReport() gen_map_global() for ", p, " error = ", err, ">\n")
          tryCatch({
            map.full.global[[p]][[err]] <- gen_map_global(map.seq, error = as.numeric(err), verbose = verbose)
          }, error = function(e) {
            message("Full-marker (", p, ", err=", err, ") gen_map_global() failed: ", conditionMessage(e))
          })
          cat("</MapPolyReport() gen_map_global() for ", p, " error = ", err, ">\n")
        }
        cat("<MapPolyReport() gen_map_prob() for ", p, ">\n")
        tryCatch({
          map.full.prob[[p]] <- gen_map_prob(map.seq, dat, verbose = verbose)
        }, error = function(e) {
          message("Full-marker (", p, ") map gen_map_prob() failed: ", conditionMessage(e))
        })
        cat("</MapPolyReport() gen_map_prob() for ", p, ">\n")
      }, error = function(e) {
        message("Full-marker (", p, ") map gen_map_sequential() failed: ", conditionMessage(e))
      })
      cat("</MapPolyReport() gen_map_sequential() for ", p, ">\n")
    }
    cat("</MapPolyReport() Building full map for ", p, ", n = ", n_mrk, ">\n")
  }
  full_map <- maps_init
  for (p in map_parentals) {
    full_map[[p]][["map.global"]] <- map.full.global[[p]]   # named list keyed by error
    full_map[[p]][["map.prob"]]   <- map.full.prob[[p]]
    full_map[[p]][["n_mrk"]]      <- map.full.n_mrk[[p]]
  }
  list(full_map = full_map)
}

# ── Export helper ─────────────────────────────────────────────────────────────
#' Export bootstrap/preprocess outputs and the two-point matrix to a tar.gz
#'
#' @description
#' Sits at the join of all three modules: it computes the two-point
#' recombination matrix from the preprocessing-derived `seq.init`, then writes
#' six RDS files (summaries, info, dat, mat, boot.refactored, full_map), moves
#' them into `results/`, and tars the directory. Lives outside the producing
#' modules because the bundled outputs span more than one module's data.
#'
#' Note: when `bootstrap_result$ran_bootstrap` is FALSE, `summaries`/`info`/
#' `dat` here will be the sentinel value `1` propagated by `mappoly_bootstrap`,
#' matching the original script's RDS contents exactly.
#'
#' @param seq.init             mappoly.sequence from preprocessing.
#' @param max_cores            Cores for `est_pairwise_rf`.
#' @param data_label,chr_suffix  Naming components.
#' @param bootstrap_result     The list returned by `mappoly_bootstrap`.
#' @param full_map             The list returned by `mappoly_full_map` (may be NULL).
#' @param output_dir           Directory in which `results/` and the tar live.
#'
#' @return Invisibly, the path to the tar.gz.
export_mappoly_results <- function(seq.init,
                                   max_cores,
                                   data_label,
                                   chr_suffix,
                                   bootstrap_result,
                                   full_map = NULL,
                                   output_dir = ".") {
  stopifnot(is.list(bootstrap_result),
            all(c("boot.refactored", "summaries", "info", "dat") %in% names(bootstrap_result)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prev_wd <- setwd(output_dir)
  on.exit(setwd(prev_wd), add = TRUE)

  # ── Two-point recombination matrix ───────────────────────────────────────
  # Estimate two-point recombination fraction of filtered sequences — saving
  # matrix to *.rds for downstream analysis and diagnostics.
  tpt <- est_pairwise_rf(input.seq = seq.init, ncpus = max_cores)
  mat <- rf_list_to_matrix(input.twopt = tpt)

  summaries       <- bootstrap_result[["summaries"]]
  info            <- bootstrap_result[["info"]]
  dat             <- bootstrap_result[["dat"]]
  boot.refactored <- bootstrap_result[["boot.refactored"]]

  saveRDS(summaries,        file = paste0(data_label, chr_suffix, "_summaries.rds"))
  saveRDS(info,             file = paste0(data_label, chr_suffix, "_info.rds"))
  saveRDS(dat,              file = paste0(data_label, chr_suffix, "_dat.rds"))
  saveRDS(mat,              file = paste0(data_label, chr_suffix, "_mat2.rds"))
  saveRDS(boot.refactored,  file = paste0(data_label, chr_suffix, "_maps.rds"))
  if (!is.null(full_map))
    saveRDS(full_map,       file = paste0(data_label, chr_suffix, "_full_map.rds"))

  dir.create("results", showWarnings = FALSE)
  rds_files <- list.files(pattern = "\\.rds$", full.names = FALSE)
  file.rename(rds_files, file.path("results", rds_files))

  tar_name <- paste0(data_label, chr_suffix, "_poly_results.tar.gz")
  system(paste0("tar -czvf ", shQuote(tar_name), " results"))
  invisible(tar_name)
}

# ── Orchestrator ──────────────────────────────────────────────────────────────
#' Drives the four modules in sequence. 
#' @description
#' Mirrors the original top-to-bottom
#' script flow exactly so `Rscript scripts/mappoly_report.R 'name=value' ...`
#' (as invoked from tasks/utilsR.wdl in the cristaniguti/reads2map:0.0.8
#' container) continues to produce the same `<data_label><chr_suffix>_poly_results.tar.gz`.
#'
#' @param vcf_file                Character string path of the variant call file.
#' @param SNPCall_program         Character string name of the program used to generate the SNP variant calls.
#' @param GenotypeCall_program    Character string: If regenotyping was performed, the name of the program/framework used to generate genotypes.
#' @param CountsFrom              Character string: Counts-source label ['vcf','bam']
#' @param parent1                 Character string identifier for parent 1.
#' @param parent2                 Character string identifier for parent 2.
#' @param prob_thres              Probability threshold for genotype calls. Set to 0 for accepting all probabilities.
#' @param repetitions             Integer number of bootstraps.
#' @param sample_size             Integer number of random markers to subset from all markers when bootstrapping linkage maps.
#' @param max_cores               Integer number of cores to run parallel bootstrapping process.
#' @param ploidy                  Integer indicating the ploidy level (2 for diploids, 4 for hexaploids, etc.)
#' @param filt_segr               Logical indicating whether to filter markers demonstrating significant segregation distortion.
#' @param global_errors           Character vector of error rates (converted to floating point numerics).
#' @param chromosome              String representing chromosome identifier to perform linkage map generation against.
#' @param build_full_map          Logical gate. When FALSE, `mappoly_full_map`
#'                                returns `list(full_map = NULL)` and no
#'                                `_full_map.rds` is written by the export step.
#' @param output_dir              Directory in which `results/` and the tar live.
#'                                Defaults to `"."` (current working directory).
#' 
# ── Defaults (mirror WDL task defaults) ───────────────────────────────────────
mappoly_report_pipeline <- function(vcf_file,
                                    SNPCall_program,
                                    GenotypeCall_program,
                                    CountsFrom,
                                    parent1,
                                    parent2,
                                    prob_thres,
                                    repetitions,
                                    sample_size,
                                    max_cores,
                                    ploidy,
                                    filt_segr,
                                    global_errors,
                                    chromosome,
                                    build_full_map,
                                    output_dir,
                                    verbose = FALSE) {
  pre <- mappoly_preprocess(vcf_file             = vcf_file,
                            SNPCall_program      = SNPCall_program,
                            GenotypeCall_program = GenotypeCall_program,
                            CountsFrom           = CountsFrom,
                            parent1              = parent1,
                            parent2              = parent2,
                            prob_thres           = prob_thres,
                            ploidy               = ploidy,
                            filt_segr            = filt_segr,
                            chromosome           = chromosome,
                            global_errors        = global_errors)

  boot <- mappoly_bootstrap(seq.init          = pre[["seq.init"]],
                            dat               = pre[["dat"]],
                            info              = pre[["info"]],
                            data_label        = pre[["data_label"]],
                            chr_suffix        = pre[["chr_suffix"]],
                            sample_size       = sample_size,
                            repetitions       = repetitions,
                            max_cores         = max_cores,
                            map_parentals     = pre[["map_parentals"]],
                            maps_init         = pre[["maps_init"]],
                            maps_init_global  = pre[["maps_init_global"]],
                            global_errors     = global_errors,
                            verbose           = verbose)

  # The full-map module reads from preprocessing's outputs (not the bootstrap
  # sentinel), so we pass `pre$dat` rather than `boot$dat`. This matches the
  # original script: `build_full_map` ran before the bootstrap-sentinel
  # overwrites of `dat`/`info` could affect anything outside the bootstrap
  # block, and operates on the live mappoly.data object.
  stopifnot(is.logical(build_full_map), length(build_full_map) == 1)
  if (!build_full_map) {
    full <- list(full_map = NULL)
  }
  full <- mappoly_full_map(seq.init          = pre[["seq.init"]],
                          dat               = pre[["dat"]],
                          data_label        = pre[["data_label"]],
                          chr_suffix        = pre[["chr_suffix"]],
                          map_parentals     = pre[["map_parentals"]],
                          maps_init         = pre[["maps_init"]],
                          maps_init_global  = pre[["maps_init_global"]],
                          global_errors     = global_errors,
                          verbose           = verbose)

  export_mappoly_results(seq.init         = pre[["seq.init"]],
                        max_cores        = max_cores,
                        data_label       = pre[["data_label"]],
                        chr_suffix       = pre[["chr_suffix"]],
                        bootstrap_result = boot,
                        full_map         = full[["full_map"]],
                        output_dir       = ".")
}

if (sys.nframe() == 0) { #The mappoly_report.R file is being called directly from an RScript. Run as 'main'
  vcf_file             <- get0("vcf_file",             ifnotfound = NULL)
  SNPCall_program      <- get0("SNPCall_program",      ifnotfound = NULL)
  GenotypeCall_program <- get0("GenotypeCall_program", ifnotfound = NULL)
  CountsFrom           <- get0("CountsFrom",           ifnotfound = NULL)
  parent1              <- get0("parent1",              ifnotfound = NULL)
  parent2              <- get0("parent2",              ifnotfound = NULL)
  prob_thres           <- get0("prob_thres",           ifnotfound = 0.8)
  repetitions          <- get0("repetitions",          ifnotfound = 100)
  sample_size          <- get0("sample_size",          ifnotfound = 30)
  max_cores            <- get0("max_cores",            ifnotfound = NULL)
  ploidy               <- get0("ploidy",               ifnotfound = NULL)
  filt_segr            <- get0("filt_segr",            ifnotfound = "TRUE")
  global_errors        <- get0("global_errors",        ifnotfound = "0.05")
  chromosome           <- get0("chromosome",           ifnotfound = "")
  build_full_map       <- get0("build_full_map",       ifnotfound = TRUE)
  output_dir           <- get0("output_dir",           ifnotfound = ".")
  verbose              <- get0("verbose",              ifnotfound = TRUE)

  # ── Parse name=value CLI args ─────────────────────────────────────────────────
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    cat("No arguments supplied; relying on defaults / pre-set globals.\n")
  } else {
    for (a in args) eval(parse(text = a))
  }

  # ── Validate & coerce ─────────────────────────────────────────────────────────
  required <- c("vcf_file", "SNPCall_program", "GenotypeCall_program",
                "CountsFrom", "parent1", "parent2", "max_cores", "ploidy")
  missing_args <- required[vapply(required, function(n) is.null(get(n)), logical(1))]
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  }

  prob_thres  <- as.numeric(prob_thres)
  repetitions <- as.integer(repetitions)
  sample_size <- as.integer(sample_size)
  max_cores   <- as.integer(max_cores)
  ploidy      <- as.integer(ploidy)
  if (is.character(build_full_map)) build_full_map <- as.logical(build_full_map)

  verbose = as.logical(verbose) #Convert to logical if not already

  # global_errors → character vector, drop empties
  if (is.character(global_errors) && length(global_errors) == 1) {
    global_errors <- unlist(strsplit(global_errors, ","))
  }
  global_errors <- as.character(global_errors)
  global_errors <- global_errors[nzchar(global_errors)]
  if (length(global_errors) == 0) {
    stop("global_errors must contain at least one rate")
  }

  # Run all output writes from output_dir so relative file names land there.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prev_wd <- setwd(output_dir)
  on.exit(setwd(prev_wd), add = TRUE)

  #Run the pipeline
  mappoly_report_pipeline(vcf_file,
                          SNPCall_program,
                          GenotypeCall_program,
                          CountsFrom,
                          parent1,
                          parent2,
                          prob_thres,
                          repetitions,
                          sample_size,
                          max_cores,
                          ploidy,
                          filt_segr,
                          global_errors,
                          chromosome,
                          build_full_map,
                          output_dir)
}
