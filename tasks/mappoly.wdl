version 1.0

task MappolyReport {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Float prob_thres = 0.8
    Int repetitions = 100
    Int sample_size = 30
    Int max_cores
    Int ploidy
    String filt_segr = "TRUE"
    Array[String] global_errors = ["0.05"]
    String? chromosome
    Boolean build_full_map = true
    Boolean run_bootstraps = false
    File mappoly_report_script
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(size(vcf_file, "MiB") * max_cores + 18000*max_cores + 4000*length(global_errors))

  command <<<
    set -euo pipefail
    cat > /tmp/mappoly_run.R << 'REOF'
source("~{mappoly_report_script}")

vcf_file             <- "~{vcf_file}"
SNPCall_program      <- "~{SNPCall_program}"
GenotypeCall_program <- "~{GenotypeCall_program}"
CountsFrom           <- "~{CountsFrom}"
parent1              <- "~{parent1}"
parent2              <- "~{parent2}"
prob_thres           <- ~{prob_thres}
repetitions          <- ~{repetitions}L
sample_size          <- ~{sample_size}L
max_cores            <- ~{max_cores}L
ploidy               <- ~{ploidy}L
filt_segr            <- "~{filt_segr}"
global_errors        <- unlist(strsplit("~{sep=',' global_errors}", ","))
chromosome           <- "~{default="" chromosome}"
build_full_map       <- ~{if build_full_map then "TRUE" else "FALSE"}
run_bootstraps       <- ~{if run_bootstraps then "TRUE" else "FALSE"}

pre <- mappoly_preprocess(
  vcf_file             = vcf_file,
  SNPCall_program      = SNPCall_program,
  GenotypeCall_program = GenotypeCall_program,
  CountsFrom           = CountsFrom,
  parent1              = parent1,
  parent2              = parent2,
  prob_thres           = prob_thres,
  ploidy               = ploidy,
  filt_segr            = filt_segr,
  chromosome           = chromosome,
  global_errors        = global_errors
)

if (run_bootstraps) {
  boot <- mappoly_bootstrap(
    seq.init         = pre[["seq.init"]],
    dat              = pre[["dat"]],
    info             = pre[["info"]],
    data_label       = pre[["data_label"]],
    chr_suffix       = pre[["chr_suffix"]],
    sample_size      = sample_size,
    repetitions      = repetitions,
    max_cores        = max_cores,
    map_parentals    = pre[["map_parentals"]],
    maps_init        = pre[["maps_init"]],
    maps_init_global = pre[["maps_init_global"]],
    global_errors    = global_errors
  )
} else {
  boot <- list(
    boot.refactored = pre[["maps_init"]],
    summaries       = 1,
    info            = pre[["info"]],
    dat             = pre[["dat"]],
    mat             = NULL,
    ran_bootstrap   = FALSE
  )
}

full <- mappoly_full_map(
  build_full_map   = build_full_map,
  seq.init         = pre[["seq.init"]],
  dat              = pre[["dat"]],
  data_label       = pre[["data_label"]],
  chr_suffix       = pre[["chr_suffix"]],
  map_parentals    = pre[["map_parentals"]],
  maps_init        = pre[["maps_init"]],
  maps_init_global = pre[["maps_init_global"]],
  global_errors    = global_errors
)

export_mappoly_results(
  seq.init         = pre[["seq.init"]],
  max_cores        = max_cores,
  data_label       = pre[["data_label"]],
  chr_suffix       = pre[["chr_suffix"]],
  bootstrap_result = boot,
  full_map         = full[["full_map"]],
  output_dir       = "."
)
REOF
    Rscript /tmp/mappoly_run.R
  >>>

  runtime {
    docker:"bliptrip/reads2map:0.1.1"
    singularity: "docker://bliptrip/reads2map:0.1.1"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MappolyReport"
    mem:"~{memory_size}G"
    time: 24
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Build linkage map using MAPpoly"
  }

  output {
    File results = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}~{if defined(chromosome) then '_' + select_first([chromosome]) else ''}_poly_results.tar.gz"
  }
}