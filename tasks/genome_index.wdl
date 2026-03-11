version 1.0

task SamtoolsFaidx {

  input {
    File ref_fasta
  }

  Int disk_size = ceil(size(ref_fasta, "GiB") * 2 + 5)
  Int memory_size = 4000

  command <<<
    set -euo pipefail
    cp ~{ref_fasta} ref.fasta
    samtools faidx ref.fasta
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    singularity: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "SamtoolsFaidx"
    mem: "~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Creates a FASTA index (.fai) for the reference genome using [samtools faidx](http://www.htslib.org/doc/samtools-faidx.html)."
  }

  output {
    File ref_fasta_index = "ref.fasta.fai"
  }
}

task PicardCreateDict {

  input {
    File ref_fasta
  }

  Int disk_size = ceil(size(ref_fasta, "GiB") * 2 + 5)
  Int memory_size = 4000

  command <<<
    set -euo pipefail
    cp ~{ref_fasta} ref.fasta
    java -jar /usr/gitc/picard.jar CreateSequenceDictionary \
      R=ref.fasta \
      O=ref.dict
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    singularity: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "PicardCreateDict"
    mem: "~{memory_size}M"
    time: 1
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Creates a sequence dictionary (.dict) for the reference genome using [Picard CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard)."
  }

  output {
    File ref_dict = "ref.dict"
  }
}

task BwaIndex {

  input {
    File ref_fasta
  }

  Int disk_size = ceil(size(ref_fasta, "GiB") * 4 + 10)
  Int memory_size = 8000

  command <<<
    set -euo pipefail
    cp ~{ref_fasta} ref.fasta
    /usr/gitc/./bwa index ref.fasta
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    singularity: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    cpu: 1
    # Cloud
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "BwaIndex"
    mem: "~{memory_size}M"
    time: 2
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Creates BWA index files (.amb, .ann, .bwt, .pac, .sa) for the reference genome using [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml)."
  }

  output {
    File ref_amb = "ref.fasta.amb"
    File ref_ann = "ref.fasta.ann"
    File ref_bwt = "ref.fasta.bwt"
    File ref_pac = "ref.fasta.pac"
    File ref_sa  = "ref.fasta.sa"
  }
}
