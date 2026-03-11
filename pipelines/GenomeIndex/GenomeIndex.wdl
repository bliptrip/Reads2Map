version 1.0

import "../../tasks/genome_index.wdl" as genome_index


workflow GenomeIndex {

  input {
    File ref_fasta
  }

  call genome_index.SamtoolsFaidx {
    input:
      ref_fasta = ref_fasta
  }

  call genome_index.PicardCreateDict {
    input:
      ref_fasta = ref_fasta
  }

  call genome_index.BwaIndex {
    input:
      ref_fasta = ref_fasta
  }

  output {
    File ref_fasta_index = SamtoolsFaidx.ref_fasta_index
    File ref_dict        = PicardCreateDict.ref_dict
    File ref_amb         = BwaIndex.ref_amb
    File ref_ann         = BwaIndex.ref_ann
    File ref_bwt         = BwaIndex.ref_bwt
    File ref_pac         = BwaIndex.ref_pac
    File ref_sa          = BwaIndex.ref_sa
  }
}
