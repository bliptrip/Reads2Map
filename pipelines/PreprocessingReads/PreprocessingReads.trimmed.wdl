version 1.0


import "structs/preprocessing_reads_structs.wdl"

import "tasks/utils.wdl"


workflow PreprocessingReads{
    input {
      Specifications spec
      Array[File] trim_seq
    }

    call utils.GenerateSamplesInfo {
      input:
        trimmed_fastqs = trim_seq,
        key_files      = spec.barcode_key_files
    }

    output {
      File samples_info = GenerateSamplesInfo.samples_info
    }
}
