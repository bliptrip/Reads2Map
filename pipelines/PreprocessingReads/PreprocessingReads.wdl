version 1.0


import "../../structs/preprocessing_reads_structs.wdl"

import "../../tasks/stacks.wdl"
import "../../tasks/cutadapt.wdl"
import "../../tasks/utils.wdl"


workflow PreprocessingReads{
    input {
      Specifications spec
    }

    scatter (i in range(length(spec.barcode_key_files))) {
      call utils.GenerateBarcodes {
        input:
          key_file = spec.barcode_key_files[i]
      }

      call stacks.ProcessRadTags {
        input:
          enzyme = spec.enzyme,
          enzyme2 = spec.enzyme2,
          fq_files = [spec.fastq_files[i]],
          barcodes = GenerateBarcodes.barcodes
      }
    }

    scatter (sequence in flatten(ProcessRadTags.seq_results)) {
      call cutadapt.RemoveAdapt {
        input:
          sequence = sequence,
          adapter = spec.adapter,
          sequence_name = basename(sequence)
      }
    }

    call utils.TarFiles {
      input:
        sequences = RemoveAdapt.trim_seq
    }

    call utils.GenerateSamplesInfo {
      input:
        trimmed_fastqs = RemoveAdapt.trim_seq,
        key_files      = spec.barcode_key_files
    }

    output {
      File         results         = TarFiles.results
      File         samples_info    = GenerateSamplesInfo.samples_info
      Array[File]  trimmed_fastqs  = RemoveAdapt.trim_seq
    }
}
