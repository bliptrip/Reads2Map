version 1.0

import "../tasks/BWA.wdl" as alg
import "../tasks/chunk_lists.wdl"
import "../tasks/utils.wdl" as utils

workflow CreateAlignmentFromFamilies {
    input {
        File         families_info
        Array[File]  trimmed_fastqs
        ReferenceFasta references
        Int max_cores
        Boolean rm_dupli
        Boolean gatk_mchap
        Boolean pair_end
        Int chunk_size
        Boolean run_merge_bams
    }

    call chunk_lists.SepareChunksFastqString {
        input:
            families_info=families_info,
            chunk_size = chunk_size
    }

    scatter (chunk in SepareChunksFastqString.chunks) {

        Array[Array[String]] sample_file = read_tsv(chunk)

        # Resolve basenames from TSV column 0 to proper WDL File references
        call utils.ResolveFastqsByBasename {
            input:
                basenames  = sample_file[0],
                all_fastqs = trimmed_fastqs
        }

        if(pair_end) {
            Array[File] pair = sample_file[3]
        }

        call alg.RunBwaAlignment {
            input:
                sampleName  = sample_file[1],
                reads1      = ResolveFastqsByBasename.resolved,
                reads2      = pair,
                libraries   = sample_file[2],
                pair_end    = pair_end,
                references  = references,
                max_cores   = max_cores,
                rm_dupli    = rm_dupli
        }
    }

    # Store for MCHap
    if (run_merge_bams) {
        call utils.MergeBams {
            input:
                bam_files = flatten(RunBwaAlignment.bam)
        }
    }

    output {
        Array[File] bam = flatten(RunBwaAlignment.bam)
        Array[File] bai = flatten(RunBwaAlignment.bai)
        Array[Array[File]] dup_metrics = RunBwaAlignment.dup_metrics
        File? merged_bam = MergeBams.merged_bam
    }
}
