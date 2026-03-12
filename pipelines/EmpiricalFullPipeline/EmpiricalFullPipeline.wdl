version 1.0

# Combined pipeline: demultiplex + adapter-trim raw reads, index the reference
# genome, then align, call SNPs, and build genetic maps.
# PreprocessingReads produces the samples_info file that EmpiricalReads2Map
# consumes, and GenomeIndex produces all BWA/GATK index files from the raw
# FASTA, so the user only needs to supply the FASTA itself.

import "structs/preprocessing_reads_structs.wdl"
import "structs/empirical_maps_structs.wdl"
import "structs/dna_seq_structs.wdl"

import "pipelines/PreprocessingReads/PreprocessingReads.wdl" as preprocessing
import "pipelines/GenomeIndex/GenomeIndex.wdl"               as indexing
import "pipelines/EmpiricalReads2Map/EmpiricalReads2Map.wdl" as mapping

workflow EmpiricalFullPipeline {

    input {
        # ---------------------------------------------------------------
        # Sequencing library specifications (demultiplexing + trimming)
        # ---------------------------------------------------------------
        Specifications spec

        # ---------------------------------------------------------------
        # Reference genome — FASTA only; all index files are generated
        # by the GenomeIndex step below.
        # ---------------------------------------------------------------
        File ref_fasta

        # ---------------------------------------------------------------
        # Population / dataset metadata.  One entry per mapping population;
        # SNP calling runs on all samples, VCF subsetting and map building
        # are scattered over this array.
        # ---------------------------------------------------------------
        Array[PopulationSpec] populations

        # ---------------------------------------------------------------
        # Compute resources
        # ---------------------------------------------------------------
        Int max_cores
        Int max_ram
        Int chunk_size

        # ---------------------------------------------------------------
        # Alignment & SNP-calling options
        # ---------------------------------------------------------------
        Boolean rm_dupli      = false
        Boolean pair_end      = false
        Boolean gatk_mchap    = false
        Boolean hardfilters   = true
        Boolean replaceAD     = true
        Boolean run_gatk      = true
        Boolean run_freebayes = true
        Boolean run_tassel    = true
        Boolean run_stacks    = true
        Int     ploidy        = 2
        Int     n_chrom

        # ---------------------------------------------------------------
        # Mapping options
        # ---------------------------------------------------------------
        String  replaceADbyMissing = "TRUE"   # passed as string into R
        String? filters
        Float?  prob_thres
        String? filt_segr
        Boolean filter_noninfo         = false
        Boolean run_updog              = true
        Boolean run_supermassa         = false
        Boolean run_polyrad            = true
        Boolean run_gusmap             = false
        Array[String] global_errors    = ["0.05"]
        Boolean genoprob_error         = true
        Array[String] genoprob_global_errors = ["0.05"]
    }

    # Step 1 — demultiplex raw FASTQs and trim adapters (runs in parallel
    #          with Step 2 since neither depends on the other)
    call preprocessing.PreprocessingReads {
        input:
            spec = spec
    }

    # Step 2 — build samtools .fai, Picard .dict, and BWA index files
    call indexing.GenomeIndex {
        input:
            ref_fasta = ref_fasta
    }

    # Assemble the ReferenceFasta struct from the raw FASTA and the
    # index files produced above.
    ReferenceFasta references = object {
        ref_fasta:       ref_fasta,
        ref_fasta_index: GenomeIndex.ref_fasta_index,
        ref_dict:        GenomeIndex.ref_dict,
        ref_amb:         GenomeIndex.ref_amb,
        ref_ann:         GenomeIndex.ref_ann,
        ref_bwt:         GenomeIndex.ref_bwt,
        ref_pac:         GenomeIndex.ref_pac,
        ref_sa:          GenomeIndex.ref_sa
    }

    # Step 3 — align, call SNPs on all samples, then build per-population maps;
    #          enzyme and key files are shared with preprocessing
    call mapping.EmpiricalReads {
        input:
            samples_info       = PreprocessingReads.samples_info,
            references         = references,
            populations        = populations,
            key_files          = spec.barcode_key_files,
            max_cores          = max_cores,
            max_ram            = max_ram,
            chunk_size         = chunk_size,
            rm_dupli           = rm_dupli,
            pair_end           = pair_end,
            gatk_mchap         = gatk_mchap,
            hardfilters        = hardfilters,
            replaceAD          = replaceAD,
            run_gatk           = run_gatk,
            run_freebayes      = run_freebayes,
            run_tassel         = run_tassel,
            run_stacks         = run_stacks,
            ploidy             = ploidy,
            n_chrom            = n_chrom,
            enzyme             = spec.enzyme,
            replaceADbyMissing     = replaceADbyMissing,
            filters                = filters,
            prob_thres             = prob_thres,
            filt_segr              = filt_segr,
            filter_noninfo         = filter_noninfo,
            run_updog              = run_updog,
            run_supermassa         = run_supermassa,
            run_polyrad            = run_polyrad,
            run_gusmap             = run_gusmap,
            global_errors          = global_errors,
            genoprob_error         = genoprob_error,
            genoprob_global_errors = genoprob_global_errors
    }

    output {
        File         preprocessing_results = PreprocessingReads.results
        File         samples_info          = PreprocessingReads.samples_info
        Array[File]  mapping_results       = EmpiricalReads.EmpiricalReads_results
    }
}
