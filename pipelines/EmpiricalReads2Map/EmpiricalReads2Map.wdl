version 1.0

import "../../structs/empirical_maps_structs.wdl"
import "../../structs/dna_seq_structs.wdl"

import "../../tasks/utils.wdl" as utils

import "../../pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl" as snpcalling
import "../../pipelines/EmpiricalMaps/EmpiricalMaps.wdl" as maps

workflow EmpiricalReads {

    input {
        File samples_info
        ReferenceFasta references

        # One entry per mapping population.  SNP calling runs on all samples
        # combined; VCFs are subsetted per population before map building.
        Array[PopulationSpec] populations

        # GBS key files used to identify which samples belong to each pedigree.
        Array[File] key_files

        Int max_cores
        Int max_ram
        Int chunk_size
        Boolean rm_dupli = false
        Boolean gatk_mchap = false
        Boolean hardfilters = true
        Boolean replaceAD = true
        String replaceADbyMissing = "TRUE" # Boolean inside R
        Boolean run_gatk = true
        Boolean run_freebayes = true
        Boolean run_tassel = true
        Boolean pair_end = false
        String? enzyme
        Int ploidy = 2
        Int n_chrom
        String? filters
        Float? prob_thres
        String? filt_segr
        Boolean filter_noninfo = false
    }

    # Step 1 — call SNPs on all samples combined
    call snpcalling.SNPCalling {
        input:
            samples_info  = samples_info,
            references    = references,
            max_cores     = max_cores,
            rm_dupli      = rm_dupli,
            gatk_mchap    = gatk_mchap,
            hardfilters   = hardfilters,
            replaceAD     = replaceAD,
            run_gatk      = run_gatk,
            run_freebayes = run_freebayes,
            run_tassel    = run_tassel,
            ploidy        = ploidy,
            n_chrom       = n_chrom,
            enzyme        = enzyme,
            pair_end      = pair_end
    }

    # Step 2 — for each population: subset VCFs, then build genetic maps
    scatter (pop in populations) {

        # Subset every caller's VCF to this population's samples
        scatter (vcf_idx in range(length(SNPCalling.vcfs))) {
            call utils.SubsetVcfToPopulation {
                input:
                    vcf_file         = SNPCalling.vcfs[vcf_idx],
                    key_files        = key_files,
                    pedigree         = pop.pedigree,
                    pedigree_aliases = pop.pedigree_aliases,
                    parent1_seedlot  = pop.dataset.parent1,
                    parent2_seedlot  = pop.dataset.parent2,
                    population_name  = pop.dataset.name
            }
        }

        # Also subset the optional GATK multi-allelic VCF when present
        if (defined(SNPCalling.gatk_multi_vcf)) {
            call utils.SubsetVcfToPopulation as SubsetMultiVcf {
                input:
                    vcf_file         = select_first([SNPCalling.gatk_multi_vcf]),
                    key_files        = key_files,
                    pedigree         = pop.pedigree,
                    pedigree_aliases = pop.pedigree_aliases,
                    parent1_seedlot  = pop.dataset.parent1,
                    parent2_seedlot  = pop.dataset.parent2,
                    population_name  = pop.dataset.name + "_multi"
            }
        }

        # Build genetic maps for this population using its subsetted VCFs
        call maps.Maps {
            input:
                dataset            = pop.dataset,
                vcfs               = SubsetVcfToPopulation.subset_vcf,
                vcfs_software      = SNPCalling.vcfs_software,
                vcfs_counts_source = SNPCalling.vcfs_counts_source,
                gatk_vcf_multi     = SubsetMultiVcf.subset_vcf,
                gatk_mchap         = gatk_mchap,
                filter_noninfo     = filter_noninfo,
                filters            = filters,
                max_cores          = max_cores,
                replaceADbyMissing = replaceADbyMissing,
                ploidy             = ploidy,
                prob_thres         = prob_thres,
                filt_segr          = filt_segr
        }
    }

    output {
        # One result archive per population
        Array[File] EmpiricalReads_results = Maps.EmpiricalReads_results
    }
}
