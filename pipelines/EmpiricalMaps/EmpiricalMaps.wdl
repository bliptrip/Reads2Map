version 1.0

import "../../structs/empirical_maps_structs.wdl"
import "../../structs/population_structs.wdl"

import "../../tasks/utils.wdl" as utils
import "../../tasks/utilsR.wdl" as utilsR
import "../../tasks/JointReports.wdl" as reports
import "../../tasks/mappoly.wdl" as mappoly_task

import "../../subworkflows/genotyping_empirical.wdl" as genotyping
import "../../subworkflows/snpcaller_maps_empirical.wdl" as snpcaller
import "../../subworkflows/gusmap_maps_empirical.wdl" as gusmap
import "../../subworkflows/mappoly_maps_empirical.wdl" as mappoly_sub

workflow Maps {
    input {
        Dataset dataset
        Array[File] vcfs
        Array[String] vcfs_software
        Array[String] vcfs_counts_source
        Boolean run_updog = true
        Boolean run_supermassa = false
        Boolean run_polyrad = true
        Boolean run_gusmap = false
        Boolean filter_noninfo
        Boolean replaceADbyMissing 
        File? gatk_vcf_multi
        String gatk_mchap
        String? genotype_dp_filter
        String? filters_include
        Int info_dp_sd_multiplier = 2
        Int max_cores
        Int ploidy
        Float prob_thres = 0.8
        String? filt_segr
        Array[String] global_errors = ["0.05"]
        Boolean genoprob_error = true
        Array[String] genoprob_global_errors = ["0.05"]
        Int repetitions = 100
        Int sample_size = 30
    }

    if (defined(genotype_dp_filter) || defined(filters_include)) {
        call utils.ApplyRandomFiltersArray {
            input:
                vcfs = vcfs,
                vcfs_SNPCall_software = vcfs_software,
                vcfs_Counts_source = vcfs_counts_source,
                vcfs_GenoCall_software = range(length(vcfs_software)),
                genotype_dp_filter    = genotype_dp_filter,
                filters_include      = filters_include,
                info_dp_sd_multiplier = info_dp_sd_multiplier
        }
    }

    Array[File] filtered_vcfs = select_first([ApplyRandomFiltersArray.vcfs_filt, vcfs])

    # Auto-detect chromosomes from the first VCF when not supplied in the dataset
    if (!defined(dataset.chromosomes)) {
        call utils.GetChromosomesFromVcf {
            input:
                vcf_file = filtered_vcfs[0]
        }
    }

    Array[String] chromosomes_to_use = select_first([dataset.chromosomes, GetChromosomesFromVcf.chromosomes])

    # Re-Genotyping with updog, supermassa and polyrad; and building maps with onemap
    scatter (idx in range(length(filtered_vcfs))) {

        call utils.SplitMarkers as splitgeno {
             input:
                vcf_file = filtered_vcfs[idx]
        }

        # Suggestion to improve performance of SuperMASSA, polyRAD and updog
        if(filter_noninfo){
            call utilsR.RemoveNonInformative {
            input:
                vcf_file = splitgeno.biallelics,
                parent1 = dataset.parent1,
                parent2 = dataset.parent2,
                replaceADbyMissing = replaceADbyMissing
            }
        }

        File vcf_up = select_first([RemoveNonInformative.vcf_filtered, splitgeno.biallelics])

        scatter (chrom in chromosomes_to_use) {
            if(run_updog) {
                call mappoly_sub.MappolyMapsEmp as updogPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "updog",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr,
                        global_errors = global_errors,
                        repetitions = repetitions,
                        sample_size = sample_size,
                        chromosome = chrom
                }
            }

            if(run_polyrad){
                call mappoly_sub.MappolyMapsEmp as polyradPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "polyrad",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr,
                        global_errors = global_errors,
                        repetitions = repetitions,
                        sample_size = sample_size,
                        chromosome = chrom
                }
            }

            if(run_supermassa){
                call mappoly_sub.MappolyMapsEmp as supermassaPolyMaps {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "supermassa",
                        CountsFrom = vcfs_counts_source[idx],
                        cross = "F1",
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr,
                        global_errors = global_errors,
                        repetitions = repetitions,
                        sample_size = sample_size,
                        chromosome = chrom
                }
            }

            if(vcfs_counts_source[idx] != "bam" && vcfs_software[idx] != "stacks" && vcfs_software[idx] != "tassel"){
                call mappoly_task.MappolyReport {
                    input:
                        vcf_file = vcf_up,
                        SNPCall_program = vcfs_software[idx],
                        GenotypeCall_program = "SNPCaller",
                        CountsFrom = vcfs_counts_source[idx],
                        parent1 = dataset.parent1,
                        parent2 = dataset.parent2,
                        max_cores = max_cores,
                        ploidy = ploidy,
                        prob_thres = prob_thres,
                        filt_segr = filt_segr,
                        global_errors = global_errors,
                        repetitions = repetitions,
                        sample_size = sample_size,
                        chromosome = chrom
                }
            }
        }
    }

    # poly callers are inside scatter(chrom), so their outputs are Array[Array[File?]?] —
    # strip outer optionals, flatten, then strip inner optionals to get a flat Array[File].
    Array[File] snpcaller_poly_flat = select_all(flatten(select_all(MappolyReport.results)))
    if (length(snpcaller_poly_flat) > 0) {
        Array[File] snpcaller_poly_results = snpcaller_poly_flat
    }
    Array[File] updog_poly_flat = select_all(flatten(select_all(updogPolyMaps.tar_gz_report)))
    if (length(updog_poly_flat) > 0) {
        Array[File] updog_poly_results = updog_poly_flat
    }
    Array[File] polyrad_poly_flat = select_all(flatten(select_all(polyradPolyMaps.tar_gz_report)))
    if (length(polyrad_poly_flat) > 0) {
        Array[File] polyrad_poly_results = polyrad_poly_flat
    }
    Array[File] supermassa_poly_flat = select_all(flatten(select_all(supermassaPolyMaps.tar_gz_report)))
    if (length(supermassa_poly_flat) > 0) {
        Array[File] supermassa_poly_results = supermassa_poly_flat
    }
    # Compress files
    call reports.JointAllReports {
        input:
            SNPCallerPolyMapsEmp = snpcaller_poly_results,
            updogPolyMaps = updog_poly_results,
            polyradPolyMaps = polyrad_poly_results,
            supermassaPolyMaps = supermassa_poly_results,
            max_cores = max_cores,
            ploidy = ploidy
    }

    output {
        File EmpiricalReads_results = JointAllReports.EmpiricalReads_results
    }
}
