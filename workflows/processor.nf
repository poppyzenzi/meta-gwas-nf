/*
========================================================================================
    PROCESSOR SUB-WORKFLOW: COHORT-LEVEL QC & MUNGING
========================================================================================
*/

// import modules
include { FORMAT_SUMSTATS; PREP_CORR } from '../modules/format'
include { SPLIT_SUMSTATS; MATCH_DBSNP_COHORT; MERGE_SUMSTATS } from '../modules/matching'
include { MUNGE_SUMSTATS; LDSC_H2_COHORT; LDSC_RG_PAIRWISE } from '../modules/ldsc'
include { CORR_EAF; SUMMARISE_COHORT_LDSC; QC_PLOTS; COMBINE_LAMBDAS; PLOT_RG_MATRIX } from '../modules/reporting'

// workflows/processor.nf

workflow COHORT_PROCESSOR {
    take:
        cohort_ch
        dbSNP_dir

    main:
        
        FORMAT_SUMSTATS(cohort_ch)
        SPLIT_SUMSTATS(FORMAT_SUMSTATS.out.formatted_stats)
        
        // transpose sends chrs separately
        MATCH_DBSNP_COHORT(SPLIT_SUMSTATS.out.chr_files.transpose(), dbSNP_dir)
        
        // regroup by [cohort, ancestry, build, ascertainment]
        merge_input = MATCH_DBSNP_COHORT.out.matched_chr.groupTuple(by: [0,1,2,3])
        MERGE_SUMSTATS(merge_input)

        // EAF correlation checks
        PREP_CORR(MERGE_SUMSTATS.out.merged_file)
        all_prepped_ch = PREP_CORR.out
            .map { cohort, ancestry, file -> file }
            .collect()
            .map { all_files -> tuple("Full_EAF_Comparison", all_files) }
        
        CORR_EAF(all_prepped_ch)

        // cohort munge and h2
        MUNGE_SUMSTATS(MERGE_SUMSTATS.out.merged_file)
        LDSC_H2_COHORT(MUNGE_SUMSTATS.out.munged_file)
        
        // collect log files for a summary report
        SUMMARISE_COHORT_LDSC(LDSC_H2_COHORT.out.log_file.collect())

        // pairwise rG between EUR cohorts
        ch_munged_for_rg = MUNGE_SUMSTATS.out.munged_file
            .filter { it[1] == "EUR" }
            .map { it -> [ name: it[0], file: it[2] ] }
            .collect()

        ch_pairs = ch_munged_for_rg.flatMap { list ->
            def result = []
            for (int i = 0; i < list.size(); i++) {
                for (int j = i + 1; j < list.size(); j++) {
                    def c1 = list[i]
                    def c2 = list[j]
                    result << ["${c1.name}_vs_${c2.name}", c1.file, c2.file]
                }
            }
            return result
        }

        LDSC_RG_PAIRWISE(ch_pairs)
        PLOT_RG_MATRIX(LDSC_RG_PAIRWISE.out.log.collect())

        // cohort QQ, manhattan and lambda
        QC_PLOTS(MERGE_SUMSTATS.out.merged_file)
        COMBINE_LAMBDAS(QC_PLOTS.out.lambda_val.collect())

    emit:
        merged_sumstats = MERGE_SUMSTATS.out.merged_file
        munged_sumstats = MUNGE_SUMSTATS.out.munged_file
}

