/*
========================================================================================
    ANALYST SUB-WORKFLOW: META-ANALYSIS & POST-HOC
========================================================================================
*/

// import modules
include { RUN_METAL; RUN_MR_MEGA } from '../modules/meta_analysis'
include { MATCH_DBSNP_META } from '../modules/matching'
include { PREP_MR_MEGA; FINAL_QC_REPORT_METAL; FINAL_QC_REPORT_MRMEGA } from '../modules/format'
include { MUNGE_META_RESULTS; LDSC_H2_META; LDSC_RG_EXTERNAL } from '../modules/ldsc'
include { SUMMARISE_META; TABULATE_EXTERNAL_RG; GENERATE_META_LDSC_REPORT } from '../modules/reporting'

// workflows/analyst.nf

workflow META_ANALYST {
    take:
        clean_sumstats_ch // from Processor: [cohort, ancestry, build, asc, file]
        munged_cohort_ch  // from Processor: [cohort, ancestry, file]
        dbSNP_dir

    main:
        // MR-MEGA branch
        PREP_MR_MEGA(clean_sumstats_ch)
        
        mrmega_list = PREP_MR_MEGA.out.mrmega_ready
            .map { c, a, b, asc, f -> "input_files/${f.name}" }
            .collectFile(name: 'mrmega_input_list.txt', newLine: true)

        RUN_MR_MEGA(mrmega_list, PREP_MR_MEGA.out.mrmega_ready.map{it[4]}.collect())

        // logic for MR-MEGA cleanup
        ch_mrmega_for_qc = RUN_MR_MEGA.out.results
            .map { file_list ->
                def resultFile = file_list.find { it.name.endsWith('.result') }
                def clean_label = resultFile.name.toString()
                                        .replaceFirst(/_results\.result\$/, "")
                                        .replaceFirst(/\.result\$/, "")
                return [ clean_label, resultFile ]
            }

        FINAL_QC_REPORT_MRMEGA(ch_mrmega_for_qc)

        // Metal branch
        ch_metal_pairs = clean_sumstats_ch.flatMap { cohort, ancestry, build, asc, file ->
            def pairs = []; def loo = ['ALSPAC', 'ABCD_US', 'UKB', 'NORDIC', 'ADDH']
            if (ancestry != 'EUR') pairs << tuple(ancestry, file)
            if (ancestry == 'EUR') {
                if (asc == 'prospective') pairs << tuple("EUR_prospective", file)
                if (asc == 'retrospective') pairs << tuple("EUR_retrospective", file)
                pairs << tuple("EUR_full_meta", file)
            }
            loo.each { t -> 
                if (!cohort.contains(t)) { 
                    pairs << tuple("PRS_LOO_ALL_ANC_${t}", file)
                    if (ancestry == 'EUR') pairs << tuple("EUR_LOO_${t}", file) 
                } 
            }
            return pairs
        }

        all_recipes = ch_metal_pairs.groupTuple().filter { label, files -> files.size() > 1 }
        RUN_METAL(all_recipes)

        // post MA QC and ref matching
        ch_metal_match_input = RUN_METAL.out.results
            .map { f -> tuple(f.name.replaceAll(/1\.txt\$/, ""), f) }
            .combine(Channel.of(1..23))
        
        MATCH_DBSNP_META(ch_metal_match_input, dbSNP_dir)
        
        // wait for all 23 chrs per meta-analysis
        FINAL_QC_REPORT_METAL(MATCH_DBSNP_META.out.matched_chunks.groupTuple(size: 23))

        // LDSC and RG for meta-analyses
        ch_meta_to_munge = FINAL_QC_REPORT_METAL.out.clean_stats
            .map { file ->
                def label = file.name.toString().replaceFirst(/_QC_passed.*/, "")
                return [ label, file ]
            }
            .filter { label, file -> label.contains("EUR") }

        MUNGE_META_RESULTS(ch_meta_to_munge)

        LDSC_H2_META(MUNGE_META_RESULTS.out.munged)

        // prep focus traits for external RG
        ch_meta_munged = MUNGE_META_RESULTS.out.munged
        ch_cohort_munged = munged_cohort_ch.map { c, a, f -> ["${c}_${a}", f] }

        // external MDD3 trait
        def mdd3_path = "${params.external_munged_dir}/pgc_mdd3_no23andMe_eur_neff.sumstats.gz"
        ch_mdd3_focus = Channel.of([ "pgc_mdd3_no23andMe_eur_neff", file(mdd3_path) ])

        final_focus_ch = ch_meta_munged
            .mix(ch_cohort_munged)
            .filter { label, file ->
                label.contains("EUR_full_meta") ||
                label.contains("EUR_prospective") ||
                label.contains("EUR_retrospective") ||
                label.contains("NORDIC_LOO_EUR")
            }
            .mix(ch_mdd3_focus)

        LDSC_RG_EXTERNAL(final_focus_ch, file(params.external_munged_dir))
        
        // meta-analyses and genetic correlation reports
        SUMMARISE_META(RUN_METAL.out.info_logs.collect())
        
        GENERATE_META_LDSC_REPORT(
            LDSC_H2_META.out.log.collect(),
            FINAL_QC_REPORT_MRMEGA.out.clean_stats.collect()
        )
        
        TABULATE_EXTERNAL_RG(LDSC_RG_EXTERNAL.out.log.collect())
}
