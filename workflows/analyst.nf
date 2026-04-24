/*
========================================================================================
    ANALYST SUB-WORKFLOW: META-ANALYSIS & POST-HOC
========================================================================================
*/

// import modules
include { RUN_METAL; RUN_MR_MEGA } from '../modules/meta_analysis'
include { MATCH_DBSNP_META } from '../modules/matching'
include { PREP_MR_MEGA; FINAL_QC_REPORT_METAL; FINAL_QC_REPORT_MRMEGA } from '../modules/format'
include { MUNGE_META_RESULTS; LDSC_H2_META; LDSC_H2_META_NONEUR; LDSC_RG_EXTERNAL; LDSC_RG_ASCERTAINMENT } from '../modules/ldsc'
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
            def pairs = []; def loo = ['ALSPAC', 'ABCD_US', 'UKB', 'NORDIC_LOO', 'ADDH', 'JANSSEN']
            if (ancestry != 'EUR') pairs << tuple(ancestry, file)
            if (ancestry == 'EUR') {
                if (asc == 'prospective') pairs << tuple("EUR_prospective", file)
                if (asc == 'retrospective') pairs << tuple("EUR_retrospective", file)
                pairs << tuple("EUR_full_meta", file)
            }
            loo.each { t -> 
                if (!cohort.contains(t)) { 
                    pairs << tuple("ALL_LOO_${t}", file)
                    if (ancestry == 'EUR') pairs << tuple("EUR_LOO_${t}", file) 
                } 
            }
            return pairs
        }

        all_recipes = ch_metal_pairs.groupTuple().filter { label, files -> files.size() > 1 }
        RUN_METAL(all_recipes)

        // post MA QC and ref matching
	ch_metal_match_input = RUN_METAL.out.results
    		.map { f -> tuple(f.simpleName.replaceAll(/1$/, ""), f) }
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
            .filter { label, file ->
		label.contains("EUR") || label.contains("AFR") ||
		label.contains("AMR") || label.contains("EAS") || label.contains("ALL_LOO")
	    }

        MUNGE_META_RESULTS(ch_meta_to_munge)

	// EUR h2 (chr split LD scores)
        LDSC_H2_META(MUNGE_META_RESULTS.out.munged)

	// non-EUR h2 (single LD file per ancestry)
        ch_noneur_h2 = MUNGE_META_RESULTS.out.munged
            .filter { label, f -> !label.contains("EUR") }
            .map { label, f ->
                def ld_code = label.contains("AMR") ? "AMR" :
                            label.contains("EAS") ? "EAS" :
                            label.contains("SAS") ? "CSA" : "AFR"
                return [ label, f, ld_code ]
            }

        LDSC_H2_META_NONEUR(ch_noneur_h2)
       
 
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

        // rG between ascertainment types
        // Pull the two EUR ascertainment types from meta munged channel
        ch_eur_types = MUNGE_META_RESULTS.out.munged
            .filter { label, file -> label in ["EUR_prospective", "EUR_retrospective"] }

        // Pull the register-based Nordic file from the cohort-level munged directory
        ch_nordic = Channel.fromPath("${params.resultsDir}/munged/NORDIC_LOO_EUR_munged.sumstats.gz")
            .map { file -> tuple("NORDIC_LOO_EUR", file) }

        // Combine into one channel and get all pairs
        ch_all_ascertainment = ch_eur_types.mix(ch_nordic)

        // Generate all pairwise combinations
        ch_rg_pairs = ch_all_ascertainment
            .combine(ch_all_ascertainment)
            .filter { label_a, file_a, label_b, file_b -> label_a < label_b } // avoid duplicates and self-pairs

        LDSC_RG_ASCERTAINMENT(ch_rg_pairs)


    emit:
        clean_stats  = FINAL_QC_REPORT_METAL.out.clean_stats
        mrmega_stats = FINAL_QC_REPORT_MRMEGA.out.clean_stats
	noneur_h2_logs = LDSC_H2_META_NONEUR.out.log

}
