/*
========================================================================================
    DOWNSTREAM SUB-WORKFLOW: GCTB / SBAYES
========================================================================================
*/

include { FORMAT_MA } from '../modules/format'
include { RUN_SBAYESS; IMPUTE_SBAYESRC; RUN_SBAYESRC; SCORE_PRS } from '../modules/sbayes' 
include { REPORT_LIABILITY_H2; REPORT_LIABILITY_H2_NONEUR; MERGE_LIABILITY_H2 } from '../modules/reporting'


workflow DOWNSTREAM_ANALYSIS {
    take:
        ma_ready_ch // tuple(label, file)
        ld_matrices // the ready flag/path from STAGE_LD_MATRICES
	noneur_h2_logs 

    main:
        // format to .ma file and run SBayesS for h2
        FORMAT_MA(ma_ready_ch)

	// filter cohorts first
        ch_sbayes_input = FORMAT_MA.out
            .filter { label, file -> label =~ /full_meta|prospective|retrospective/ }

        RUN_SBAYESS(ch_sbayes_input, ld_matrices)

        // collect all parRes files and send to the reporter (EUR sbayesS)
        REPORT_LIABILITY_H2(RUN_SBAYESS.out.sbayess_results.collect())

        // non-EUR LDSC liability h2
        REPORT_LIABILITY_H2_NONEUR(noneur_h2_logs.collect())

        // merge EUR sbayesS and non-EUR LDSC liability h2 into one file
        MERGE_LIABILITY_H2(
            REPORT_LIABILITY_H2.out.h2_report
                .mix(REPORT_LIABILITY_H2_NONEUR.out.h2_report)
                .collect()
        )


	// SBayesRC for PRS SNP predictions
	// Reference panel mapping: ancestry -> LD panel name
        def ref_map = [
            EUR: "ukbEUR_HM3",
            AFR: "ukbAFR_HM3",
            EAS: "ukbEAS_HM3"
        ]

        // Map: label pattern -> target ancestries for PRS
        // These match the FORMAT_MA output labels e.g. EUR_LOO_ABCD_US, LOO_ALL_ANC_ALSPAC
        def sbrc_targets = [
            [~/EUR_LOO_ABCD/,       ["EUR", "AFR"]],
            [~/EUR_LOO_ALSPAC/,     ["EUR"]],
            [~/EUR_LOO_ADDH/,       ["EUR", "AFR", "EAS"]],
            [~/ALL_LOO_ABCD/,   ["EUR", "AFR"]],
            [~/ALL_LOO_ALSPAC/, ["EUR"]],
            [~/ALL_LOO_ADDH/,   ["EUR", "AFR", "EAS"]]
        ]

        // Filter to SBayesRC-relevant labels, then fan out by ancestry
        ch_sbrc_input = FORMAT_MA.out.ma_file
            .filter { label, file ->
                sbrc_targets.any { pattern, _ -> label =~ pattern }
            }
            .flatMap { label, file ->
                def match  = sbrc_targets.find { pattern, _ -> label =~ pattern }
                def ancestries = match ? match[1] : []
                ancestries.collect { anc ->
                    tuple(label, anc, ref_map[anc], file)
                }
            }

        // Step 1: QC + impute
        IMPUTE_SBAYESRC(ch_sbrc_input)

        // Step 2: run SBayesRC → SNP weights
        RUN_SBAYESRC(IMPUTE_SBAYESRC.out.imputed_ma)


	// PRS scores plink
        // AMR uses EUR SNP predictions as proxy (no AMR LD panel available)
        // Target scoring map: [label_pattern, cohort, [scoring ancestries]]
        // scoring ancestry != LD ancestry for AMR
        def prs_scoring_map = [
            [~/EUR_LOO_ABCD/,       "ABCD",   ["EUR", "AMR", "AFR"]],
            [~/ALL_LOO_ABCD/,   "ABCD",   ["EUR", "AMR", "AFR"]],
            [~/EUR_LOO_ALSPAC/,     "ALSPAC", ["EUR"]],
            [~/ALL_LOO_ALSPAC/, "ALSPAC", ["EUR"]]
        ]

        // LD panel used for AMR scoring is EUR (proxy) — map scoring anc -> snpRes anc
        def scoring_to_ld_anc = [EUR: "EUR", AMR: "EUR", AFR: "AFR"]

        ch_prs_input = RUN_SBAYESRC.out.sbrc_results
            // only keep .snpRes files
            .map { label, ld_anc, ref_panel, files ->
                def snp_res = files instanceof List
                    ? files.find { it.name.endsWith('.snpRes') }
                    : (files.name.endsWith('.snpRes') ? files : null)
                snp_res ? tuple(label, ld_anc, ref_panel, snp_res) : null
            }
            .filter { it != null }
            // expand into scoring jobs
            .flatMap { label, ld_anc, ref_panel, snp_res ->
                def match = prs_scoring_map.find { pat, _, __ -> label =~ pat }
                if (!match) return []
                def (pat, cohort, scoring_ancs) = match
                // only emit if this snpRes ld_anc matches what we need for scoring anc
                scoring_ancs
                    .findAll { scoring_to_ld_anc[it] == ld_anc }
                    .collect { scoring_anc -> tuple(label, cohort, scoring_anc, snp_res) }
            }

        SCORE_PRS(ch_prs_input)

}
