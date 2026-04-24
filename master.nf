// master.nf
nextflow.enable.dsl=2

// import subworkflows
include { COHORT_PROCESSOR } from './workflows/processor'
include { META_ANALYST } from './workflows/analyst'
include { DOWNSTREAM_ANALYSIS } from './workflows/downstream'
include { STAGE_LD_MATRICES } from './modules/sbayes'

workflow {
    // reference setup
    dbSNP_dir = file(params.dbSNP_dir)

    // input setup
    cohort_ch = Channel.fromPath("${projectDir}/data/cohorts.csv")
    .splitCsv(header: true)
    .map { row -> 
        // assing build by cohort name
        def build = (row.cohort =~ /NORDIC|TRAILS|RAINE|GLAD|ALLOFUS/) ? "b38" : "b37"
        
        return tuple(
            row.cohort, 
            row.ancestry, 
            build, 
            row.ascertainment, 
            file("${params.cohorts_dir}/${row.rel_path}")
        ) 
    }

    // execute stage 1
    COHORT_PROCESSOR(cohort_ch, dbSNP_dir)

    // execute stage 2 using outputs from stage 1
    META_ANALYST(COHORT_PROCESSOR.out.merged_sumstats, COHORT_PROCESSOR.out.munged_sumstats, dbSNP_dir)

    // Define specific patterns for downstream analyses
    def target_patterns = ~/^(EUR_full_meta|EUR_prospective|EUR_retrospective|EUR_LOO_|ALL_LOO_)/

    // prep LD matrices (stage from datastore to scratch)
    STAGE_LD_MATRICES()

    // filter the channel for downstream sumstats
    ch_for_downstream = META_ANALYST.out.clean_stats
        .filter { file -> file.name =~ target_patterns }
        .map { file ->
            def label = file.name.toString().replaceFirst(/1\.txt_QC_passed.*/, "")
            return [ label, file ]
        }

    // execute stage 3: pass both the data and staged LD matrices into the subworkflow
    DOWNSTREAM_ANALYSIS(ch_for_downstream, STAGE_LD_MATRICES.out.ready, META_ANALYST.out.noneur_h2_logs)
    
}
