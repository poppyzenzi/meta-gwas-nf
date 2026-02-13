// master.nf
nextflow.enable.dsl=2

// import subworkflows
include { COHORT_PROCESSOR } from './workflows/processor'
include { META_ANALYST } from './workflows/analyst'

workflow {
    // reference setup
    dbSNP_dir = file(params.dbSNP_dir)

    // input setup
    cohort_ch = Channel.fromPath("${projectDir}/data/test_cohorts.csv")
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
    //META_ANALYST(COHORT_PROCESSOR.out.merged_sumstats, COHORT_PROCESSOR.out.munged_sumstats, dbSNP_dir)
}
