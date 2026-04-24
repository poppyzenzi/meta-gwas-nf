# GWAS meta-analysis of adolescent-onset depression

This repository contains a Nextflow pipeline to QC and meta-analyse multi-ancestry GWAS sumstats. The pipeline is set up to perform:
- cross-ancestry meta-analysis with `MR-MEGA`
- ancestry-stratified and phenotype ascertainment-stratiefied meta-analyses with `METAL`
- a series of leave-one-out meta-analyses for downstream methods like PRS and MR

The pipeline is configured to run on the `Eddie` HPCC system which uses `SGE` job scheduling. After starting an interactive session in `tmux`, load the Nexflow module: `module load roslin/nextflow/25.10.2`

The pipeline is controlled by the `master.nf` workflow which calls two subworkflows:
1. `processor.nf` which performs all pre-meta-analysis QC formatting and diagnostics
2. `analyst.nf` which performs the meta-analyses and post-hoc QC

The processes are arranged into the following modules:
```
// modules/format.nf
process FORMAT_SUMSTATS {} // format initial sumstats
process PREP_CORR {} // prep for EAF corr
process PREP_MR_MEGA {} // prep for MR-MEGA
process FINAL_QC_REPORT_METAL {} //final metal sumstats QC
process FINAL_QC_REPORT_MRMEGA {} //final mrmega sumstats QC
process FORMAT_MA{} // .ma format for gctb

// modules/matching.nf
process SPLIT_SUMSTATS {} // split by chr for faster matching
process MATCH_DBSNP_COHORT {} // match by chr to dbsnp ref
process MERGE_SUMSTATS {} // merge back together
process MATCH_DBSNP_META {} // rematch the meta-analysed sumstats with dbSNP ref

// modules/reporting.nf
process CORR_EAF {} // run EAF corr check
process SUMMARISE_COHORT_LDSC {} // cohort LDSC checks
process PLOT_RG_MATRIX {} // rgs between cohorts
process QC_PLOTS {} // cohort QQ and manhattan plots
process COMBINE_LAMBDAS {} // cohort lambda report
process SUMMARISE_META {} // meta-analysis METAL report
process GENERATE_META_LDSC_REPORT {} // ldsc results for meta analyses
process TABULATE_EXTERNAL_RG {} // table for rG external trait results
process REPORT_LIABILITY_H2 {} // liability scale sbayess h2
process REPORT_LIABILITY_H2_NONEUR {} // liability scale ldsc h2
process MERGE_LIABILITY_H2 {}

// modules/ldsc.nf
process MUNGE_SUMSTATS {} // munge cohort level sumstats
process LDSC_H2_COHORT {} // run ldsc h2 on cohort level
process LDSC_RG_PAIRWISE {} // run rg between cohorts
process MUNGE_META_RESULTS {} // munge meta analysis results
process LDSC_H2_META {}
process LDSC_H2_META_NONEU R{} 
process LDSC_RG_EXTERNAL {} // rG with external traits
process LDSC_RG_ASCERTAINMENT {} // rG btw ascertainment types

// modules/meta_analysis.nf
process RUN_MR_MEGA {}
process RUN_METAL {}

// modules.sbayesnf
process STAGE_LD_MATRICES {}
process RUN_SBAYESS {}
process IMPUTE_SBAYESRC {}
process RUN_SBAYESRC {}
process SCORE_PRS {}

```
