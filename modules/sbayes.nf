// run sbayess for h2 and sbayesrc for SNP predictions

process STAGE_LD_MATRICES {
    // Use the staging queue specifically for dataStore transfers
    clusterOptions = "-q staging -l h_rt=12:00:00"

    output:
    val true, emit: ready

    script:
    """
    SOURCE="/exports/igmm/datastore/GenScotDepression/data/resources/SBayesR_matrices/ukb_50k_bigset_2.8M"
    DESTINATION="/exports/eddie/scratch/s2421111/"
    
    # Check if a key file exists to skip unnecessary rsync
    if [ ! -f "\${DESTINATION}/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.new.mldmlist" ]; then
        echo "Data missing or incomplete. Starting rsync..."
        rsync -rl \${SOURCE} \${DESTINATION}
    else
        echo "LD matrices already present in scratch. Skipping copy."
    fi
    """
}

process RUN_SBAYESS {
    tag "${trait_label}"
    publishDir "${params.resultsDir}/sbayess", mode: 'copy'

    input:
    tuple val(trait_label), path(ma_file)
    val ready // Dependency on staging process

    output:
    path("sbayess_${trait_label}.*"), emit: sbayess_results

    script:
    """
    # capture abs path of input file as when 'cd' the relative path to the work dir breaks
    GWAS_ABS=\$(readlink -f ${ma_file})
    
    # capture current work dir
    OUT_DIR=\$(pwd)

    # cd into scratch so GCTB can see relative LD paths
    cd /exports/eddie/scratch/s2421111

    # run GCTB
    # point --out back to the Nextflow directory so the files aren't lost in scratch
    ${params.gctb} \
        --sbayes S \
        --mldm ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.new.mldmlist \
        --gwas-summary \$GWAS_ABS \
        --burn-in 2000 \
        --chain-length 10000 \
        --out \$OUT_DIR/sbayess_${trait_label}

    """
}

process IMPUTE_SBAYESRC {
    tag "${label}_${ancestry}"

    input:
    tuple val(label), val(ancestry), val(ref_panel), path(ma_file)

    output:
    tuple val(label), val(ancestry), val(ref_panel), path("${label}_${ancestry}_${ref_panel}.imputed.ma"), emit: imputed_ma

    script:
    """
    ${params.gctb} \
        --ldm-eigen ${params.gctb_ldm_sbrc}/${ref_panel} \
        --gwas-summary ${ma_file} \
        --impute-summary \
        --out ${label}_${ancestry}_${ref_panel} \
        --thread ${task.cpus}
    """
}


process RUN_SBAYESRC {
    tag "${label}_${ancestry}"
    publishDir "${params.resultsDir}/sbayesrc/${label}", mode: 'copy'

    input:
    tuple val(label), val(ancestry), val(ref_panel), path(imputed_ma)

    output:
    tuple val(label), val(ancestry), val(ref_panel), path("${label}_${ancestry}_${ref_panel}_sbrc.*"), emit: sbrc_results

    script:
    """
    ${params.gctb} \
        --ldm-eigen ${params.gctb_ldm_sbrc}/${ref_panel} \
        --gwas-summary ${imputed_ma} \
        --sbayes RC \
        --annot ${params.gctb_annot} \
        --out ${label}_${ancestry}_${ref_panel}_sbrc \
        --thread ${task.cpus}
    """
}


process SCORE_PRS {
    tag "${label}_${target_cohort}_${target_ancestry}"
    publishDir "${params.resultsDir}/prs/${target_cohort}", mode: 'copy'

    input:
    tuple val(label), val(target_cohort), val(target_ancestry), path(snp_res)

    output:
    tuple val(label), val(target_cohort), val(target_ancestry), path("${label}_${target_cohort}_${target_ancestry}_sbayesrc.profile"), emit: prs_scores
    path "${label}_${target_cohort}_${target_ancestry}_sbayesrc.log"

    script:
    // resolve which bfile to use
    def bfile = (target_cohort == "ALSPAC") ? params.plink_alspac :
                (target_ancestry == "EUR")  ? params.plink_abcd_eur :
                (target_ancestry == "AMR")  ? params.plink_abcd_amr :
                (target_ancestry == "AFR")  ? params.plink_abcd_afr : null

    if (!bfile) error "No bfile found for cohort=${target_cohort} ancestry=${target_ancestry}"

    """
    plink \\
        --bfile ${bfile} \\
        --score ${snp_res} 2 5 8 header sum center \\
        --out ${label}_${target_cohort}_${target_ancestry}_sbayesrc
    """
}
