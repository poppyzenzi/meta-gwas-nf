// all LDSC-related processes

process MUNGE_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/munged", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(matched_file)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_munged.sumstats.gz"), emit: munged_file

script:
    """
    # 1. Process the file
    # We use a math-based check for N: adding 0 to a string in awk 
    # results in 0 if it is not a number.
    tr -d '\\r' < ${matched_file} | awk '
        NR > 1 {
            # Check if column 11 is numeric and greater than 0
            # and check that the last column is not missing
            is_numeric = (\$11 + 0 == \$11)
            if (is_numeric && \$11 > 0 && \$NF != "missing") {
                print \$0
            }
        }' > filtered.txt
    
    # 2. Re-attach the header
    head -n 1 ${matched_file} > cleaned.txt
    cat filtered.txt >> cleaned.txt
    
    # 3. Determine if OR or BETA
    if grep -q -w "OR" cleaned.txt; then SIGNED="OR,1"; else SIGNED="BETA,0"; fi
    
    # 4. Run LDSC
    python ${params.ldsc_path}/munge_sumstats.py \
        --sumstats cleaned.txt \
        --out ${cohort}_${ancestry}_munged \
        --snp rsID_build37 \
        --p P \
        --a1 EA \
        --a2 NEA \
        --N-col N \
        --signed-sumstats \$SIGNED \
        --frq EAF \
        --ignore RSID,SNP,MARKER_build37,SNPrs,INFO \
        --chunksize 50000 \
        --merge-alleles ${params.hm_snplist}
    """

}

process LDSC_H2_COHORT {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/ldsc_results", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), path(munged_file)

    output:
    path "${cohort}_${ancestry}_ldsc_h2.log", emit: log_file

    script:
    """
    python ${params.ldsc_path}/ldsc.py \
        --h2 ${munged_file} \
        --ref-ld-chr ${params.ld_ref_dir}/ \
        --w-ld-chr ${params.ld_ref_dir}/ \
        --out ${cohort}_${ancestry}_ldsc_h2
    """
}

process LDSC_RG_PAIRWISE {
    tag "${label}"
    errorStrategy 'ignore' // Skip pairs that fail (e.g. low SNP overlap)

    input:
    tuple val(label), path(file1), path(file2)

    output:
    path "${label}.log", emit: log

    script:
    """
    python ${params.ldsc_path}/ldsc.py \
        --rg ${file1},${file2} \
        --ref-ld-chr ${params.ld_ref_dir}/ \
        --w-ld-chr ${params.ld_ref_dir}/ \
        --out ${label}
    """
}

process MUNGE_META_RESULTS {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/munged", mode: 'copy'

    input:
    tuple val(meta_label), path(sumstats)

    output:
    tuple val(meta_label), path("${meta_label}_munged.sumstats.gz"), emit: munged

    script:
    """
    python ${params.ldsc_path}/munge_sumstats.py \
        --sumstats ${sumstats} \
        --out ${meta_label}_munged \
        --snp rsID \
        --p P \
        --a1 EA \
        --a2 NEA \
        --signed-sumstats BETA,0 \
        --frq EAF \
        --N-col N \
	--ignore varID,MAF \
        --merge-alleles ${params.hm_snplist} \
        --chunksize 50000
    """
}

process LDSC_H2_META {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/ldsc", mode: 'copy'

    input:
    tuple val(meta_label), path(munged_file)

    output:
    path "${meta_label}_ldsc.log", emit: log

    script:
    """
    python ${params.ldsc_path}/ldsc.py \
        --h2 ${munged_file} \
        --ref-ld-chr ${params.ld_ref_dir}/ \
        --w-ld-chr ${params.ld_ref_dir}/ \
        --out ${meta_label}_ldsc
    """
}

process LDSC_RG_EXTERNAL {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/external_rg", mode: 'copy'

    input:
    tuple val(meta_label), path(focus_sumstats)
    path external_dir // where external trait sumstats are saved

    output:
    path "${meta_label}_vs_external.log", emit: log

    script:
    """
    EXT_FILES=\$(ls ${external_dir}/*.sumstats.gz | grep -v "${focus_sumstats.name}" | paste -sd "," -)

    python ${params.ldsc_path}/ldsc.py \
        --rg ${focus_sumstats},\${EXT_FILES} \
        --ref-ld-chr ${params.ld_ref_dir}/ \
        --w-ld-chr ${params.ld_ref_dir}/ \
        --out ${meta_label}_vs_external
    """
}

