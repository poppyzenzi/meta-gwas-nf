// meta-analysis processes (metal and mr-mega)

process RUN_MR_MEGA {
    publishDir "${params.resultsDir}/meta_analysis/mrmega", mode: 'copy'
    
    input:
    path file_list
    path "input_files/*"

    output:
    path "cross_ancestry_results*", emit: results

    script:
    """
    ${params.mrmega_path}/MR-MEGA --pc 4 --filelist ${file_list} --out cross_ancestry_results
    """
}

process RUN_METAL {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/metal/${meta_label}", mode: 'copy'

    input:
    tuple val(meta_label), path(files)

    output:
    path "${meta_label}1.txt", emit: results
    path "${meta_label}1.txt.info", emit: info_logs

    script:
    """
    . /etc/profile.d/modules.sh
    module load igmm/apps/metal/2020-05-05
    
    metal << EOF
    SCHEME STDERR
    MARKER MARKER_build37
    ALLELE EA NEA
    WEIGHT N
    EFFECT BETA
    STDERRLABEL SE
    PVALUELABEL P
    FREQLABEL EAF
    AVERAGEFREQ ON
    CUSTOMVARIABLE N
    LABEL N as N
    ${files.collect { "PROCESS $it" }.join('\n    ')}
    OUTFILE ${meta_label} .txt
    ANALYZE
    QUIT
    EOF
    """
}

