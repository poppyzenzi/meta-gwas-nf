nextflow.enable.dsl = 2

// Parameters
params.input_csv   = "${projectDir}/data/cohorts.csv"
params.cohorts_dir = '/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/cohorts'
params.outdir      = "${projectDir}/results/formatted_sumstats"
params.resultsDir  = "${projectDir}/results"
params.dbSNP_dir   = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/checks/ref_panels/split_dbSNP"
params.mrmega_path = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/MR-MEGA"
params.external_munged_dir = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/ldsc_rg/munged"
params.ld_ref_dir = "/exports/igmm/eddie/GenScotDepression/users/poppy/gsem/munging/eur_w_ld_chr"

// Processes

process FORMAT_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(ascertainment), path(raw_sumstats)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("${cohort}_${ancestry}_${build}_formatted.txt"), emit: formatted_stats

    script:
    build = (cohort =~ /NORDIC|TRAILS|RAINE|GLAD|ALLOFUS/) ? "b38" : "b37"
    """
    LC_ALL=C gzip -dc -f ${raw_sumstats} | awk -v cohort="${cohort}" '
    BEGIN {
        map["CHR"] = "CHROM CHR"; map["POS"] = "GENPOS POS BP"; map["SNP"] = "SNP ID";
        map["EA"] = "EA A1 ALLELE1"; map["NEA"] = "NEA A2 ALLELE0";
        map["EAF"] = "EAF A1FREQ FREQ";
        map["EAFCA"] = "EAFCA A1FREQ_CASES";
        map["EAFCO"] = "EAFCO A1FREQ_CONTROLS";
        map["BETA"] = "BETA LOG(OR) LOG_OR LOG-OR";
        map["OR"] = "OR";
        map["SE"] = "SE"; map["P"] = "P PVAL"; map["LOG10P"] = "LOG10P";
        map["N"] = "N N_TOTAL"; map["NCA"] = "NCA N_CASES"; map["NCO"] = "NCO N_CONTROLS";
        map["RSID"] = "RSID SNPRS RSID_HG38"; map["INFO"] = "INFO";
    }
    NR == 1 {
        for (i=1; i<=NF; i++) { header[toupper(\$i)] = i }
        c["CHR"]=f("CHR"); c["POS"]=f("POS"); c["SNP"]=f("SNP");
        c["EA"]=f("EA"); c["NEA"]=f("NEA"); c["SE"]=f("SE");
        c["BETA"]=f("BETA"); c["OR"]=f("OR"); c["P"]=f("P");
        c["LOG10P"]=f("LOG10P"); c["N"]=f("N"); c["NCA"]=f("NCA"); c["NCO"]=f("NCO");
        c["EAF"]=f("EAF"); c["EAFCA"]=f("EAFCA"); c["EAFCO"]=f("EAFCO");
        c["RSID"]=f("RSID"); c["INFO"]=f("INFO");
        print "CHR", "POS", "SNP", "RSID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N", "INFO"
        next
    }
    {
        chr=\$c["CHR"];
        gsub(/^chr/,"",chr); if(chr=="X") chr="23";
        pos=\$c["POS"]; snp=\$c["SNP"]; ea=\$c["EA"]; nea=\$c["NEA"]; se=\$c["SE"];
        
        if(c["BETA"]) beta=\$c["BETA"]; else beta=log(\$c["OR"]);
        if(c["P"]) pval=\$c["P"]; else pval=exp(-\$c["LOG10P"] * log(10));
        if(c["N"]) n=\$c["N"]; else n=\$c["NCA"] + \$c["NCO"];

        if (c["EAF"] && \$c["EAF"] != "NA" && \$c["EAF"] != "") eaf_val = \$c["EAF"];
        else if (c["EAFCO"] && \$c["EAFCO"] != "NA" && \$c["EAFCO"] != "") eaf_val = \$c["EAFCO"];
        else if (c["EAFCA"] && \$c["EAFCA"] != "NA" && \$c["EAFCA"] != "") eaf_val = \$c["EAFCA"];
        else eaf_val = "NA";

        rsid = (c["RSID"] ? \$c["RSID"] : "NA");
        info = (c["INFO"] ? \$c["INFO"] : "NA");
        
        print chr, pos, snp, rsid, ea, nea, eaf_val, beta, se, pval, n, info
    }
    function f(s) {
        split(map[s], a, " ");
        for (j in a) { if (header[a[j]]) return header[a[j]] }
        return 0
    }
    ' > ${cohort}_${ancestry}_${build}_formatted.txt
    """
}


process SPLIT_SUMSTATS {
    tag "${cohort}_${ancestry}"
    
    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(formatted_txt)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("split_chr_*.txt"), emit: chr_files

    script:
    """
    tail -n +2 ${formatted_txt} | awk '{print > "split_chr_"\$1".txt"}'
    """
}

process MATCH_DBSNP_PARALLEL {
    tag "${cohort}_${ancestry}_chr${chr_num}"

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(chr_file)
    path dbSNP_dir

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("matched_chr_${chr_num}.txt"), emit: matched_chr

    script:
    chr_num = chr_file.baseName.replaceAll(/split_chr_/, "")
    def ref_file = "${dbSNP_dir}/dbSNP_chr${chr_num}.txt"
    if ( build == "b37" )
        """
        awk 'NR==FNR { marker[\$5,\$6,\$2,\$3]=\$4; rsid[\$5,\$6,\$2,\$3]=\$1; next }
        { if ((\$1,\$2,\$5,\$6) in marker) print \$0,marker[\$1,\$2,\$5,\$6],rsid[\$1,\$2,\$5,\$6];
          else if ((\$1,\$2,\$6,\$5) in marker) print \$0,marker[\$1,\$2,\$6,\$5],rsid[\$1,\$2,\$6,\$5];
          else print \$0,"missing","missing"; }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
        """
    else
        """
        awk 'NR==FNR { marker[\$11,\$12,\$8,\$9]=\$4; rsid[\$11,\$12,\$8,\$9]=\$1; next }
        { if ((\$1,\$2,\$5,\$6) in marker) print \$0,marker[\$1,\$2,\$5,\$6],rsid[\$1,\$2,\$5,\$6];
          else if ((\$1,\$2,\$6,\$5) in marker) print \$0,marker[\$1,\$2,\$6,\$5],rsid[\$1,\$2,\$6,\$5];
          else print \$0,"missing","missing"; }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
        """
}

process MERGE_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/matched", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(matched_files)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("${cohort}_${ancestry}_matched_full.txt"), emit: merged_file

    script:
    """
    echo "CHR POS SNP RSID EA NEA EAF BETA SE P N INFO MARKER_build37 rsID_build37" > ${cohort}_${ancestry}_matched_full.txt
    cat ${matched_files} >> ${cohort}_${ancestry}_matched_full.txt
    """
}

process PREP_CORR {
    tag "${cohort}_${ancestry}"
    
    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(sumstats)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_reduced.txt")

    script:
    """
    # Column mapping for MERGE_SUMSTATS output:
    # \$13 is MARKER_build37 (CHR:POS:REF:ALT or similar)
    # \$7 is EAF
    awk 'BEGIN {FS=" "; OFS="\\t"} 
         NR==1 {next} 
         \$13 != "missing" {print \$13, \$7}' ${sumstats} | LC_ALL=C sort -k1,1 > ${cohort}_${ancestry}_reduced.txt
    """
}

process CORR_EAF {
    tag "${group_label}"
    publishDir "${params.resultsDir}/qc/eaf_correlations", mode: 'copy'

    input:
    tuple val(group_label), path(files)

    output:
    path "${group_label}_eaf_corr.csv"

    script:
    """
    export LC_ALL=C
    
    # create header and clean filenames to just the cohort_ancestry labels
    labels=\$(ls *_reduced.txt | sed 's/_reduced.txt//g')
    echo "ID \$(echo \$labels | tr '\\n' ' ')" > headers.txt

    # iterative join
    files=( *_reduced.txt )
    cp "\${files[0]}" combined_temp.txt

    for ((i=1; i<\${#files[@]}; i++)); do
        join -t\$'\\t' combined_temp.txt "\${files[i]}" > combined_new.txt
        mv combined_new.txt combined_temp.txt
    done

    # combine headers and data
    cat headers.txt combined_temp.txt > final_matrix.txt

    Rscript - <<EOF
    library(data.table)
    dt <- fread("final_matrix.txt", header = TRUE, fill = TRUE)
    
    if (nrow(dt) < 10) {
        stop(paste("Insufficient overlapping SNPs found. Only", nrow(dt), "SNPs overlap across all cohorts."))
    }
    
    # Convert to numeric, handling any potential 'NA' strings
    mat_cols <- names(dt)[-1]
    dt[, (mat_cols) := lapply(.SD, as.numeric), .SDcols = mat_cols]
    
    mat <- as.matrix(dt[, -1, with = FALSE])
    corr_matrix <- cor(mat, use = "pairwise.complete.obs")
    
    write.csv(corr_matrix, "${group_label}_eaf_corr.csv")
EOF
    """
}

process MUNGE_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/munged", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(matched_file)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_munged.sumstats.gz"), emit: munged_file

    script:
    """
    cat ${matched_file} | tr -d '\\r' | awk '\$1=\$1' > cleaned.txt
    if grep -q -w "OR" cleaned.txt; then SIGNED="OR,1"; else SIGNED="BETA,0"; fi

    python ${params.ldsc_path}/munge_sumstats.py \
        --sumstats cleaned.txt \
        --out ${cohort}_${ancestry}_munged \
        --snp rsID_build37 \
        --p P \
        --a1 EA \
        --a2 NEA \
        --signed-sumstats \$SIGNED \
        --frq EAF \
        --ignore RSID,SNP,MARKER_build37,SNPrs,INFO \
        --chunksize 50000 \
        --merge-alleles ${params.hm_snplist}
    """
}

process LDSC_H2 {
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
        --ref-ld-chr ${params.hm_refs}/ \
        --w-ld-chr ${params.hm_refs}/ \
        --out ${cohort}_${ancestry}_ldsc_h2
    """
}

process GENERATE_REPORT {
    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    path logs 

    output:
    path "LDSCintercept_summary.txt"

    script:
    """
    echo "Cohort_Ancestry: Intercept = Intercept, SD = SD, 95% CI = [Lower, Upper]" > LDSCintercept_summary.txt
    for log_file in ${logs}; do
        intercept=\$(grep "Intercept:" \$log_file | awk '{print \$2}' | head -1)
        sd=\$(grep "Intercept:" \$log_file | awk '{print \$3}' | sed 's/[()]//g' | head -1)
        echo | awk -v ival="\$intercept" -v s="\$sd" -v name="\$log_file" 'BEGIN { OFS=" " } { if (ival == "" || s == "") { printf "%s: Intercept data missing in log\\n", name } else { low = ival - (1.96 * s); high = ival + (1.96 * s); printf "%s: Intercept = %f, SD = %f, 95%% CI = [%.6f, %.6f]\\n", name, ival, s, low, high } }' >> LDSCintercept_summary.txt
    done
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
        --ref-ld-chr ${params.hm_refs}/ \
        --w-ld-chr ${params.hm_refs}/ \
        --out ${label}
    """
}

process PLOT_RG_MATRIX {
    publishDir "${params.resultsDir}/genetic_correlation", mode: 'copy'

    input:
    path "logs/*"

    output:
    path "RG_matrix_full_table.csv"
    path "RG_matrix_wide.csv"
    path "RG_matrix_heatmap.png"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(ggplot2)
    library(tidyr)

    # load log files
    log_files <- list.files("logs", pattern = "*.log", full.names = TRUE)
    results <- lapply(log_files, function(f) {
        lines <- readLines(f)
        start_idx <- grep("^p1", lines)
        if(length(start_idx) == 0) return(NULL)
        dt <- fread(text = paste(lines[start_idx:(start_idx+1)], collapse="\\n"))
        return(dt)
    })

    master_dt <- rbindlist(results, fill = TRUE)
    
    # clean names: remove suffixes and paths
    clean_name <- function(x) gsub("_munged.sumstats.gz", "", basename(x))
    master_dt[, p1 := clean_name(p1)]
    master_dt[, p2 := clean_name(p2)]

    # save long table
    fwrite(master_dt, "RG_matrix_full_table.csv")

    # build symmetric matrix
    # We need A-B, B-A, and the diagonal A-A
    ids <- unique(c(master_dt\$p1, master_dt\$p2))
    
    # create 'mirror' images
    mirror_dt <- master_dt[, .(p1 = p2, p2 = p1, rg)]
    # Create the diagonal
    diag_dt <- data.table(p1 = ids, p2 = ids, rg = 1)
    
    # combine
    plot_dt <- rbind(master_dt[, .(p1, p2, rg)], mirror_dt, diag_dt, fill = TRUE)
    
    # remove any potential dups and ensure rg is numeric
    plot_dt <- unique(plot_dt, by = c("p1", "p2"))
    plot_dt[, rg := as.numeric(rg)]

    # save wide version
    wide_dt <- dcast(plot_dt, p1 ~ p2, value.var = "rg")
    fwrite(wide_dt, "RG_matrix_wide.csv")

    # plot heatmap
    p <- ggplot(plot_dt, aes(x = p1, y = p2, fill = rg)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_gradient2(
            low = "#377eb8", 
            high = "#e41a1c", 
            mid = "white", 
            midpoint = 0, 
            limit = c(-1, 1), 
            na.value = "grey95",
            name = "rg"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            aspect.ratio = 1
        ) +
        labs(
            title = "Genetic Correlation (rg) Matrix",
            subtitle = "Grey cells indicate failed/non-significant heritability estimation",
            x = "", y = ""
        )

    # dynamic sizing based on number of cohorts
    plot_size <- max(8, length(ids) * 0.4)
    ggsave("RG_matrix_heatmap.png", plot = p, width = plot_size, height = plot_size - 1, dpi = 300)
    """
}

process QC_PLOTS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.resultsDir}/lambda", mode: 'copy', pattern: "*.lambda.txt"

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(matched_file)

    output:
    path "*.png"
    path "${cohort}_${ancestry}_lambda.txt", emit: lambda_val

    script:
    """
    #!/usr/bin/env Rscript
    library(qqman)
    library(data.table)
    library(dplyr)

    df <- fread("${matched_file}", data.table = FALSE)
    df\$P <- as.numeric(df\$P)
    df <- df[!is.na(df\$P) & df\$P > 0 & df\$P <= 1, ]
    df\$CHR <- gsub("chr", "", as.character(df\$CHR), ignore.case = TRUE)
    df\$CHR <- gsub("X", "23", df\$CHR, ignore.case = TRUE)
    df\$CHR <- as.numeric(df\$CHR)
    df <- df[!is.na(df\$CHR) & df\$CHR >= 1 & df\$CHR <= 23, ]

    if (nrow(df) > 0) {
        # cohort lambda
        chisq <- qchisq(1 - df\$P, df = 1)
        lambda <- median(chisq) / qchisq(0.5, df = 1)
        write.table(data.frame(Cohort="${cohort}", Ancestry="${ancestry}", Lambda=round(lambda, 3)), 
                    file="${cohort}_${ancestry}_lambda.txt", sep="\t", quote=FALSE, row.names=FALSE)
        # cohort QQ plot
        png("${cohort}_${ancestry}_QQ.png", width=800, height=600)
        qqman::qq(df\$P, main=paste("QQ Plot:", "${cohort}", "${ancestry}"))
        dev.off()
        # cohort manhattan plot
        png("${cohort}_${ancestry}_Manhattan.png", width=1200, height=600)
        result <- try(
            qqman::manhattan(df, chr="CHR", bp="POS", p="P", snp="SNP", 
                main=paste("Manhattan:", "${cohort}", "${ancestry}"),
                suggestiveline = 1e-05, genomewideline = 5e-08, col = c("blue4", "skyblue"))
        )
        dev.off()
    }
    """    
}

process COMBINE_LAMBDAS {
    publishDir "${params.resultsDir}", mode: 'copy'
    input: path all_lambdas 
    output: path "all_cohorts_lambda_summary.txt"
    script:
    """
    echo -e "Cohort\\tAncestry\\tLambda" > all_cohorts_lambda_summary.txt
    tail -q -n +2 ${all_lambdas} >> all_cohorts_lambda_summary.txt
    """
}

process PREP_MR_MEGA {
    tag "${cohort}_${ancestry}"

    input:
    tuple val(cohort), val(ancestry), val(ascertainment), val(build), path(matched_file)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_mrmega_ready.txt.gz"), emit: mrmega_ready

    script:
        """
        awk '
        BEGIN { FS="[[:space:]\\t]+"; OFS="\\t" }
        NR==1 {
            for (i=1; i<=NF; i++) col[\$i] = i
            eaf_idx = ("EAFCO" in col) ? col["EAFCO"] : (("EAF" in col) ? col["EAF"] : 0);
            print "MARKERNAME", "EA", "NEA", "OR", "OR_95L", "OR_95U", "EAF", "N", "CHROMOSOME", "POSITION"
            next
        }
        {
            rsid = \$(col["rsID_build37"]); marker_b37 = \$(col["MARKER_build37"])
            split(marker_b37, a, ":"); chr = a[1]; pos = a[2]
            if (match(chr, /^[0-9]+\$/) == 0 || chr < 1 || chr > 22) next
            if (match(pos, /^[0-9]+\$/) == 0 || rsid == "" || rsid == "NA" || rsid == "missing") next 
            if (col["OR"] && \$(col["OR"]) != "NA" && \$(col["OR"]) != 0) { or_val = \$(col["OR"]) }
            else if (col["BETA"] && \$(col["BETA"]) != "NA") { or_val = exp(\$(col["BETA"])) }
            else { next }
            se = \$(col["SE"]); if (se == "" || se == "NA" || se == 0) next
            or_l = or_val * exp(-1.96 * se); or_u = or_val * exp(1.96 * se)
            ea = \$(col["EA"]); nea = \$(col["NEA"]); n_size = \$(col["N"])
            eaf = (eaf_idx > 0) ? \$(eaf_idx) : "NA"
            if (ea == "NA" || nea == "NA" || eaf == "NA" || eaf == "") next
            print rsid, ea, nea, or_val, or_l, or_u, eaf, n_size, chr, pos
        }
        ' ${matched_file} | gzip > ${cohort}_${ancestry}_mrmega_ready.txt.gz
        """
}

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

process SUMMARISE_META {
    publishDir "${params.resultsDir}/meta_analysis", mode: 'copy'
    
    input:
    path "info_files/*"

    output:
    path "meta_analysis_summary.tsv"

    script:
    """
    echo -e "analysis_label\\tstudy_count\\ttotal_N" > meta_analysis_summary.tsv
    for f in info_files/*.info; do
        label=\$(basename \$f 1.txt.info)
        count=\$(grep -c "Input File" \$f)
        total_n=\$(grep "Max N" \$f | awk '{print \$NF}' || echo "N/A")
        echo -e "\$label\\t\$count\\t\$total_n" >> meta_analysis_summary.tsv
    done
    """
}

process MATCH_QC_PARALLEL {
    tag "${meta_label}_chr${chr}"

    input:
    tuple val(meta_label), path(metal_out), val(chr)
    path dbSNP_dir

    output:
    tuple val(meta_label), path("${meta_label}_chr${chr}_matched.txt"), emit: matched_chunks

    script:
    def ref_file = "${dbSNP_dir}/dbSNP_chr${chr}.txt"
    """
    awk -v c="${chr}" 'BEGIN {FS="\\t"; OFS="\\t"}
    NR==1 {print \$0; next}
    {
        # Catch X and 23 and convert to X to match dbSNP ref
        if (c == "23") {
            gsub(/^23:/, "X:", \$1);
            gsub(/^X:/, "X:", \$1);
        }
        
        # filter for current chromosome
        if (\$1 ~ "^"c":" || (c == "23" && \$1 ~ "^X:")) {
            \$1 = toupper(\$1)
            print \$0
        }
    }' ${metal_out} > temp_chr.txt

    # match with dbSNP
    awk 'BEGIN {FS="[\\t ]+"; OFS="\\t"}
    # Load Ref: Normalize any "23:" to "X:" so it matches our temp_chr.txt
    NR==FNR { 
        ref_key = toupper(\$4);
        gsub(/^23:/, "X:", ref_key);
        rsid[ref_key] = \$1; 
        next 
    }
    # Headers
    FNR==1 { print \$0, "RSID"; next }
    # Perform Match
    { 
        print \$0, ((\$1 in rsid) ? rsid[\$1] : \$1) 
    }' ${ref_file} temp_chr.txt > matched_with_x.txt

    # convert X back to 23 for final sumstats
    awk 'BEGIN {FS="\\t"; OFS="\\t"}
    {
        gsub(/^X:/, "23:", \$1);
        print \$0
    }' matched_with_x.txt > ${meta_label}_chr${chr}_matched.txt
    """
}

process FINAL_QC_REPORT {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/final_sumstats", mode: 'copy'

    input:
    tuple val(meta_label), path(matched_files)

    output:
    path "${meta_label}_QC_passed.txt.gz", emit: clean_stats
    path "${meta_label}_QC_report.txt", emit: report

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    files <- list.files(pattern = "matched.txt")
    dt <- rbindlist(lapply(files, fread))
    orig_count <- nrow(dt)

    # clean temporary METAL headers
    setnames(dt, old = names(dt), new = gsub("[[:space:].]+", "_", names(dt)))

    # numeric conversions
    num_cols <- c("P-value", "Effect", "StdErr", "Freq1", "N")
    for (col in num_cols) {
        if (col %in% names(dt)) dt[[col]] <- as.numeric(as.character(dt[[col]]))
    }

    # tidy CHR and POS
    dt[, c("CHR", "POS") := tstrsplit(MarkerName, ":", fixed=TRUE, keep=1:2)]
    dt[CHR == "X", CHR := "23"]
    dt <- dt[CHR %in% as.character(1:23)]
    dt[, CHR := as.integer(CHR)]
    dt[, POS := as.integer(POS)]

    # QC filters: unmatched, missing, invalid P, Zscore outliers MAF
    unmatched_markers <- nrow(dt[grepl("missing|no_match", MarkerName)])
    dt <- dt[!grepl("missing|no_match", MarkerName)]

    na_outliers <- nrow(dt[!complete.cases(dt[, .(`P-value`, Effect, StdErr, Freq1)])])
    dt <- dt[complete.cases(dt[, .(`P-value`, Effect, StdErr, Freq1)])]

    invalid_ps <- nrow(dt[`P-value` <= 0 | `P-value` > 1])
    dt <- dt[`P-value` > 0 & `P-value` <= 1]

    dt[, Z := Effect / StdErr]
    z_outliers <- nrow(dt[abs(Z) > 30])
    dt <- dt[abs(Z) <= 30]

    dt[, MAF := ifelse(Freq1 <= 0.5, Freq1, 1 - Freq1)]
    maf_outliers <- nrow(dt[MAF < 0.005])
    dt <- dt[MAF >= 0.005]

    # QC filter: N eff 0.5 on autosomes only
    n_col <- grep("^N", names(dt), value = TRUE)[1]
    # Calculate max_n using only autosomes (1-22) and filter
    max_n_autosomes <- dt[CHR <= 22, max(get(n_col), na.rm = TRUE)]    
    threshold <- 0.5 * max_n_autosomes
    low_n <- nrow(dt[CHR <= 22 & get(n_col) < threshold])
    # keep if (ch23) OR (sample size meets the autosome-derived threshold)
    dt <- dt[CHR == 23 | get(n_col) >= threshold]

    mapping <- c(
        "MarkerName" = "varID",
        "Allele1"    = "EA",
        "Allele2"    = "NEA",
        "Freq1"      = "EAF",
        "FreqSE"     = "EAF_SE",
        "Effect"     = "BETA",
        "StdErr"     = "SE",
        "P-value"    = "P",
        "RSID"       = "rsID"
    )

    # rename only the columns that exist
    existing_cols_to_rename <- names(mapping)[names(mapping) %in% names(dt)]
    setnames(dt, old = existing_cols_to_rename, new = mapping[existing_cols_to_rename])

    # capitalise alleles
    if ("EA" %in% names(dt)) dt[, EA := toupper(EA)]
    if ("NEA" %in% names(dt)) dt[, NEA := toupper(NEA)]

    # set col order
    desired_lead_cols <- c("varID", "rsID", "CHR", "POS", "EA", "NEA")
    final_order <- intersect(desired_lead_cols, names(dt))
    setcolorder(dt, final_order)

    report_content <- paste0("=== QC Report for: ${meta_label} ===\\n",
                             "Original SNPs: ", orig_count, "\\n",
                             "Missing/unmatched marker names removed: ", unmatched_markers, "\\n",
                             "SNPs missing statistics removed: ", na_outliers, "\\n",
                             "Invalid P-values removed: ", invalid_ps, "\\n",
                             "Z-score outliers removed: ", z_outliers, "\\n",
                             "Low MAF SNPs removed: ", maf_outliers, "\\n",
                             "Low N SNPs removed: ", low_n, "\\n",
                             "Final SNP count: ", nrow(dt), "\\n")

    writeLines(report_content, "${meta_label}_QC_report.txt")
    fwrite(dt, "${meta_label}_QC_passed.txt.gz", sep="\\t")
    """
}

process FINAL_QC_REPORT_MRMEGA {
    tag "${meta_label}"
    publishDir "${params.resultsDir}/meta_analysis/final_sumstats", mode: 'copy'

    input:
    tuple val(meta_label), path(mrmega_results)

    output:
    path "${meta_label}_QC_passed.txt.gz", emit: clean_stats
    path "${meta_label}_QC_report.txt", emit: report

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    dt <- fread("${mrmega_results}")
    orig_count <- nrow(dt)

    # standardisation mapping to catch MR-MEGA strings
    setnames(dt, "Chromosome", "CHR", skip_absent = TRUE)
    setnames(dt, "Position", "BP", skip_absent = TRUE)
    setnames(dt, "P-value_association", "P", skip_absent = TRUE)
    setnames(dt, "MarkerName", "rsID", skip_absent = TRUE)
    setnames(dt, "Nsample", "Nsample", skip_absent = TRUE)
    setnames(dt, "N_sample", "Nsample", skip_absent = TRUE)
    setnames(dt, "Ncohort", "Ncohort", skip_absent = TRUE)
    setnames(dt, "N_cohort", "Ncohort", skip_absent = TRUE)

    # clean remaining headers (removing dots/spaces)
    setnames(dt, old = names(dt), new = gsub("[[:space:].]+", "_", names(dt)))

    # tidy column format
    dt[, CHR := as.character(CHR)]
    dt[CHR == "X", CHR := "23"]
    dt <- dt[CHR %in% as.character(1:23)]
    dt[, CHR := as.integer(CHR)]
    dt[, BP := as.integer(BP)]

    num_cols <- c("P", "EAF", "chisq_association")
    for (col in num_cols) {
        if (col %in% names(dt)) dt[[col]] <- as.numeric(as.character(dt[[col]]))
    }

    # QC filters: missing, invalid P, MAF
    dt <- dt[!is.na(P) & !is.na(EAF)]
    invalid_ps <- nrow(dt[P <= 0 | P > 1])
    dt <- dt[P > 0 & P <= 1]
    dt[, MAF := ifelse(EAF <= 0.5, EAF, 1 - EAF)]
    maf_outliers <- nrow(dt[MAF < 0.005])
    dt <- dt[MAF >= 0.005]

    # uppercase alleles
    if ("EA" %in% names(dt)) dt[, EA := toupper(as.character(EA))]
    if ("NEA" %in% names(dt)) dt[, NEA := toupper(as.character(NEA))]

    # make varID
    dt[, varID := paste(CHR, BP, EA, NEA, sep=":")]

    # arrange cols
    keep_cols <- c("rsID", "CHR", "BP", "EA", "NEA", "EAF", "Nsample", 
                   "Ncohort", "Effects", "chisq_association", "P", "MAF", "varID")
    
    existing_cols <- intersect(keep_cols, names(dt))
    dt <- dt[, ..existing_cols]
    setorder(dt, CHR, BP)

    report_content <- paste0("=== MR-MEGA QC Report for: ${meta_label} ===\\n",
                             "Original SNPs: ", orig_count, "\\n",
                             "Invalid P-values removed: ", invalid_ps, "\\n",
                             "Low MAF SNPs removed: ", maf_outliers, "\\n",
                             "Final SNP count: ", nrow(dt), "\\n")

    writeLines(report_content, "${meta_label}_QC_report.txt")
    fwrite(dt, "${meta_label}_QC_passed.txt.gz", sep="\\t")
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
        --ref-ld-chr ${params.hm_refs}/ \
        --w-ld-chr ${params.hm_refs}/ \
        --out ${meta_label}_ldsc
    """
}

process GENERATE_META_LDSC_REPORT {
    publishDir "${params.resultsDir}/meta_analysis", mode: 'copy'

    input:
    path logs            // From LDSC_H2_META
    path mrmega_results  // From FINAL_QC_REPORT_MRMEGA

    output:
    path "Meta_Analysis_LDSC_Metrics.csv"

    script:
    """
    # file header
    echo "Analysis,Lambda_GC,Mean_Chi2,Intercept,Polygenicity_Percent" > Meta_Analysis_LDSC_Metrics.csv

    # process log files from LDSC
    for log_file in ${logs}; do
        [ ! -f "\$log_file" ] && continue
        label=\$(basename \$log_file _ldsc.log)
        intercept=\$(grep "Intercept:" \$log_file | awk '{print \$2}')
        mean_chi2=\$(grep "Mean Chi^2:" \$log_file | awk '{print \$3}')
        lambda=\$(grep "Lambda GC:" \$log_file | awk '{print \$3}' || echo "NA")

        poly=\$(awk -v i="\$intercept" -v m="\$mean_chi2" 'BEGIN {
            if (m > 1 && i != "") {
                p = 1 - ((i - 1) / (m - 1));
                printf "%.1f%%", p * 100
            } else {
                print "NA"
            }
        }')

        echo "\$label,\$lambda,\$mean_chi2,\$intercept,\$poly" >> Meta_Analysis_LDSC_Metrics.csv
    done

    # MR-MEGA results
    Rscript - <<EOF
    library(data.table)

    mrmega_files <- list.files(pattern = "_QC_passed.txt.gz")

    for (f in mrmega_files) {
        label <- gsub("_QC_passed.txt.gz", "", f)
        dt <- fread(f)

        p_vals <- as.numeric(dt\\\$P)
        p_vals <- p_vals[!is.na(p_vals) & p_vals > 0 & p_vals <= 1]

        chisq_values <- qchisq(1 - p_vals, df = 1)
        lambda_gc <- round(median(chisq_values) / qchisq(0.5, df = 1), 3)
        mean_chi2 <- round(mean(chisq_values), 3)

        # these are NA for MR-MEGA analysis
        intercept <- "NA"
        poly <- "NA"

        entry <- data.table(
            Analysis = label, 
            Lambda_GC = lambda_gc, 
            Mean_Chi2 = mean_chi2, 
            Intercept = intercept, 
            Polygenicity_Percent = poly
        )
        
        fwrite(entry, "Meta_Analysis_LDSC_Metrics.csv", append = TRUE, col.names = FALSE)
    }
EOF
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
        --ref-ld-chr ${params.hm_refs}/ \
        --w-ld-chr ${params.hm_refs}/ \
        --out ${meta_label}_vs_external
    """
}

process TABULATE_EXTERNAL_RG {
    publishDir "${params.resultsDir}/meta_analysis/external_rg", mode: 'copy'

    input:
    path "logs/*" // collects all logs from the LDSC_RG_EXTERNAL process

    output:
    path "ldsc_correlation_results_table.txt", emit: table
    path "master_RG_results_for_plotting.csv", emit: csv

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    library(readr)

    extract_ldsc_results <- function(log_file) {
      lines <- readLines(log_file)
      start_line <- grep("Summary of Genetic Correlation Results", lines)
      
      if (length(start_line) == 0) return(NULL)

      results_lines <- lines[(start_line + 2):length(lines)]
      
      # tryCatch as some logs might end quickly
      results_df <- tryCatch({
          read.table(text = paste(results_lines, collapse = "\n"), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
      }, error = function(e) return(NULL))

      if (is.null(results_df) || ncol(results_df) < 10) return(NULL)
      
      results_clean <- results_df %>%
        select(V1, V2, V3, V4, V6) %>%
        rename(sumstats1 = V1, sumstats2 = V2, rg = V3, se = V4, p = V6) %>%
        mutate(
          sumstats1 = basename(sumstats1),
          sumstats2 = basename(sumstats2),
          sumstats1 = str_replace(sumstats1, ".sumstats.gz", ""),
          sumstats2 = str_replace(sumstats2, ".sumstats.gz", "")
        ) %>%
        filter(!is.na(sumstats2),
               !grepl("Analysis|Total", sumstats1))

      return(results_clean)
    }

    # stage files into the 'logs/' directory
    log_files <- list.files("logs", pattern = ".log", full.names = TRUE)
    results <- bind_rows(lapply(log_files, extract_ldsc_results))

    # sumstats1 is focus trait; sumstats2 is external trait
    results <- results %>%
      mutate(
        Category = case_when(
          str_detect(sumstats1, "full_meta") ~ "Meta-analysis",
          str_detect(sumstats1, "prospective") ~ "Prospective",
          str_detect(sumstats1, "retrospective") ~ "Retrospective",
          str_detect(sumstats1, "NORDIC") ~ "Register-based",
          str_detect(sumstats1, "mdd3") ~ "MDD3",
          TRUE ~ "Other"
        )
      )

    write.table(results, "ldsc_correlation_results_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    write_csv(results, "master_RG_results_for_plotting.csv")
    """
}

// WORKFLOW

workflow {

    // reference data setup to track directory
    dbSNP_dir = file(params.dbSNP_dir)

    // locate and parse master input file > map into tuple
    cohort_ch = Channel.fromPath(params.input_csv).splitCsv(header: true)
        .map { row -> tuple(row.cohort, row.ancestry, row.ascertainment, file("${params.cohorts_dir}/${row.rel_path}")) }

    FORMAT_SUMSTATS(cohort_ch)
    SPLIT_SUMSTATS(FORMAT_SUMSTATS.out.formatted_stats)
    MATCH_DBSNP_PARALLEL(SPLIT_SUMSTATS.out.chr_files.transpose(), dbSNP_dir)
    merge_input = MATCH_DBSNP_PARALLEL.out.matched_chr.groupTuple(by: [0,1,2,3])
    MERGE_SUMSTATS(merge_input)

 
    PREP_CORR(MERGE_SUMSTATS.out.merged_file)

    // Group by ancestry (or any other logic) to run correlations
    // To compare all cohorts together:
    all_prepped_ch = PREP_CORR.out
        .map { cohort, ancestry, file -> file }
        .collect()
        .map { all_files -> tuple("Full_EAF_Comparison", all_files) }

    CORR_EAF(all_prepped_ch)

    MUNGE_SUMSTATS(MERGE_SUMSTATS.out.merged_file)
    LDSC_H2(MUNGE_SUMSTATS.out.munged_file)
    GENERATE_REPORT(LDSC_H2.out.log_file.collect())

    // collect all EUR munged files > create pairs > run pairwise LDSC jobs
    ch_munged_for_rg = MUNGE_SUMSTATS.out.munged_file
        .filter { it[1] == "EUR" }
        .map { it -> [ name: it[0], file: it[2] ] }
        .collect()

    ch_pairs = ch_munged_for_rg.flatMap { list ->
        def result = []
        // Standard double-loop to get all unique pairs
        for (int i = 0; i < list.size(); i++) {
            for (int j = i + 1; j < list.size(); j++) {
                def c1 = list[i]
                def c2 = list[j]
                // We create a simple list of 3 items
                result << ["${c1.name}_vs_${c2.name}", c1.file, c2.file]
            }
        }
        return result
    }

    LDSC_RG_PAIRWISE(ch_pairs)
    PLOT_RG_MATRIX(LDSC_RG_PAIRWISE.out.log.collect())

    QC_PLOTS(MERGE_SUMSTATS.out.merged_file)
    COMBINE_LAMBDAS(QC_PLOTS.out.lambda_val.collect())

    mrmega_prep_ch = PREP_MR_MEGA(MERGE_SUMSTATS.out.merged_file)
    mrmega_list = mrmega_prep_ch.mrmega_ready.map { c, a, f -> "input_files/${f.name}" }.collectFile(name: 'mrmega_input_list.txt', newLine: true)
    RUN_MR_MEGA(mrmega_list, mrmega_prep_ch.mrmega_ready.map{it[2]}.collect())

    // set up METAL recipes for ascertainment type and LOO
    ch_raw_data = MERGE_SUMSTATS.out.merged_file.map { row ->
        def fileObj = row.find { it.toString().endsWith('.txt') }; return [cohort: row[0], ancestry: row[1], asc: row[3], file: fileObj]
    }

    ch_metal_pairs = ch_raw_data.flatMap { d ->
        def pairs = []; def loo = ['ALSPAC', 'ABCD_US', 'UKB', 'NORDIC', 'ADDH']
        if (d.ancestry != 'EUR') pairs << tuple(d.ancestry, d.file)
        if (d.ancestry == 'EUR') {
            if (d.asc == 'prospective') pairs << tuple("EUR_prospective", d.file)
            if (d.asc == 'retrospective') pairs << tuple("EUR_retrospective", d.file)
            pairs << tuple("EUR_full_meta", d.file)
        }
        loo.each { t -> if (!d.cohort.contains(t)) { pairs << tuple("PRS_LOO_ALL_ANC_${t}", d.file); if (d.ancestry == 'EUR') pairs << tuple("EUR_LOO_${t}", d.file) } }
        return pairs
    }

    all_recipes = ch_metal_pairs.groupTuple().filter { label, files -> files.size() > 1 }
    RUN_METAL(all_recipes)
    SUMMARISE_META(RUN_METAL.out.info_logs.collect())


    // match METAL output to dbSNP ref file and QC final sumstats
    ch_matching_input = RUN_METAL.out.results
        .map { f -> tuple(f.name.replaceAll(/1\.txt\$/, ""), f) }
        .combine(Channel.of(1..23))

    MATCH_QC_PARALLEL(ch_matching_input, dbSNP_dir)

    ch_metal_for_qc = MATCH_QC_PARALLEL.out.matched_chunks
        .groupTuple(size: 23)

    FINAL_QC_REPORT(ch_metal_for_qc)

    // match MRMEGA output to dbSNP ref file and QC final sumstats
    ch_mrmega_for_qc = RUN_MR_MEGA.out.results
        .map { file_list ->
            def resultFile = file_list.find { it.name.endsWith('.result') }
            def clean_label = resultFile.name.toString()
                                    .replaceFirst(/_results\.result\$/, "")
                                    .replaceFirst(/\.result\$/, "")
            return [ clean_label, resultFile ] 
        }

    FINAL_QC_REPORT_MRMEGA(ch_mrmega_for_qc)


    // run LDSC on METAL EUR sumstats only
    ch_meta_to_munge = FINAL_QC_REPORT.out.clean_stats
        .map { file ->
            // Gets label from filename: "EUR_prospective1_QC_passed.txt.gz" -> "EUR_prospective1"
            def label = file.name.toString().replaceFirst(/_QC_passed.*/, "")
            return [ label, file ]
        }
        .filter { label, file ->
            // Strict filter: must be EUR and we ignore anything else here
            label.contains("EUR")
        }

    MUNGE_META_RESULTS(ch_meta_to_munge)
    LDSC_H2_META(MUNGE_META_RESULTS.out.munged)

    // collect all MR-MEGA QC results into a list
    ch_mrmega_for_report = FINAL_QC_REPORT_MRMEGA.out.clean_stats.collect()

    // lambda etc report
    GENERATE_META_LDSC_REPORT(
        LDSC_H2_META.out.log.collect(), 
        ch_mrmega_for_report
    )

    // LDSC RG with external traits
    // normalise meta-analysis and cohort-level channels to [label, file]
    ch_meta_munged = MUNGE_META_RESULTS.out.munged
    ch_cohort_munged = MUNGE_SUMSTATS.out.munged_file
        .map { cohort, ancestry, file -> 
            def label = "${cohort}_${ancestry}"
            return [ label, file ] 
        }

    // mix and filter for focus traits
    ch_focus_traits = ch_meta_munged
        .mix(ch_cohort_munged)
        .filter { label, file ->
            // Use logical OR (||) or regex to catch exactly the ones you want
            label.contains("EUR_full_meta") || 
            label.contains("EUR_prospective") || 
            label.contains("EUR_retrospective") || 
            label.contains("NORDIC_LOO_EUR")
        }

    // add MDD3 as external focus trait
    def mdd3_path = "${params.external_munged_dir}/pgc_mdd3_no23andMe_eur_neff.sumstats.gz"
    ch_mdd3_focus = Channel.of([ "pgc_mdd3_no23andMe_eur_neff", file(mdd3_path) ])
    
    final_focus_ch = ch_focus_traits.mix(ch_mdd3_focus)

    // run correlations and make table
    LDSC_RG_EXTERNAL(
        final_focus_ch, 
        file(params.external_munged_dir))
	TABULATE_EXTERNAL_RG(LDSC_RG_EXTERNAL.out.log.collect())


}
