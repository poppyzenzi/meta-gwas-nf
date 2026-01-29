nextflow.enable.dsl = 2

// --- Parameters ---
params.input_csv   = "${projectDir}/data/cohorts.csv"
params.cohorts_dir = '/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/cohorts'
params.outdir      = "${projectDir}/results/formatted_sumstats"
params.resultsDir  = "${projectDir}/results"
params.dbSNP_dir   = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/checks/ref_panels/split_dbSNP"
params.mrmega_path = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/MR-MEGA"


// --- Processes ---

process FORMAT_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(ascertainment), path(raw_sumstats)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("${cohort}_${ancestry}_${build}_formatted.txt"), emit: formatted_stats

    script:
    build = (cohort =~ /NORDIC|TRAILS|RAINE|GLAD/) ? "b38" : "b37"
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
        map["RSID"] = "RSID SNPRS"; map["INFO"] = "INFO";
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
    memory '8 GB'

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

    # 1. Load data
    df <- fread("${matched_file}", data.table = FALSE)
    
    # 2. Clean and Format (Crucial for Manhattan plots)
    df\$P <- as.numeric(df\$P)
    df <- df[!is.na(df\$P) & df\$P > 0 & df\$P <= 1, ]

    # Convert CHR to numeric (Remove 'chr' prefix if present, change X to 23)
    df\$CHR <- gsub("chr", "", as.character(df\$CHR), ignore.case = TRUE)
    df\$CHR <- gsub("X", "23", df\$CHR, ignore.case = TRUE)
    df\$CHR <- as.numeric(df\$CHR)
    
    # Filter to keep only standard chromosomes to avoid plotting errors
    df <- df[!is.na(df\$CHR) & df\$CHR >= 1 & df\$CHR <= 23, ]

    if (nrow(df) > 0) {
        # 3. Lambda Calculation
        chisq <- qchisq(1 - df\$P, df = 1)
        lambda <- median(chisq) / qchisq(0.5, df = 1)
        write.table(data.frame(Cohort="${cohort}", Ancestry="${ancestry}", Lambda=round(lambda, 3)), 
                    file="${cohort}_${ancestry}_lambda.txt", sep="\t", quote=FALSE, row.names=FALSE)

        # 4. QQ Plot
        png("${cohort}_${ancestry}_QQ.png", width=800, height=600)
        qqman::qq(df\$P, main=paste("QQ Plot:", "${cohort}", "${ancestry}"))
        dev.off()

        # 5. Manhattan Plot (Using qqman)
        # We use a try() block to catch errors so the rest of the pipeline can finish
        png("${cohort}_${ancestry}_Manhattan.png", width=1200, height=600)
        result <- try(
            qqman::manhattan(df, 
                chr="CHR", 
                bp="POS", 
                p="P", 
                snp="SNP", 
                main=paste("Manhattan:", "${cohort}", "${ancestry}"),
                suggestiveline = 1e-05, 
                genomewideline = 5e-08,
                col = c("blue4", "skyblue")
            )
        )
        if(class(result) == "try-error") print("Manhattan plot failed - check log")
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
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_mrmega_ready.txt"), emit: mrmega_ready

    script:
    """
    awk '
    BEGIN { FS="[[:space:]\\t]+"; OFS="\\t" }

    NR==1 {
        for (i=1; i<=NF; i++) col[\$i] = i
        print "MARKERNAME", "MARKER_build37", "EA", "NEA", "EAFCO", "N", "OR", "SE", "CHR", "POS", "OR_95L", "OR_95U"
        next
    }

    {
        rsid = \$(col["rsID_build37"])
        # 1. Stricter RSID filtering: Skip "no_rsid" and "no_match..."
        if (rsid == "" || rsid == "NA" || rsid == "no_rsid" || rsid == "no_match_in_build37") next
        
        marker_b37 = \$(col["MARKER_build37"])
        split(marker_b37, a, ":")
        chr = a[1]; pos = a[2]

        # 2. Skip non-autosomal chromosomes (X, Y, MT) to avoid "mismatch" errors
        if (chr == "X" || chr == "Y" || chr == "MT" || chr == "23") next

        ea = \$(col["EA"]); nea = \$(col["NEA"])
        eaf = \$(col["EAF"]); n_size = \$(col["N"])
        beta = \$(col["BETA"]); se = \$(col["SE"])
        
        or_val = exp(beta)
        or_l   = or_val * exp(-1.96 * se)
        or_u   = or_val * exp(1.96 * se)
        
        print rsid, marker_b37, ea, nea, eaf, n_size, or_val, se, chr, pos, or_l, or_u
    }
    ' ${matched_file} > ${cohort}_${ancestry}_mrmega_ready.txt
    """
}


process RUN_MR_MEGA {
    publishDir "${params.resultsDir}/meta_analysis/mrmega", mode: 'copy'
    
    input:
    path file_list
    path "input_files/*" // Collects all files into a subfolder for MR-MEGA to find

    output:
    path "cross_ancestry_results*", emit: results

    script:
    """
    ${params.mrmega_path}/MR-MEGA --pc 4 \
        --filelist ${file_list} \
        --name_eaf EAFCO \
        --name_chr CHR \
        --name_pos POS \
        --out cross_ancestry_results
    """
}



workflow {
    dbSNP_dir = file(params.dbSNP_dir)

    cohort_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.cohort, row.ancestry, row.ascertainment, file("${params.cohorts_dir}/${row.rel_path}")) }

    FORMAT_SUMSTATS(cohort_ch)
    SPLIT_SUMSTATS(FORMAT_SUMSTATS.out.formatted_stats)
    MATCH_DBSNP_PARALLEL(SPLIT_SUMSTATS.out.chr_files.transpose(), dbSNP_dir)
    
    merge_input = MATCH_DBSNP_PARALLEL.out.matched_chr.groupTuple(by: [0,1,2,3])
    MERGE_SUMSTATS(merge_input)

    MUNGE_SUMSTATS(MERGE_SUMSTATS.out.merged_file)
    LDSC_H2(MUNGE_SUMSTATS.out.munged_file)
    GENERATE_REPORT(LDSC_H2.out.log_file.collect())
	
    QC_PLOTS(MERGE_SUMSTATS.out.merged_file)
    COMBINE_LAMBDAS(QC_PLOTS.out.lambda_val.collect())

    mrmega_prep_ch = PREP_MR_MEGA(MERGE_SUMSTATS.out.merged_file)
    
    // Create the file list for MR-MEGA: [Path, CohortName]
    mrmega_list = mrmega_prep_ch.mrmega_ready
    .map { cohort, ancestry, file -> "input_files/${file.name}\t${cohort}_${ancestry}" }
    .collectFile(name: 'mrmega_input_list.txt', newLine: true)


    RUN_MR_MEGA(mrmega_list, mrmega_prep_ch.mrmega_ready.map{it[2]}.collect())

}
