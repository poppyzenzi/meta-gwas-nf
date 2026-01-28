nextflow.enable.dsl = 2

// --- Parameters ---
params.input_csv   = "${projectDir}/data/cohorts.csv"
params.cohorts_dir = '/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/cohorts'
params.outdir      = "${projectDir}/results/formatted_sumstats"
params.resultsDir  = "${projectDir}/results"
params.dbSNP_dir   = "/exports/igmm/eddie/GenScotDepression/users/poppy/aGWAS/checks/ref_panels/split_dbSNP"

// --- Processes ---

process FORMAT_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), path(raw_sumstats)

    output:
    tuple val(cohort), val(ancestry), val(build), path("${cohort}_${ancestry}_${build}_formatted.txt"), emit: formatted_stats

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

        # Priority logic for EAF: checks EAF, then EAFCO, then EAFCA
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
    tuple val(cohort), val(ancestry), val(build), path(formatted_txt)

    output:
    tuple val(cohort), val(ancestry), val(build), path("split_chr_*.txt"), emit: chr_files

    script:
    """
    tail -n +2 ${formatted_txt} | awk '{print > "split_chr_"\$1".txt"}'
    """
}

process MATCH_DBSNP_PARALLEL {
    tag "${cohort}_${ancestry}_chr${chr_num}"
    memory '8 GB'
    cpus 1

    input:
    tuple val(cohort), val(ancestry), val(build), path(chr_file)
    path dbSNP_dir

    output:
    tuple val(cohort), val(ancestry), path("matched_chr_${chr_num}.txt"), emit: matched_chr

    script:
    // This allows both the tag and the ref_file path to use the chr number
    chr_num = chr_file.baseName.replaceAll(/split_chr_/, "")
    def ref_file = "${dbSNP_dir}/dbSNP_chr${chr_num}.txt"
    
    if ( build == "b37" )
        """
        awk 'NR==FNR {
            marker[\$5,\$6,\$2,\$3]=\$4; rsid[\$5,\$6,\$2,\$3]=\$1; next
        }
        {
            if ((\$1,\$2,\$5,\$6) in marker) 
                print \$0,marker[\$1,\$2,\$5,\$6],rsid[\$1,\$2,\$5,\$6];
            else if ((\$1,\$2,\$6,\$5) in marker) 
                print \$0,marker[\$1,\$2,\$6,\$5],rsid[\$1,\$2,\$6,\$5];
            else 
                print \$0,"missing","missing";
        }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
        """
    else
        """
        awk 'NR==FNR {
            marker[\$11,\$12,\$8,\$9]=\$4; rsid[\$11,\$12,\$8,\$9]=\$1; next
        }
        {
            if ((\$1,\$2,\$5,\$6) in marker) 
                print \$0,marker[\$1,\$2,\$5,\$6],rsid[\$1,\$2,\$5,\$6];
            else if ((\$1,\$2,\$6,\$5) in marker) 
                print \$0,marker[\$1,\$2,\$6,\$5],rsid[\$1,\$2,\$6,\$5];
            else 
                print \$0,"missing","missing";
        }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
        """
}

process MERGE_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/matched", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), path(matched_files)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_matched_full.txt")

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
    tuple val(cohort), val(ancestry), path(matched_file)

    output:
    tuple val(cohort), val(ancestry), path("${cohort}_${ancestry}_munged.sumstats.gz"), emit: munged_file

    script:
    """
    # 1. Sanitize the file: Remove carriage returns and ensure standard spacing
    # This fixes the "No objects to concatenate" error caused by hidden characters
    cat ${matched_file} | tr -d '\\r' | awk '\$1=\$1' > cleaned.txt

    # Detect if OR column exists to set the signed-sumstats flag
    # In Step 1 we converted OR to BETA, but if you want to keep the logic:
    if grep -q -w "OR" cleaned.txt; then
        SIGNED="OR,1"
    else
        SIGNED="BETA,0"
    fi

    # Standardize the environment further
    # Explicitly use the python version inside the ldsc conda environment
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
        # Extract Intercept and SD
        intercept=\$(grep "Intercept:" \$log_file | awk '{print \$2}' | head -1)
        sd=\$(grep "Intercept:" \$log_file | awk '{print \$3}' | sed 's/[()]//g' | head -1)
        
        # Calculate 95% CI using awk
        echo | awk -v ival="\$intercept" -v s="\$sd" -v name="\$log_file" 'BEGIN { OFS=" " } { if (ival == "" || s == "") { printf "%s: Intercept data missing in log\\n", name } else { low = ival - (1.96 * s); high = ival + (1.96 * s); printf "%s: Intercept = %f, SD = %f, 95%% CI = [%.6f, %.6f]\\n", name, ival, s, low, high } }' >> LDSCintercept_summary.txt
    done
    """
}


// --- Workflow ---
workflow {
    dbSNP_dir = file(params.dbSNP_dir)

    // Setup input channel from CSV
    cohort_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.cohort, row.ancestry, file("${params.cohorts_dir}/${row.rel_path}")) }
        // .take(1)  // Uncomment this to test on the first cohort only

    // 1. Standardize formatting
    formatted      = FORMAT_SUMSTATS(cohort_ch)
    
    // 2. Split into chromosomes for parallel matching
    split_ch       = SPLIT_SUMSTATS(formatted.formatted_stats)
    
    // 3. Match against dbSNP reference by chromosome
    matching_input = split_ch.chr_files.transpose()
    matched_ch     = MATCH_DBSNP_PARALLEL(matching_input, dbSNP_dir)
    
    // 4. Group chromosome files back together and merge into one file per cohort
    merge_input    = matched_ch.matched_chr.groupTuple(by: [0,1])
    merged_output  = MERGE_SUMSTATS(merge_input)

    // 5. Munge sumstats for LDSC (Step 1 of your original LDSC script)
    // We pass the tuple [cohort, ancestry, merged_file]
    munged_ch      = MUNGE_SUMSTATS(merged_output)

    // 6. Run LDSC Heritability/Intercept (Step 2 of your original LDSC script)
    ldsc_logs      = LDSC_H2(munged_ch.munged_file)

    // 7. Generate final summary report
    // .collect() waits for all cohorts to finish and passes all logs as a list
    GENERATE_REPORT(ldsc_logs.collect())
}

