// matching to dbSNP ref file processes

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

process MATCH_DBSNP_COHORT {
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
        awk 'BEGIN {OFS=" "}
        # 1. Load dbSNP b37: Key is [Chr, Pos, Ref, Alt] -> [5, 6, 2, 3]
        NR==FNR { marker[\$5,\$6,\$2,\$3]=\$4; rsid[\$5,\$6,\$2,\$3]=\$1; next }

        # 2. Process Sumstats: 1:CHR, 2:POS, 5:EA, 6:NEA, 7:EAF, 8:BETA
        {
            # CASE 1: Your NEA (6) is Ref and EA (5) is Alt -> Already Aligned
            if ((\$1,\$2,\$6,\$5) in marker) {
                print \$0, marker[\$1,\$2,\$6,\$5], rsid[\$1,\$2,\$6,\$5]
            }
            # CASE 2: Your EA (5) is Ref and NEA (6) is Alt -> Needs Flip
            else if ((\$1,\$2,\$5,\$6) in marker) {
                orig_ea=\$5; orig_nea=\$6;
                \$5=orig_nea; \$6=orig_ea;            # Swap so EA is Alt, NEA is Ref
                if (\$7 != "NA") \$7=1-\$7;           # Flip Frequency
                if (\$8 != "NA") \$8=-1*\$8;          # Flip Beta
                
                # Fetch IDs using the key as it exists in dbSNP [Chr, Pos, Ref, Alt]
                # which is now [1, 2, 6, 5] after our swap
                print \$0, marker[\$1,\$2,\$6,\$5], rsid[\$1,\$2,\$6,\$5]
            }
            # CASE 3: Not in dbSNP
            else {
                print \$0, "missing", "missing"
            }
        }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
        """
    else
        """
        awk 'BEGIN {OFS=" "}
        # 1. Load dbSNP b38: Key is [Chr, Pos, Ref, Alt] -> [11, 12, 8, 9]
        NR==FNR { marker[\$11,\$12,\$8,\$9]=\$10; rsid[\$11,\$12,\$8,\$9]=\$7; next }

        # 2. Process Sumstats
        {
            if ((\$1,\$2,\$6,\$5) in marker) {
                print \$0, marker[\$1,\$2,\$6,\$5], rsid[\$1,\$2,\$6,\$5]
            }
            else if ((\$1,\$2,\$5,\$6) in marker) {
                orig_ea=\$5; orig_nea=\$6;
                \$5=orig_nea; \$6=orig_ea; 
                if (\$7 != "NA") \$7=1-\$7;
                if (\$8 != "NA") \$8=-1*\$8;
                
                print \$0, marker[\$1,\$2,\$6,\$5], rsid[\$1,\$2,\$6,\$5]
            }
            else {
                print \$0, "missing", "missing"
            }
        }' ${ref_file} ${chr_file} > matched_chr_${chr_num}.txt
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


process MATCH_DBSNP_META {
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

