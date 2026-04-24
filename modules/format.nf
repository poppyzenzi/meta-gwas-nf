// formatting processes

process FORMAT_SUMSTATS {
    tag "${cohort}_${ancestry}"
    publishDir "${params.resultsDir}/formatted_sumstats", mode: 'copy'

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(raw_sumstats)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("${cohort}_${ancestry}_${build}_formatted.txt"), emit: formatted_stats

    script:
        def build_val = (cohort =~ /NORDIC_LOO|TRAILS|RAINE|GLAD|ALLOFUS/) ? "b38" : "b37"
        """
        LC_ALL=C gzip -dc -f ${raw_sumstats} | awk -v cohort="${cohort}" '
        BEGIN {
            map["CHR"] = "CHROM CHR"; map["POS"] = "GENPOS POS BP POS_38"; map["SNP"] = "SNP ID";
            map["EA"] = "EA A1 ALLELE1"; map["NEA"] = "NEA A2 ALLELE0";
            map["EAF"] = "EAF A1FREQ FREQ";
            map["EAFCA"] = "EAFCA A1FREQ_CASES";
            map["EAFCO"] = "EAFCO A1FREQ_CONTROLS";
            map["BETA"] = "BETA LOG(OR) LOG_OR LOG-OR";
            map["OR"] = "OR";
            map["SE"] = "SE"; map["P"] = "P PVAL"; map["LOG10P"] = "LOG10P";
            map["N"] = "N N_TOTAL"; map["NCA"] = "NCA N_CASES N_CAS"; map["NCO"] = "NCO N_CONTROLS N_CON";
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
            pos=\$c["POS"]; snp = (c["SNP"] ? \$c["SNP"] : chr":"pos); ea=\$c["EA"]; nea=\$c["NEA"]; se=\$c["SE"];

            if(c["BETA"]) beta=\$c["BETA"]; else if(c["OR"] && \$c["OR"] > 0) beta=log(\$c["OR"]); else beta="NA";
	    if(c["P"]) pval=\$c["P"]; else pval=exp(-\$c["LOG10P"] * log(10));
            if(c["N"]) n=\$c["N"]; else n=\$c["NCA"] + \$c["NCO"];

            if (c["EAF"] && \$c["EAF"] != "NA" && \$c["EAF"] != "") eaf_val = \$c["EAF"];
            else if (c["EAFCO"] && \$c["EAFCO"] != "NA" && \$c["EAFCO"] != "") eaf_val = \$c["EAFCO"];
            else if (c["EAFCA"] && \$c["EAFCA"] != "NA" && \$c["EAFCA"] != "") eaf_val = \$c["EAFCA"];
            else eaf_val = "NA";

            rsid = (c["RSID"] ? \$c["RSID"] : "NA");
            info = (c["INFO"] ? \$c["INFO"] : "NA");

            if(beta == "NA" || se == "NA" || se == 0) next;
print chr, pos, snp, rsid, ea, nea, eaf_val, beta, se, pval, n, info
        }

        function f(s) {
            split(map[s], a, " ");
            for (j in a) { if (header[a[j]]) return header[a[j]] }
            return 0
        }
        ' > ${cohort}_${ancestry}_${build_val}_formatted.txt
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

process PREP_MR_MEGA {
    tag "${cohort}_${ancestry}"

    input:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path(matched_file)

    output:
    tuple val(cohort), val(ancestry), val(build), val(ascertainment), path("${cohort}_${ancestry}_mrmega_ready.txt.gz"), emit: mrmega_ready

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


process FINAL_QC_REPORT_METAL {
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
        "POS"  = "BP",
        "Allele1" = "EA",
        "Allele2" = "NEA",
        "Freq1"   = "EAF",
        "FreqSE"  = "EAF_SE",
        "Effect"  = "BETA",
        "StdErr"  = "SE",
        "P-value"  = "P",
	"RSID"	  = "rsID"
    )

    # rename only the columns that exist
    existing_cols_to_rename <- names(mapping)[names(mapping) %in% names(dt)]
    setnames(dt, old = existing_cols_to_rename, new = mapping[existing_cols_to_rename])

    # capitalise alleles
    if ("EA" %in% names(dt)) dt[, EA := toupper(EA)]
    if ("NEA" %in% names(dt)) dt[, NEA := toupper(NEA)]

    # sort by chr and bp
    setorder(dt, CHR, BP)

    # set col order
    canonical_order <- c(
        "CHR","BP","varID","rsID","EA","NEA",
        "EAF","EAF_SE","BETA","SE","P",
        "N","Z","MAF",
        "INFO","Direction"
    )

    existing <- canonical_order[canonical_order %in% names(dt)]
    setcolorder(dt, c(existing, setdiff(names(dt), existing)))


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

process FORMAT_MA {
    tag "${label}"
    publishDir "${params.resultsDir}/gctb_input", mode: 'copy'

    input:
    tuple val(label), path(sumstats_gz)

    output:
    tuple val(label), path("${label}.ma"), emit: ma_file

    script:
    """
    gunzip -c ${sumstats_gz} | awk '
    BEGIN { FS="[[:space:]\\t]+"; OFS="\\t" }
    NR == 1 {
        # Define the header for GCTB
        print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"
        next
    }
    {
        print \$4, \$5, \$6, \$7, \$9, \$10, \$11, \$12
    }' > ${label}.ma
    """
}
