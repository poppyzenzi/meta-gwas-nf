// processes to generate QC reports and plots

process CORR_EAF {
    tag "${group_label}"
    publishDir "${params.resultsDir}/qc/eaf_correlations", mode: 'copy'

    input:
    tuple val(group_label), path(files)

    output:
    path "${group_label}_eaf_corr.csv"
    path "${group_label}_overlap_counts.csv"

    script:
    """
    export LC_ALL=C

    # header as tabs>  clean the labels and join them with tabs to match the join command output
    labels=\$(ls *_reduced.txt | sed 's/_reduced.txt//g' | tr '\\n' '\\t' | sed 's/\\t\$/\\n/')
    echo -e "ID\\t\$labels" > headers.txt

    # iterative full outer join
    files=( *_reduced.txt )
    cp "\${files[0]}" combined_temp.txt

    for ((i=1; i<\${#files[@]}; i++)); do
        # use -o auto to ensure the 'join' respects the tab-delim structure
        join -t\$'\\t' -a1 -a2 -e "NA" -o auto combined_temp.txt "\${files[i]}" > combined_new.txt
        mv combined_new.txt combined_temp.txt
    done

    # combine tab sep files
    cat headers.txt combined_temp.txt > final_matrix.txt

    Rscript - <<EOF
    library(data.table)
    dt <- fread("final_matrix.txt", header = TRUE, na.strings="NA", sep="\\t")

    if (nrow(dt) == 0) stop("The matrix is empty!")

    # make cols numeric
    cols <- names(dt)[-1]
    dt[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

    mat <- as.matrix(dt[, ..cols])

    # correlation matrix
    corr_matrix <- cor(mat, use = "pairwise.complete.obs")
    write.csv(corr_matrix, "${group_label}_eaf_corr.csv")

    # overlap counts 
    count_overlap <- function(x, y) sum(!is.na(x) & !is.na(y))
    overlap_mat <- matrix(NA, nrow=length(cols), ncol=length(cols), dimnames=list(cols, cols))
    for(i in seq_along(cols)) {
        for(j in seq_along(cols)) {
            overlap_mat[i,j] <- count_overlap(mat[,i], mat[,j])
        }
    }
    write.csv(overlap_mat, "${group_label}_overlap_counts.csv")
    EOF
    """
}

process SUMMARISE_COHORT_LDSC {
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
