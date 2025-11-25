library(DESeq2)
library(tidyverse)

#' Perform DESeq2 Differential Expression Analysis
#'
#' @param count_matrix A numeric matrix of raw counts (genes x samples)
#' @param sample_info A data.frame with sample metadata, rownames must match colnames of count_matrix
#' @param group_col Character, column name in sample_info for grouping (e.g., "group", "genotype")
#' @param contrast Character vector of length 3: c(group_col, "numerator", "denominator")
#'                For example: c("genotype", "ko", "wt") means ko vs wt (fold change = ko/wt)
#' @param design Formula for DESeq2 design (default: ~group_col)
#' @param filter_low_expr Logical, whether to filter low expression genes (default: TRUE)
#' @param min_count Minimum count threshold for filtering (default: 10)
#' @param min_samples Minimum number of samples passing min_count (default: 3)
#' @param padj_cutoff Adjusted p-value cutoff for significance (default: 0.05)
#' @param lfc_cutoff Log2 fold change cutoff for significance (default: 1)
#' @param shrink_lfc Logical, whether to shrink log2 fold changes (default: TRUE)
#' @param output_dir Character, directory to save results (default: NULL, no output)
#' @param prefix Character, prefix for output file names (default: "DEG")
#' 
#' @return A list containing:
#'   - dds: DESeqDataSet object
#'   - res: DESeq2 results table
#'   - res_sig: Significant DEGs only
#'   - summary: Summary statistics
#'   
#' @export
run_deseq2_analysis <- function(count_matrix,
                                sample_info,
                                group_col,
                                contrast,
                                design = NULL,
                                filter_low_expr = TRUE,
                                min_count = 10,
                                min_samples = 3,
                                padj_cutoff = 0.05,
                                lfc_cutoff = 1,
                                shrink_lfc = TRUE,
                                output_dir = NULL,
                                prefix = "DEG") {
    
    # ========== 1. Input validation ==========
    cat("========== Starting DESeq2 Analysis ==========\n")
    
    if (!all(colnames(count_matrix) %in% rownames(sample_info))) {
        stop("Error: Not all samples in count_matrix are found in sample_info rownames")
    }
    
    if (!group_col %in% colnames(sample_info)) {
        stop("Error: group_col '", group_col, "' not found in sample_info")
    }
    
    # Reorder sample_info to match count_matrix
    sample_info <- sample_info[colnames(count_matrix), , drop = FALSE]
    
    # Check contrast validity
    if (length(contrast) != 3) {
        stop("Error: contrast must be a vector of length 3: c(variable, numerator, denominator)")
    }
    
    group_levels <- unique(sample_info[[group_col]])
    if (!all(contrast[2:3] %in% group_levels)) {
        stop("Error: contrast levels not found in ", group_col, 
             "\nAvailable levels: ", paste(group_levels, collapse = ", "))
    }
    
    cat("Input validation passed\n")
    cat("Samples:", ncol(count_matrix), "\n")
    cat("Genes:", nrow(count_matrix), "\n")
    cat("Contrast:", contrast[2], "vs", contrast[3], "\n\n")
    
    
    # ========== 2. Create DESeqDataSet ==========
    if (is.null(design)) {
        design <- as.formula(paste0("~ ", group_col))
    }
    
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = sample_info,
        design = design
    )
    
    cat("DESeqDataSet created with design:", deparse(design), "\n")
    
    
    # ========== 3. Filter low expression genes ==========
    if (filter_low_expr) {
        keep <- rowSums(counts(dds) >= min_count) >= min_samples
        dds <- dds[keep, ]
        cat("Low expression filtering: kept", sum(keep), "genes (removed", 
            sum(!keep), "genes)\n")
        cat("Filter criteria: count >=", min_count, "in at least", min_samples, "samples\n\n")
    }
    
    
    # ========== 4. Run DESeq2 ==========
    cat("Running DESeq2 normalization and differential expression...\n")
    dds <- DESeq(dds, quiet = TRUE)
    cat("DESeq2 analysis completed\n\n")
    
    
    # ========== 5. Extract results ==========
    cat("Extracting results for contrast:", paste(contrast, collapse = " "), "\n")
    
    if (shrink_lfc) {
        # Get coefficient name for shrinkage
        res_names <- resultsNames(dds)
        coef_name <- res_names[grep(paste0(group_col, "_", contrast[2], "_vs_", contrast[3]), 
                                     res_names)]
        
        if (length(coef_name) == 0) {
            warning("Could not find coefficient for shrinkage, using lfcShrink with contrast")
            res <- lfcShrink(dds, contrast = contrast, type = "ashr", quiet = TRUE)
        } else {
            res <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
        }
        cat("Log2 fold change shrinkage applied\n")
    } else {
        res <- results(dds, contrast = contrast)
    }
    
    res <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj, desc(abs(log2FoldChange)))
    
    
    # ========== 6. Identify significant DEGs ==========
    res_sig <- res %>%
        filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
        mutate(regulation = ifelse(log2FoldChange > 0, "Up", "Down"))
    
    n_up <- sum(res_sig$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res_sig$regulation == "Down", na.rm = TRUE)
    
    cat("\n========== DEG Summary ==========\n")
    cat("Total genes tested:", nrow(res), "\n")
    cat("Significant DEGs (padj <", padj_cutoff, "& |log2FC| >", lfc_cutoff, "):", nrow(res_sig), "\n")
    cat("  - Up-regulated:", n_up, "\n")
    cat("  - Down-regulated:", n_down, "\n\n")
    
    
    # ========== 7. Create summary object ==========
    summary_stats <- list(
        total_genes = nrow(res),
        sig_genes = nrow(res_sig),
        up_genes = n_up,
        down_genes = n_down,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        contrast = contrast
    )
    
    
    # ========== 8. Save results ==========
    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Save all results
        write.table(res,
                    file.path(output_dir, paste0(prefix, "_all_results.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Save significant DEGs
        write.table(res_sig,
                    file.path(output_dir, paste0(prefix, "_significant.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Save summary
        sink(file.path(output_dir, paste0(prefix, "_summary.txt")))
        cat("DESeq2 Analysis Summary\n")
        cat("======================\n\n")
        cat("Contrast:", paste(contrast, collapse = " "), "\n")
        cat("Total genes tested:", summary_stats$total_genes, "\n")
        cat("Significant DEGs:", summary_stats$sig_genes, "\n")
        cat("  Up-regulated:", summary_stats$up_genes, "\n")
        cat("  Down-regulated:", summary_stats$down_genes, "\n")
        cat("\nCutoffs:\n")
        cat("  Adjusted p-value:", summary_stats$padj_cutoff, "\n")
        cat("  Log2 fold change:", summary_stats$lfc_cutoff, "\n")
        sink()
        
        cat("Results saved to:", output_dir, "\n\n")
    }
    
    
    # ========== 9. Return results ==========
    return(list(
        dds = dds,
        res = res,
        res_sig = res_sig,
        summary = summary_stats
    ))
}


#' Batch DEG Analysis for Multiple Contrasts
#'
#' @param count_matrix A numeric matrix of raw counts (genes x samples)
#' @param sample_info A data.frame with sample metadata
#' @param group_col Character, column name for grouping
#' @param contrast_list List of contrasts, each as c(group_col, numerator, denominator)
#' @param output_dir Character, base directory to save results (default: NULL, no output)
#' @param create_subdir Logical, whether to create subdirectories for each contrast (default: TRUE)
#' @param ... Additional arguments passed to run_deseq2_analysis (except output_dir and prefix)
#' 
#' @return A named list of results for each contrast
#' 
#' @export
batch_deseq2_analysis <- function(count_matrix,
                                  sample_info,
                                  group_col,
                                  contrast_list,
                                  output_dir = NULL,
                                  create_subdir = TRUE,
                                  ...) {
    
    cat("========== Batch DEG Analysis ==========\n")
    cat("Total contrasts:", length(contrast_list), "\n\n")
    
    results <- list()
    
    # Remove output_dir and prefix from additional arguments if present
    extra_args <- list(...)
    extra_args$output_dir <- NULL
    extra_args$prefix <- NULL
    
    for (i in seq_along(contrast_list)) {
        contrast <- contrast_list[[i]]
        contrast_name <- paste(contrast[2], "vs", contrast[3], sep = "_")
        
        cat(">>> Processing contrast", i, "/", length(contrast_list), ":", contrast_name, "\n")
        
        # Set up output directory and prefix for this contrast
        if (!is.null(output_dir)) {
            if (create_subdir) {
                # Create subdirectory for each contrast
                contrast_dir <- file.path(output_dir, contrast_name)
                contrast_prefix <- "DEG"
            } else {
                # Use base directory with unique prefix
                contrast_dir <- output_dir
                contrast_prefix <- contrast_name
            }
        } else {
            contrast_dir <- NULL
            contrast_prefix <- "DEG"
        }
        
        tryCatch({
            results[[contrast_name]] <- do.call(
                run_deseq2_analysis,
                c(
                    list(
                        count_matrix = count_matrix,
                        sample_info = sample_info,
                        group_col = group_col,
                        contrast = contrast,
                        output_dir = contrast_dir,
                        prefix = contrast_prefix
                    ),
                    extra_args
                )
            )
        }, error = function(e) {
            cat("Error in contrast", contrast_name, ":", e$message, "\n\n")
            results[[contrast_name]] <- NULL
        })
        
        cat("\n")
    }
    
    cat("========== Batch Analysis Completed ==========\n")
    
    return(results)
}


#' Quick Visualization of DEG Results
#'
#' @param deg_result Result object from run_deseq2_analysis
#' @param plot_type Character, one of "volcano", "ma" (default: "volcano")
#' @param title Character, plot title
#' 
#' @return A ggplot object
#' 
#' @export
plot_deg_results <- function(deg_result, 
                             plot_type = "volcano",
                             title = NULL) {
    
    res <- deg_result$res
    padj_cutoff <- deg_result$summary$padj_cutoff
    lfc_cutoff <- deg_result$summary$lfc_cutoff
    
    res <- res %>%
        mutate(
            significance = case_when(
                is.na(padj) ~ "NS",
                padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Up",
                padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
                TRUE ~ "NS"
            )
        )
    
    color_values <- c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey70")
    
    if (plot_type == "volcano") {
        p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
            geom_point(alpha = 0.6, size = 1.5) +
            scale_color_manual(values = color_values) +
            geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey30") +
            geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey30") +
            labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
                 title = title %||% "Volcano Plot") +
            theme_classic() +
            theme(legend.position = "top")
        
    } else if (plot_type == "ma") {
        p <- ggplot(res, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
            geom_point(alpha = 0.6, size = 1.5) +
            scale_color_manual(values = color_values) +
            geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey30") +
            geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
            labs(x = "Log10 Mean Expression", y = "Log2 Fold Change",
                 title = title %||% "MA Plot") +
            theme_classic() +
            theme(legend.position = "top")
    }
    
    return(p)
}