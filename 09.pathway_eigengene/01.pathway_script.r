################################################################################
# Driver Genes GSEA Analysis
# Description: Simplified GSEA pipeline for pre-annotated driver gene data
# Based on: original comprehensive_gsea_functions.R
################################################################################

# Source the original GSEA functions
source("/disk2/cai113/data/stateTrans/07.pathway_analysis/01.pathway_GSEAandORA_script.r")

library(tidyverse)

#' Run GSEA on Driver Genes Data
#'
#' @param driver_genes_df Data frame with driver genes (must contain gene_name, log2FoldChange)
#' @param methods Vector of enrichment methods: "GO", "KEGG", "MSigDB" (default: all)
#' @param go_ontologies GO ontologies to test (default: c("BP", "CC", "MF"))
#' @param msigdb_categories MSigDB categories (default: c("H", "C2", "C5"))
#' @param gsea_method "GSEA" or "ORA" (default: "GSEA")
#' @param enrich_pvalue P-value cutoff (default: 0.05)
#' @param min_geneset_size Minimum gene set size (default: 10)
#' @param max_geneset_size Maximum gene set size (default: 500)
#' @param output_dir Output directory (default: NULL, no save)
#' @param prefix Output file prefix (default: "DriverGenes_GSEA")
#' 
#' @return List of enrichment results
#' 
#' @export
run_driver_genes_gsea <- function(driver_genes_df,
                                  methods = c("GO", "KEGG", "MSigDB"),
                                  go_ontologies = c("BP", "CC", "MF"),
                                  msigdb_categories = c("H", "C2", "C5"),
                                  gsea_method = "GSEA",
                                  enrich_pvalue = 0.05,
                                  min_geneset_size = 5,
                                  max_geneset_size = 100,
                                  output_dir = NULL,
                                  prefix = "DriverGenes_GSEA") {
    
    cat("========== Driver Genes GSEA Analysis ==========\n\n")
    
    # ========== 1. Data validation ==========
    required_cols <- c("gene_name", "log2FoldChange")
    if (!all(required_cols %in% colnames(driver_genes_df))) {
        stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
    }
    
    cat("Input genes:", nrow(driver_genes_df), "\n")
    cat("Columns detected:", paste(colnames(driver_genes_df), collapse = ", "), "\n\n")
    
    
    # ========== 2. Convert gene symbols to ENTREZ IDs ==========
    cat("Converting gene symbols to ENTREZ IDs...\n")
    
    # Remove any genes without gene_name
    driver_genes_clean <- driver_genes_df %>%
        filter(!is.na(gene_name), gene_name != "")
    
    cat("Genes with valid symbols:", nrow(driver_genes_clean), "\n")
    
    # Convert using bitr
    id_conversion <- bitr(driver_genes_clean$gene_name,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Mm.eg.db)
    
    cat("Successfully converted:", nrow(id_conversion), "genes\n")
    cat("Conversion rate:", 
        round(nrow(id_conversion) / nrow(driver_genes_clean) * 100, 2), "%\n\n")
    
    # Merge back with original data
    driver_genes_with_entrez <- driver_genes_clean %>%
        left_join(id_conversion, by = c("gene_name" = "SYMBOL")) %>%
        filter(!is.na(ENTREZID))
    
    cat("Final gene count with ENTREZ IDs:", nrow(driver_genes_with_entrez), "\n\n")
    
    
    # ========== 3. Prepare gene lists ==========
    # For GSEA: ranked list by log2FC or contribution_score
    if ("contribution_score" %in% colnames(driver_genes_with_entrez)) {
        cat("Using contribution_score for ranking\n")
        rank_metric <- driver_genes_with_entrez$contribution_score
    } else {
        cat("Using log2FoldChange for ranking\n")
        rank_metric <- driver_genes_with_entrez$log2FoldChange
    }
    
    gene_list_entrez <- setNames(rank_metric, driver_genes_with_entrez$ENTREZID)
    gene_list_symbol <- setNames(rank_metric, driver_genes_with_entrez$gene_name)
    
    # For ORA
    gene_set_entrez <- driver_genes_with_entrez$ENTREZID
    gene_set_symbol <- driver_genes_with_entrez$gene_name
    
    cat("\n")
    
    
    # ========== 4. Run enrichment analyses ==========
    results <- list()
    results$input_genes <- driver_genes_with_entrez
    
    
    # ========== 4.1 GO Enrichment ==========
    if ("GO" %in% methods) {
        cat(">>>>>>>>>> GO Enrichment Analysis <<<<<<<<<<\n")
        
        for (ont in go_ontologies) {
            result_name <- paste0("GO_", ont)
            
            if (gsea_method == "GSEA") {
                results[[result_name]] <- run_go_enrichment(
                    gene_list = gene_list_entrez,
                    method = "GSEA",
                    ont = ont,
                    pvalueCutoff = enrich_pvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            } else {
                results[[result_name]] <- run_go_enrichment(
                    gene_set = gene_set_entrez,
                    method = "ORA",
                    ont = ont,
                    pvalueCutoff = enrich_pvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            }
        }
    }
    
    
    # ========== 4.2 KEGG Enrichment ==========
    if ("KEGG" %in% methods) {
        cat(">>>>>>>>>> KEGG Pathway Analysis <<<<<<<<<<\n")
        
        if (gsea_method == "GSEA") {
            results$KEGG <- run_kegg_enrichment(
                gene_list = gene_list_entrez,
                method = "GSEA",
                pvalueCutoff = enrich_pvalue,
                minGSSize = min_geneset_size,
                maxGSSize = max_geneset_size
            )
        } else {
            results$KEGG <- run_kegg_enrichment(
                gene_set = gene_set_entrez,
                method = "ORA",
                pvalueCutoff = enrich_pvalue,
                minGSSize = min_geneset_size,
                maxGSSize = max_geneset_size
            )
        }
    }
    
    
    # ========== 4.3 MSigDB Enrichment ==========
    if ("MSigDB" %in% methods) {
        cat(">>>>>>>>>> MSigDB Gene Set Analysis <<<<<<<<<<\n")
        
        for (cat_name in msigdb_categories) {
            result_name <- paste0("MSigDB_", cat_name)
            
            if (gsea_method == "GSEA") {
                results[[result_name]] <- run_msigdb_enrichment(
                    gene_list = gene_list_symbol,
                    method = "GSEA",
                    category = cat_name,
                    pvalueCutoff = enrich_pvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            } else {
                results[[result_name]] <- run_msigdb_enrichment(
                    gene_set = gene_set_symbol,
                    method = "ORA",
                    category = cat_name,
                    pvalueCutoff = enrich_pvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            }
        }
    }
    
    
    # ========== 5. Save results ==========
    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        cat("\n>>>>>>>>>> Saving Results <<<<<<<<<<\n")
        
        # Save whole results as RData
        save(results, file = file.path(output_dir, paste0(prefix, "_results.RData")))
        cat("Saved RData file\n")
        
        # Save input genes with ENTREZ IDs
        write.table(driver_genes_with_entrez,
                   file.path(output_dir, paste0(prefix, "_genes_with_entrez.txt")),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        cat("Saved gene list with ENTREZ IDs\n")
        
        # Save each enrichment result
        for (result_name in names(results)) {
            if (result_name == "input_genes") next
            
            result <- results[[result_name]]
            if (!is.null(result) && nrow(result) > 0) {
                # Save full result
                result_df <- as.data.frame(result)
                write.table(result_df,
                           file.path(output_dir, paste0(prefix, "_", result_name, ".txt")),
                           sep = "\t", row.names = FALSE, quote = FALSE)
                
                # Save top 20
                top_result <- head(result_df, 20)
                write.table(top_result,
                           file.path(output_dir, paste0(prefix, "_", result_name, "_top20.txt")),
                           sep = "\t", row.names = FALSE, quote = FALSE)
                
                cat("Saved", result_name, "(", nrow(result_df), "terms )\n")
                
                # Generate plots
                tryCatch({
                    if (nrow(result) >= 2) {
                        p <- dotplot(result, showCategory = min(20, nrow(result))) +
                            ggtitle(paste(result_name, "Enrichment"))
                        ggsave(file.path(output_dir, paste0(prefix, "_", result_name, "_dotplot.pdf")),
                              p, width = 10, height = 8)
                    }
                    
                    if (gsea_method == "GSEA" && nrow(result) >= 1) {
                        p2 <- gseaplot2(result, geneSetID = 1:min(3, nrow(result)))
                        ggsave(file.path(output_dir, paste0(prefix, "_", result_name, "_gseaplot.pdf")),
                              p2, width = 10, height = 6)
                    }
                }, error = function(e) {
                    cat("  Warning: Could not generate plot for", result_name, "\n")
                })
            }
        }
        
        cat("\nAll results saved to:", output_dir, "\n")
    }
    
    
    # ========== 6. Print summary ==========
    cat("\n========== Analysis Summary ==========\n")
    cat("Input genes:", nrow(driver_genes_df), "\n")
    cat("Genes with ENTREZ IDs:", nrow(driver_genes_with_entrez), "\n")
    cat("Method:", gsea_method, "\n\n")
    
    cat("Enrichment results:\n")
    for (result_name in names(results)) {
        if (result_name == "input_genes") next
        result <- results[[result_name]]
        if (!is.null(result)) {
            cat("  ", result_name, ":", nrow(result), "terms\n")
        }
    }
    
    cat("\n========== Analysis Completed ==========\n\n")
    
    return(results)
}


#' Batch GSEA for Multiple Driver Gene Sets
#'
#' @param driver_genes_list Named list of driver gene data frames
#' @param output_base_dir Base output directory
#' @param ... Additional arguments for run_driver_genes_gsea
#' 
#' @return Named list of all GSEA results
#' 
#' @export
batch_driver_genes_gsea <- function(driver_genes_list,
                                    output_base_dir = NULL,
                                    ...) {
    
    cat("========== Batch Driver Genes GSEA ==========\n")
    cat("Total contrasts:", length(driver_genes_list), "\n\n")
    
    all_results <- list()
    
    for (i in seq_along(driver_genes_list)) {
        contrast_name <- names(driver_genes_list)[i]
        driver_genes <- driver_genes_list[[i]]
        
        cat("\n")
        cat("=", rep("=", 70), "=\n", sep = "")
        cat("Processing", i, "/", length(driver_genes_list), ":", contrast_name, "\n")
        cat("=", rep("=", 70), "=\n\n", sep = "")
        
        # Set output directory
        if (!is.null(output_base_dir)) {
            contrast_outdir <- file.path(output_base_dir, contrast_name)
        } else {
            contrast_outdir <- NULL
        }
        
        # Run GSEA
        tryCatch({
            all_results[[contrast_name]] <- run_driver_genes_gsea(
                driver_genes_df = driver_genes,
                output_dir = contrast_outdir,
                prefix = "DriverGenes_GSEA",
                ...
            )
        }, error = function(e) {
            cat("Error in", contrast_name, ":", e$message, "\n\n")
            all_results[[contrast_name]] <- NULL
        })
    }
    
    cat("\n========== Batch Analysis Completed ==========\n")
    
    return(all_results)
}