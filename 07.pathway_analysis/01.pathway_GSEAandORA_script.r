################################################################################
# Gene Set Enrichment Analysis (GSEA) Pipeline
# Organism: Mouse (Mus musculus, GRCm39)
# Description: Comprehensive GSEA analysis with GO, KEGG, and MSigDB
# please use conda env r43 to do next analysis
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
library(rtracklayer)
library(msigdbr)  # For MSigDB gene sets

#' Load and Parse GTF File to Create Gene Annotation
#'
#' @param gtf_file Path to GTF file
#' 
#' @return A data.frame with gene annotations
#' 
#' @export
load_gene_annotation <- function(gtf_file) {
    
    cat("Loading GTF file:", gtf_file, "\n")
    
    # Read GTF file
    gtf <- import(gtf_file)
    gtf_df <- as.data.frame(gtf)
    
    # Extract gene-level information
    gene_info <- gtf_df %>%
        filter(type == "gene") %>%
        dplyr::select(gene_id, gene_name, gene_biotype) %>%
        distinct()
    
    # Clean gene_id (remove version numbers if present)
    gene_info$gene_id_clean <- gsub("\\..*", "", gene_info$gene_id)
    
    cat("Loaded", nrow(gene_info), "genes from GTF\n\n")
    
    return(gene_info)
}


#' Convert ENSEMBL IDs to ENTREZ IDs
#'
#' @param ensembl_ids Character vector of ENSEMBL gene IDs
#' @param org_db Organism database (default: org.Mm.eg.db for mouse)
#' 
#' @return A data.frame with ENSEMBL to ENTREZ mapping
#' 
#' @export
convert_ensembl_to_entrez <- function(ensembl_ids, org_db = org.Mm.eg.db) {
    
    # Remove version numbers if present
    ensembl_clean <- gsub("\\..*", "", ensembl_ids)
    
    # Convert to ENTREZ
    conversion <- bitr(ensembl_clean,
                      fromType = "ENSEMBL",
                      toType = c("ENTREZID", "SYMBOL"),
                      OrgDb = org_db)
    
    cat("Converted", nrow(conversion), "genes out of", length(unique(ensembl_clean)), "\n")
    cat("Conversion rate:", round(nrow(conversion) / length(unique(ensembl_clean)) * 100, 2), "%\n\n")
    
    return(conversion)
}


#' Perform GO Enrichment Analysis
#'
#' @param gene_list Named numeric vector (ENTREZ IDs as names, log2FC as values)
#' @param gene_set Character vector of gene ENTREZ IDs for ORA (default: NULL for GSEA)
#' @param method Character, "GSEA" or "ORA" (default: "GSEA")
#' @param ont GO ontology: "BP", "CC", "MF", or "ALL" (default: "BP")
#' @param pvalueCutoff P-value cutoff (default: 0.05)
#' @param qvalueCutoff Q-value cutoff (default: 0.2)
#' @param minGSSize Minimum gene set size (default: 10)
#' @param maxGSSize Maximum gene set size (default: 500)
#' 
#' @return enrichResult object
#' 
#' @export
run_go_enrichment <- function(gene_list,
                              gene_set = NULL,
                              method = "GSEA",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2,
                              minGSSize = 10,
                              maxGSSize = 500) {
    
    cat("Running GO", ont, "enrichment (", method, ")...\n")
    
    if (method == "GSEA") {
        # GSEA requires ranked gene list
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        result <- gseGO(geneList = gene_list,
                       OrgDb = org.Mm.eg.db,
                       ont = ont,
                       minGSSize = minGSSize,
                       maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
    } else {
        # ORA requires gene set only
        result <- enrichGO(gene = gene_set,
                          OrgDb = org.Mm.eg.db,
                          ont = ont,
                          pvalueCutoff = pvalueCutoff,
                          qvalueCutoff = qvalueCutoff,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pAdjustMethod = "BH",
                          readable = TRUE)
    }
    
    if (!is.null(result) && nrow(result) > 0) {
        cat("  Found", nrow(result), "enriched terms\n\n")
    } else {
        cat("  No enriched terms found\n\n")
    }
    
    return(result)
}


#' Perform KEGG Enrichment Analysis
#'
#' @param gene_list Named numeric vector (ENTREZ IDs as names, log2FC as values)
#' @param gene_set Character vector of gene ENTREZ IDs for ORA
#' @param method Character, "GSEA" or "ORA" (default: "GSEA")
#' @param organism KEGG organism code (default: "mmu" for mouse)
#' @param pvalueCutoff P-value cutoff (default: 0.05)
#' @param minGSSize Minimum gene set size (default: 10)
#' @param maxGSSize Maximum gene set size (default: 500)
#' 
#' @return enrichResult object
#' 
#' @export
run_kegg_enrichment <- function(gene_list,
                                gene_set = NULL,
                                method = "GSEA",
                                organism = "mmu",
                                pvalueCutoff = 0.05,
                                minGSSize = 10,
                                maxGSSize = 500) {
    
    cat("Running KEGG enrichment (", method, ")...\n")
    
    if (method == "GSEA") {
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        result <- gseKEGG(geneList = gene_list,
                         organism = organism,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = "BH",
                         verbose = FALSE)
    } else {
        result <- enrichKEGG(gene = gene_set,
                            organism = organism,
                            pvalueCutoff = pvalueCutoff,
                            minGSSize = minGSSize,
                            maxGSSize = maxGSSize,
                            pAdjustMethod = "BH")
    }
    
    if (!is.null(result) && nrow(result) > 0) {
        cat("  Found", nrow(result), "enriched pathways\n\n")
    } else {
        cat("  No enriched pathways found\n\n")
    }
    
    return(result)
}


#' Perform MSigDB Enrichment Analysis
#'
#' @param gene_list Named numeric vector (gene symbols as names, log2FC as values)
#' @param gene_set Character vector of gene symbols for ORA
#' @param method Character, "GSEA" or "ORA" (default: "GSEA")
#' @param category MSigDB category: "H", "C1"-"C8" (default: "H" for hallmark)
#' @param subcategory MSigDB subcategory (default: NULL)
#' @param species Species name (default: "Mus musculus")
#' @param pvalueCutoff P-value cutoff (default: 0.05)
#' @param minGSSize Minimum gene set size (default: 10)
#' @param maxGSSize Maximum gene set size (default: 500)
#' 
#' @return enrichResult object
#' 
#' @export
run_msigdb_enrichment <- function(gene_list,
                                  gene_set = NULL,
                                  method = "GSEA",
                                  category = "H",
                                  subcategory = NULL,
                                  species = "Mus musculus",
                                  pvalueCutoff = 0.05,
                                  minGSSize = 10,
                                  maxGSSize = 500) {
    
    cat("Running MSigDB", category, "enrichment (", method, ")...\n")
    
    # Get MSigDB gene sets
    msigdb_sets <- msigdbr(species = species, 
                          category = category,
                          subcategory = subcategory)
    
    msigdb_t2g <- msigdb_sets %>%
        dplyr::select(gs_name, gene_symbol) %>%
        as.data.frame()
    
    if (method == "GSEA") {
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        result <- GSEA(geneList = gene_list,
                      TERM2GENE = msigdb_t2g,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = "BH",
                      verbose = FALSE)
    } else {
        result <- enricher(gene = gene_set,
                          TERM2GENE = msigdb_t2g,
                          pvalueCutoff = pvalueCutoff,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pAdjustMethod = "BH")
    }
    
    if (!is.null(result) && nrow(result) > 0) {
        cat("  Found", nrow(result), "enriched gene sets\n\n")
    } else {
        cat("  No enriched gene sets found\n\n")
    }
    
    return(result)
}


#' Comprehensive GSEA Analysis on DEG Results
#'
#' @param deg_file Path to DEG_significant.txt or DEG_all_results.txt
#' @param gtf_file Path to GTF annotation file
#' @param gene_filter Character, "all", "up", or "down" (default: "all")
#' @param padj_cutoff Adjusted p-value cutoff for filtering (default: 0.05)
#' @param lfc_cutoff Log2 fold change cutoff for filtering (default: 1)
#' @param methods Vector of enrichment methods to run (default: c("GO", "KEGG", "MSigDB"))
#' @param go_ontologies Vector of GO ontologies (default: c("BP", "CC", "MF"))
#' @param msigdb_categories Vector of MSigDB categories (default: c("H", "C2", "C5"))
#' @param gsea_method Character, "GSEA" or "ORA" (default: "GSEA")
#' @param enrich_pvalue P-value cutoff for enrichment (default: 0.05)
#' @param enrich_qvalue Q-value cutoff for enrichment (default: 0.2)
#' @param min_geneset_size Minimum gene set size (default: 10)
#' @param max_geneset_size Maximum gene set size (default: 500)
#' @param output_dir Directory to save results (default: NULL, no output)
#' @param prefix Prefix for output files (default: "GSEA")
#' 
#' @return A list containing all enrichment results
#' 
#' @export
run_comprehensive_gsea <- function(deg_file,
                                   gtf_file,
                                   gene_filter = "all",
                                   padj_cutoff = 0.05,
                                   lfc_cutoff = 0,
                                   methods = c("GO", "KEGG", "MSigDB"),
                                   go_ontologies = c("BP", "CC", "MF"),
                                   msigdb_categories = c("H", "C2", "C5"),
                                   gsea_method = "GSEA",
                                   enrich_pvalue = 0.05,
                                   enrich_qvalue = 0.2,
                                   min_geneset_size = 10,
                                   max_geneset_size = 500,
                                   output_dir = NULL,
                                   prefix = "GSEA") {
    
    cat("========== Starting Comprehensive GSEA Analysis ==========\n\n")
    
    # ========== 1. Load gene annotation ==========
    gene_annot <- load_gene_annotation(gtf_file)
    
    
    # ========== 2. Load DEG results ==========
    cat("Loading DEG results:", deg_file, "\n")
    deg_data <- read.table(deg_file, header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    cat("Total genes in DEG file:", nrow(deg_data), "\n")
    
    # Filter genes based on criteria
    if (gene_filter == "up") {
        deg_data <- deg_data %>% filter(regulation == "Up")
        cat("Filtered to", nrow(deg_data), "up-regulated genes\n")
    } else if (gene_filter == "down") {
        deg_data <- deg_data %>% filter(regulation == "Down")
        cat("Filtered to", nrow(deg_data), "down-regulated genes\n")
    }
    
    # Apply additional cutoffs if needed
    if (padj_cutoff < 1) {
        deg_data <- deg_data %>% filter(padj < padj_cutoff)
        cat("Applied padj cutoff <", padj_cutoff, ":", nrow(deg_data), "genes\n")
    }
    
    if (lfc_cutoff > 0) {
        deg_data <- deg_data %>% filter(abs(log2FoldChange) > lfc_cutoff)
        cat("Applied |log2FC| cutoff >", lfc_cutoff, ":", nrow(deg_data), "genes\n")
    }
    
    cat("\n")
    
    
    # ========== 3. Annotate with gene names and biotypes ==========
    cat("Annotating genes with GTF information...\n")
    deg_data$gene_id_clean <- gsub("\\..*", "", deg_data$gene_id)
    
    deg_annotated <- deg_data %>%
        left_join(gene_annot %>% select(-gene_id), by = c("gene_id_clean" = "gene_id_clean")) %>%
        dplyr::select(gene_id, gene_name, gene_biotype, everything(), -gene_id_clean)
    
    cat("Annotated", sum(!is.na(deg_annotated$gene_name)), "genes with gene names\n\n")
    
    
    # ========== 4. Convert to ENTREZ IDs ==========
    cat("Converting ENSEMBL IDs to ENTREZ IDs...\n")
    id_conversion <- convert_ensembl_to_entrez(deg_annotated$gene_id)
    
    # Merge with DEG data
    deg_with_entrez <- deg_annotated %>%
        mutate(gene_id_clean = gsub("\\..*", "", gene_id)) %>%
        left_join(id_conversion, by = c("gene_id_clean" = "ENSEMBL")) %>%
        filter(!is.na(ENTREZID))
    
    cat("Final gene count with ENTREZ IDs:", nrow(deg_with_entrez), "\n\n")
    
    
    # ========== 5. Prepare gene lists ==========
    # For GSEA: ranked list by log2FC
    gene_list_entrez <- setNames(deg_with_entrez$log2FoldChange, 
                                 deg_with_entrez$ENTREZID)
    gene_list_symbol <- setNames(deg_with_entrez$log2FoldChange, 
                                 deg_with_entrez$SYMBOL)
    
    # For ORA: gene set only
    gene_set_entrez <- deg_with_entrez$ENTREZID
    gene_set_symbol <- deg_with_entrez$SYMBOL
    
    
    # ========== 6. Run enrichment analyses ==========
    results <- list()
    results$annotated_genes <- deg_with_entrez
    
    
    # ========== 6.1 GO Enrichment ==========
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
                    qvalueCutoff = enrich_qvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            } else {
                results[[result_name]] <- run_go_enrichment(
                    gene_set = gene_set_entrez,
                    method = "ORA",
                    ont = ont,
                    pvalueCutoff = enrich_pvalue,
                    qvalueCutoff = enrich_qvalue,
                    minGSSize = min_geneset_size,
                    maxGSSize = max_geneset_size
                )
            }
        }
    }
    
    
    # ========== 6.2 KEGG Enrichment ==========
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
    
    
    # ========== 6.3 MSigDB Enrichment ==========
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
    
    # ========== 7. Save results ==========
    if (!is.null(output_dir)) {
        # Create subdirectory based on gene_filter to avoid overwriting
        if (gene_filter != "all") {
            final_output_dir <- file.path(output_dir, gene_filter)
        } else {
            final_output_dir <- output_dir
        }
        
        if (!dir.exists(final_output_dir)) {
            dir.create(final_output_dir, recursive = TRUE)
        }
        
        cat("\n>>>>>>>>>> Saving Results <<<<<<<<<<\n")
        
        # Save the whole results as RData
        save(results, file = file.path(final_output_dir, paste0(prefix, "_whole_results.RData")))
        cat("Saved RData file of all results\n")
        
        # Save annotated gene list
        write.table(deg_with_entrez,
                   file.path(final_output_dir, paste0(prefix, "_annotated_genes.txt")),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        cat("Saved annotated gene list\n")
        
        # Save each enrichment result
        for (result_name in names(results)) {
            if (result_name == "annotated_genes") next
            
            result <- results[[result_name]]
            if (!is.null(result) && nrow(result) > 0) {
                # Save full result
                result_df <- as.data.frame(result)
                write.table(result_df,
                           file.path(final_output_dir, paste0(prefix, "_", result_name, ".txt")),
                           sep = "\t", row.names = FALSE, quote = FALSE)
                
                # Save top 20 for quick view
                top_result <- head(result_df, 20)
                write.table(top_result,
                           file.path(final_output_dir, paste0(prefix, "_", result_name, "_top20.txt")),
                           sep = "\t", row.names = FALSE, quote = FALSE)
                
                cat("Saved", result_name, "results (", nrow(result_df), "terms )\n")
                
                # Generate plots
                tryCatch({
                    # Dotplot
                    if (nrow(result) >= 2) {
                        p <- dotplot(result, showCategory = min(20, nrow(result))) +
                            ggtitle(paste(result_name, "Enrichment"))
                        ggsave(file.path(final_output_dir, paste0(prefix, "_", result_name, "_dotplot.pdf")),
                              p, width = 10, height = 8)
                    }
                    
                    # For GSEA results, generate GSEA plot
                    if (gsea_method == "GSEA" && nrow(result) >= 1) {
                        p2 <- gseaplot2(result, geneSetID = 1:min(3, nrow(result)))
                        ggsave(file.path(final_output_dir, paste0(prefix, "_", result_name, "_gseaplot.pdf")),
                              p2, width = 10, height = 6)
                    }
                }, error = function(e) {
                    cat("  Warning: Could not generate plot for", result_name, "\n")
                })
            }
        }
        
        cat("\nAll results saved to:", final_output_dir, "\n")
    }
    
    # ========== 8. Print summary ==========
    cat("\n========== Analysis Summary ==========\n")
    cat("Total genes analyzed:", nrow(deg_with_entrez), "\n")
    cat("Gene filter:", gene_filter, "\n")
    cat("Method:", gsea_method, "\n\n")
    
    cat("Enrichment results:\n")
    for (result_name in names(results)) {
        if (result_name == "annotated_genes") next
        result <- results[[result_name]]
        if (!is.null(result)) {
            cat("  ", result_name, ":", nrow(result), "terms\n")
        }
    }
    
    cat("\n========== Analysis Completed ==========\n\n")
    
    return(results)
}


#' Batch GSEA Analysis for Multiple Contrasts
#'
#' @param deg_base_dir Base directory containing DEG results
#' @param contrast_list List of contrasts, each as c(group_col, numerator, denominator)
#' @param gtf_file Path to GTF annotation file
#' @param output_base_dir Base output directory for GSEA results
#' @param create_subdir Logical, create subdirectory for each contrast (default: TRUE)
#' @param ... Additional arguments passed to run_comprehensive_gsea
#' 
#' @return A named list of GSEA results for each contrast
#' 
#' @export
batch_gsea_analysis <- function(deg_base_dir,
                                contrast_list,
                                gtf_file,
                                output_base_dir = NULL,
                                create_subdir = TRUE,
                                ...) {
    
    cat("========== Batch GSEA Analysis ==========\n")
    cat("Total contrasts:", length(contrast_list), "\n\n")
    
    all_results <- list()
    
    for (i in seq_along(contrast_list)) {
        contrast <- contrast_list[[i]]
        contrast_name <- paste(contrast[2], "vs", contrast[3], sep = "_")
        
        cat("\n")
        cat("=" , rep("=", 70), "=\n", sep = "")
        cat("Processing contrast", i, "/", length(contrast_list), ":", contrast_name, "\n")
        cat("=" , rep("=", 70), "=\n\n", sep = "")
        
        # Construct DEG file path
        deg_file <- file.path(deg_base_dir, contrast_name, "DEG_significant.txt")
        
        if (!file.exists(deg_file)) {
            cat("Warning: DEG file not found:", deg_file, "\n")
            cat("Skipping this contrast...\n\n")
            next
        }
        
        # Set up output directory
        if (!is.null(output_base_dir)) {
            if (create_subdir) {
                contrast_outdir <- file.path(output_base_dir, contrast_name)
            } else {
                contrast_outdir <- output_base_dir
                prefix <- contrast_name
            }
        } else {
            contrast_outdir <- NULL
        }
        
        # Run GSEA analysis
        tryCatch({
            all_results[[contrast_name]] <- run_comprehensive_gsea(
                deg_file = deg_file,
                gtf_file = gtf_file,
                output_dir = contrast_outdir,
                prefix = ifelse(create_subdir, "GSEA", contrast_name),
                ...
            )
        }, error = function(e) {
            cat("Error in contrast", contrast_name, ":", e$message, "\n\n")
            all_results[[contrast_name]] <- NULL
        })
    }
    
    cat("\n========== Batch GSEA Analysis Completed ==========\n")
    
    return(all_results)
}