# subsequent analysis of DEG results, only fit for ./00.DEG_DESeq2_script.r
library(UpSetR)
library(ggplot2)
library(tidyverse)
library(openxlsx)

#' Plot UpSet diagram for differential expression analysis results
#' 
#' @param deg_dir Directory containing DEG results
#' @param contrast_list List of comparison pairs
#' @param output_dir Output directory for plots and results
#' @param regulation_type Type of regulation: "Up", "Down", "All"
#' @param padj_cutoff Adjusted p-value cutoff
#' @param lfc_cutoff Log2 fold change cutoff
#' @param min_set_size Minimum set size to display
#' @param top_n_sets Number of top intersections to show
#' @return List containing UpSet plot object and intersection data
plot_deg_upset <- function(deg_dir, contrast_list, output_dir = "./upset_plots",
                          regulation_type = "All", padj_cutoff = 0.05, 
                          lfc_cutoff = 1, min_set_size = 1, top_n_sets = 20) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("Loading differential expression data...\n")
  deg_data <- read_all_deg_files(deg_dir, contrast_list, regulation_type, padj_cutoff, lfc_cutoff)
  
  if (is.null(deg_data) || length(deg_data) == 0) {
    stop("No valid DEG data found. Please check file paths and parameters.")
  }
  
  cat("Successfully loaded DEG data from", length(deg_data), "comparisons\n")
  gene_sets <- create_gene_sets(deg_data)
  
  if (length(gene_sets) == 0) {
    stop("No differentially expressed genes found meeting the criteria.")
  }
  
  upset_plot <- create_upset_plot(gene_sets, min_set_size, top_n_sets, output_dir, regulation_type)
  intersection_data <- save_intersection_results(gene_sets, output_dir, regulation_type)
  generate_upset_summary(gene_sets, output_dir, regulation_type)
  
  cat("UpSet analysis completed! Results saved to:", output_dir, "\n")
  
  return(list(
    upset_plot = upset_plot,
    intersection_data = intersection_data,
    gene_sets = gene_sets,
    parameters = list(
      regulation_type = regulation_type,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff,
      contrast_count = length(contrast_list)
    )
  ))
}

#' Read all DEG result files from specified comparisons
read_all_deg_files <- function(deg_dir, contrast_list, regulation_type, padj_cutoff, lfc_cutoff) {
  deg_data <- list()
  
  for (contrast in contrast_list) {
    contrast_name <- paste(contrast[2], "vs", contrast[3], sep = "_")
    deg_file <- file.path(deg_dir, contrast_name, "DEG_significant.txt")
    
    cat("Reading:", deg_file, "\n")
    
    if (file.exists(deg_file)) {
      deg_df <- read.delim(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
      filtered_genes <- filter_deg_genes(deg_df, regulation_type, padj_cutoff, lfc_cutoff)
      
      if (length(filtered_genes) > 0) {
        deg_data[[contrast_name]] <- filtered_genes
        cat("  -", contrast_name, ":", length(filtered_genes), "genes\n")
      } else {
        cat("  -", contrast_name, ": No genes meeting criteria\n")
      }
    } else {
      warning("File not found: ", deg_file)
    }
  }
  return(deg_data)
}

#' Filter DEGs based on regulation type and statistical cutoffs
filter_deg_genes <- function(deg_df, regulation_type, padj_cutoff, lfc_cutoff) {
  # Handle column name variations
  if (!"padj" %in% colnames(deg_df) && "FDR" %in% colnames(deg_df)) {
    deg_df$padj <- deg_df$FDR
  }
  if (!"log2FoldChange" %in% colnames(deg_df) && "logFC" %in% colnames(deg_df)) {
    deg_df$log2FoldChange <- deg_df$logFC
  }
  
  if (!"padj" %in% colnames(deg_df) || !"log2FoldChange" %in% colnames(deg_df)) {
    warning("Required columns (padj/FDR and log2FoldChange/logFC) not found")
    return(character(0))
  }
  
  # Apply filters based on regulation type
  if (regulation_type == "Up") {
    filtered <- deg_df %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff)
  } else if (regulation_type == "Down") {
    filtered <- deg_df %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff)
  } else {
    filtered <- deg_df %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  }
  
  # Extract gene IDs
  if (nrow(filtered) > 0) {
    if ("gene_id" %in% colnames(filtered)) return(unique(filtered$gene_id))
    if ("GeneID" %in% colnames(filtered)) return(unique(filtered$GeneID))
    return(unique(rownames(filtered)))
  }
  return(character(0))
}

#' Create gene sets from filtered DEG data
create_gene_sets <- function(deg_data) {
  non_empty_sets <- deg_data[sapply(deg_data, length) > 0]
  return(non_empty_sets)
}

#' Generate UpSet plot from gene sets
create_upset_plot <- function(gene_sets, min_set_size, top_n_sets, output_dir, regulation_type) {
  all_genes <- unique(unlist(gene_sets))
  binary_matrix <- matrix(0, nrow = length(all_genes), ncol = length(gene_sets))
  colnames(binary_matrix) <- names(gene_sets)
  rownames(binary_matrix) <- all_genes
  
  for (i in 1:length(gene_sets)) {
    set_name <- names(gene_sets)[i]
    binary_matrix[all_genes %in% gene_sets[[set_name]], i] <- 1
  }
  
  upset_df <- as.data.frame(binary_matrix)
  output_file <- file.path(output_dir, paste0("upset_plot_", regulation_type, ".pdf"))
  
  pdf(output_file, width = 8, height = 6)
  p <- upset(upset_df,
            nsets = length(gene_sets),
            nintersects = top_n_sets,
            sets = names(gene_sets),
            mb.ratio = c(0.6, 0.4),
            order.by = "freq",
            decreasing = TRUE,
            number.angles = 30,
            point.size = 3.5,
            line.size = 1,
            mainbar.y.label = "Gene Intersections",
            sets.x.label = paste("DEGs per Comparison (", regulation_type, ")"),
            text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
  )
  print(p)
  dev.off()
  
  cat("UpSet plot saved:", output_file, "\n")
  return(p)
}

#' Save detailed intersection results to Excel file
save_intersection_results <- function(gene_sets, output_dir, regulation_type) {
  all_genes <- unique(unlist(gene_sets))
  set_names <- names(gene_sets)
  
  intersection_df <- data.frame(gene_id = all_genes)
  for (set_name in set_names) {
    intersection_df[[set_name]] <- ifelse(all_genes %in% gene_sets[[set_name]], 1, 0)
  }
  
  intersection_df$set_count <- rowSums(intersection_df[, set_names])
  intersection_df$intersection_pattern <- apply(intersection_df[, set_names], 1, 
    function(x) paste(set_names[which(x == 1)], collapse = "|"))
  
  intersection_df <- intersection_df %>% arrange(desc(set_count), intersection_pattern)
  
  output_file <- file.path(output_dir, paste0("intersection_details_", regulation_type, ".xlsx"))
  wb <- createWorkbook()
  
  addWorksheet(wb, "All Intersection Genes")
  writeData(wb, "All Intersection Genes", intersection_df)
  
  max_count <- max(intersection_df$set_count)
  for (i in max_count:1) {
    subset_data <- intersection_df %>% filter(set_count == i)
    if (nrow(subset_data) > 0) {
      sheet_name <- paste0("In_", i, "_Sets")
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, subset_data)
    }
  }
  
  core_genes <- intersection_df %>% filter(set_count == length(set_names))
  if (nrow(core_genes) > 0) {
    addWorksheet(wb, "Core_Genes")
    writeData(wb, "Core_Genes", core_genes)
  }
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Intersection details saved:", output_file, "\n")
  return(intersection_df)
}

#' Generate summary statistics for UpSet analysis
generate_upset_summary <- function(gene_sets, output_dir, regulation_type) {
  set_names <- names(gene_sets)
  set_sizes <- sapply(gene_sets, length)
  total_genes <- length(unique(unlist(gene_sets)))
  core_genes <- Reduce(intersect, gene_sets)
  
  summary_file <- file.path(output_dir, paste0("upset_summary_", regulation_type, ".txt"))
  sink(summary_file)
  cat("DEG UpSet Analysis Summary\n")
  cat("===========================\n\n")
  cat("Regulation type:", regulation_type, "\n")
  cat("Number of comparisons:", length(gene_sets), "\n")
  cat("Total unique DEGs:", total_genes, "\n")
  cat("Core genes (all comparisons):", length(core_genes), "\n\n")
  
  cat("DEGs per comparison:\n")
  for (i in 1:length(set_sizes)) {
    cat("  -", set_names[i], ":", set_sizes[i], "\n")
  }
  
  cat("\nCore gene list:\n")
  if (length(core_genes) > 0) {
    for (gene in core_genes) {
      cat("  -", gene, "\n")
    }
  } else {
    cat("  No core genes found\n")
  }
  sink()
  cat("Analysis summary saved:", summary_file, "\n")
}


#-------------------------- run analysis example
base_dir <- "/disk2/cai113/data/stateTrans"
DEG_dir <- file.path(base_dir, "06.DEG_pathway", "DEG_results")
output_dir <- file.path(base_dir, "06.DEG_pathway", "Upset_analysis")

contrast_list <- list(
  c("cp", "c1", "c1_star"),
  c("cp", "c2", "c1_star"), 
  c("cp", "c3", "c1_star"),
  c("cp", "c2", "c1"),
  c("cp", "c3", "c1"),
  c("cp", "c3", "c2")
)

result_all <- plot_deg_upset(
  deg_dir = DEG_dir,
  contrast_list = contrast_list,
  output_dir = output_dir,
  regulation_type = "All",
  padj_cutoff = 0.05,
  lfc_cutoff = 1.5,
  top_n_sets = 20
)
