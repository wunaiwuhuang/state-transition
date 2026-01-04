# =============================================================================
# Driver Gene Differential Expression Analysis in Single-Cell RNA-seq Data
# Author: wuguojia
# Date: 2026-01-02
# Purpose: Identify and visualize driver genes differentially expressed between 
#          WT and Ncf2_KO conditions, both globally and within cell types.
# =============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

# =============================================================================
# Configuration & Constants
# =============================================================================
BASE_DIR <- "/disk2/cai113/data/stateTrans/10.scRNAseq_analysis/"
SCDATA_PATH <- file.path(BASE_DIR, "yhz_data/mian/sce_celltype.RData")
DRIVER_PATH <- "/disk2/cai113/data/stateTrans/08.eigengene/driver_gene_selection_results.rdata"

# Output file paths (for clarity and consistency)
OUTPUT_FILES <- list(
  genes_info_csv       = file.path(BASE_DIR, "00.genes_infomation.csv"),
  overall_sig_csv      = file.path(BASE_DIR, "01.Overall_sig_KO_vs_WT.csv"),
  volcano_pdf          = file.path(BASE_DIR, "01.volcano_driver_genes.pdf"),
  violin_pdf_overall   = file.path(BASE_DIR, "01.violin_overall_genes.pdf"),
  celltype_sig_csv     = file.path(BASE_DIR, "02.Celltype_sig_KO_vs_WT.csv"),
  summary_csv          = file.path(BASE_DIR, "02.celltype_deg_summary.csv"),
  gene_celltype_map    = file.path(BASE_DIR, "02.gene_to_celltype_mapping.csv"),
  heatmap_pdf          = file.path(BASE_DIR, "02.heatmap_driver_genes_celltype.pdf"),
  dotplot_pdf          = file.path(BASE_DIR, "02.dotplot_driver_genes_celltype.pdf"),
  violin_pdf_ct        = file.path(BASE_DIR, "02.violin_celltype_genes.pdf"),
  violin_pdf_ct_fc     = file.path(BASE_DIR, "02.violin_celltype_genes_fc.pdf"),
  analysis_summary_txt = file.path(BASE_DIR, "00.analysis_summary.txt")
)

# =============================================================================
# Utility Functions
# =============================================================================

#' Save a data frame to CSV with consistent formatting
#'
#' @param data Data frame to save
#' @param path File path to save the CSV
save_csv <- function(data, path) {
  write.csv(data, file = path, row.names = FALSE)
}


#' Perform differential expression analysis using Seurat's FindMarkers
#' @param object Seurat object
#' @param ident1 Experimental group (e.g., "Ncf2_KO")
#' @param ident2 Control group (e.g., "WT")
#' @param features Genes to test
#' @param min_pct Minimum fraction of cells expressing the gene in either group
#' @return Data frame of DEG results with metadata
run_de_analysis <- function(object, ident1 = "Ncf2_KO", ident2 = "WT", features, min_pct = 0.01) {
  deg <- FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    features = features,
    logfc.threshold = 0,
    min.pct = min_pct,
    test.use = "wilcox",
    verbose = FALSE
  )
  
  if (nrow(deg) == 0) return(NULL)
  
  deg$gene <- rownames(deg)
  deg$max_pct <- pmax(deg$pct.1, deg$pct.2)
  deg$significant <- ifelse(
    deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 1,
    "Significant",
    "Not Significant"
  )
  
  deg$expression_category <- case_when(
    deg$max_pct < 0.1 ~ "Rare (<10% cells)",
    deg$max_pct < 0.3 ~ "Low (10-30% cells)",
    deg$max_pct < 0.5 ~ "Medium (30-50% cells)",
    TRUE ~ "High (>50% cells)"
  )
  
  return(deg)
}


#' Filter DEG results to keep only genes with consistent direction 
#' compared to bulk-derived driver gene definition.
#' @param deg_df DEG result from Seurat (must have 'gene' column)
#' @param direction_ref Reference data frame with columns: gene_name, is_upregulated, is_downregulated
#' @return Filtered deg_df with consistent direction and significance
filter_consistent_direction <- function(deg_df, direction_ref) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  # Ensure deg_df has 'gene' column (as created in run_de_analysis)
  deg_df <- deg_df %>%
    left_join(direction_ref, by = c("gene" = "gene_name"))
  
  # Check for missing direction info (should not happen, but safe)
  if (any(is.na(deg_df$is_upregulated))) {
    warning("Some genes missing direction info; they will be excluded.")
    deg_df <- deg_df[!is.na(deg_df$is_upregulated), ]
  }
  
  # Define expected direction:
  # - If bulk says upregulated → require avg_log2FC > 0
  # - If bulk says downregulated → require avg_log2FC < 0
  consistent <- with(deg_df, {
    (is_upregulated & avg_log2FC > 0) |
      (is_downregulated & avg_log2FC < 0)
  })
  
  # Keep only significant AND direction-consistent genes
  filtered <- deg_df %>%
    filter(
      p_val_adj < 0.05,
      abs(avg_log2FC) > 1,
      consistent
    )
  
  return(filtered)
}


#' Generate a volcano plot for DEG results
#'
#' @param deg_df Data frame of DEG results (must have avg_log2FC and p_val_adj columns)
#' @param title Title of the plot
#' @return ggplot object of the volcano plot
create_volcano_plot <- function(deg_df, title = "Driver Genes: Ncf2_KO vs WT") {
  ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Significant" = "#DC143C", "Not Significant" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
    labs(
      title = title,
      x = "log2 Fold Change (KO / WT)",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_classic() +
    theme(legend.position = "right")
}


#' Plot average expression heatmap for single-cell data
#'
#' @param seurat_obj A Seurat object
#' @param genes A vector of genes to display
#' @param group.by Primary grouping variable (typically cell type), default "celltype"
#' @param split.by Secondary grouping variable for comparison, default "orig.ident"
#' @param assay Assay to use, default "RNA"
#' @param slot Data slot to use, default "data" (log-normalized)
#' @param scale_method Scaling method: "zscore" (row-wise Z-score) or "none", default "zscore"
#' @param cluster_rows Whether to cluster rows (genes), default TRUE
#' @param cluster_cols Whether to cluster columns (cell type × condition), default TRUE
#' @param show_rownames Whether to display gene names, default TRUE
#' @param show_colnames Whether to display column names, default TRUE
#' @param fontsize Base font size, default 10
#' @param fontsize_row Font size for row labels, default 9
#' @param fontsize_col Font size for column labels, default 9
#' @param color_low Color for low expression, default "#2A52BE" (blue)
#' @param color_mid Color for mid expression, default "white"
#' @param color_high Color for high expression, default "#D73027" (red)
#' @param color_breaks Color breakpoints, default c(-2, 2); values outside are saturated
#' @param annotation_colors A list of annotation colors, default NULL (auto-generated)
#' @param plot_width Width of the output PDF, default 14
#' @param plot_height Height of the output PDF, default 10
#' @param output_file Output file path, default NULL (no file saved)
#' @param return_data Whether to return the heatmap matrix data, default FALSE
#' @param main Heatmap title, default "Average Expression Heatmap"
#' 
#' @return A pheatmap object or a list (if return_data = TRUE)
plot_heatmap_sc <- function(
    seurat_obj,
    genes,
    group.by = "celltype",
    split.by = "orig.ident",
    assay = "RNA",
    slot = "data",
    scale_method = "zscore",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 9,
    color_low = "#2A52BE",
    color_mid = "white",
    color_high = "#D73027",
    color_breaks = c(-2, 2),
    annotation_colors = NULL,
    plot_width = 14,
    plot_height = 10,
    output_file = NULL,
    return_data = FALSE,
    main = "Average Expression Heatmap"
) {
  
  # ========== Load required packages ==========
  require(pheatmap)
  require(dplyr)
  require(tidyr)
  
  # ========== Argument validation ==========
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  
  if (length(genes) == 0) {
    stop("`genes` must not be empty.")
  }
  
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "does not exist in the Seurat object."))
  }
  
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Grouping variable", group.by, "was not found in meta.data."))
  }
  
  if (!split.by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Splitting variable", split.by, "was not found in meta.data."))
  }
  
  # ========== Data extraction ==========
  cat("Extracting expression data...\n")
  
  if (slot == "data") {
    expr_mat <- seurat_obj[[assay]]@data
  } else if (slot == "counts") {
    expr_mat <- seurat_obj[[assay]]@counts
  } else {
    stop("`slot` must be either 'data' or 'counts'.")
  }
  
  # Validate genes
  valid_genes <- intersect(genes, rownames(expr_mat))
  missing_genes <- setdiff(genes, rownames(expr_mat))
  
  if (length(missing_genes) > 0) {
    warning(paste(
      "The following genes were not found:",
      paste(head(missing_genes, 10), collapse = ", ")
    ))
  }
  
  if (length(valid_genes) == 0) {
    stop("No valid genes were found in the expression matrix.")
  }
  
  cat(paste(
    "Using", length(valid_genes),
    "genes (out of", length(genes), "input genes)\n"
  ))
  
  expr_mat <- expr_mat[valid_genes, , drop = FALSE]
  
  # Retrieve metadata
  meta <- seurat_obj@meta.data
  meta$cell_id <- rownames(meta)
  
  # ========== Compute average expression ==========
  cat("Computing average expression per group...\n")
  
  expr_df <- t(expr_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id")
  
  expr_with_meta <- expr_df %>%
    inner_join(
      meta %>% select(cell_id, all_of(c(group.by, split.by))),
      by = "cell_id"
    )
  
  avg_expr <- expr_with_meta %>%
    group_by(across(all_of(c(group.by, split.by)))) %>%
    summarise(
      across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  avg_long <- avg_expr %>%
    pivot_longer(
      cols = all_of(valid_genes),
      names_to = "gene",
      values_to = "expr"
    )
  
  avg_long <- avg_long %>%
    unite(
      col = "col_group",
      all_of(c(group.by, split.by)),
      sep = "_",
      remove = FALSE
    )
  
  avg_wide <- avg_long %>%
    select(gene, col_group, expr) %>%
    pivot_wider(names_from = col_group, values_from = expr)
  
  mat <- avg_wide %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  mat[is.na(mat)] <- 0
  
  cat(paste(
    "Heatmap matrix dimensions:",
    nrow(mat), "genes ×", ncol(mat), "groups\n"
  ))
  
  # ========== Scaling ==========
  if (scale_method == "zscore") {
    cat("Applying row-wise Z-score normalization...\n")
    mat_scaled <- t(scale(t(mat)))
    mat_scaled[is.na(mat_scaled)] <- 0
  } else {
    mat_scaled <- mat
  }
  
  # ========== Column annotations ==========
  cat("Preparing column annotations...\n")
  
  col_groups <- colnames(mat_scaled)
  split_pattern <- strsplit(col_groups, "_")
  
  conditions <- sapply(split_pattern, function(x) x[length(x)])
  groups <- sapply(split_pattern, function(x) {
    if (length(x) > 1) {
      paste(x[1:(length(x) - 1)], collapse = "_")
    } else {
      x[1]
    }
  })
  
  ann_df <- data.frame(
    row.names = col_groups,
    stringsAsFactors = FALSE
  )
  
  ann_df[[split.by]] <- factor(conditions)
  
  # ========== Annotation colors ==========
  if (is.null(annotation_colors)) {
    split_levels <- levels(ann_df[[split.by]])
    
    if (length(split_levels) == 2) {
      split_colors <- c("#4E79A7", "#F28E2B")
      names(split_colors) <- split_levels
    } else {
      split_colors <- scales::hue_pal()(length(split_levels))
      names(split_colors) <- split_levels
    }
    
    annotation_colors <- list()
    annotation_colors[[split.by]] <- split_colors
  }
  
  # ========== Color scale ==========
  color_palette <- colorRampPalette(
    c(color_low, color_mid, color_high)
  )(50)
  
  if (scale_method == "zscore") {
    breaks <- c(
      -Inf,
      seq(color_breaks[1], color_breaks[2], length.out = 49),
      Inf
    )
  } else {
    breaks <- NULL
  }
  
  # ========== Draw heatmap ==========
  cat("Drawing heatmap...\n")
  
  p <- pheatmap(
    mat = mat_scaled,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    color = color_palette,
    breaks = breaks,
    annotation_col = ann_df,
    annotation_colors = annotation_colors,
    annotation_names_row = FALSE,
    annotation_names_col = TRUE,
    border_color = NA,
    main = main,
    silent = TRUE
  )
  
  # ========== Save output ==========
  if (!is.null(output_file)) {
    cat(paste("Saving to:", output_file, "\n"))
    
    pdf(output_file, width = plot_width, height = plot_height)
    print(p)
    dev.off()
    
    # matrix_file <- sub("\\.pdf$", "_matrix.csv", output_file)
    # write.csv(mat_scaled, matrix_file)
    # cat(paste("Matrix data saved to:", matrix_file, "\n"))
  }
  
  # ========== Return results ==========
  if (return_data) {
    return(list(
      plot = p,
      matrix_raw = mat,
      matrix_scaled = mat_scaled,
      annotation = ann_df,
      genes_used = valid_genes,
      genes_missing = missing_genes
    ))
  } else {
    return(p)
  }
}


#' Draw violin plots for single-cell expression data
#'
#' @param seurat_obj A Seurat object
#' @param genes A vector of gene names (character vector)
#' @param group.by Grouping variable, typically cell type (default: "celltype")
#' @param split.by Splitting variable for group comparison (default: "orig.ident")
#' @param colors A vector of colors corresponding to levels in split.by
#' @param ncol Number of columns in the combined plot (default: 3)
#' @param assay Assay to use (default: "RNA")
#' @param slot Data slot to use: "data" (normalized) or "counts" (raw) (default: "data")
#' @param add_stats Whether to add statistical significance annotations (default: FALSE)
#' @param show_points Whether to display individual data points (default: FALSE)
#' @param point_size Size of individual data points (default: 0.1)
#' @param plot_width Width of the output PDF (default: 18)
#' @param plot_height Height of the output PDF (default: NULL, auto-calculated)
#' @param output_file Output file path (default: NULL, no file saved)
#' @param return_data Whether to return the plotting data (default: FALSE)
#' 
#' @return A ggplot object or a list containing plots and data
plot_violin_sc <- function(
    seurat_obj,
    genes,
    group.by = "celltype",
    split.by = "orig.ident",
    colors = NULL,
    ncol = 3,
    assay = "RNA",
    slot = "data",
    add_stats = FALSE,
    show_points = FALSE,
    point_size = 0.1,
    plot_width = 18,
    plot_height = NULL,
    output_file = NULL,
    return_data = FALSE
) {
  
  # Load required packages
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # ========== Argument validation ==========
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  
  if (length(genes) == 0) {
    stop("`genes` must not be empty.")
  }
  
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "does not exist in the Seurat object."))
  }
  
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Grouping variable", group.by, "was not found in meta.data."))
  }
  
  if (!split.by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Splitting variable", split.by, "was not found in meta.data."))
  }
  
  # ========== Data extraction ==========
  cat("Extracting expression data...\n")
  
  # Retrieve expression matrix
  if (slot == "data") {
    expr_matrix <- seurat_obj[[assay]]@data
  } else if (slot == "counts") {
    expr_matrix <- seurat_obj[[assay]]@counts
  } else {
    stop("`slot` must be either 'data' or 'counts'.")
  }
  
  # Retrieve metadata
  meta_data <- seurat_obj@meta.data
  
  # Check gene availability
  genes_exist <- genes %in% rownames(expr_matrix)
  if (!all(genes_exist)) {
    missing <- genes[!genes_exist]
    warning(paste("The following genes were not found and will be skipped:",
                  paste(missing, collapse = ", ")))
    genes <- genes[genes_exist]
  }
  
  if (length(genes) == 0) {
    stop("No valid genes available for plotting.")
  }
  
  cat(paste("Preparing to plot", length(genes), "genes.\n"))
  
  # ========== Color setup ==========
  split_levels <- unique(meta_data[[split.by]])
  
  if (is.null(colors)) {
    # Default color scheme
    if (length(split_levels) == 2) {
      colors <- c("#4DBBD5", "#E64B35")
      names(colors) <- split_levels
    } else {
      colors <- scales::hue_pal()(length(split_levels))
      names(colors) <- split_levels
    }
  } else if (is.null(names(colors))) {
    # Automatically assign names if not provided
    names(colors) <- split_levels[1:length(colors)]
  }
  
  # ========== Generate plots for each gene ==========
  plot_list <- list()
  data_list <- list()
  
  for (gene in genes) {
    cat(paste("Plotting:", gene, "\n"))
    
    # Prepare data frame
    plot_df <- data.frame(
      expression = as.numeric(expr_matrix[gene, ]),
      group = factor(meta_data[[group.by]]),
      condition = factor(meta_data[[split.by]])
    )
    
    # Remove missing values
    plot_df <- plot_df[complete.cases(plot_df), ]
    
    # Base plot
    p <- ggplot(plot_df, aes(x = group, y = expression, fill = condition)) +
      geom_violin(
        position = position_dodge(0.9),
        scale = "width",
        trim = TRUE,
        alpha = 0.8
      ) +
      scale_fill_manual(values = colors) +
      labs(
        title = gene,
        x = "",
        y = "Expression Level",
        fill = split.by
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey90", size = 0.3)
      )
    
    # Add jittered points
    if (show_points) {
      p <- p + geom_jitter(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
        size = point_size,
        alpha = 0.3
      )
    }
    
    # Add statistical significance annotations
    if (add_stats) {
      if (requireNamespace("ggpubr", quietly = TRUE)) {
        p <- p + ggpubr::stat_compare_means(
          aes(group = condition),
          method = "wilcox.test",
          label = "p.signif",
          hide.ns = FALSE
        )
      } else {
        warning("Package 'ggpubr' is required to add statistical annotations.")
      }
    }
    
    plot_list[[gene]] <- p
    
    if (return_data) {
      data_list[[gene]] <- plot_df
    }
  }
  
  # ========== Combine plots ==========
  cat("Combining plots...\n")
  
  n_genes <- length(plot_list)
  nrow_plots <- ceiling(n_genes / ncol)
  
  combined_plot <- wrap_plots(plot_list, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # ========== Save output ==========
  if (!is.null(output_file)) {
    # Automatically calculate plot height
    if (is.null(plot_height)) {
      plot_height <- 4 * nrow_plots + 1
    }
    
    cat(paste("Saving plot to:", output_file, "\n"))
    ggsave(
      filename = output_file,
      plot = combined_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300
    )
  }
  
  # ========== Return results ==========
  if (return_data) {
    return(list(
      plot = combined_plot,
      data = data_list,
      individual_plots = plot_list
    ))
  } else {
    return(combined_plot)
  }
}


#' Faceted violin plot by cell type
#' 
#' @param seurat_obj A Seurat object
#' @param gene A single gene name
#' @param group.by Grouping variable (typically cell type)
#' @param split.by Splitting variable for comparison
#' @param colors A named vector of colors
#' @param add_stats Whether to add statistical significance annotations (default: FALSE)
#' @param output_file Output file path
#' 
#' @return A ggplot object
plot_violin_facet <- function(
    seurat_obj,
    gene,
    group.by = "celltype",
    split.by = "orig.ident",
    colors = c("WT" = "#4DBBD5", "Ncf2_KO" = "#E64B35"),
    add_stats = FALSE,
    output_file = NULL
) {
  
  require(ggplot2)
  require(dplyr)

  # Extract data
  expr_data <- seurat_obj[["RNA"]]@data
  meta_data <- seurat_obj@meta.data
  
  # Check gene availability
  if (!gene %in% rownames(expr_data)) {
    stop(paste("Gene", gene, "was not found in the expression matrix."))
  }
  
  # Prepare data
  plot_df <- data.frame(
    expression = as.numeric(expr_data[gene, ]),
    celltype = factor(meta_data[[group.by]]),
    condition = factor(meta_data[[split.by]])
  )
  plot_df <- plot_df[complete.cases(plot_df), ] # remove NAs

  # Base plot
  p <- ggplot(plot_df, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(add = c(0, 0.6))) +
    labs(
      title = gene,
      x = "",
      y = "Expression Level"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
    
  # Add statistical significance annotations (per facet)
  if (add_stats) {
    if (requireNamespace("ggpubr", quietly = TRUE)) {
      conditions <- levels(plot_df$condition) #conditions = c("WT", "Ncf2_KO")
      p <- p + ggpubr::stat_compare_means(
        method = "wilcox.test",
        comparisons = list(conditions), #conditions = c("WT", "Ncf2_KO")
        label = "p.signif",
        hide.ns = FALSE      
      )
    } else {
      warning("Package 'ggpubr' is required to add statistical annotations.")
    }
  }
  
  # Save output
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 16, height = 12, dpi = 300)
  }
  
  return(p)
}


# violin plot overall, no cell type separation
#' 
#' @param seurat_obj A Seurat object
#' @param genes A vector of gene names
#' @param group.by Grouping variable for comparison (default: "orig.ident")
#' @param colors A named vector of colors
#' @param ncol Number of columns in the combined plot
#' @param add_stats Whether to perform and display statistical tests
#' @param show_points Whether to display individual data points
#' @param output_file Output file path
#' 
#' @return A ggplot object
plot_violin_overall <- function(
    seurat_obj,
    genes,
    group.by = "orig.ident",
    colors = c("WT" = "#4DBBD5", "Ncf2_KO" = "#E64B35"),
    ncol = 3,
    assay = "RNA",
    slot = "data",
    add_stats = TRUE,
    show_points = FALSE,
    point_size = 0.1,
    plot_width = 12,
    plot_height = NULL,
    output_file = NULL
) {
  
  require(ggplot2)
  require(patchwork)
  
  # Extract expression data
  if (slot == "data") {
    expr_matrix <- seurat_obj[[assay]]@data
  } else {
    expr_matrix <- seurat_obj[[assay]]@counts
  }
  
  meta_data <- seurat_obj@meta.data
  
  # Check gene availability
  genes_exist <- genes %in% rownames(expr_matrix)
  if (!all(genes_exist)) {
    warning(paste(
      "The following genes were not found and will be skipped:",
      paste(genes[!genes_exist], collapse = ", ")
    ))
    genes <- genes[genes_exist]
  }
  
  if (length(genes) == 0) {
    stop("No valid genes available for plotting.")
  }
  
  # Generate plots for each gene
  plot_list <- list()
  
  for (gene in genes) {
    cat(paste("Plotting:", gene, "\n"))
    
    plot_df <- data.frame(
      expression = as.numeric(expr_matrix[gene, ]),
      group = factor(meta_data[[group.by]])
    )
    
    # Compute statistics
    groups <- levels(plot_df$group)
    if (length(groups) == 2) {
      expr1 <- plot_df$expression[plot_df$group == groups[1]]
      expr2 <- plot_df$expression[plot_df$group == groups[2]]
      
      # Wilcoxon rank-sum test
      test_result <- wilcox.test(expr1, expr2)
      p_value <- test_result$p.value
      
      # Compute log2 fold change
      mean1 <- mean(expr1[expr1 > 0])
      mean2 <- mean(expr2[expr2 > 0])
      log2fc <- log2((mean2 + 0.01) / (mean1 + 0.01))
      
      # Construct subtitle
      if (p_value < 0.001) {
        p_text <- "p < 0.001"
      } else {
        p_text <- sprintf("p = %.3f", p_value)
      }
      
      subtitle <- sprintf("log2FC = %.2f, %s", log2fc, p_text)
    } else {
      subtitle <- ""
    }
    
    # Base plot
    p <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = colors) +
      labs(
        title = gene,
        subtitle = subtitle,
        x = "",
        y = "Expression Level"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90", size = 0.3)
      )
    
    # Add jittered points
    if (show_points) {
      p <- p + geom_jitter(width = 0.2, size = point_size, alpha = 0.3)
    }
    
    # Add statistical significance annotation
    if (add_stats && requireNamespace("ggpubr", quietly = TRUE)) {
      p <- p + ggpubr::stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = list(groups),
        hide.ns = FALSE
      )
    }
    
    plot_list[[gene]] <- p
  }
  
  # Combine plots
  n_genes <- length(plot_list)
  nrow_plots <- ceiling(n_genes / ncol)
  
  combined_plot <- wrap_plots(plot_list, ncol = ncol)
  
  # Save output
  if (!is.null(output_file)) {
    if (is.null(plot_height)) {
      plot_height <- 4 * nrow_plots + 1
    }
    
    ggsave(
      output_file,
      combined_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300
    )
    cat(paste("Saved to:", output_file, "\n"))
  }
  
  return(combined_plot)
}


# =============================================================================
# 1. Load scRNA seq Data
# =============================================================================
message("Loading single-cell dataset...")
load(SCDATA_PATH)  # Loads LUNG tissue data
cat("Dataset overview:\n")
cat("  Total cells:", ncol(sce), "\n")
cat("  Sample groups:\n")
print(table(sce$orig.ident))
cat("  Cell types:\n")
print(table(sce$celltype))

# =============================================================================
# 2. Load Driver Genes
# =============================================================================
message("Loading driver gene list...")
load(DRIVER_PATH)  # Loads 'results' object

gene_records <- list()
for (group_name in names(results$driver_genes)) {
  df <- results$driver_genes[[group_name]]
  subset_df <- data.frame(
    gene_name = df$gene_name,
    is_upregulated = df$is_upregulated,
    is_downregulated = df$is_downregulated,
    group = group_name,
    stringsAsFactors = FALSE
  )
  gene_records[[group_name]] <- subset_df
}
all_genes_df <- do.call(rbind, gene_records) 

summary_df <- all_genes_df %>%
  group_by(gene_name) %>%
  summarise(
    is_upregulated = first(is_upregulated),
    is_downregulated = first(is_downregulated),
    occurance = paste(group, collapse = ", "),
    .groups = "drop"
  )

driver_genes <- unique(unlist(lapply(results$driver_genes, function(df) df$gene_name)))
genes_in_data <- driver_genes[driver_genes %in% rownames(sce)]
summary_df$in_scRNA <- ifelse(summary_df$gene_name %in% genes_in_data, "TRUE", "FALSE")

# Create a direction lookup table for fast matching
direction_lookup <- summary_df %>%
  select(gene_name, is_upregulated, is_downregulated)

cat("Driver gene coverage:\n")
cat("  Total driver genes:", length(driver_genes), "\n")
cat("  Found in scRNA-seq data:", length(genes_in_data), "\n")
cat("  Missing:", length(driver_genes) - length(genes_in_data), "\n")

save_csv(summary_df, OUTPUT_FILES$genes_info_csv)

# =============================================================================
# 3. Subset to WT and Ncf2_KO Samples
# =============================================================================
sce_subset <- subset(sce, subset = orig.ident %in% c("WT", "Ncf2_KO"))
cat("After subsetting to WT/Ncf2_KO:\n")
print(table(sce_subset$orig.ident))

# =============================================================================
# 4. Global Differential Expression Analysis
# =============================================================================
message("Running global DEG analysis (across all cell types)...")

Idents(sce_subset) <- sce_subset$orig.ident
DefaultAssay(sce_subset) <- "RNA"

deg_overall <- run_de_analysis(
  object = sce_subset,
  features = genes_in_data,
  min_pct = 0.01  # Very permissive to capture cell-type-specific signals
)

if (!is.null(deg_overall)) {

  sig_overall <- filter_consistent_direction(
      deg_df = deg_overall,
      direction_ref = direction_lookup
      )

  save_csv(sig_overall, OUTPUT_FILES$overall_sig_csv)
  
  cat("Global DEG results:\n")
  cat("  Significant genes:", nrow(sig_overall), "\n")
  cat("    Up in KO:", sum(sig_overall$avg_log2FC > 0), "\n")
  cat("    Down in KO:", sum(sig_overall$avg_log2FC < 0), "\n")
} else {
  stop("No genes passed filtering in global DEG analysis.")
}

# Volcano plot
p_volcano <- create_volcano_plot(deg_overall)
ggsave(OUTPUT_FILES$volcano_pdf, p_volcano, width = 8, height = 6)

# =============================================================================
# 5. Cell-Type-Specific Differential Expression
# =============================================================================
message("Running cell-type-specific DEG analysis...")

Idents(sce_subset) <- sce_subset$celltype
cell_types <- levels(Idents(sce_subset))
celltype_deg_list <- list()

for (ct in cell_types) {
  cells_ct <- subset(sce_subset, idents = ct)
  group_counts <- table(cells_ct$orig.ident)
  group_counts <- group_counts[group_counts > 0]
  
  if (length(group_counts) < 2 || min(group_counts) < 10) {
    message("Skipping ", ct, ": insufficient cells (", paste(names(group_counts), group_counts, collapse = ", "), ")")
    next
  }
  
  message("Analyzing ", ct, " (WT=", group_counts["WT"], ", KO=", group_counts["Ncf2_KO"], ")")
  
  Idents(cells_ct) <- cells_ct$orig.ident
  
  tryCatch({
    deg_ct <- run_de_analysis(
      object = cells_ct,
      features = genes_in_data,
      min_pct = 0.5  # Require expression in ≥50% of cells within the type
    )
    
    if (!is.null(deg_ct)) {
      deg_ct_consistent <- filter_consistent_direction(
        deg_df = deg_ct,
        direction_ref = direction_lookup
      )
      
      if (!is.null(deg_ct_consistent) && nrow(deg_ct_consistent) > 0) {
        deg_ct_consistent$celltype <- ct
        deg_ct_consistent$n_cells_WT <- group_counts["WT"]
        deg_ct_consistent$n_cells_KO <- group_counts["Ncf2_KO"]
        
        celltype_deg_list[[ct]] <- deg_ct_consistent
        
        n_sig <- nrow(deg_ct_consistent)
        message("  Completed: ", nrow(deg_ct), " genes tested, ", n_sig, " significant & direction-consistent")
      } else {
        message("  No direction-consistent significant genes in ", ct)
      }
    }
  }, error = function(e) {
    message("  Error in ", ct, ": ", e$message)
  })
}

# Combine and summarize cell-type results
if (length(celltype_deg_list) > 0) {
  # significant genes by cell type
  sig_by_celltype <- bind_rows(celltype_deg_list)
  save_csv(sig_by_celltype, OUTPUT_FILES$celltype_sig_csv)
  
  # summary table
  sig_summary <- sig_by_celltype %>%
    group_by(celltype) %>%
    summarise(
      n_sig_genes = n(),
      n_upregulated = sum(avg_log2FC > 0),
      n_downregulated = sum(avg_log2FC < 0),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_genes))
  save_csv(sig_summary, OUTPUT_FILES$summary_csv)
  
  cat("\nCell-type-specific DEG summary:\n")
  print(sig_summary)
  
  # gene to cell type mapping
  gene_map <- sig_by_celltype %>%
  group_by(gene) %>%
  summarise(
    n_celltypes_sig = n(),
    celltypes = paste(celltype, collapse = "; "),
    max_log2FC = max(abs(avg_log2FC)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_celltypes_sig), desc(max_log2FC))
  save_csv(gene_map, OUTPUT_FILES$gene_celltype_map)
} else {
  sig_by_celltype <- NULL
  cat("No cell types met the minimum cell count criteria.\n")
}


# =============================================================================
# 6. Heatmap: Average Expression by Cell Type and Group
# =============================================================================
message("Generating heatmap of top driver genes...")

# Get top 50 direction-consistent significant genes by |log2FC|
top_genes <- if (nrow(sig_overall) > 0) {
  sig_overall %>% arrange(desc(abs(avg_log2FC))) %>% pull(gene) %>% head(50)
} else {
  character(0)
}
if (length(top_genes) == 0) stop("No significant driver genes to plot.")
# Ensure genes exist in RNA assay
valid_genes <- intersect(top_genes, rownames(sce_subset[["RNA"]]))
if (length(valid_genes) == 0) stop("None of top genes found in RNA assay.")

plot_heatmap_sc(
  seurat_obj = sce_subset,
  genes = top_genes,
  main = "Overall Top Driver Genes top50",
  output_file = OUTPUT_FILES$heatmap_pdf,
)

# =============================================================================
# 7. Dot Plot: Top Genes by Cell Type
# =============================================================================
top_dot_genes <- if (nrow(sig_overall) > 0) {
  sig_overall %>% arrange(desc(abs(avg_log2FC))) %>% pull(gene) %>% head(50)
} else {
  character(0)
}
Idents(sce_subset) <- sce_subset$celltype

p_dot <- DotPlot(
  sce_subset,
  features = top_dot_genes,
  split.by = "orig.ident",
  cols = c("WT" = "#4DBBD5", "Ncf2_KO" = "#E64B35")
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top Driver Genes Expression by Cell Type")

ggsave(OUTPUT_FILES$dotplot_pdf, p_dot, width = 12, height = 12)

# =============================================================================
# 8. Violin Plots: Top Cell-Type-Specific or Global Genes
# =============================================================================
target_genes <- sig_overall %>%
  arrange(desc(abs(avg_log2FC))) %>%
  pull(gene) %>%
  unique() %>%
  head(9)

# target_genes <- sig_by_celltype %>% filter(celltype == "Macrophage") %>% arrange(desc(abs(avg_log2FC))) %>% pull(gene) %>% head(6)

p1 <- plot_violin_overall(
  sce_subset,
  genes = target_genes,
  ncol = 3,
  add_stats = TRUE,
  output_file = OUTPUT_FILES$violin_pdf_overall,
)

p2 <- plot_violin_sc(
  sce_subset,
  genes = target_genes,
  ncol = 3,
  add_stats = TRUE,
  output_file = OUTPUT_FILES$violin_pdf_ct,
)

p3 <- plot_violin_facet(
  sce_subset,
  gene = target_genes[1],
  add_stats = TRUE,
  output_file = OUTPUT_FILES$violin_pdf_ct_fc
)

# =============================================================================
# 9. Generate Summary Report
# =============================================================================
report_header <- c(
  "==============================================",
  "        Driver Gene Analysis Summary",
  "==============================================",
  "",
  "1. Dataset Overview:",
  paste("   - Total cells:", ncol(sce_subset)),
  paste("   - WT cells:", sum(sce_subset$orig.ident == "WT")),
  paste("   - Ncf2_KO cells:", sum(sce_subset$orig.ident == "Ncf2_KO")),
  paste("   - Cell types analyzed:", length(unique(sce_subset$celltype))),
  "",
  "2. Gene Coverage:",
  paste("   - Submitted driver genes:", length(driver_genes)),
  paste("   - Detected in data:", length(genes_in_data)),
  paste("   - Missing:", length(driver_genes) - length(genes_in_data)),
  "",
  "3. Global DEG Results (|log2FC| > 1, adj.P.Val < 0.05):",
  paste("   - Significant genes:", nrow(sig_overall)),
  paste("     - Up in KO:", sum(sig_overall$avg_log2FC > 0)),
  paste("     - Down in KO:", sum(sig_overall$avg_log2FC < 0)),
  "",
  "4. Cell-Type-Specific Analysis:",
  paste("   - Cell types analyzed:", length(celltype_deg_list)),
  paste("   - Total cell-type-specific significant DEGs:", if (!is.null(sig_by_celltype)) nrow(sig_by_celltype) else 0),
  "   - Overview of cell types with significant DEGs:"
)
# sig_summary
sig_table_lines <- paste0("     ", capture.output(print(sig_summary)))
report_footer <- c(
  "",
  "5. Key Outputs:",
  "   Global:",
  "     - 01.Overall_sig_KO_vs_WT.csv",
  "     - volcano_plot_driver_genes.pdf",
  "   By Cell Type:",
  "     - 02.Celltype_all_KO_vs_WT.csv",
  "     - 02.Celltype_sig_KO_vs_WT.csv",
  "     - celltype_deg_summary.csv",
  "   Visualizations:",
  "     - heatmap_driver_genes_by_celltype.pdf",
  "     - dotplot_top_driver_genes.pdf",
  "     - violin_plot_*.pdf",
  "   Other:",
  "     - missing_genes.txt",
  "",
  "6. Notes:",
  "   - log2FC > 0: higher in Ncf2_KO; log2FC < 0: higher in WT",
  "   - Global analysis uses min.pct=0.01 to retain rare signals",
  "   - Cell-type analysis requires ≥50% expression within type",
  "",
  "7. Gene–Cell Type Mapping:",
  paste("   - Genes significant in multiple cell types:", sum(gene_map$n_celltypes_sig > 1)),
  paste("   - Genes specific to one cell type:", sum(gene_map$n_celltypes_sig == 1)),
  "   - See: gene_to_celltype_mapping.csv",
  "=============================================="
)
report_lines <- c(report_header, sig_table_lines, report_footer)
writeLines(report_lines, con = OUTPUT_FILES$analysis_summary_txt)

message("Analysis completed successfully! See output files in: ", BASE_DIR)