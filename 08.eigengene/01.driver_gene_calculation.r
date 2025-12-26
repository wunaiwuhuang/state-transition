# ============================================================================
# Driver Gene Selection for State-Transition Model
# Based on Eigengene Geometry and Disease Contribution
# ============================================================================

library(tidyverse)
library(rtracklayer)

BASE_DIR <- "/disk2/cai113/data/stateTrans"

# ============================================================================
# PART 1: Data Loading and Preparation
# ============================================================================

#' Load and Prepare All Required Data
#' 
#' @param comparisons Character vector of DEG comparisons to load
#' @return List containing all necessary data objects
load_project_data <- function(comparisons = c("c2_vs_c1", "c3_vs_c1", "c3_vs_c2")) {
  
  # 1. Load PCA results and sample info
  load(file.path(BASE_DIR, "05.model/sample_info_with_PCs.rdata"))
  
  # 2. Load VST normalized expression matrix
  vst_mat <- as.matrix(read.delim(
    file.path(BASE_DIR, "04.standardization/03.vst_matrix.txt"),
    row.names = 1, check.names = FALSE
  ))
  
  # 3. Load PCA object to get loadings (eigengenes)
  # Re-run PCA to get rotation matrix
  X <- t(vst_mat)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  pca_result <- prcomp(X_centered, center = FALSE, scale. = FALSE)
  
  # 4. Load gene annotations from GTF
  gtf_file <- file.path(BASE_DIR, "02.alignment/STAR_ref/mouse_GRCm39/Mus_musculus.GRCm39.110.gtf")
  gene_annotations <- parse_gtf_annotations(gtf_file)
  
  # 5. Load DEG results
  deg_results <- load_deg_results(comparisons = comparisons)
  
  # 6. Load model results (critical points)
  model_file <- file.path(BASE_DIR, "05.model/model_results.rdata")
  if (file.exists(model_file)) {
    load(model_file)  # Should contain critical points c1, c2, c3, etc.
    critical_points <- model_results$critical_points[c("c1","c2","c3","c1_star")]
  } else {
    warning("Model results not found. Please ensure model is built.")
    critical_points <- NULL
  }
  
  return(list(
    vst_mat = vst_mat,
    sample_info_pcs = sample_info_pcs,
    x_cgd = x_cgd,
    pca_result = pca_result,
    gene_annotations = gene_annotations,
    deg_results = deg_results,
    critical_points = critical_points
  ))
}


#' Parse GTF File to Extract Gene Annotations
#' 
#' @param gtf_file Path to GTF file
#' @return Data frame with gene annotations
parse_gtf_annotations <- function(gtf_file) {
  
  cat("Parsing GTF annotations...\n")
  gtf <- import(gtf_file)
  
  # Extract gene-level information
  gene_info <- as.data.frame(gtf) %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_biotype, 
                  seqnames, start, end, strand) %>%
    distinct()
  
  cat("Parsed", nrow(gene_info), "genes\n")
  
  return(gene_info)
}

#' Load DEG Results from Multiple Comparisons
#'
#' @param deg_dir Root directory of the DEG project
#' @param comparisons Character vector of comparisons
#' @param filename DEG result file name
#'
#' @return Named list of DEG data frames
load_deg_results <- function(deg_dir = file.path(BASE_DIR, "06.DEG_analysis/DEG_results"), comparisons = c("c3_vs_c1", "c2_vs_c1", "c3_vs_c2"), filename = "DEG_significant.txt") {
  cat("Loading DEG results...\n")
  deg_list <- lapply(comparisons, function(comp) {
    file_path <- file.path(deg_dir, comp, filename)
    
    if (!file.exists(file_path)) {
      stop("DEG file not found: ", file_path)
    }
    
    read.delim(
      file_path,
      row.names   = 1,
      check.names = FALSE
    )
  })

  names(deg_list) <- comparisons
  return(deg_list)
}


# ============================================================================
# PART 2: Eigengene Analysis
# ============================================================================

#' Calculate Eigengene Vectors for All Genes
#' 
#' @param pca_result PCA object from prcomp
#' @return Data frame with eigengene coordinates
calculate_eigengenes <- function(pca_result) {
  
  # Get loadings (rotation matrix): genes x PCs
  loadings <- pca_result$rotation
  
  # Extract first two PCs (disease axis and orthogonal axis)
  # Note: PC1 was negated in your analysis, so adjust accordingly
  eigengenes <- data.frame(
    gene_id = rownames(loadings),
    v1 = loadings[, "PC1"],      # PC1 loading
    v2 = loadings[, "PC2"],      # PC2 loading
    cgd_loading = -loadings[, "PC1"]  # Disease axis (negated PC1)
  )
  
  # Calculate polar coordinates
  eigengenes <- eigengenes %>%
    mutate(
      radius = sqrt(v1^2 + v2^2),
      angle_rad = atan2(v2, -v1),  # Use -v1 because x_cgd = -PC1
      angle_deg = angle_rad * 180 / pi
    )
  #        PC2 (v2)
  #         ↑ 90°
  #         |
  #-180° ←--+--→ 0° (-v1, 疾病轴方向)
  #         |
  #         ↓ -90°
  return(eigengenes)
}


#' Define Disease Cone Angle
#' 
#' @param x_cgd Disease axis values for samples
#' @param sample_info_pcs Sample information with PC scores
#' @param critical_points List with c2 value
#' @return Angle range defining disease-associated genes
define_disease_cone <- function(x_cgd, sample_info_pcs, critical_points) {
  
  # --- Step 1: Select disease samples ---
  ko_samples <- sample_info_pcs$genotype == "Ncf2"
  if (!is.null(critical_points)) {
    disease_samples <- ko_samples & (x_cgd < critical_points$c2)
  } else {
    disease_samples <- ko_samples & (x_cgd < quantile(x_cgd[ko_samples], 0.25))
  }
  
  if (sum(disease_samples) == 0) stop("No disease samples selected.")
  
  # --- Step 2: Compute angles in your coordinate system ---
  # x-axis = -PC1 (disease direction), y-axis = PC2
  # So angle = atan2(y, x) = atan2(PC2, -PC1)
  pc1_neg <- -sample_info_pcs$PC1[disease_samples]
  pc2_val <-  sample_info_pcs$PC2[disease_samples]
  angles_raw <- atan2(pc2_val, pc1_neg) * 180 / pi  # in [-180, 180]
  
  # --- Step 3: Find minimal covering arc on circle ---
  angles_360 <- ifelse(angles_raw < 0, angles_raw + 360, angles_raw)
  angles_sorted <- sort(angles_360)
  n <- length(angles_sorted)
  
  if (n == 1) {
    start_360 <- angles_sorted[1]
    end_360   <- angles_sorted[1]
  } else {
    # Compute gaps between consecutive points (including wrap-around)
    gaps <- c(diff(angles_sorted), angles_sorted[1] + 360 - angles_sorted[n])
    max_gap_idx <- which.max(gaps)
    
    if (max_gap_idx == n) {
      # No wrap: arc from first to last
      start_360 <- angles_sorted[1]
      end_360   <- angles_sorted[n]
    } else {
      # Wrap case: arc starts after the largest gap
      start_360 <- angles_sorted[max_gap_idx + 1]
      end_360   <- angles_sorted[max_gap_idx] + 360
    }
  }
  
  span_deg <- end_360 - start_360
  if (span_deg > 180) warning("Disease cone span > 180° — may not be biologically meaningful.")
  
  # Map back to [-180, 180] for final output
  normalize_angle <- function(a) ((a + 180) %% 360) - 180
  
  angle_start <- normalize_angle(start_360)
  angle_end   <- normalize_angle(end_360)
  
  # --- Step 4: Compute symmetric cone via y-axis reflection (θ → 180° − θ) ---
  angle_start_sym_raw <- 180 - angle_end     # note: use end for start_sym!
  angle_end_sym_raw   <- 180 - angle_start   # and start for end_sym!
  
  angle_start_sym <- normalize_angle(angle_start_sym_raw)
  angle_end_sym   <- normalize_angle(angle_end_sym_raw)
  
  # --- Step 5: Output ---
  cat("\nCoordinate system:\n")
  cat("         ↑ 90° (PC2)\n")
  cat("         |\n")
  cat("-180° ←--+--→ 0° (-PC1, disease axis)\n")
  cat("         |\n")
  cat("         ↓ -90°\n\n")
  
  cat("Disease cone: from ", angle_start, "° to ", angle_end, "° (逆时针)\n", sep = "")
  cat("Symmetric cone: from ", angle_start_sym, "° to ", angle_end_sym, "° (逆时针)\n", sep = "")
  cat("All disease angles: ", paste(round(angles_raw, 3), collapse = ", "), "\n\n")
  
  # --- Return ---
  list(
    angle_start       = angle_start,
    angle_end         = angle_end,
    angle_start_sym   = angle_start_sym,
    angle_end_sym     = angle_end_sym,
    angles_all_samples = angles_raw,
    cone_span_degrees = span_deg
  )
}


#' Filter Genes by Disease Cone Angle
#' 
#' @param eigengenes Data frame from calculate_eigengenes
#' @param cone_angles List from define_disease_cone
#' @return Filtered eigengenes
filter_by_disease_cone <- function(eigengenes, cone_info) {
  # --- Step 1: helper function to check if angle(s) fall in cone ---
  in_cone <- function(angle, start, end) {
    # Normalize all to [-180, 180)
    norm <- function(a) ((a + 180) %% 360) - 180
    angle <- norm(angle)
    start <- norm(start)
    end   <- norm(end)
    
    if (start <= end) {
      # No wrap: e.g., 10° to 50°
      return(angle >= start & angle <= end)
    } else {
      # Wrap case: e.g., 156° to -167° (i.e., 156° → 180° → -180° → -167°)
      return(angle >= start | angle <= end)
    }
  }
  # --- Step 2: Apply to eigengenes ---
  # Extract cone boundaries
  start1 <- cone_info$angle_start
  end1   <- cone_info$angle_end
  start2 <- cone_info$angle_start_sym
  end2   <- cone_info$angle_end_sym
  
  # Check membership in either cone
  in_original <- in_cone(eigengenes$angle_deg, start1, end1)
  in_symmetric <- in_cone(eigengenes$angle_deg, start2, end2)
  
  # Union: in original OR symmetric cone
  selected <- in_original | in_symmetric
  
  eigengenes_filtered <- eigengenes[selected, , drop = FALSE]
  
  cat("Genes in disease cone (original or symmetric):", 
      nrow(eigengenes_filtered), "out of", nrow(eigengenes), "\n")
  
  return(eigengenes_filtered)
}


# ============================================================================
# PART 3: Disease Contribution Analysis
# ============================================================================
#' Calculate Disease Contribution for Each Gene Across Multiple Contrasts
#'
#' Logic from Reference Paper 2:
#' - Pro-CGD: (Up-regulated AND cis-loading) OR (Down-regulated AND trans-loading)
#' - Anti-CGD: (Up-regulated AND trans-loading) OR (Down-regulated AND cis-loading)
#'
#' In our case:
#' - cgd_loading < 0 means gene expression increases as disease worsens (pro-disease, "cis")
#' - cgd_loading > 0 means gene expression decreases as disease worsens (anti-disease, "trans")
#'
#' @param eigengenes Eigengene data frame (must contain 'gene_id' and 'cgd_loading')
#' @param deg_results A list of DEG results, where each element corresponds to a contrast (e.g., c3_vs_c1)
#' @param contrasts A character vector or list of contrast names to process (e.g., c("c2_vs_c1", "c3_vs_c1"))
#' @return A named list of data frames, each with disease contribution for the respective contrast
calculate_disease_contribution <- function(eigengenes, deg_results, contrasts = c("c2_vs_c1", "c3_vs_c1", "c3_vs_c2")) {
  
  # Validate inputs
  if (!all(contrasts %in% names(deg_results))) {
    stop("Some contrasts not found in deg_results. Available: ", paste(names(deg_results), collapse = ", "))
  }
  
  # Process each contrast
  result_list <- lapply(contrasts, function(contrast_name) {
    
    # Extract DEG results for this contrast
    deg_table <- deg_results[[contrast_name]]
    
    # Merge with eigengenes
    contrib_data <- eigengenes %>%
      left_join(
        deg_table %>%
          rownames_to_column("gene_id") %>%
          dplyr::select(gene_id, log2FoldChange, padj, baseMean),
        by = "gene_id"
      )
    
    # Define cis/trans based on cgd_loading sign
    contrib_data <- contrib_data %>%
      mutate(
        is_cis = cgd_loading < 0,
        is_trans = cgd_loading > 0,
        is_upregulated = log2FoldChange > 0,
        is_downregulated = log2FoldChange < 0,
        is_significant = !is.na(padj) & padj < 0.05
      )
    
    # Calculate contribution
    contrib_data <- contrib_data %>%
      mutate(
        disease_contribution = case_when(
          # Pro-disease combinations
          (is_upregulated & is_cis) | (is_downregulated & is_trans) ~ "Pro-disease",
          # Anti-disease combinations
          (is_upregulated & is_trans) | (is_downregulated & is_cis) ~ "Anti-disease",
          TRUE ~ "Neutral"
        ),
        contribution_score = ifelse(
          disease_contribution == "Pro-disease",
          abs(cgd_loading) * abs(log2FoldChange),
          -abs(cgd_loading) * abs(log2FoldChange)
        )
      )
    
    # Optional: print summary per contrast
    cat("\nDisease contribution summary for contrast:", contrast_name, "\n")
    print(table(contrib_data$disease_contribution, 
                contrib_data$is_significant))
    
    return(contrib_data)
  })
  
  # Name the list elements by contrast
  names(result_list) <- contrasts
  
  return(result_list)
}
# ============================================================================
# PCA LOADINGS tell us:
# - The DIRECTION a gene points in the disease state-space
# - Derived from expression patterns across ALL samples (global trend)
# - Example: loading = -0.05 means "if this gene increases, samples 
#   move toward disease state"
# - BUT: Does NOT tell us if the gene actually changes in our disease model

# DEG ANALYSIS tells us:
# - Whether a gene ACTUALLY changes between disease vs health states
# - Derived from statistical comparison (e.g., c3 vs c1)
# - Example: log2FC = 2, padj < 0.05 means "gene is significantly 
#   upregulated in disease"
# - BUT: Does NOT tell us if this change is DIRECTIONALLY CONSISTENT 
#   with disease progression

# COMBINING BOTH identifies TRUE DRIVER GENES:
# 
#   Gene Loading (Direction) + Expression Change (DEG) = Contribution
#   ----------------------------------------------------------------
#   Cis (pro-disease)        + Up-regulated          = PRO-disease ✓
#   Cis (pro-disease)        + Down-regulated        = ANTI-disease
#   Trans (anti-disease)     + Up-regulated          = ANTI-disease
#   Trans (anti-disease)     + Down-regulated        = PRO-disease ✓
#
# Key insight: A gene can have strong loading but no expression change
# (not involved in THIS disease), or strong DEG but wrong direction
# (compensatory response, not driver).
#
# Only genes with BOTH:
#   1) Significant loading (contribute to disease axis)
#   2) Significant DEG (actually change in disease)
#   3) Consistent direction (loading + DEG agree)
# → are TRUE MECHANISTIC DRIVERS of state transition

# Analogy:
# - Loading = "This road LEADS TO the destination" (potential path)
# - DEG = "People ARE WALKING on this road" (actual traffic)
# - Combined = "This is the MAIN ROUTE people take" (driver pathway)
# ============================================================================



# ============================================================================
# PART 4: Multi-Layer Filtering for Driver Gene Selection
# ============================================================================

#' Apply Comprehensive Filters to Identify Driver Genes
#' 
#' @param contrib_data List with disease contribution dataframes
#' @param gene_annotations Gene annotations from GTF
#' @param vst_mat VST expression matrix
#' @return Filtered driver genes
select_driver_genes <- function(contrib_data_list, gene_annotations, vst_mat,
                                min_baseMean = 50,
                                min_abs_log2FC = 1,
                                max_padj = 0.05,
                                min_abs_cgd_loading = 0.01,
                                only_protein_coding = TRUE) {
  
  cat("\n=== Applying multi-layer filters ===\n")
  
  # Calculate mean expression across all samples (once for all contrasts)
  mean_expr <- rowMeans(vst_mat)
  expr_df <- data.frame(gene_id = names(mean_expr), mean_vst = mean_expr)
  
  # Process each contrast
  result_list <- lapply(names(contrib_data_list), function(contrast_name) {
    
    cat("\n--- Processing contrast:", contrast_name, "---\n")
    
    contrib_data <- contrib_data_list[[contrast_name]]
    
    # Merge all information
    driver_candidates <- contrib_data %>%
      left_join(gene_annotations, by = "gene_id") %>%
      left_join(expr_df, by = "gene_id")
    
    initial_n <- nrow(driver_candidates)
    
    # Filter 1: Pro-disease contribution
    driver_candidates <- driver_candidates %>%
      filter(disease_contribution == "Pro-disease")
    cat("After pro-disease filter:", nrow(driver_candidates), "/", initial_n, "\n")
    
    # Filter 2: Significant DEG
    driver_candidates <- driver_candidates %>%
      filter(is_significant)
    cat("After significance filter:", nrow(driver_candidates), "\n")
    
    # Filter 3: Expression level (avoid lowly expressed genes)
    driver_candidates <- driver_candidates %>%
      filter(baseMean >= min_baseMean | mean_vst >= 5)
    cat("After expression filter:", nrow(driver_candidates), "\n")
    
    # Filter 4: Strong fold change
    driver_candidates <- driver_candidates %>%
      filter(abs(log2FoldChange) >= min_abs_log2FC)
    cat("After fold-change filter:", nrow(driver_candidates), "\n")
    
    # Filter 5: Substantial loading on disease axis
    driver_candidates <- driver_candidates %>%
      filter(abs(cgd_loading) >= min_abs_cgd_loading)
    cat("After loading filter:", nrow(driver_candidates), "\n")
    
    # Filter 6: Protein-coding genes only (optional)
    if (only_protein_coding) {
      driver_candidates <- driver_candidates %>%
        filter(gene_biotype == "protein_coding")
      cat("After protein-coding filter:", nrow(driver_candidates), "\n")
    }
    # Filter 7: Remove sex chromosome genes
    driver_candidates <- driver_candidates %>%
      filter(!seqnames %in% c("X", "Y", "chrX", "chrY"))
    cat("After sex chromosome filter:", nrow(driver_candidates), "\n")

    # Rank by contribution score
    driver_candidates <- driver_candidates %>%
      arrange(desc(abs(contribution_score))) %>%
      mutate(rank = row_number())
    
    cat("=== Final driver genes for", contrast_name, ":", nrow(driver_candidates), "===\n")
    
    return(driver_candidates)
  })
  
  # Name the list elements by contrast
  names(result_list) <- names(contrib_data_list)
  
  # Print overall summary
  cat("\n=== Overall Summary ===\n")
  for (contrast_name in names(result_list)) {
    cat(sprintf("%-12s: %4d driver genes\n", 
                contrast_name, nrow(result_list[[contrast_name]])))
  }
  
  return(result_list)
}

#' Classify Driver Genes by Expression Pattern
#' 
#' @param driver_genes_list Filtered driver genes list
#' @return Driver genes list with pattern classification
classify_gene_patterns <- function(driver_genes_list) {
  
  # Process each contrast
  result_list <- lapply(names(driver_genes_list), function(contrast_name) {
    
    cat("\n=== Classifying patterns for:", contrast_name, "===\n")
    
    driver_genes <- driver_genes_list[[contrast_name]]
    
    driver_genes <- driver_genes %>%
      mutate(
        expression_pattern = case_when(
          is_upregulated & is_cis ~ "Upregulated-Cis (Strong Pro-disease)",
          is_downregulated & is_trans ~ "Downregulated-Trans (Pro-disease)",
          TRUE ~ "Other"
        ),
        loading_strength = case_when(
          abs(cgd_loading) >= 0.04 ~ "Strong",
          abs(cgd_loading) >= 0.02 ~ "Moderate",
          TRUE ~ "Weak"
        ),
        fc_strength = case_when(
          abs(log2FoldChange) >= 2 ~ "Strong (>4-fold)",
          abs(log2FoldChange) >= 1.5 ~ "Moderate (>2.8-fold)",
          TRUE ~ "Weak (>2-fold)"
        )
      )
    
    # Summary by pattern
    cat("\nGene patterns:\n")
    print(table(driver_genes$expression_pattern))
    
    cat("\nLoading strength:\n")
    print(table(driver_genes$loading_strength))
    
    cat("\nFold-change strength:\n")
    print(table(driver_genes$fc_strength))
    
    return(driver_genes)
  })
  
  # Name the list elements by contrast
  names(result_list) <- names(driver_genes_list)
  
  return(result_list)
}


# ============================================================================
# PART 5: Main Execution Pipeline
# ============================================================================

#' Main Function to Execute Driver Gene Selection
#' 
#' @return List with all results
run_driver_gene_analysis <- function() {
  
  cat("==============================================\n")
  cat("Starting Driver Gene Selection Pipeline\n")
  cat("==============================================\n\n")
  
  # Step 1: Load data
  cat("Step 1: Loading project data...\n")
  contrasts <- c("c1_vs_c1_star", "c2_vs_c1_star", "c3_vs_c1_star", "c2_vs_c1", "c3_vs_c1", "c3_vs_c2")
  cat("Contrasts to analyze:", paste(contrasts, collapse = ", "), "\n")
  project_data <- load_project_data(comparisons = contrasts)
  
  # Step 2: Calculate eigengenes
  cat("\nStep 2: Calculating eigengenes...\n")
  eigengenes <- calculate_eigengenes(project_data$pca_result)
  
  # Step 3: Define disease cone
  cat("\nStep 3: Defining disease cone...\n")
  cone_angles <- define_disease_cone(
    project_data$x_cgd,
    project_data$sample_info_pcs,
    project_data$critical_points
  )
  
  # Step 4: Filter by cone (optional, can skip this filter)
  eigengenes_in_cone <- filter_by_disease_cone(eigengenes, cone_angles)
  
  # Step 5: Calculate disease contribution
  cat("\nStep 5: Calculating disease contribution...\n")
  contrib_data <- calculate_disease_contribution(
    eigengenes_in_cone,
    project_data$deg_results,
    contrasts = contrasts
  )
  
  # Step 6: Select driver genes
  cat("\nStep 6: Selecting driver genes with filters...\n")
  driver_genes <- select_driver_genes(
    contrib_data,
    project_data$gene_annotations,
    project_data$vst_mat,
    min_baseMean = 50,
    min_abs_log2FC = 1,
    max_padj = 0.05,
    min_abs_cgd_loading = 0.01,
    only_protein_coding = TRUE
  )
  
  # Step 7: Classify patterns
  cat("\nStep 7: Classifying gene patterns...\n")
  driver_genes <- classify_gene_patterns(driver_genes)

  return(list(
    project_data = project_data,
    eigengenes = eigengenes,
    cone_angles = cone_angles,
    eigengenes_in_cone = eigengenes_in_cone,
    contrib_data = contrib_data,
    driver_genes = driver_genes
  ))
}


# ============================================================================
# USAGE
# ============================================================================

# Run the complete analysis
results <- run_driver_gene_analysis()

# Save results
cat("\nSaving results...\n")
output_dir <- file.path(BASE_DIR, "08.eigengene")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
save(results, file = file.path(output_dir, "driver_gene_selection_results.rdata"))
cat("\n==============================================\n")
cat("Analysis Complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("==============================================\n")