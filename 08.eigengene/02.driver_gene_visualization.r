# ============================================================================
# Driver Gene visualization
# Based on Driver Gene Selection Pipeline
# ============================================================================

library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(ComplexUpset)
library(dplyr)

BASE_DIR <- "/disk2/cai113/data/stateTrans"
# ============================================================================
# PART 1: Gene Distribution in PC1-PC2 Space
# ============================================================================

#' Plot Gene Distribution in PC1-PC2 Space
#' 
#' Visualizes all genes as points in the principal component space,
#' with PC1 (disease axis) on x-axis and PC2 (orthogonal axis) on y-axis
#' 
#' @param pca_result PCA result object from prcomp()
#' @param highlight_genes Optional vector of gene IDs to highlight
#' @param highlight_color Color for highlighted genes (default: red)
#' @param point_size Size of points (default: 1)
#' @param point_alpha Transparency of points (default: 0.3)
#' @return ggplot object
plot_gene_pc_space <- function(pca_result, 
                               highlight_genes = NULL,
                               highlight_color = "#E41A1C",
                               point_size = 1,
                               point_alpha = 0.3) {
  
  # Extract gene loadings (rotation matrix)
  loadings <- pca_result$rotation
  
  # Create data frame with PC1 and PC2 loadings
  gene_df <- data.frame(
    gene_id = rownames(loadings),
    PC1 = loadings[, "PC1"],
    PC2 = loadings[, "PC2"],
    neg_PC1 = -loadings[, "PC1"],  # Disease axis (negated PC1)
    stringsAsFactors = FALSE
  )
  
  # Add highlight status if genes provided
  if (!is.null(highlight_genes)) {
    gene_df$is_highlighted <- gene_df$gene_id %in% highlight_genes
  } else {
    gene_df$is_highlighted <- FALSE
  }
  
  # Calculate variance explained for axis labels
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  pc1_var <- round(var_explained[1], 2)
  pc2_var <- round(var_explained[2], 2)
  
  # Create plot
  p <- ggplot(gene_df, aes(x = neg_PC1, y = PC2)) +
    # Background genes
    geom_point(data = gene_df %>% filter(!is_highlighted),
               color = "gray70", size = point_size, alpha = point_alpha) +
    # Highlighted genes (if any)
    {if(any(gene_df$is_highlighted)) {
      geom_point(data = gene_df %>% filter(is_highlighted),
                 color = highlight_color, size = point_size * 2, alpha = 0.8)
    }} +
    # Add origin lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Labels and theme
    labs(
      x = sprintf("Disease Axis (-PC1, %.2f%% variance)", pc1_var),
      y = sprintf("PC2 (%.2f%% variance)", pc2_var),
      title = "Gene Distribution in Principal Component Space",
      subtitle = sprintf("Total genes: %s%s", 
                        format(nrow(gene_df), big.mark = ","),
                        ifelse(any(gene_df$is_highlighted), 
                              sprintf(" | Highlighted: %d", sum(gene_df$is_highlighted)),
                              ""))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_fixed(ratio = 1)  # Equal aspect ratio for proper geometry
  
  return(p)
}


# ============================================================================
# PART 2: Sample Distribution in PC Space with Disease Cone
# ============================================================================

#' Plot Sample Distribution in PC Space with Disease Cone
#' 
#' Visualizes all samples in PC1-PC2 space with disease cone boundaries,
#' distinguishing WT, healthy KO, and disease KO samples
#' 
#' @param sample_info_pcs Data frame with sample information and PC scores
#' @param cone_angles List from define_disease_cone containing cone boundaries
#' @param critical_points List with c1, c2, c3, c1_star values
#' @param cone_alpha Transparency of cone shading (default: 0.15)
#' @param point_size Size of sample points (default: 3)
#' @return ggplot object
plot_sample_cone_space <- function(sample_info_pcs, 
                                   cone_angles, 
                                   critical_points,
                                   cone_alpha = 0.15,
                                   point_size = 3) {
  
  # --- Step 1: Prepare sample data ---
  sample_df <- sample_info_pcs %>%
    mutate(
      neg_PC1 = -PC1,  # Disease axis
      sample_type = case_when(
        genotype == "WT" ~ "WT (Control)",
        genotype == "Ncf2" & x_cgd < critical_points$c2 ~ "KO (Disease, <c2)",
        genotype == "Ncf2" & x_cgd >= critical_points$c2 ~ "KO (Health-like, ≥c2)",
        TRUE ~ "Other"
      )
    )
  
  # --- Step 2: Create cone polygon data ---
  # Main disease cone
  cone_main <- create_cone_polygon(
    angle_start = cone_angles$angle_start,
    angle_end = cone_angles$angle_end,
    radius = get_plot_radius(sample_df)
  )
  cone_main$cone_type <- "Disease Cone"
  
  # Symmetric cone (y-axis reflection)
  cone_sym <- create_cone_polygon(
    angle_start = cone_angles$angle_start_sym,
    angle_end = cone_angles$angle_end_sym,
    radius = get_plot_radius(sample_df)
  )
  cone_sym$cone_type <- "Symmetric Cone"
  
  # Combine polygons
  cone_polygons <- rbind(cone_main, cone_sym)
  
  # --- Step 3: Create boundary lines with angle labels ---
  boundaries <- data.frame(
    angle = c(cone_angles$angle_start, cone_angles$angle_end,
              cone_angles$angle_start_sym, cone_angles$angle_end_sym),
    cone = c("Main", "Main", "Symmetric", "Symmetric"),
    label = c(
      sprintf("%.1f°", cone_angles$angle_start),
      sprintf("%.1f°", cone_angles$angle_end),
      sprintf("%.1f°", cone_angles$angle_start_sym),
      sprintf("%.1f°", cone_angles$angle_end_sym)
    )
  )
  
  # Convert angles to line endpoints
  max_radius <- get_plot_radius(sample_df)
  boundaries <- boundaries %>%
    mutate(
      x_end = max_radius * cos(angle * pi / 180),
      y_end = max_radius * sin(angle * pi / 180),
      label_x = max_radius * 1.15 * cos(angle * pi / 180),
      label_y = max_radius * 1.15 * sin(angle * pi / 180)
    )
  
  # --- Step 4: Create plot ---
  p <- ggplot() +
    # Cone shading
    geom_polygon(data = cone_polygons, 
                 aes(x = x, y = y, fill = cone_type),
                 alpha = cone_alpha) +
    # Cone boundary lines
    geom_segment(data = boundaries,
                 aes(x = 0, y = 0, xend = x_end, yend = y_end, 
                     color = cone, linetype = cone),
                 linewidth = 0.8) +
    # Angle labels
    geom_text(data = boundaries,
              aes(x = label_x, y = label_y, label = label),
              size = 3.5, fontface = "bold", color = "gray20") +
    # Origin axes
    geom_hline(yintercept = 0, linetype = "solid", 
               color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", 
               color = "gray50", linewidth = 0.5) +
    # Sample points
    geom_point(data = sample_df,
               aes(x = neg_PC1, y = PC2, 
                   color = sample_type, shape = sample_type),
               size = point_size, alpha = 0.8, stroke = 1.2) +
    # Sample labels (optional, for small datasets)
    {if(nrow(sample_df) <= 40) {
      geom_text(data = sample_df,
                aes(x = neg_PC1, y = PC2, label = sample),
                size = 2.5, vjust = -1.2, hjust = 0.5, color = "gray30")
    }} +
    # Color and shape scales
    scale_color_manual(
      name = "Sample Type",
      values = c(
        "WT (Control)" = "#377EB8",
        "KO (Health-like, ≥c2)" = "#4DAF4A",
        "KO (Disease, <c2)" = "#E41A1C",
        "Main" = "#984EA3",
        "Symmetric" = "#FF7F00"
      ),
      breaks = c("WT (Control)", "KO (Health-like, ≥c2)", "KO (Disease, <c2)")
    ) +
    scale_shape_manual(
      name = "Sample Type",
      values = c(
        "WT (Control)" = 17,
        "KO (Health-like, ≥c2)" = 16,
        "KO (Disease, <c2)" = 15
      )
    ) +
    scale_fill_manual(
      name = "Cone Region",
      values = c("Disease Cone" = "#E41A1C", "Symmetric Cone" = "#FF7F00")
    ) +
    scale_linetype_manual(
      name = "Boundary",
      values = c("Main" = "dashed", "Symmetric" = "dotted")
    ) +
    # Labels and theme
    labs(
      x = "Disease Axis (-PC1)",
      y = "PC2 (Orthogonal Axis)",
      title = "Sample Distribution with Disease Cone",
      subtitle = sprintf(
        "Cone span: %.1f° | Disease samples: %d | c2 threshold: %.2f",
        cone_angles$cone_span_degrees,
        sum(sample_df$sample_type == "KO (Disease, <c2)"),
        critical_points$c2
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_fixed(ratio = 1)  # Equal aspect ratio
  
  return(p)
}

#' Create Cone Polygon for Shading
#' 
#' @param angle_start Starting angle in degrees
#' @param angle_end Ending angle in degrees
#' @param radius Maximum radius for cone
#' @param n_points Number of points for smooth arc (default: 100)
#' @return Data frame with x, y coordinates for polygon
create_cone_polygon <- function(angle_start, angle_end, radius, n_points = 100) {
  
  # Normalize angles to [0, 360)
  normalize_360 <- function(a) ((a %% 360) + 360) %% 360
  start_360 <- normalize_360(angle_start)
  end_360 <- normalize_360(angle_end)
  
  # Determine arc direction (counter-clockwise from start to end)
  if (end_360 > start_360) {
    angles_arc <- seq(start_360, end_360, length.out = n_points)
  } else {
    # Wrap case: go through 0°/360°
    angles_arc <- c(
      seq(start_360, 360, length.out = n_points/2),
      seq(0, end_360, length.out = n_points/2)
    )
  }
  
  # Convert back to [-180, 180] for consistency
  angles_arc <- ((angles_arc + 180) %% 360) - 180
  
  # Create polygon: origin → arc → back to origin
  x_coords <- c(0, radius * cos(angles_arc * pi / 180), 0)
  y_coords <- c(0, radius * sin(angles_arc * pi / 180), 0)
  
  data.frame(x = x_coords, y = y_coords)
}


#' Get Appropriate Plot Radius Based on Sample Distribution
#' 
#' @param sample_df Sample data frame with neg_PC1 and PC2
#' @return Radius value for cone boundary lines
get_plot_radius <- function(sample_df) {
  max_dist <- max(sqrt(sample_df$neg_PC1^2 + sample_df$PC2^2))
  return(max_dist * 1.3)  # 30% beyond furthest sample
}


# ============================================================================
# PART 3: Gene Distribution in PC Space with Disease Cone
# ============================================================================

#' Plot Gene Distribution in PC Space with Disease Cone
#' 
#' Visualizes all genes in PC1-PC2 space with disease cone boundaries,
#' highlighting genes within the disease-associated cone
#' 
#' @param pca_result PCA result object from prcomp()
#' @param cone_angles List from define_disease_cone containing cone boundaries
#' @param cone_alpha Transparency of cone shading (default: 0.15)
#' @param point_size Size of gene points (default: 0.8)
#' @param point_alpha Transparency of gene points (default: 0.3)
#' @return ggplot object
plot_gene_cone_space <- function(pca_result,
                                 cone_angles,
                                 cone_alpha = 0.15,
                                 point_size = 0.8,
                                 point_alpha = 0.3) {
  
  # --- Step 1: Extract gene loadings ---
  loadings <- pca_result$rotation
  
  gene_df <- data.frame(
    gene_id = rownames(loadings),
    PC1 = loadings[, "PC1"],
    PC2 = loadings[, "PC2"],
    neg_PC1 = -loadings[, "PC1"],  # Disease axis
    stringsAsFactors = FALSE
  )
  
  # Calculate angles for each gene
  gene_df <- gene_df %>%
    mutate(
      angle = atan2(PC2, neg_PC1) * 180 / pi
    )
  
  # --- Step 2: Determine which genes are in cone ---
  gene_df$in_cone <- is_in_cone(
    gene_df$angle,
    cone_angles$angle_start,
    cone_angles$angle_end
  ) | is_in_cone(
    gene_df$angle,
    cone_angles$angle_start_sym,
    cone_angles$angle_end_sym
  )
  
  # --- Step 3: Create cone polygon data ---
  max_radius <- max(sqrt(gene_df$neg_PC1^2 + gene_df$PC2^2)) * 1.2
  
  # Main disease cone
  cone_main <- create_cone_polygon(
    angle_start = cone_angles$angle_start,
    angle_end = cone_angles$angle_end,
    radius = max_radius
  )
  cone_main$cone_type <- "Disease Cone"
  
  # Symmetric cone
  cone_sym <- create_cone_polygon(
    angle_start = cone_angles$angle_start_sym,
    angle_end = cone_angles$angle_end_sym,
    radius = max_radius
  )
  cone_sym$cone_type <- "Symmetric Cone"
  
  # Combine polygons
  cone_polygons <- rbind(cone_main, cone_sym)
  
  # --- Step 4: Create boundary lines ---
  boundaries <- data.frame(
    angle = c(cone_angles$angle_start, cone_angles$angle_end,
              cone_angles$angle_start_sym, cone_angles$angle_end_sym),
    cone = c("Main", "Main", "Symmetric", "Symmetric"),
    label = c(
      sprintf("%.1f°", cone_angles$angle_start),
      sprintf("%.1f°", cone_angles$angle_end),
      sprintf("%.1f°", cone_angles$angle_start_sym),
      sprintf("%.1f°", cone_angles$angle_end_sym)
    )
  )
  
  boundaries <- boundaries %>%
    mutate(
      x_end = max_radius * cos(angle * pi / 180),
      y_end = max_radius * sin(angle * pi / 180),
      label_x = max_radius * 1.08 * cos(angle * pi / 180),
      label_y = max_radius * 1.08 * sin(angle * pi / 180)
    )
  
  # --- Step 5: Calculate variance explained ---
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  pc1_var <- round(var_explained[1], 2)
  pc2_var <- round(var_explained[2], 2)
  
  # --- Step 6: Create plot ---
  n_in_cone <- sum(gene_df$in_cone)
  n_total <- nrow(gene_df)
  pct_in_cone <- round(n_in_cone / n_total * 100, 1)
  
  subtitle_text <- sprintf(
    "Genes in cone: %s / %s (%.1f%%) | Cone span: %.1f°",
    format(n_in_cone, big.mark = ","),
    format(n_total, big.mark = ","),
    pct_in_cone,
    cone_angles$cone_span_degrees
  )
  
  p <- ggplot() +
    # Cone shading
    geom_polygon(data = cone_polygons,
                 aes(x = x, y = y, fill = cone_type),
                 alpha = cone_alpha) +
    # Cone boundary lines
    geom_segment(data = boundaries,
                 aes(x = 0, y = 0, xend = x_end, yend = y_end,
                     color = cone, linetype = cone),
                 linewidth = 0.8) +
    # Angle labels
    geom_text(data = boundaries,
              aes(x = label_x, y = label_y, label = label),
              size = 3, fontface = "bold", color = "gray20") +
    # Origin axes
    geom_hline(yintercept = 0, linetype = "solid",
               color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid",
               color = "gray50", linewidth = 0.5) +
    # Gene points - outside cone first (background)
    geom_point(data = gene_df %>% filter(!in_cone),
               aes(x = neg_PC1, y = PC2),
               color = "gray70", size = point_size, alpha = point_alpha) +
    # Gene points - inside cone (highlighted)
    geom_point(data = gene_df %>% filter(in_cone),
               aes(x = neg_PC1, y = PC2),
               color = "#E41A1C", size = point_size * 1.5, alpha = 0.6) +
    # Color scales
    scale_fill_manual(
      name = "Cone Region",
      values = c("Disease Cone" = "#E41A1C", "Symmetric Cone" = "#FF7F00")
    ) +
    scale_color_manual(
      name = "Boundary",
      values = c("Main" = "#984EA3", "Symmetric" = "#FF7F00")
    ) +
    scale_linetype_manual(
      name = "Boundary",
      values = c("Main" = "dashed", "Symmetric" = "dotted")
    ) +
    # Labels and theme
    labs(
      x = sprintf("Disease Axis (-PC1, %.2f%% variance)", pc1_var),
      y = sprintf("PC2 (%.2f%% variance)", pc2_var),
      title = "Gene Distribution with Disease Cone",
      subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_fixed(ratio = 1)
  
  return(p)
}


#' Check if an Angle Falls Within a Cone Range
#' 
#' @param angle Angle(s) to check (in degrees, [-180, 180])
#' @param angle_start Cone start angle
#' @param angle_end Cone end angle
#' @return Logical vector indicating if angle is in cone
is_in_cone <- function(angle, angle_start, angle_end) {
  
  # Normalize all angles to [0, 360)
  normalize_360 <- function(a) ((a %% 360) + 360) %% 360
  
  angle_norm <- normalize_360(angle)
  start_norm <- normalize_360(angle_start)
  end_norm <- normalize_360(angle_end)
  
  # Check if in cone (counter-clockwise from start to end)
  if (end_norm > start_norm) {
    # Normal case: no wrap
    in_cone <- (angle_norm >= start_norm) & (angle_norm <= end_norm)
  } else {
    # Wrap case: cone crosses 0°/360°
    in_cone <- (angle_norm >= start_norm) | (angle_norm <= end_norm)
  }
  
  return(in_cone)
}


# ============================================================================
# PART 4: Gene Distribution in Polar Coordinates with Disease Cone
# ============================================================================

#' Plot Gene Distribution in Polar Coordinates with Disease Cone
#' 
#' Visualizes genes in polar coordinate system (radius = loading magnitude,
#' angle = direction in PC space), clearly showing disease cone boundaries
#' 
#' @param eigengenes Data frame with gene eigengene information
#' @param cone_angles List from define_disease_cone containing cone boundaries
#' @param max_radius_quantile Quantile for maximum radius display (default: 0.99)
#' @param point_size Size of gene points (default: 1)
#' @param point_alpha Transparency of gene points (default: 0.4)
#' @return ggplot object
plot_gene_polar_cone <- function(eigengenes,
                                 cone_angles,
                                 max_radius_quantile = 0.99,
                                 point_size = 1,
                                 point_alpha = 0.4) {
  
  # --- Step 1: Prepare gene data ---
  gene_polar <- eigengenes %>%
    mutate(
      angle = angle_deg,
      r = radius
    )
  
  # --- Step 2: Determine which genes are in cone ---
  gene_polar$in_cone <- is_in_cone(
    gene_polar$angle,
    cone_angles$angle_start,
    cone_angles$angle_end
  ) | is_in_cone(
    gene_polar$angle,
    cone_angles$angle_start_sym,
    cone_angles$angle_end_sym
  )
  
  # --- Step 3: Set radius limits for better visualization ---
  max_r <- quantile(gene_polar$r, probs = max_radius_quantile, na.rm = TRUE)
  gene_polar <- gene_polar %>%
    mutate(r_capped = pmin(r, max_r))
  
  # --- Step 4: Statistics ---
  n_in_cone <- sum(gene_polar$in_cone)
  n_total <- nrow(gene_polar)
  pct_in_cone <- round(n_in_cone / n_total * 100, 1)
  
  # --- Step 5: Create polar plot ---
  p <- ggplot() +
    # Gene points - outside cone (background)
    geom_point(data = gene_polar %>% filter(!in_cone),
               aes(x = angle, y = r_capped),
               color = "gray70", size = point_size, alpha = point_alpha) +
    # Gene points - inside cone (highlighted)
    geom_point(data = gene_polar %>% filter(in_cone),
               aes(x = angle, y = r_capped),
               color = "#E41A1C", size = point_size * 1.5, alpha = 0.7) +
    # Cone boundary lines
    geom_vline(xintercept = cone_angles$angle_start,
               linetype = "dashed", color = "#984EA3", linewidth = 1.2) +
    geom_vline(xintercept = cone_angles$angle_end,
               linetype = "dashed", color = "#984EA3", linewidth = 1.2) +
    geom_vline(xintercept = cone_angles$angle_start_sym,
               linetype = "dotted", color = "#FF7F00", linewidth = 1.2) +
    geom_vline(xintercept = cone_angles$angle_end_sym,
               linetype = "dotted", color = "#FF7F00", linewidth = 1.2) +
    # Angle labels on boundaries
    annotate("text", x = cone_angles$angle_start, y = max_r * 0.8,
             label = sprintf("%.1f°", cone_angles$angle_start),
             size = 3.5, fontface = "bold", color = "#984EA3") +
    annotate("text", x = cone_angles$angle_end, y = max_r * 0.8,
             label = sprintf("%.1f°", cone_angles$angle_end),
             size = 3.5, fontface = "bold", color = "#984EA3") +
    annotate("text", x = cone_angles$angle_start_sym, y = max_r * 0.8,
             label = sprintf("%.1f°", cone_angles$angle_start_sym),
             size = 3.5, fontface = "bold", color = "#FF7F00") +
    annotate("text", x = cone_angles$angle_end_sym, y = max_r * 0.8,
             label = sprintf("%.1f°", cone_angles$angle_end_sym),
             size = 3.5, fontface = "bold", color = "#FF7F00") +
    # Scales
    scale_y_continuous(
      name = "Loading Magnitude",
      limits = c(0, max_r * 1.05),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = seq(-180, 180, 30),
      limits = c(-180, 180)
    ) +
    # Convert to polar coordinates with correct rotation
    # start = pi/2 makes 0° at top, and positive angles go clockwise (correct!)
    coord_polar(theta = "x", start = pi/2, direction = -1) +
    # Labels and theme
    labs(
      title = "Gene Distribution in Polar Coordinates",
      subtitle = sprintf(
        "Genes in cone: %s / %s (%.1f%%) | Cone span: %.1f°",
        format(n_in_cone, big.mark = ","),
        format(n_total, big.mark = ","),
        pct_in_cone,
        cone_angles$cone_span_degrees
      ),
      x = "Angle (degrees)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10),
      legend.position = "none",  # Remove legend since no fill/color legend needed
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.2)
    )
  
  return(p)
}


# ============================================================================
# PART 5: Driver Gene Visualization in PC Space with Cone
# ============================================================================

#' Plot Driver Genes Distribution in PC Space with Cone
#' 
#' Visualizes driver genes in PC1-PC2 space, showing their differential
#' expression status (up/down-regulated) against background of cone genes
#' 
#' @param eigengenes_in_cone Data frame with genes in disease cone
#' @param driver_genes_list List of driver genes data frames by contrast
#' @param cone_angles List from define_disease_cone containing cone boundaries
#' @param contrasts Vector of contrast names to display, or "ALL" for union
#' @param point_size Size of gene points (default: 1.5)
#' @param bg_alpha Transparency of background genes (default: 0.2)
#' @return ggplot object
plot_driver_genes_in_cone <- function(eigengenes_in_cone,
                                      driver_genes_list,
                                      cone_angles,
                                      contrasts = "ALL",
                                      point_size = 1.5,
                                      bg_alpha = 0.2) {
  
  # --- Step 1: Prepare background (all genes in cone) ---
  bg_genes <- eigengenes_in_cone %>%
    mutate(neg_PC1 = -v1)
  
  # --- Step 2: Extract driver genes based on contrasts ---
  if (length(contrasts) == 1 && contrasts[1] == "ALL") {
    # Union of all driver genes
    all_driver_ids <- unique(unlist(lapply(driver_genes_list, function(df) df$gene_id)))
    
    # Get driver gene info (just IDs, no expression info)
    driver_display <- eigengenes_in_cone %>%
      filter(gene_id %in% all_driver_ids) %>%
      mutate(
        neg_PC1 = -v1,
        regulation = "Driver Gene",
        contrast = "ALL"
      )
    
    plot_title <- "Driver Genes Distribution (All Contrasts)"
    use_color <- FALSE
    
  } else {
    # Specific contrasts with up/down regulation
    driver_list <- lapply(contrasts, function(contrast_name) {
      if (!contrast_name %in% names(driver_genes_list)) {
        warning(paste("Contrast", contrast_name, "not found. Skipping."))
        return(NULL)
      }
      driver_genes_list[[contrast_name]] %>%
        dplyr::mutate(
          neg_PC1 = -v1,
          regulation = ifelse(is_upregulated, "Up-regulated", "Down-regulated"),
          contrast = contrast_name
        ) %>%
        dplyr::select(gene_id, gene_name, neg_PC1, v2, cgd_loading, 
                      log2FoldChange, regulation, contrast)      
    })
    
    driver_display <- do.call(rbind, driver_list[!sapply(driver_list, is.null)])
    
    if (length(contrasts) == 1) {
      plot_title <- sprintf("Driver Genes Distribution (%s)", contrasts[1])
    } else {
      plot_title <- sprintf("Driver Genes Distribution (%d contrasts)", length(contrasts))
    }
    use_color <- TRUE
  }
  
  # --- Step 3: Calculate plot boundaries ---
  max_radius <- max(sqrt(bg_genes$neg_PC1^2 + bg_genes$v2^2)) * 1.2
  
  # --- Step 4: Create cone boundary lines ---
  boundaries <- data.frame(
    angle = c(cone_angles$angle_start, cone_angles$angle_end,
              cone_angles$angle_start_sym, cone_angles$angle_end_sym),
    cone = c("Main", "Main", "Symmetric", "Symmetric"),
    label = c(
      sprintf("%.1f°", cone_angles$angle_start),
      sprintf("%.1f°", cone_angles$angle_end),
      sprintf("%.1f°", cone_angles$angle_start_sym),
      sprintf("%.1f°", cone_angles$angle_end_sym)
    )
  )
  
  boundaries <- boundaries %>%
    mutate(
      x_end = max_radius * cos(angle * pi / 180),
      y_end = max_radius * sin(angle * pi / 180),
      label_x = max_radius * 1.05 * cos(angle * pi / 180),
      label_y = max_radius * 1.05 * sin(angle * pi / 180)
    )
  
  # --- Step 5: Create plot ---
  n_drivers <- nrow(driver_display)
  n_bg <- nrow(bg_genes)
  
  subtitle_text <- sprintf(
    "Driver genes: %d / %d genes in cone (%.1f%%)",
    n_drivers, n_bg, n_drivers/n_bg * 100
  )
  
  if (use_color) {
    n_up <- sum(driver_display$regulation == "Up-regulated")
    n_down <- sum(driver_display$regulation == "Down-regulated")
    subtitle_text <- paste0(subtitle_text, 
                           sprintf(" | Up: %d, Down: %d", n_up, n_down))
  }
  
  p <- ggplot() +
    # Background: all genes in cone
    geom_point(data = bg_genes,
               aes(x = neg_PC1, y = v2),
               color = "gray70", size = point_size * 0.6, alpha = bg_alpha) +
    # Cone boundary lines
    geom_segment(data = boundaries,
                 aes(x = 0, y = 0, xend = x_end, yend = y_end,
                     color = cone, linetype = cone),
                 linewidth = 0.8, show.legend = FALSE) +
    # Angle labels
    geom_text(data = boundaries,
              aes(x = label_x, y = label_y, label = label),
              size = 3, fontface = "bold", color = "gray20") +
    # Origin axes
    geom_hline(yintercept = 0, linetype = "solid",
               color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid",
               color = "gray50", linewidth = 0.5)
  
  # Add driver genes - different styling for ALL vs specific contrasts
  if (use_color) {
    # Colored by regulation status
    p <- p +
      geom_point(data = driver_display,
                 aes(x = neg_PC1, y = v2, color = regulation, shape = regulation),
                 size = point_size, alpha = 0.8, stroke = 1) +
      scale_color_manual(
        name = "Regulation",
        values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#4DAF4A")
      ) +
      scale_shape_manual(
        name = "Regulation",
        values = c("Up-regulated" = 16, "Down-regulated" = 17)
      )
  } else {
    # All black for union
    p <- p +
      geom_point(data = driver_display,
                 aes(x = neg_PC1, y = v2),
                 color = "black", size = point_size, alpha = 0.7)
  }
  
  # Boundary line colors (not in legend)
  p <- p +
    scale_linetype_manual(
      values = c("Main" = "dashed", "Symmetric" = "dotted")
    ) +
    # Labels and theme
    labs(
      x = "Disease Axis (-PC1)",
      y = "PC2 (Orthogonal Axis)",
      title = plot_title,
      subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_fixed(ratio = 1)
  
  return(p)
}


# ============================================================================
# PART 6: UpSet Plot of Driver Genes Across Contrasts
# ============================================================================

#' Plot UpSet Diagram for Driver Genes Across Contrasts (Simple Version)
#'
#' This function visualizes the overlap of driver genes across multiple contrasts
#' using an UpSet plot. All advanced features (e.g., gene labeling) are removed
#' for maximum stability and compatibility.
#'
#' @param driver_genes_list A named list of data frames, each with at least:
#'   - `gene_id`: unique gene identifier
#' @param contrasts Character vector of contrast names to include.
#'   If NULL, all contrasts in `driver_genes_list` are used.
#' @param filter_pattern Expression pattern to retain (e.g., "Upregulated-Cis...").
#'   If NULL, no filtering.
#' @param filter_loading_strength Loading strength filter (e.g., "Strong").
#' @param filter_fc_strength Fold-change strength filter (e.g., "Strong (>4-fold)").
#' @param top_n_intersections Number of top intersections to show. Default: 20.
#'
#' @return A `ggplot` object.
plot_driver_genes_upset <- function(driver_genes_list,
                                    contrasts = NULL,
                                    filter_pattern = NULL,
                                    filter_loading_strength = NULL,
                                    filter_fc_strength = NULL,
                                    top_n_intersections = 20) {
  
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' is required. Install via: install.packages('ComplexUpset')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  library(ComplexUpset)
  library(dplyr)
  
  # Step 1: Select contrasts
  if (is.null(contrasts)) {
    contrasts <- names(driver_genes_list)
  } else {
    invalid <- contrasts[!contrasts %in% names(driver_genes_list)]
    if (length(invalid) > 0) {
      stop("Invalid contrast names: ", paste(invalid, collapse = ", "))
    }
  }
  
  cat("=== Creating UpSet plot ===\n")
  cat("Contrasts:", paste(contrasts, collapse = ", "), "\n")
  
  # Step 2: Apply filters and extract gene sets
  gene_sets <- lapply(contrasts, function(contrast_name) {
    df <- driver_genes_list[[contrast_name]]
    if (!is.null(filter_pattern)) {
      df <- df %>% filter(expression_pattern == filter_pattern)
    }
    if (!is.null(filter_loading_strength)) {
      df <- df %>% filter(loading_strength == filter_loading_strength)
    }
    if (!is.null(filter_fc_strength)) {
      df <- df %>% filter(fc_strength == filter_fc_strength)
    }
    return(df$gene_id)
  })
  names(gene_sets) <- contrasts
  
  # Print filters and gene counts
  if (!is.null(filter_pattern) || !is.null(filter_loading_strength) || !is.null(filter_fc_strength)) {
    cat("\nFilters applied:\n")
    if (!is.null(filter_pattern)) cat("  Pattern:", filter_pattern, "\n")
    if (!is.null(filter_loading_strength)) cat("  Loading:", filter_loading_strength, "\n")
    if (!is.null(filter_fc_strength)) cat("  FC:", filter_fc_strength, "\n")
  }
  
  cat("\nGenes per contrast:\n")
  for (i in seq_along(gene_sets)) {
    cat(sprintf("  %-15s: %4d\n", names(gene_sets)[i], length(gene_sets[[i]])))
  }
  
  # Step 3: Build universe
  all_genes <- unique(unlist(gene_sets))
  if (length(all_genes) == 0) stop("No genes after filtering!")
  cat("\nTotal unique genes:", length(all_genes), "\n")
  
  # Step 4: Create binary matrix
  upset_data <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  for (contrast in contrasts) {
    upset_data[[contrast]] <- upset_data$gene_id %in% gene_sets[[contrast]]
  }
  
  # Step 5: Build title
  title_parts <- "Driver Genes Overlap"
  if (!is.null(filter_pattern) || !is.null(filter_loading_strength) || !is.null(filter_fc_strength)) {
    filters <- c()
    if (!is.null(filter_pattern)) filters <- c(filters, gsub(" \\(.*\\)", "", filter_pattern))
    if (!is.null(filter_loading_strength)) filters <- c(filters, filter_loading_strength)
    if (!is.null(filter_fc_strength)) filters <- c(filters, gsub(" \\(.*\\)", "", filter_fc_strength))
    title_parts <- paste0(title_parts, " (", paste(filters, collapse = ", "), ")")
  }
  
  # Step 6: Generate plot — NO LABELS, NO EXTRA LAYERS
  p <- upset(
    data = upset_data,
    intersect = contrasts,
    width_ratio = 0.4,
    sort_intersections = "descending",
    n_intersections = top_n_intersections,
    themes = upset_default_themes(
      text = element_text(size = 12),
      strip_text = element_text(size = 11, face = "bold")
    )
  ) +
    labs(title = title_parts) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  cat("\nUpSet plot created successfully!\n")
  return(p)
}




# ============================================================================
# PART 7: Focused Gene selection and Visualization Functions
# ============================================================================

#' Get Intersection Gene Lists
#'
#' Extract specific intersection gene lists from UpSet analysis,
#' supporting any identifier column (e.g., gene_id, gene_name).
#'
#' @param driver_genes_list List of driver genes data frames.
#' @param contrasts Vector of contrast names (must be names in driver_genes_list).
#' @param intersection_type Type of intersection to extract.
#'   Options: 
#'   - "union": all genes appearing in any contrast
#'   - "intersection": genes common to all contrasts
#'   - "unique_to_X": genes only in contrast X (replace X with actual name)
#' @param column Name of the column to extract (e.g., "gene_id", "gene_name"). Default: "gene_id".
#' @return A character vector of values from the specified column.
#' @examples
#' # Get union of gene names
#' get_intersection_genes(results$driver_genes, c("c2_vs_c1", "c3_vs_c2"), 
#'                        intersection_type = "union", column = "gene_name")
get_intersection_genes <- function(driver_genes_list, 
                                  contrasts,
                                  intersection_type = "union",
                                  column = "gene_id") {
  
  # Validate inputs
  if (is.null(contrasts) || length(contrasts) == 0) {
    stop("At least one contrast must be provided.")
  }
  missing_contrasts <- contrasts[!contrasts %in% names(driver_genes_list)]
  if (length(missing_contrasts) > 0) {
    stop("Contrast(s) not found in driver_genes_list: ", paste(missing_contrasts, collapse = ", "))
  }
  
  # Check that all data frames have the requested column
  for (cn in contrasts) {
    df <- driver_genes_list[[cn]]
    if (!column %in% names(df)) {
      stop("Column '", column, "' not found in contrast '", cn, "'. Available columns: ", 
           paste(names(df), collapse = ", "))
    }
  }
  
  # Extract the specified column from each contrast
  gene_sets <- lapply(contrasts, function(x) {
    vals <- driver_genes_list[[x]][[column]]
    # Ensure it's a character vector (handles factors, etc.)
    as.character(vals)
  })
  names(gene_sets) <- contrasts
  
  # Compute intersection/union/etc.
  if (intersection_type == "union") {
    result <- unique(unlist(gene_sets))
    
  } else if (intersection_type == "intersection") {
    if (length(gene_sets) == 1) {
      result <- unique(gene_sets)
    } else {
      result <- Reduce(intersect, gene_sets)
    }
    
  } else if (grepl("^unique_to_", intersection_type)) {
    target <- gsub("unique_to_", "", intersection_type)
    if (!target %in% contrasts) {
      stop("Target contrast '", target, "' not in provided contrasts: ", paste(contrasts, collapse = ", "))
    }
    others <- setdiff(contrasts, target)
    target_vals <- gene_sets[[target]]
    other_vals <- unlist(gene_sets[others])
    result <- setdiff(target_vals, other_vals)
    
  } else {
    stop("Unknown intersection_type: '", intersection_type, "'. Supported: 'union', 'intersection', 'unique_to_X'.")
  }
  
  return(result)
}


#' Plot Gene Loading Vectors in PC Space with Pagination
#' 
#' Visualizes gene loading vectors decomposed into PC1 and PC2 components,
#' with automatic pagination if genes exceed the layout capacity
#' 
#' @param gene_list Vector of gene IDs or gene names to plot
#' @param pca_result PCA result object from prcomp()
#' @param gene_annotations Data frame with gene annotations
#' @param arrow_size Size of arrow heads (default: 0.3)
#' @param arrow_width Width of arrow lines (default: 1)
#' @param pc1_color Color for PC1 component arrow (default: "#E41A1C" red)
#' @param pc2_color Color for PC2 component arrow (default: "gray60")
#' @param total_color Color for total vector arrow (default: "black")
#' @param use_gene_names Logical, whether gene_list contains gene names (default: TRUE)
#' @param ncol Number of columns in facet layout (default: 4)
#' @param nrow Number of rows in facet layout (default: 4)
#' @return List of ggplot objects (one per page)
plot_gene_loading_vectors <- function(gene_list,
                                           pca_result,
                                           gene_annotations,
                                           arrow_size = 0.3,
                                           arrow_width = 1,
                                           pc1_color = "#E41A1C",
                                           pc2_color = "gray60",
                                           total_color = "black",
                                           use_gene_names = TRUE,
                                           ncol = 4,
                                           nrow = 4) {
  
  # --- Step 1: Match gene_list to gene_id ---
  if (use_gene_names) {
    # Convert gene names to gene IDs
    gene_map <- gene_annotations %>%
      filter(gene_name %in% gene_list) %>%
      dplyr::select(gene_id, gene_name) %>%
      distinct()
    
    if (nrow(gene_map) == 0) {
      stop("No genes found matching the provided gene names!")
    }
    
    missing_genes <- setdiff(gene_list, gene_map$gene_name)
    if (length(missing_genes) > 0) {
      warning("The following genes were not found: ", 
              paste(missing_genes, collapse = ", "))
    }
    
    # Preserve original order
    gene_map <- gene_map[match(intersect(gene_list, gene_map$gene_name), gene_map$gene_name), ]
    gene_ids <- gene_map$gene_id
    
  } else {
    # Assume gene_list contains gene IDs
    gene_ids <- gene_list
    
    gene_map <- gene_annotations %>%
      filter(gene_id %in% gene_ids) %>%
      dplyr::select(gene_id, gene_name) %>%
      distinct()
    
    if (nrow(gene_map) == 0) {
      stop("No genes found matching the provided gene IDs!")
    }
    
    # Preserve original order
    gene_map <- gene_map[match(intersect(gene_ids, gene_map$gene_id), gene_map$gene_id), ]
  }
  
  # --- Step 2: Calculate pagination ---
  genes_per_page <- ncol * nrow
  n_genes <- nrow(gene_map)
  n_pages <- ceiling(n_genes / genes_per_page)
  
  cat("Total genes:", n_genes, "| Layout:", ncol, "x", nrow, 
      "| Genes per page:", genes_per_page, "| Pages:", n_pages, "\n")
  
  # --- Step 3: Extract all gene loadings ---
  loadings <- pca_result$rotation
  
  all_gene_loadings <- data.frame(
    gene_id = rownames(loadings),
    PC1 = loadings[, "PC1"],
    PC2 = loadings[, "PC2"],
    stringsAsFactors = FALSE
  ) %>%
    filter(gene_id %in% gene_map$gene_id) %>%
    left_join(gene_map, by = "gene_id") %>%
    mutate(
      neg_PC1 = -PC1,
      magnitude = sqrt(neg_PC1^2 + PC2^2)
    )
  
  # Calculate global max for consistent scales across pages
  global_max_abs <- max(abs(c(all_gene_loadings$neg_PC1, all_gene_loadings$PC2))) * 1.2
  
  # --- Step 4: Create plots for each page ---
  plot_list <- list()
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * genes_per_page + 1
    end_idx <- min(page * genes_per_page, n_genes)
    
    genes_this_page <- gene_map$gene_id[start_idx:end_idx]
    
    cat("  Page", page, ": genes", start_idx, "to", end_idx, 
        "(", length(genes_this_page), "genes )\n")
    
    # Filter loadings for this page
    gene_loadings_page <- all_gene_loadings %>%
      filter(gene_id %in% genes_this_page)
    
    # Create arrow data for this page
    arrow_data_list <- lapply(1:nrow(gene_loadings_page), function(i) {
      gene_row <- gene_loadings_page[i, ]
      
      data.frame(
        gene_name = rep(gene_row$gene_name, 3),
        component = c("PC1 (-PC1)", "PC2", "Total"),
        x_start = c(0, 0, 0),
        y_start = c(0, 0, 0),
        x_end = c(gene_row$neg_PC1, 0, gene_row$neg_PC1),
        y_end = c(0, gene_row$PC2, gene_row$PC2),
        stringsAsFactors = FALSE
      )
    })
    
    arrow_data <- do.call(rbind, arrow_data_list)
    
    # Set factor levels for proper ordering
    arrow_data$component <- factor(
      arrow_data$component,
      levels = c("Total", "PC1 (-PC1)", "PC2")
    )
    
    # Ensure gene_name is factor with correct order for faceting
    gene_loadings_page$gene_name <- factor(
      gene_loadings_page$gene_name,
      levels = gene_loadings_page$gene_name
    )
    arrow_data$gene_name <- factor(
      arrow_data$gene_name,
      levels = levels(gene_loadings_page$gene_name)
    )
    
    # Create plot
    p <- ggplot() +
      # Origin axes
      geom_hline(yintercept = 0, linetype = "solid", 
                 color = "gray80", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "solid", 
                 color = "gray80", linewidth = 0.5) +
      # Arrow vectors
      geom_segment(data = arrow_data,
                   aes(x = x_start, y = y_start, 
                       xend = x_end, yend = y_end,
                       color = component, 
                       linetype = component,
                       linewidth = component),
                   arrow = arrow(length = unit(arrow_size, "cm"), 
                                type = "closed")) +
      # Component labels
      geom_text(data = gene_loadings_page,
                aes(x = neg_PC1/2, y = -global_max_abs * 0.15, 
                    label = sprintf("%.4f", neg_PC1)),
                size = 2.5, color = pc1_color, fontface = "bold") +
      geom_text(data = gene_loadings_page,
                aes(x = -global_max_abs * 0.15, y = PC2/2, 
                    label = sprintf("%.4f", PC2)),
                size = 2.5, color = pc2_color, fontface = "bold", angle = 90) +
      # Scales
      scale_color_manual(
        name = "Component",
        values = c("Total" = total_color, 
                   "PC1 (-PC1)" = pc1_color, 
                   "PC2" = pc2_color)
      ) +
      scale_linetype_manual(
        name = "Component",
        values = c("Total" = "solid", 
                   "PC1 (-PC1)" = "solid", 
                   "PC2" = "solid")
      ) +
      scale_linewidth_manual(
        name = "Component",
        values = c("Total" = arrow_width * 1.2, 
                   "PC1 (-PC1)" = arrow_width, 
                   "PC2" = arrow_width)
      ) +
      # Facet by gene
      facet_wrap(~ gene_name, ncol = ncol, nrow = nrow) +
      # Labels and theme
      labs(
        x = "Disease Axis (-PC1)",
        y = "PC2 (Orthogonal Axis)",
        title = "Gene Loading Vector Decomposition",
        subtitle = sprintf(
          "Page %d/%d | Genes %d-%d of %d | Red: PC1, Gray: PC2, Black: Total",
          page, n_pages, start_idx, end_idx, n_genes
        )
      ) +
      coord_fixed(ratio = 1, xlim = c(-global_max_abs, global_max_abs), 
                  ylim = c(-global_max_abs, global_max_abs)) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "gray95", color = "gray80"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )
    
    plot_list[[page]] <- p
  }
  
  cat("Created", length(plot_list), "plot(s)\n")
  
  return(plot_list)
}


#' Plot Multiple Gene Lists with Automatic Pagination
#' 
#' Process multiple gene lists and combine all plots
#' 
#' @param gene_lists Named list of gene vectors
#' @param ... Other parameters passed to plot_gene_loading_vectors_auto
#' @return Named list of plot lists
plot_multiple_gene_lists <- function(gene_lists, ...) {
  
  all_plots <- list()
  
  for (list_name in names(gene_lists)) {
    cat("\n=== Processing gene list:", list_name, "===\n")
    
    plots <- plot_gene_loading_vectors(
      gene_list = gene_lists[[list_name]],
      ...
    )
    
    all_plots[[list_name]] <- plots
  }
  
  # Flatten to single list
  flat_plots <- unlist(all_plots, recursive = FALSE)
  
  cat("\n=== Total plots created:", length(flat_plots), "===\n")
  
  return(flat_plots)
}


# ============================================================================
# USAGE
# ============================================================================
# load results from driver gene calculation pipeline

BASE_DIR <- "/disk2/cai113/data/stateTrans"
load(file.path(BASE_DIR, "08.eigengene", "driver_gene_selection_results.rdata"))
# draw plots
p1 <- plot_gene_pc_space(
  pca_result = results$project_data$pca_result,
  highlight_genes = NULL,
  highlight_color = "#E41A1C",
  point_size = 1,
  point_alpha = 0.3
)

p2 <- plot_sample_cone_space(
  sample_info_pcs = results$project_data$sample_info_pcs,
  cone_angles = results$cone_angles,
  critical_points = results$project_data$critical_points,
  cone_alpha = 0.15,
  point_size = 3
)

p3 <- plot_gene_cone_space(
  pca_result = results$project_data$pca_result,
  cone_angles = results$cone_angles,
  cone_alpha = 0.15,
  point_size = 0.8,
  point_alpha = 0.3
)

p4 <- plot_gene_polar_cone(
  eigengenes = results$eigengenes,
  cone_angles = results$cone_angles,
  max_radius_quantile = 1,
  point_size = 1,
  point_alpha = 0.4
)

p5.1 <- plot_driver_genes_in_cone(
  eigengenes_in_cone = results$eigengenes_in_cone,
  driver_genes_list = results$driver_genes,
  cone_angles = results$cone_angles,
  contrasts = "ALL",
  point_size = 1.5,
  bg_alpha = 0.2
)
p5.2 <- plot_driver_genes_in_cone(
  eigengenes_in_cone = results$eigengenes_in_cone,
  driver_genes_list = results$driver_genes,
  cone_angles = results$cone_angles,
  contrasts = c("c1_vs_c1_star"),
  point_size = 1.5,
  bg_alpha = 0.2
)
p5.3 <- plot_driver_genes_in_cone(
  eigengenes_in_cone = results$eigengenes_in_cone,
  driver_genes_list = results$driver_genes,
  cone_angles = results$cone_angles,
  contrasts = c("c2_vs_c1"),
  point_size = 1.5,
  bg_alpha = 0.2
)
p5.4 <- plot_driver_genes_in_cone(
  eigengenes_in_cone = results$eigengenes_in_cone,
  driver_genes_list = results$driver_genes,
  cone_angles = results$cone_angles,
  contrasts = c("c3_vs_c2"),
  point_size = 1.5,
  bg_alpha = 0.2
)

p6.1 <- plot_driver_genes_upset(
  driver_genes_list = results$driver_genes,
  contrasts = c("c1_vs_c1_star", "c2_vs_c1", "c3_vs_c2"),
  filter_pattern = NULL,
  filter_loading_strength = NULL,
  filter_fc_strength = NULL,
  top_n_intersections = 15
)
p6.2 <- plot_driver_genes_upset(
  driver_genes_list = results$driver_genes,
  contrasts = c("c1_vs_c1_star", "c2_vs_c1", "c3_vs_c2"),
  filter_pattern = "Upregulated-Cis (Strong Pro-disease)",
  filter_loading_strength = NULL,
  filter_fc_strength = NULL,
  top_n_intersections = 15
)
p6.3 <- plot_driver_genes_upset(
  driver_genes_list = results$driver_genes,
  contrasts = c("c1_vs_c1_star", "c2_vs_c1", "c3_vs_c2"),
  filter_pattern = "Downregulated-Trans (Pro-disease)",
  filter_loading_strength = NULL,
  filter_fc_strength = NULL,
  top_n_intersections = 15
)

gene_list.1 <- get_intersection_genes(
                driver_genes_list = results$driver_genes,
                contrasts = c("c1_vs_c1_star", "c2_vs_c1"),
                intersection_type = "intersection",
                column = "gene_name"
              )
gene_list.2 <- get_intersection_genes(
                driver_genes_list = results$driver_genes,
                contrasts = c("c2_vs_c1", "c3_vs_c2"),
                intersection_type = "intersection",
                column = "gene_name"
              )
p7 <- plot_multiple_gene_lists(
        gene_lists = list(
          "c1_vs_c1_star & c2_vs_c1" = gene_list.1,
          "c2_vs_c1 & c3_vs_c2" = gene_list.2
        ),
        pca_result = results$project_data$pca_result,
        gene_annotations = results$project_data$gene_annotations,
        arrow_size = 0.25,
        arrow_width = 0.8,
        pc1_color = "#E41A1C",
        pc2_color = "gray60",
        total_color = "black",
        use_gene_names = TRUE,
        ncol = 4,
        nrow = 4
      )

# Save results
cat("\nSaving results...\n")
output_dir <- file.path(BASE_DIR, "08.eigengene/figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(filename = file.path(output_dir, "01.gene_pc_space.pdf"), plot = p1, width = 6, height = 6)

ggsave(filename = file.path(output_dir, "02.sample_cone_space.pdf"), plot = p2, width = 12, height = 5)

ggsave(filename = file.path(output_dir, "03.gene_cone_space.pdf"), plot = p3, width = 6, height = 6)

ggsave(filename = file.path(output_dir, "04.gene_polar_cone.pdf"), plot = p4, width = 6, height = 6)

plots_5 <- list(p5.1, p5.2, p5.3, p5.4)
pdf(file.path(output_dir, "05.driver_genes_cone_plots.pdf"), width = 6, height = 6)
walk(plots_5, print)
dev.off()

plots_6 <- list(p6.1, p6.2, p6.3)
pdf(file.path(output_dir, "06.driver_genes_upset_plots.pdf"), width = 8, height = 6)
walk(plots_6, print)
dev.off()

pdf(file.path(output_dir, "07.driver_genes_multiple_gene_lists.pdf"), width = 14, height = 14)
walk(p7, print)
dev.off()

cat("Plots saved to:", output_dir, "\n")