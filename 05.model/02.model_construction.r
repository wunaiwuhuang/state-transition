# ============================================================================
# State-Transition Model with Double-Well Quasi-Potential
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

# ============================================================================
# PART 1: Critical Point Estimation and Model Construction
# ============================================================================

#' Estimate Critical Points Using K-means Clustering
#' 
#' @param x_wt Numeric vector of PC values for wild-type samples
#' @param x_ko Numeric vector of PC values for knockout samples
#' @param n_clusters Number of clusters for KO samples (default: 3)
#' @return List containing critical points and clustering results
estimate_critical_points <- function(x_wt, x_ko, n_clusters = 3) {
  
  # c1* (control/healthy reference): centroid of all WT samples
  c1_star <- mean(x_wt)
  
  # K-means clustering on KO samples to identify 3 clusters
  set.seed(20030622)  # for reproducibility
  kmeans_result <- kmeans(x_ko, centers = n_clusters, nstart = 25)
  
  # Sort cluster centers (higher PC = healthier, lower PC = more diseased)
  sorted_centers <- sort(kmeans_result$centers, decreasing = TRUE)
  
  # Assign critical points
  # c1: centroid of cluster closest to health (highest PC value)
  # c3: centroid of cluster in disease state (lowest PC value)
  c1 <- sorted_centers[1]
  c3 <- sorted_centers[3]
  
  # Get cluster membership sorted by center values
  cluster_labels <- kmeans_result$cluster
  center_order <- order(kmeans_result$centers, decreasing = TRUE)
  c2_cluster_samples <- x_ko[cluster_labels == center_order[2]]
  
  return(list(
    c1_star = c1_star,
    c1 = c1,
    c3 = c3,
    c2_cluster_samples = c2_cluster_samples,
    kmeans_result = kmeans_result,
    sorted_centers = sorted_centers
  ))
}


#' Calculate Quasi-Potential Derivative
#' 
#' @param x Position in state-space
#' @param alpha Scaling parameter
#' @param c1 Critical point 1 (near health)
#' @param c2 Critical point 2 (unstable)
#' @param c3 Critical point 3 (disease)
#' @return Derivative value
potential_derivative <- function(x, alpha, c1, c2, c3) {
  alpha * (x - c1) * (x - c2) * (x - c3)
}


#' Calculate Quasi-Potential by Integration
#' 
#' @param x Position in state-space
#' @param alpha Scaling parameter
#' @param c1 Critical point 1
#' @param c2 Critical point 2
#' @param c3 Critical point 3
#' @return Potential value
calculate_potential <- function(x, alpha, c1, c2, c3) {
  # U_p = integral of U_p' dx
  # Expanding: U_p' = alpha * (x - c1)(x - c2)(x - c3)
  # After expansion and integration:
  a4 <- alpha / 4
  a3 <- -alpha * (c1 + c2 + c3) / 3
  a2 <- alpha * (c1*c2 + c1*c3 + c2*c3) / 2
  a1 <- -alpha * c1 * c2 * c3
  
  potential <- a4 * x^4 + a3 * x^3 + a2 * x^2 + a1 * x
  return(potential)
}


#' Calculate Boltzmann Ratio
#' 
#' @param U_c1 Potential at c1
#' @param U_c3 Potential at c3
#' @param kBT Temperature parameter (default: 1)
#' @return Boltzmann ratio
calculate_boltzmann_ratio <- function(U_c1, U_c3, kBT = 1) {
  exp(-(U_c3 - U_c1) / kBT)
}


#' Find Optimal c2 by Matching Observed Disease/Health Ratio
#' 
#' @param critical_points List from estimate_critical_points
#' @param sample_info Data frame with sample information
#' @param alpha Scaling parameter
#' @param kBT Temperature parameter
#' @param n_grid Number of grid points for c2 search
#' @return List with optimal c2 and search results
find_optimal_c2 <- function(critical_points, sample_info, alpha, kBT, n_grid = 200) {
  #------------ find KO group optimal c2
  c1 <- critical_points$c1
  c3 <- critical_points$c3
  # Define search range for c2 (from c1 to c3)
  c2_range <- seq(from = max(c1, c3), to = min(c1, c3), length.out = n_grid)
  # Calculate observed disease ratio in KO samples. Less than or equal to week 6 is healthy.
  # ko_indices <- which(sample_info$genotype == "Ncf2")
  # healthy_samples <- sample_info$week[ko_indices] <= 6
  # disease_samples <- sample_info$week[ko_indices] > 6
  # observed_ratio <- sum(disease_samples) / length(ko_indices)
  observed_ratio <- 0.99  # NCF2 KO mice is 100% diseased, but there i supposed 99%
  target_br <- observed_ratio / (1 - observed_ratio)
  # Initialize storage
  br_values <- numeric(n_grid)
  U_c1_values <- numeric(n_grid)
  U_c3_values <- numeric(n_grid)
  # Calculate BR for each c2
  for (i in seq_along(c2_range)) {
    c2_test <- c2_range[i]
    U_c1 <- calculate_potential(c1, alpha, c1, c2_test, c3)
    U_c3 <- calculate_potential(c3, alpha, c1, c2_test, c3)
    br_values[i] <- calculate_boltzmann_ratio(U_c1, U_c3, kBT)
    U_c1_values[i] <- U_c1
    U_c3_values[i] <- U_c3
  }
  # Find c2 that best matches target BR in KO
  br_diff <- abs(br_values - target_br)
  best_idx <- which.min(br_diff)
  c2_optimal <- c2_range[best_idx]

  #------------ find WT group optimal c2
  c1_star <- critical_points$c1_star
  c3 <- critical_points$c3
  c2_wt_range <- seq(from = max(c1_star, c3), to = min(c1_star, c3), length.out = n_grid)
  target_br_wt <- (1 - observed_ratio) / observed_ratio # for WT
  br_values_wt <- numeric(n_grid)
  U_c1_values_wt <- numeric(n_grid)
  U_c3_values_wt <- numeric(n_grid)
  for (i in seq_along(c2_wt_range)) {
    c2_wt_test <- c2_wt_range[i]

    U_c1_wt <- calculate_potential(c1_star, alpha, c1_star, c2_wt_test, c3)
    U_c3_wt <- calculate_potential(c3, alpha, c1_star, c2_wt_test, c3)

    br_values_wt[i] <- calculate_boltzmann_ratio(U_c1_wt, U_c3_wt, kBT)

    U_c1_values_wt[i] <- U_c1_wt
    U_c3_values_wt[i] <- U_c3_wt
  }
  # Find c2 that best matches target BR in WT
  br_diff_wt <- abs(br_values_wt - target_br_wt)
  best_idx_wt <- which.min(br_diff_wt)
  c2_optimal_wt <- c2_wt_range[best_idx_wt]
  
  # Return results
  return(list(
    c2_optimal = c2_optimal,
    c2_optimal_wt = c2_optimal_wt,
    c2_range = c2_range,
    c2_wt_range = c2_wt_range,
    br_values = br_values,
    br_values_wt = br_values_wt,
    target_br = target_br,
    target_br_wt = target_br_wt,
    best_idx = best_idx,
    best_idx_wt = best_idx_wt,
    observed_ratio = observed_ratio,
    U_c1_values = U_c1_values,
    U_c1_values_wt = U_c1_values_wt,
    U_c3_values = U_c3_values,
    U_c3_values_wt = U_c3_values_wt
  ))
}


#' Main Function to Build State-Transition Model
#' 
#' @param x_wt WT sample PC values
#' @param x_ko KO sample PC values
#' @param sample_info Sample information data frame
#' @param alpha Scaling parameter
#' @param kBT Temperature parameter
#' @return Complete model results
build_state_transition_model <- function(x_wt, x_ko, sample_info, alpha = 4.85e-8, kBT = 1) {
  
  # Step 1: Estimate critical points
  critical_points <- estimate_critical_points(x_wt, x_ko)
  
  # Step 2: Find optimal c2
  c2_results <- find_optimal_c2(critical_points, sample_info, alpha, kBT)
  
  # Step 3: Calculate final potentials
  # KO group
  c1 <- critical_points$c1
  c2 <- c2_results$c2_optimal
  c3 <- critical_points$c3  
  U_c1_ko <- calculate_potential(c1, alpha, c1, c2, c3)
  U_c2_ko <- calculate_potential(c2, alpha, c1, c2, c3)
  U_c3_ko <- calculate_potential(c3, alpha, c1, c2, c3)
  energy_barrier_ko <- U_c2_ko - U_c1_ko # energy barriers
  # WT group
  c1_star <- critical_points$c1_star
  c2_wt <- c2_results$c2_optimal_wt
  U_c1_wt <- calculate_potential(c1_star, alpha, c1_star, c2_wt, c3)
  U_c2_wt <- calculate_potential(c2_wt, alpha, c1_star, c2_wt, c3)
  U_c3_wt <- calculate_potential(c3, alpha, c1_star, c2_wt, c3)
  energy_barrier_wt <- U_c2_wt - U_c1_wt # energy barriers

  # step 4: calculate barrier ratio
  barrier_down_ratio <- (energy_barrier_wt - energy_barrier_ko) / energy_barrier_wt

  # Return complete model results
  return(list(
    critical_points = list(
      c1_star = c1_star,
      c2_wt = c2_wt,
      c1 = c1,
      c2 = c2,
      c3 = c3
    ),
    c2_search = c2_results,
    energy_barrier_ko = energy_barrier_ko,
    energy_barrier_wt = energy_barrier_wt,
    barrier_down_ratio = barrier_down_ratio,
    alpha = alpha,
    kBT = kBT,
    x_wt = x_wt,
    x_ko = x_ko
  ))
}

# ============================================================================
# PART 2: Visualization Functions
# ============================================================================

#' Plot 1: Final Potential Curves with Sample Points
#' 
#' @param model_results Results from build_state_transition_model
#' @return ggplot object
plot_potential_curves <- function(model_results) {
  
  cp <- model_results$critical_points
  alpha <- model_results$alpha
  x_wt <- model_results$x_wt
  x_ko <- model_results$x_ko
  
  # Generate x range (from high PC to low PC, left to right)
  x_range <- seq(from = max(c(x_wt, x_ko)) + 5,
                 to = min(c(x_wt, x_ko)) - 5,
                 length.out = 500)
  
  # Calculate KO potential with lower barrier (disease-prone)
  U_ko <- calculate_potential(x_range, alpha, cp$c1, cp$c2, cp$c3)
  
  # Calculate WT potential with higher barrier (resistant to disease)
  # Use c1* as the healthy state, and assume c2_wt and c3 for potential structure
  # Scale alpha by a factor to create much higher energy barrier for WT
  alpha_wt <- alpha  # use same alpha for simplicity
  U_wt <- calculate_potential(x_range, alpha_wt, cp$c1_star, cp$c2_wt, cp$c3) # use same c3

  # Calculate sample potentials
  U_ko_samples <- calculate_potential(x_ko, alpha, cp$c1, cp$c2, cp$c3)
  U_wt_samples <- calculate_potential(x_wt, alpha_wt, cp$c1_star, cp$c2_wt, cp$c3)

  # Normalize potentials to have minimum at 0
  # U_ko_samples <- U_ko_samples - min(U_ko)
  # U_wt_samples <- U_wt_samples - min(U_wt)
  # U_ko <- U_ko - min(U_ko)
  # U_wt <- U_wt - min(U_wt)
  
  # Create data frame
  df_curve <- data.frame(
    PC_value = rep(x_range, 2),
    Potential = c(U_ko, U_wt),
    Group = rep(c("KO", "WT"), each = length(x_range))
  )
  
  df_samples <- data.frame(
    PC_value = c(x_ko, x_wt),
    Potential = c(U_ko_samples, U_wt_samples),
    Group = rep(c("KO", "WT"), c(length(x_ko), length(x_wt)))
  )
  
  df_critical <- data.frame(
    PC_value = c(cp$c1_star, cp$c1, cp$c2, cp$c3),
    Label = c("c1* (WT)", "c1 (Health-like)", "c2 (Unstable)", "c3 (Disease)")
  )
  
  # Plot
  p <- ggplot() +
    geom_line(data = df_curve, aes(x = PC_value, y = Potential, 
                                   color = Group, linetype = Group), 
              linewidth = 1) +
    geom_point(data = df_samples, aes(x = PC_value, y = Potential, 
                                      color = Group, shape = Group),
               size = 2.5, alpha = 0.6) +
    geom_vline(data = df_critical, aes(xintercept = PC_value),
               linetype = "dashed", color = "gray40", alpha = 0.5) +
    geom_text(data = df_critical, aes(x = PC_value, y = max(df_curve$Potential) * 0.95,
                                      label = Label),
              angle = 90, vjust = -0.3, hjust = 1, size = 3.5, color = "gray20") +
    scale_color_manual(values = c("KO" = "#E41A1C", "WT" = "#377EB8")) +
    scale_shape_manual(values = c("KO" = 16, "WT" = 17)) +
    scale_linetype_manual(values = c("KO" = "solid", "WT" = "solid")) +
    labs(x = "PC Value (Health ← → Disease)",
         y = "Quasi-Potential (a.u.)",
         title = "Double-Well Quasi-Potential Landscape",
         color = "Genotype", shape = "Genotype", linetype = "Genotype") +
    theme_minimal() +
    scale_x_reverse() +        # x axis from high to low PC
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}


#' Plot 2: Sample Density Distribution
#' 
#' @param model_results Results from build_state_transition_model
#' @return ggplot object
plot_sample_density <- function(model_results) {
  
  x_wt <- model_results$x_wt
  x_ko <- model_results$x_ko
  
  df_density <- data.frame(
    PC_value = c(x_ko, x_wt),
    Group = rep(c("KO", "WT"), c(length(x_ko), length(x_wt)))
  )
  
  p <- ggplot(df_density, aes(x = PC_value, fill = Group, color = Group)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    geom_rug(aes(color = Group), alpha = 0.6, length = unit(0.05, "npc")) +
    scale_fill_manual(values = c("KO" = "#E41A1C", "WT" = "#377EB8")) +
    scale_color_manual(values = c("KO" = "#E41A1C", "WT" = "#377EB8")) +
    scale_x_reverse() +  # High PC (health) on left, low PC (disease) on right
    labs(x = "PC Value (Health ← → Disease)",
         y = "Density",
         title = "Sample Distribution Along Disease Progression Axis",
         fill = "Genotype", color = "Genotype") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}


#' Plot 3: Representative Potential Curves for Different c2 Values
#' 
#' @param model_results Results from build_state_transition_model
#' @param n_curves Number of representative curves to show
#' @return ggplot object
plot_c2_exploration <- function(model_results, n_curves = 20) {
  
  alpha <- model_results$alpha
  cp <- model_results$critical_points
  c2_search <- model_results$c2_search
  c2_optimal_idx <- c2_search$best_idx
  
  # Select representative c2 values
  total_points <- length(c2_search$c2_range)
  selected_indices <- round(seq(1, total_points, length.out = n_curves))
  selected_indices <- unique(c(selected_indices, c2_optimal_idx, 1, total_points))

  x_range <- seq(from = max(c(model_results$x_wt, model_results$x_ko)) + 5, to = min(c(model_results$x_wt, model_results$x_ko)) - 5, length.out = 500)
  
  # Calculate potentials for each selected c2
  df_list <- list()
  for (i in seq_along(selected_indices)) {
    idx <- selected_indices[i]
    c2_test <- c2_search$c2_range[idx]
    U_test <- calculate_potential(x_range, alpha, cp$c1, c2_test, cp$c3)
    # U_test <- U_test - min(U_test)
    
    df_list[[i]] <- data.frame(
      PC_value = x_range,
      Potential = U_test,
      c2_value = c2_test,
      c2_label = sprintf("c2 = %.2f (BR = %.2f)", 
                        c2_test, c2_search$br_values[idx]),
      is_optimal = (idx == c2_optimal_idx)
    )
  }
  
  df_curves <- do.call(rbind, df_list)
  df_optimal <- df_curves[df_curves$is_optimal, ]
  
  p <- ggplot() +
    geom_line(data = df_curves, 
              aes(x = PC_value, y = Potential, 
                  group = c2_label, color = c2_value),
              alpha = 0.5, linewidth = 0.8) +
    geom_line(data = df_optimal,
              aes(x = PC_value, y = Potential),
              color = "red", linewidth = 1.5) +
    scale_color_gradient(low = "#4575B4", high = "#D73027",
                        name = "c2 Value") +
    labs(x = "PC Value (Health ← → Disease)",
         y = "Quasi-Potential (a.u.)",
         title = "Exploration of Potential Landscapes with Different c2 Values",
         subtitle = "Red curve: optimal c2 matching observed disease ratio") +
    theme_minimal() +
    scale_x_reverse() +        # x axis from high to low PC
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}


#' Plot 4: C2 Selection Rationale
#' 
#' @param model_results Results from build_state_transition_model
#' @return ggplot object
plot_c2_selection <- function(model_results) {
  
  c2_search <- model_results$c2_search
  
  df_br <- data.frame(
    c2_value = c2_search$c2_range,
    BR = c2_search$br_values,
    BR_diff = abs(c2_search$br_values - c2_search$target_br)
  )
  
  optimal_c2 <- c2_search$c2_optimal
  optimal_br <- c2_search$br_values[c2_search$best_idx]
  
  p <- ggplot(df_br, aes(x = c2_value, y = BR)) +
    geom_line(color = "#377EB8", linewidth = 1) +
    geom_point(color = "#377EB8", size = 2, alpha = 0.5) +
    geom_hline(yintercept = c2_search$target_br, 
               linetype = "dashed", color = "#E41A1C", linewidth = 1) +
    geom_point(aes(x = optimal_c2, y = optimal_br),
               color = "#E41A1C", size = 5, shape = 18) +
    annotate("text", x = optimal_c2, y = optimal_br * 1.2,
             label = sprintf("Optimal c2 = %.2f\nBR = %.2f", optimal_c2, optimal_br),
             color = "#E41A1C", fontface = "bold", size = 4) +
    annotate("text", x = max(df_br$c2_value) * 0.95, y = c2_search$target_br * 1.1,
             label = sprintf("Target BR = %.2f\n(Observed ratio: %.2f)", 
                           c2_search$target_br, c2_search$observed_ratio),
             color = "#E41A1C", hjust = 1, size = 3.5) +
    labs(x = "c2 Value",
         y = "Boltzmann Ratio (BR)",
         title = "Selection of Optimal c2 Value",
         subtitle = "Matching predicted BR to observed disease/health ratio") +
    theme_minimal() +
    scale_x_reverse() +        # x axis from high to low PC
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# Load your data
basedir <- "/disk2/cai113/data/stateTrans"
load(file.path(basedir, "05.model/sample_info_with_PCs.rdata"))

# Prepare data
ko_samples <- sample_info_pcs$genotype == "Ncf2"
wt_samples <- sample_info_pcs$genotype == "WT"
x_wt <- x_cgd[wt_samples]
x_ko <- x_cgd[ko_samples]

# Estimate alpha
evaluate_alpha_sweep <- function(x_wt, x_ko, sample_info_pcs, alpha_seq = seq(5e-8, 2e-7, length.out = 40), kBT = 1) {

  x_wt = x_wt
  x_ko = x_ko
  sample_info_pcs = sample_info_pcs
  alpha_seq = alpha_seq
  kBT = kBT
  results_list <- vector("list", length(alpha_seq))

  for (i in seq_along(alpha_seq)) {
    alpha_val <- alpha_seq[i]
    # Build model with current alpha
    model_res <- build_state_transition_model(x_wt, x_ko, sample_info_pcs, alpha = alpha_val, kBT = kBT)
    # Extract metrics
    eb_ko  <- model_res$energy_barrier_ko
    eb_wt  <- model_res$energy_barrier_wt
    down_ratio  <- model_res$barrier_down_ratio  # note: this may be defined as EB_KO / EB_WT or similar
    # Compute barrier ratio (WT / KO)
    barrier_ratio <- ifelse(eb_ko > 0, eb_wt / eb_ko, NA)
    # Store
    results_list[[i]] <- data.frame(
      alpha = alpha_val,
      energy_barrier_ko = eb_ko,
      energy_barrier_wt = eb_wt,
      barrier_down_ratio = down_ratio,
      eb_ratio_wt_over_ko = barrier_ratio,
      stringsAsFactors = FALSE
    )
  }
  # Combine all rows
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  return(results_df)
}
alpha_results <- evaluate_alpha_sweep(x_wt, x_ko, sample_info_pcs, alpha_seq = seq(5e-8, 2e-7, length.out = 40), kBT = 1)
alpha_optimal <- alpha_results[alpha_results$eb_ratio_wt_over_ko <= 8 & alpha_results$eb_ratio_wt_over_ko >= 5, ]
# Choose alpha in the middle of the range, alpha = 1.4e-7

# Build model
model_results <- build_state_transition_model(x_wt, x_ko, sample_info_pcs,alpha = 1.4e-7, kBT = 1)

# Generate plots
p1 <- plot_potential_curves(model_results)
p2 <- plot_sample_density(model_results)
p3 <- plot_c2_exploration(model_results, n_curves = 10)
p4 <- plot_c2_selection(model_results)

# Display or save plots
print(p1)
print(p2)
print(p3)
print(p4)

ggsave(paste0(file.path(basedir, "05.model/"), "potential_curves.pdf"), p1, width = 8, height = 6)
ggsave(paste0(file.path(basedir, "05.model/"), "sample_density.pdf"), p2, width = 8, height = 6)
ggsave(paste0(file.path(basedir, "05.model/"), "c2_exploration.pdf"), p3, width = 8, height = 6)
ggsave(paste0(file.path(basedir, "05.model/"), "c2_selection.pdf"), p4, width = 8, height = 6)