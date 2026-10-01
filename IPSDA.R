###############################################################################
# Purpose:   The pipeline for Inverse-PSO Longitudinal Swarm Dynamics Analysis (IPSDA)
#            It targets Eco-Thermodynamics, Phase Transitions & Keystone Driver Identification
#            
# Version:   1.1
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
#            Ali Zeeshan Ijaz (ali.zeeshan@isb.nu.edu.pk)
# Created:   2026-10-01
# License:   Copyright (c) 2026 Environmental'Omics Lab, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)          # GAMM derivative estimation
  library(nlme)          # lme controls for convergence
  library(compositions)  # Compositional data analysis (rCLR)
  library(quadprog)      # Constrained quadratic programming
  library(patchwork)     # Camera-ready figure layouts
  library(viridis)       # Colour palettes
  library(ggpubr)        # Camera-ready themes & stats
})

###############################################################################
# Main pipeline wrapper
###############################################################################
run_advanced_ipsda <- function(
    count_file,
    meta_file,
    output_dir         = "./ipsda_unified_results",
    sample_id_col      = "SampleID",
    time_col           = "Day",
    subject_col        = "Pen_No",   # Random effect grouping (e.g., Pen, Patient, Subject)
    target_vars        = c(),        # Vector of metadata column names (continuous or categorical)
    window_candidates  = c(3.0, 4.0, 5.0, 6.0, 7.0), # Vector of sliding window spans to optimise
    window_step        = 1.0,        # Sliding step increment
    inertia_candidates = seq(0.1, 0.8, by = 0.1),    # Vector of swarm inertia weights (omega) to optimise
    n_modules          = 4,          # gbest dynamic clusters
    gamm_k             = 8           # Spline knots for GAMM
) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("================================================================\n")
  cat("[INFO] Initialising Fully Integrated Eco-Thermodynamic IPSDA Pipeline\n")
  cat("[INFO] Output Directory:", output_dir, "\n")
  cat("================================================================\n")
  
  ###############################################################################
  # Data loading & Robust CLR pre-processing (Supports CSV and .biom)
  ###############################################################################
  metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  # Define configurable list of 16S rRNA contaminants and off-target sequences
  contaminants <- c(
    "Eukaryota",
    "Mitochondria",
    "Chloroplast",
    "Unassigned"
  )
  
  if (requireNamespace("biomformat", quietly = TRUE) && grepl("\\.biom$", count_file, ignore.case = TRUE)) {
    biom_obj <- biomformat::read_biom(count_file)
    raw_counts <- as.matrix(biomformat::biom_data(biom_obj))
    
    # Taxonomic cleaning
    obs_md <- biomformat::observation_metadata(biom_obj)
    if (!is.null(obs_md)) {
      if ("taxonomy" %in% colnames(obs_md)) {
        if (is.list(obs_md$taxonomy)) {
          tax_paths <- sapply(obs_md$taxonomy, function(x) paste(x, collapse = ";"))
        } else {
          tax_paths <- as.character(obs_md$taxonomy)
        }
      } else {
        tax_paths <- apply(obs_md, 1, function(x) paste(na.omit(x), collapse = ";"))
      }
      
      # Clean consecutive semicolons and trim trailing/leading ones
      tax_paths <- gsub(";+", ";", tax_paths)
      tax_paths <- gsub("^;|;$", "", tax_paths)
      
      # Concatenate OTU ID with its cleaned full taxonomic path
      rownames(raw_counts) <- paste(rownames(raw_counts), tax_paths, sep = "|")
    }
    counts <- t(raw_counts)
  } else {
    counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
    if (all(metadata[[sample_id_col]] %in% colnames(counts))) counts <- t(counts)
  }
  
  # Synchronise samples between metadata and count matrix
  common_samples <- intersect(metadata[[sample_id_col]], rownames(counts))
  metadata <- metadata %>% filter(.data[[sample_id_col]] %in% common_samples) %>% arrange(.data[[time_col]])
  counts <- counts[metadata[[sample_id_col]], ]
  
  # Contaminants and off-target filtering
  contam_regex <- paste(contaminants, collapse = "|")
  is_contam <- grepl(contam_regex, colnames(counts), ignore.case = TRUE)
  
  cat(sprintf("[INFO] Removed %d contaminant/off-target taxa matching: %s\n", 
              sum(is_contam), paste(contaminants, collapse = ", ")))
  
  counts <- counts[, !is_contam, drop = FALSE]
  
  # Filter low-prevalence taxa (< 10% sample occupancy)
  prev <- colSums(counts > 0) / nrow(counts)
  counts <- counts[, prev >= 0.10]
  
  # Robust CLR (rCLR): Geometric mean over non-zero features
  rclr <- function(mat) {
    apply(mat, 1, function(row) {
      pos <- row[row > 0]
      if (length(pos) == 0) return(rep(0, length(row)))
      gmean <- exp(mean(log(pos)))
      log(pmax(row, 0.5) / gmean)
    }) %>% t()
  }
  
  clr_matrix <- rclr(counts)
  cat(sprintf("[INFO] Preprocessed matrix: %d samples, %d taxa.\n", nrow(clr_matrix), ncol(clr_matrix)))
  
  ###############################################################################
  # GAMM derivative engine (velocity & acceleration with pen/subject awareness)
  ###############################################################################
  
  cat(sprintf("[INFO] Fitting GAMMs with subject random intercepts (k=%d)...\n", gamm_k))
  time_vec <- metadata[[time_col]]
  subj_vec <- as.factor(metadata[[subject_col]])
  
  vel_mat <- matrix(0, nrow = nrow(clr_matrix), ncol = ncol(clr_matrix), dimnames = dimnames(clr_matrix))
  acc_mat <- matrix(0, nrow = nrow(clr_matrix), ncol = ncol(clr_matrix), dimnames = dimnames(clr_matrix))
  
  eps <- 1e-4
  lme_ctrl <- nlme::lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim", returnObject = TRUE)
  
  for (j in 1:ncol(clr_matrix)) {
    df_j <- data.frame(y = clr_matrix[, j], Time = time_vec, Subject = subj_vec)
    
    fit_gam <- tryCatch({
      res_gamm <- gamm(
        y ~ s(Time, k = gamm_k, bs = "cr"), 
        random = list(Subject = ~1), 
        data = df_j,
        control = lme_ctrl
      )
      res_gamm$gam
    }, error = function(e) {
      suppressWarnings(gam(y ~ s(Time, k = gamm_k, bs = "cr"), data = df_j, method = "REML"))
    })
    
    if (!is.null(fit_gam)) {
      df_eps <- df_j; df_eps$Time <- df_eps$Time + eps
      df_neg <- df_j; df_neg$Time <- df_neg$Time - eps
      
      pred_0   <- predict(fit_gam, newdata = df_j)
      pred_pos <- predict(fit_gam, newdata = df_eps)
      pred_neg <- predict(fit_gam, newdata = df_neg)
      
      vel_mat[, j] <- (pred_pos - pred_neg) / (2 * eps)
      acc_mat[, j] <- (pred_pos - 2 * pred_0 + pred_neg) / (eps^2)
    }
  }
  
  ###############################################################################
  # Joint 2D sliding-window & inertia optimisation (window size x inertia weight)
  ###############################################################################
  cat("[INFO] Executing 2D Joint Optimisation Surface (Window Size x Inertia Weight)...\n")
  
  optimisation_scores <- data.frame()
  candidate_results   <- list()
  candidate_meta      <- list()
  
  for (w_sz in window_candidates) {
    for (w_inertia in inertia_candidates) {
      
      min_t <- min(time_vec)
      max_t <- max(time_vec)
      windows <- seq(min_t, max_t - w_sz, by = window_step)
      
      results <- list()
      window_meta_list <- list()
      
      for (w_start in windows) {
        w_end <- w_start + w_sz
        w_center <- round(w_start + (w_sz / 2)) 
        w_idx <- which(time_vec >= w_start & time_vec <= w_end)
        
        if (length(w_idx) < 8) next
        
        sub_meta <- metadata[w_idx, ]
        meta_row <- data.frame(Window_Center = w_center)
        for (v in target_vars) {
          if (v %in% colnames(metadata)) {
            vec <- sub_meta[[v]]
            if (is.numeric(vec)) {
              meta_row[[v]] <- mean(vec, na.rm = TRUE)
            } else {
              vec <- vec[!is.na(vec)]
              if (length(vec) > 0) {
                tbl <- table(vec)
                meta_row[[v]] <- names(tbl)[which.max(tbl)]
              } else {
                meta_row[[v]] <- NA
              }
            }
          }
        }
        window_meta_list[[length(window_meta_list) + 1]] <- meta_row
        
        w_clr <- clr_matrix[w_idx, , drop=FALSE]
        w_vel <- vel_mat[w_idx, , drop=FALSE]
        w_acc <- acc_mat[w_idx, , drop=FALSE]
        
        hc <- hclust(dist(t(w_clr)), method = "ward.D2")
        clusters <- cutree(hc, k = min(n_modules, ncol(w_clr)))
        
        gbest_mat <- matrix(0, nrow = length(w_idx), ncol = n_modules)
        for (k in 1:n_modules) {
          k_taxa <- names(clusters[clusters == k])
          gbest_mat[, k] <- if(length(k_taxa) == 1) w_clr[, k_taxa] else rowMeans(w_clr[, k_taxa, drop=FALSE])
        }
        
        for (j in 1:ncol(w_clr)) {
          x_ij <- w_clr[, j]
          v_ij <- w_vel[, j]
          a_ij <- w_acc[, j]
          
          pbest_ij <- max(x_ij)
          gbest_ij <- gbest_mat[, clusters[j]]
          
          # Dynamical equation incorporating candidate inertia weight (omega)
          Y  <- a_ij + (1 - w_inertia) * v_ij
          X1 <- pbest_ij - x_ij
          X2 <- gbest_ij - x_ij
          
          X_mat <- cbind(X1, X2)
          Dmat  <- t(X_mat) %*% X_mat + 1e-5 * diag(2)
          dvec  <- t(X_mat) %*% Y
          Amat  <- diag(2)
          
          qp_fit <- tryCatch({ solve.QP(Dmat, dvec, Amat, bvec = c(0, 0)) }, error = function(e) NULL)
          c_vals <- if (!is.null(qp_fit)) qp_fit$solution else c(0, 0)
          
          # Integrated eco-thermodynamic engine
          taxon_ke <- 0.5 * mean(v_ij^2)
          taxon_u  <- -0.5 * c_vals[2] * mean(X2^2)
          
          results[[length(results) + 1]] <- data.frame(
            Window_Center   = w_center,
            Taxon           = colnames(clr_matrix)[j],
            Module          = clusters[j],
            c1_Cognitive    = c_vals[1],
            c2_Social       = c_vals[2],
            Dominance_Ratio = log((c_vals[2] + 1e-4) / (c_vals[1] + 1e-4)),
            Kinetic_Energy  = taxon_ke,
            Potential_Depth = taxon_u,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(results) > 0) {
        temp_pso_df  <- bind_rows(results)
        temp_meta_df <- bind_rows(window_meta_list)
        
        valid_wins <- length(unique(temp_pso_df$Window_Center))
        
        # Calculate stability metric across windows
        win_summary <- temp_pso_df %>%
          group_by(Window_Center) %>%
          summarise(
            se_c1 = sd(c1_Cognitive, na.rm=TRUE) / sqrt(n()),
            se_c2 = sd(c2_Social, na.rm=TRUE) / sqrt(n()),
            .groups="drop"
          ) %>%
          summarise(
            mean_se_c1 = mean(se_c1, na.rm=TRUE),
            mean_se_c2 = mean(se_c2, na.rm=TRUE)
          )
        
        # Multi-objective score: maximises temporal resolution while minimising force estimation variance
        score <- (valid_wins / w_sz) / (win_summary$mean_se_c1 + win_summary$mean_se_c2 + 1e-5)
        
        optimisation_scores <- bind_rows(optimisation_scores, data.frame(
          Window_Size    = w_sz,
          Inertia_Weight = w_inertia,
          Valid_Windows  = valid_wins,
          Mean_SE        = win_summary$mean_se_c1 + win_summary$mean_se_c2,
          Score          = score
        ))
        
        key_str <- sprintf("w_%.1f_i_%.2f", w_sz, w_inertia)
        candidate_results[[key_str]] <- temp_pso_df
        candidate_meta[[key_str]]    <- temp_meta_df
      }
    }
  }
  
  if (nrow(optimisation_scores) == 0) stop("No valid windows found for candidate parameters. Check temporal resolution.")
  
  # Extract optimal window & inertia pair
  opt_row <- optimisation_scores %>% arrange(desc(Score)) %>% slice(1)
  opt_w_sz      <- opt_row$Window_Size
  opt_inertia_w <- opt_row$Inertia_Weight
  
  cat(sprintf("[INFO] Optimal Parameter Selection Identified:\n"))
  cat(sprintf("       -> Window Size (w): %.1f days\n", opt_w_sz))
  cat(sprintf("       -> Swarm Inertia Weight (\u03c9): %.2f\n", opt_inertia_w))
  cat(sprintf("       -> Optimisation Score: %.3f\n", opt_row$Score))
  
  opt_key <- sprintf("w_%.1f_i_%.2f", opt_w_sz, opt_inertia_w)
  pso_df <- candidate_results[[opt_key]]
  df_window_meta <- candidate_meta[[opt_key]]
  
  write.csv(optimisation_scores, file.path(output_dir, "Window_Inertia_Optimisation_Scores.csv"), row.names = FALSE)
  
  ###############################################################################
  # Top driver identification & taxon ranking
  ###############################################################################
  cat("[INFO] Identifying and ranking top social (c2) and niche (c1) drivers for optimal configuration...\n")
  
  taxon_summary <- pso_df %>%
    group_by(Taxon) %>%
    summarise(
      Mean_c1  = mean(c1_Cognitive, na.rm = TRUE),
      SD_c1    = sd(c1_Cognitive, na.rm = TRUE),
      Mean_c2  = mean(c2_Social, na.rm = TRUE),
      SD_c2    = sd(c2_Social, na.rm = TRUE),
      Mean_Dom = mean(Dominance_Ratio, na.rm = TRUE),
      Mean_KE  = mean(Kinetic_Energy, na.rm = TRUE),
      Mean_U   = mean(Potential_Depth, na.rm = TRUE),
      .groups  = "drop"
    )
  
  top_niche_drivers  <- taxon_summary %>% arrange(desc(Mean_c1)) %>% slice_head(n = 20)
  top_social_drivers <- taxon_summary %>% arrange(desc(Mean_c2)) %>% slice_head(n = 20)
  
  write.csv(top_niche_drivers, file.path(output_dir, "Top_Niche_Drivers_c1.csv"), row.names = FALSE)
  write.csv(top_social_drivers, file.path(output_dir, "Top_Social_Drivers_c2.csv"), row.names = FALSE)
  write.csv(taxon_summary, file.path(output_dir, "All_Taxa_Driver_Summary.csv"), row.names = FALSE)
  
  ###############################################################################
  # Export feature to guild (module) mapping with taxonomic paths
  ###############################################################################
  feature_guild_mapping <- pso_df %>%
    group_by(Taxon) %>%
    summarise(
      # Identify the primary (most frequent) module assigned across sliding windows
      Guild_Module = as.integer(names(which.max(table(Module)))),
      Mean_c1_Cognitive = mean(c1_Cognitive, na.rm = TRUE),
      Mean_c2_Social    = mean(c2_Social, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Feature_ID     = ifelse(grepl("\\|", Taxon), sub("\\|.*$", "", Taxon), Taxon),
      Taxonomic_Path = ifelse(grepl("\\|", Taxon), sub("^[^|]*\\|", "", Taxon), "Unassigned")
    ) %>%
    select(Feature_ID, Taxonomic_Path, Guild_Module, Mean_c1_Cognitive, Mean_c2_Social, Full_Feature_ID = Taxon) %>%
    arrange(Guild_Module, Feature_ID)
  
  write.csv(feature_guild_mapping, file.path(output_dir, "Feature_Guild_Assignments.csv"), row.names = FALSE)
  
  ###############################################################################
  # Ecological & thermodynamic visualisation
  ###############################################################################
  cat("[INFO] Generating publication-grade figure layouts...\n")
  
  theme_set(theme_pubr(base_size = 10) + theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ))
  
  ##############################################################################
  # Figure 1: 2D joint optimisation surface & inertia sensitivity
  ##############################################################################
  p_1a <- ggplot(optimisation_scores, aes(x = Window_Size, y = Inertia_Weight, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis_c(option = "viridis", name = "Optimisation\nScore") +
    geom_vline(xintercept = opt_w_sz, linetype = "dashed", color = "#D95F02", linewidth = 0.8) +
    geom_hline(yintercept = opt_inertia_w, linetype = "dashed", color = "#D95F02", linewidth = 0.8) +
    geom_point(data = opt_row, aes(x = Window_Size, y = Inertia_Weight), 
               color = "#D95F02", fill = "#FFFFBF", shape = 23, size = 5, stroke = 1.5) +
    annotate("label", x = opt_w_sz, y = opt_inertia_w, 
             label = sprintf("Optimal: w = %.1f, \u03c9 = %.2f", opt_w_sz, opt_inertia_w), 
             fill = "white", color = "#D95F02", fontface = "bold", size = 3.2, vjust = -0.45, hjust = 0.5) +
    scale_x_continuous(breaks = unique(optimisation_scores$Window_Size)) +
    scale_y_continuous(breaks = unique(optimisation_scores$Inertia_Weight)) +
    labs(x = "Candidate Window Size (Timeline Units)", y = "Swarm Inertia Weight (\u03c9)", 
         title = "(a) 2D Joint Parameter Optimisation Surface") +
    theme(panel.grid.major = element_blank())
  
  p_1b <- ggplot(optimisation_scores, aes(x = Inertia_Weight, y = Score, 
                                          color = as.factor(Window_Size), 
                                          group = as.factor(Window_Size))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_viridis_d(option = "plasma", name = "Window Size") +
    scale_x_continuous(breaks = unique(optimisation_scores$Inertia_Weight)) +
    labs(x = "Swarm Inertia Weight (\u03c9)", y = "Optimisation Score", 
         title = "(b) Inertia Sensitivity Across Window Candidates")
  
  fig1 <- p_1a / p_1b + plot_layout(heights = c(1.2, 1))
  fig1_path <- file.path(output_dir, "Figure1_Window_Inertia_Optimisation.pdf")
  ggsave(fig1_path, fig1, width = 8, height = 9, dpi = 300, device = cairo_pdf)
  
  ##############################################################################
  # Process summary for optimal window & inertia Pair
  ##############################################################################
  summary_pso <- pso_df %>%
    group_by(Window_Center) %>%
    summarise(
      Mean_c1  = mean(c1_Cognitive),
      Mean_c2  = mean(c2_Social),
      SE_c1    = sd(c1_Cognitive) / sqrt(n()),
      SE_c2    = sd(c2_Social) / sqrt(n()),
      Mean_Dom = mean(Dominance_Ratio),
      Mean_KE  = mean(Kinetic_Energy),
      Mean_U   = mean(Potential_Depth),
      .groups  = "drop"
    ) %>%
    mutate(Escape_Index = Mean_KE / (abs(Mean_U) + 1e-5))
  
  merged_perf <- left_join(summary_pso, df_window_meta, by = "Window_Center")
  time_breaks <- unique(summary_pso$Window_Center)
  
  ##############################################################################
  # Figure 2: Temporal dynamics & phase space
  ##############################################################################
  p_2a <- ggplot(summary_pso, aes(x = Window_Center)) +
    geom_line(aes(y = Mean_c1, color = "Niche Memory (c1)"), linewidth = 1.1) +
    geom_ribbon(aes(ymin = Mean_c1 - SE_c1, ymax = Mean_c1 + SE_c1, fill = "Niche Memory (c1)"), alpha = 0.2) +
    geom_line(aes(y = Mean_c2, color = "Environmental Filtering (c2)"), linewidth = 1.1) +
    geom_ribbon(aes(ymin = Mean_c2 - SE_c2, ymax = Mean_c2 + SE_c2, fill = "Environmental Filtering (c2)"), alpha = 0.2) +
    scale_x_continuous(breaks = time_breaks) +
    scale_color_manual(values = c("Niche Memory (c1)" = "#D95F02", "Environmental Filtering (c2)" = "#1B9E77")) +
    scale_fill_manual(values = c("Niche Memory (c1)" = "#D95F02", "Environmental Filtering (c2)" = "#1B9E77")) +
    labs(x = paste(time_col, "(Timeline Center)"), y = "Acceleration Coefficient", 
         title = "(a) Temporal Dynamics of Community Drivers", color = NULL, fill = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_2b <- ggplot(summary_pso, aes(x = Mean_c1, y = Mean_c2, color = Window_Center)) +
    geom_path(linewidth = 1.0, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_point(size = 3.0) +
    scale_color_viridis_c(option = "plasma", name = time_col) +
    labs(x = "Cognitive Driver (c1: Niche Memory)", y = "Social Driver (c2: Filtering)", 
         title = "(b) Swarm Phase Space Trajectory")
  
  fig2 <- p_2a / p_2b + plot_layout(heights = c(1, 1))
  fig2_path <- file.path(output_dir, "Figure2_Temporal_Dynamics_Phase_Space.pdf")
  ggsave(fig2_path, fig2, width = 8, height = 10, dpi = 300, device = cairo_pdf)
  
  ##############################################################################
  # Figure 3: Thermodynamics & escape velocity
  ##############################################################################
  p_3a <- ggplot(pso_df, aes(x = c1_Cognitive, y = c2_Social, z = Potential_Depth)) +
    stat_summary_2d(fun = mean, bins = 25, alpha = 0.9) +
    scale_fill_viridis_c(option = "mako", name = "Potential\nDepth (U)", direction = -1) +
    labs(x = "Cognitive Vector (c1)", y = "Social Vector (c2)", 
         title = "(a) Thermodynamic Energy Landscape") +
    coord_cartesian(
      xlim = c(0, quantile(pso_df$c1_Cognitive, 0.95, na.rm = TRUE)), 
      ylim = c(0, quantile(pso_df$c2_Social, 0.95, na.rm = TRUE))
    )
  
  p_3b <- ggplot(summary_pso, aes(x = Window_Center, y = Escape_Index)) +
    geom_area(fill = "#9E0142", alpha = 0.35) +
    geom_line(color = "#9E0142", linewidth = 1.2) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.9) +
    scale_x_continuous(breaks = time_breaks) +
    annotate("text", x = min(summary_pso$Window_Center), y = 1.5, 
             label = "Escape Velocity Threshold (KE > |U|)", hjust = 0, fontface = "bold", size = 3) +
    labs(x = paste(time_col, "(Timeline Center)"), y = "Dysbiosis Risk (KE / |U|)", 
         title = "(b) Swarm Escape Velocity & Dysbiosis Forecasting") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  fig3 <- p_3a / p_3b + plot_layout(heights = c(1, 1))
  fig3_path <- file.path(output_dir, "Figure3_Thermodynamics_Escape_Velocity.pdf")
  ggsave(fig3_path, fig3, width = 8, height = 10, dpi = 300, device = cairo_pdf)
  
  ##############################################################################
  # Figure 4: Metadata attractor regressions and boxplots
  ##############################################################################
  param_plots <- list()
  plot_idx <- 1
  for (v in target_vars) {
    if (!v %in% colnames(merged_perf)) next
    
    panel_letter <- letters[plot_idx]
    
    if (is.numeric(merged_perf[[v]])) {
      pv <- ggplot(merged_perf, aes(x = Mean_c2, y = .data[[v]])) +
        geom_point(aes(color = Window_Center), size = 3.0) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        # Use breaks_pretty to prevent text crowding and format the colorbar guide
        scale_color_viridis_c(
          option = "magma", 
          name = time_col, 
          breaks = scales::breaks_pretty(n = 6),
          guide = guide_colorbar(barheight = unit(3.5, "cm"), title.position = "top")
        ) +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
        labs(x = "Social Driver (c2)", y = v, title = paste("(",panel_letter, ") Swarm Attractor vs", v))
    } else {
      pv <- ggplot(merged_perf, aes(x = as.factor(.data[[v]]), y = Mean_c2, fill = as.factor(.data[[v]]))) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
        stat_compare_means(method = "kruskal.test", label.y.npc = "bottom", size = 3.5, label = "p.format") +
        scale_fill_viridis_d(option = "turbo", guide = "none") +
        labs(x = v, y = "Social Driver (c2)", title = paste("(",panel_letter, ") Social Attractor across", v))
    }
    param_plots[[v]] <- pv
    plot_idx <- plot_idx + 1
  }
  
  if (length(param_plots) > 0) {
    # Collect guides into a single non-overlapping legend sidebar and slightly widen the canvas
    fig4 <- wrap_plots(param_plots, ncol = min(2, length(param_plots))) + 
      plot_layout(guides = "collect")
    
    fig4_path <- file.path(output_dir, "Figure4_Metadata_Attractor_Regressions.pdf")
    ggsave(fig4_path, fig4, width = 11, height = max(4, ceiling(length(param_plots)/2) * 4), dpi = 300, device = cairo_pdf)
  }
  ##############################################################################
  # Extract last resolved taxonomic level for visualisation
  ##############################################################################
  get_last_resolved_taxon <- function(tax_string) {
    parts <- strsplit(tax_string, "\\|")[[1]]
    id_part <- ifelse(length(parts) > 1, parts[1], "")
    path_part <- ifelse(length(parts) > 1, paste(parts[-1], collapse = "|"), parts[1])
    
    ranks <- strsplit(path_part, ";")[[1]]
    ranks <- trimws(ranks)
    
    valid_ranks <- ranks[!grepl("^[a-zA-Z]__\\(|^uncultured.*|^unidentified.*|^\\)|Ambiguous", ranks, ignore.case = TRUE)]
    
    if (length(valid_ranks) > 0) {
      last_res <- tail(valid_ranks, 1)
      last_res <- sub("^[a-zA-Z]__", "", last_res)
    } else {
      last_res <- "Unknown"
    }
    
    if (id_part != "") {
      return(paste(id_part, last_res, sep = " | "))
    } else {
      return(last_res)
    }
  }
  
  top_combined <- bind_rows(
    top_niche_drivers %>% mutate(Driver_Type = "Top Niche Memory (c1)", Value = Mean_c1),
    top_social_drivers %>% mutate(Driver_Type = "Top Social Filtering (c2)", Value = Mean_c2)
  ) %>%
    mutate(Taxon_Clean = sapply(as.character(Taxon), get_last_resolved_taxon))
  
  ##############################################################################
  # Figure 5: Top keystone drivers
  ##############################################################################
  
  # Create a unique Taxon-Facet key to enforce independent sorting per facet
  top_combined_sorted <- top_combined %>%
    mutate(
      Taxon_Facet = reorder(paste(Taxon_Clean, Driver_Type, sep = "___"), Value)
    )
  
  p_5a <- ggplot(top_combined_sorted, aes(x = Taxon_Facet, y = Value, color = Driver_Type)) +
    geom_segment(aes(xend = Taxon_Facet, y = 0, yend = Value), linewidth = 1) +
    geom_point(size = 3) +
    coord_flip() +
    facet_wrap(~Driver_Type, scales = "free") +
    # Strip the delimiter suffix so axis labels remain clean
    scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_color_manual(values = c("Top Niche Memory (c1)" = "#D95F02", "Top Social Filtering (c2)" = "#1B9E77")) +
    labs(
      x = "Taxon (Lowest Resolved)", 
      y = "Mean Acceleration Coefficient", 
      title = "Top Keystone Autecological & Synecological Drivers"
    ) +
    theme(legend.position = "none", axis.text.y = element_text(size = 8))
  
  fig5_path <- file.path(output_dir, "Figure5_Top_Drivers_Lollipop.pdf")
  ggsave(fig5_path, p_5a, width = 10, height = 6, dpi = 300, device = cairo_pdf) 
  
  ##############################################################################
  # Figure 6: Ecological guild (module) dynamics & thermodynamic niches
  ##############################################################################
  cat("[INFO] Generating Figure 6: Module-level dynamics and thermodynamics...\n")
  
  # Prepare & Sanitise Module-Level Data
  eps <- 1e-6  # Tiny offset to prevent log10(0) = -Inf
  
  pso_df_mod <- pso_df %>% 
    mutate(
      Module = as.factor(Module),
      abs_U  = abs(Potential_Depth),
      # Pre-compute safe log transformations directly to keep all observations intact
      log_U  = log10(abs_U + eps),
      log_KE = log10(Kinetic_Energy + eps)
    ) %>%
    # Filter out any non-finite entries upfront
    filter(
      is.finite(c1_Cognitive), 
      is.finite(c2_Social), 
      is.finite(log_U), 
      is.finite(log_KE)
    )
  
  # Calculate time-series statistics per Module
  module_time_stats <- pso_df_mod %>%
    group_by(Window_Center, Module) %>%
    summarise(
      Mean_c1 = mean(c1_Cognitive, na.rm = TRUE),
      SE_c1   = sd(c1_Cognitive, na.rm = TRUE) / sqrt(n()),
      Mean_c2 = mean(c2_Social, na.rm = TRUE),
      SE_c2   = sd(c2_Social, na.rm = TRUE) / sqrt(n()),
      Mean_KE = mean(Kinetic_Energy, na.rm = TRUE),
      Mean_U  = mean(Potential_Depth, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Extract Keystone Taxa per Module
  keystone_taxa <- pso_df_mod %>%
    group_by(Module, Taxon) %>%
    summarise(Total_Force = mean(c1_Cognitive + c2_Social, na.rm = TRUE), .groups = "drop") %>%
    group_by(Module) %>%
    slice_max(order_by = Total_Force, n = 1) %>%
    mutate(Keystone_Taxon_Clean = sapply(as.character(Taxon), get_last_resolved_taxon)) %>%
    select(Module, Keystone_Taxon_Clean)
  
  # 3. Calculate Global Module Statistics & Export CSV
  module_overall_stats <- pso_df_mod %>%
    group_by(Module) %>%
    summarise(
      Total_Observations  = n(),
      Unique_Taxa_Count   = n_distinct(Taxon),
      Avg_c1_Cognitive    = mean(c1_Cognitive, na.rm = TRUE),
      Avg_c2_Social       = mean(c2_Social, na.rm = TRUE),
      Avg_Kinetic_Energy  = mean(Kinetic_Energy, na.rm = TRUE),
      Avg_Potential_Depth = mean(Potential_Depth, na.rm = TRUE),
      Stability_Index     = mean(abs_U / (Kinetic_Energy + eps), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(keystone_taxa, by = "Module") %>%
    arrange(desc(Stability_Index))
  
  mod_stats_file <- file.path(output_dir, "Figure6_Module_Thermodynamic_Stats.csv")
  write.csv(module_overall_stats, mod_stats_file, row.names = FALSE)
  cat(sprintf("[INFO] Module statistics successfully saved to: %s\n", mod_stats_file))
  
  # 4. Generate Multi-Panel Publication Figure
  
  # Panel A: Temporal Evolution of Niche Strategy
  p_6a <- ggplot(module_time_stats, aes(x = Window_Center, y = Mean_c1, color = Module, fill = Module)) +
    geom_line(linewidth = 1.1) +
    geom_ribbon(aes(ymin = Mean_c1 - SE_c1, ymax = Mean_c1 + SE_c1), alpha = 0.18, color = NA) +
    scale_color_viridis_d(option = "turbo", name = "Ecological Guild\n(Module)") +
    scale_fill_viridis_d(option = "turbo", name = "Ecological Guild\n(Module)") +
    labs(
      x = paste(time_col, "(Timeline Center)"),
      y = "Mean Niche Memory (c1 ± SE)",
      title = "(a) Temporal Niche Strategy by Ecological Guild"
    ) +
    theme(legend.position = "right")
  
  # Create dedicated subsets for active dynamics
  pso_thermo_active <- pso_df_mod %>% 
    filter(abs_U > 1e-5) # Strips the log10(1e-6) = -6 zero artifact
  
  pso_forces_active <- pso_df_mod %>% 
    filter(c1_Cognitive > 1e-4 | c2_Social > 1e-4) # Strips the (0,0) origin cluster
  
  # Panel B: Cleaned Thermodynamic Plot
  p_6b <- ggplot(pso_thermo_active, aes(x = log_U, y = log_KE, color = Module, fill = Module)) +
    geom_point(alpha = 0.35, size = 1.1, stroke = 0) +
    stat_ellipse(geom = "polygon", alpha = 0.18, type = "norm", level = 0.85, na.rm = TRUE) +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    labs(
      x = expression("Potential Depth " * log[10] * "(|U| + " * 10^-6 * ") (Stability)"),
      y = expression("Kinetic Energy " * log[10] * "(KE + " * 10^-6 * ") (Volatility)"),
      title = "(b) Thermodynamic Niche Partitioning (Active Core)"
    )
  
  # Panel C: Cleaned Evolutionary Force Trade-offs
  p_6c <- ggplot(pso_forces_active, aes(x = c1_Cognitive, y = c2_Social, color = Module, fill = Module)) +
    geom_point(alpha = 0.35, size = 1.1, stroke = 0) +
    stat_ellipse(geom = "polygon", alpha = 0.18, type = "norm", level = 0.85, na.rm = TRUE) +
    scale_color_viridis_d(option = "turbo", name = "Ecological Guild\n(Module)") +
    scale_fill_viridis_d(option = "turbo", name = "Ecological Guild\n(Module)") +
    labs(
      x = "Cognitive Vector (c1, Autecology)",
      y = "Social Vector (c2, Synecology)",
      title = "(c) Evolutionary Force Trade-offs (Active Agents)"
    )
  
  # Compile Figure 6 using patchwork
  fig6 <- (p_6a) / (p_6b | p_6c) + plot_layout(heights = c(1, 1.2), guides = "collect")
  fig6_path <- file.path(output_dir, "Figure6_Module_Eco_Thermodynamics.pdf")
  ggsave(fig6_path, fig6, width = 11, height = 9, dpi = 300, device = cairo_pdf)
  
  cat(sprintf("[INFO] Multiple camera-ready figures successfully saved to: %s\n", output_dir))
  
  return(list(
    optimal_window_size = opt_w_sz,
    optimal_inertia_w   = opt_inertia_w,
    optimisation_scores = optimisation_scores,
    pso_data            = pso_df,
    taxon_summary       = taxon_summary,
    feature_guild_mapping = feature_guild_mapping,   
    top_niche_drivers   = top_niche_drivers,
    top_social_drivers  = top_social_drivers,
    metadata_aggregated = df_window_meta,
    summary_thermo      = summary_pso,
    plot_fig1           = fig1,
    plot_fig2           = fig2,
    plot_fig3           = fig3,
    plot_fig5           = p_5a,
    plot_fig6           = fig6,                 
    module_stats        = module_overall_stats  
  ))
}

##############################################################################
# Execution example
##############################################################################
# result <- run_advanced_ipsda(
#   count_file         = "Ijaz2018/feature_w_tax.biom",
#   meta_file          = "Ijaz2018/meta_table.csv",
#   output_dir         = "./ipsda_unified_results",
#   sample_id_col      = "SampleID",
#   time_col           = "Day", # Continuous time variable
#   subject_col        = "Pen_No",  # Controls for repeated measures
#   target_vars        = c("BW_Mean", "Gain", "FI", "FCR", "Diet", "Sampling_Groups"), # Metadata mapped to attractors
#   window_candidates  = c(3.0, 4.0, 5.0, 6.0, 7.0), # Temporal resolution grid
#   window_step        = 1.0,
#   inertia_candidates = seq(0.1, 0.8, by = 0.1), # Swarm memory grid
#   n_modules          = 4, # Number of synecological guilds
#   gamm_k             = 8 # Spline knots for differentiation
# )