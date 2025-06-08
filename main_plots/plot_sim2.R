# check, install, and load necessary packages
packages <- c("data.table", "tidyverse", "patchwork")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


########## Data Preparation ##########

# load function
source("functions/get_estimates.R")

# specify input file and load results
# selected_dgp <- "DGP_1"
selected_dgp <- "DGP_2"
input_file <- paste0("results/sim2_results_", selected_dgp, ".rds")
results_sim2 <- readRDS(input_file)

# define trimming
trunc <- 0.01

# function to compute aggregated estimates for subset of cross-fitting repetitions
get_agg_estimates <- function(psi_a_array, psi_b_array) {
  results <- list()
  num_reps <- dim(psi_a_array)[2]
  
  for (i in 1:num_reps) {
    if (i %% 50 == 0) cat("Processing cumulative repetition", i, "of", num_reps, "\n")
    
    # subset psi arrays up to the current repetition
    subset_psi_a <- psi_a_array[, 1:i, 1, drop = FALSE]
    subset_psi_b <- psi_b_array[, 1:i, 1, drop = FALSE]
    
    # compute point and standard error estimates from current subset
    estimates <- get_estimates(subset_psi_a, subset_psi_b)
    
    # store estimates together with current cumulative repetition index
    results[[i]] <- cbind(data.frame(cum_reps = i), as.data.frame(estimates))
  }
  
  # combine all cumulative results into a single dataframe
  return(do.call(rbind, results))
}

# compute psi arrays for each learner
psi_arrays_all_learners <- list()

for (learner_name in unique(results_sim2$learner)) {
  cat("\nProcessing learner:", learner_name, "\n")
  row <- results_sim2[learner == learner_name]
  
  exhat_mat <- row$exhat[[1]]
  mxhat_mat <- row$mxhat[[1]]
  mwhat0_mat <- row$mwhat0[[1]]
  mwhat1_mat <- row$mwhat1[[1]]
  y <- row$y[[1]]
  w <- row$w[[1]]
  
  n_obs <- length(y)
  n_cv_rep <- ncol(exhat_mat)
  max_reps <- n_cv_rep
  
  
  # initialize empty psi arrays
  psi_a_plr_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
  psi_b_plr_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
  psi_a_irm_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
  psi_b_irm_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
  
  for (rep_idx in 1:n_cv_rep) {
    exhat <- exhat_mat[, rep_idx]
    mxhat <- mxhat_mat[, rep_idx]
    mwhat0 <- mwhat0_mat[, rep_idx]
    mwhat1 <- mwhat1_mat[, rep_idx]
    
    # apply trimming
    keep <- (exhat > trunc) & (exhat < 1 - trunc)
    if (sum(keep) < 20) {
      cat("Trimming threshold too strict at repetition", rep_idx, "(", sum(keep), "remaining observations); skipping.\n")
      next
    }
    y_trim <- y[keep]
    w_trim <- w[keep]
    exhat_trim <- exhat[keep]
    mxhat_trim <- mxhat[keep]
    mwhat0_trim <- mwhat0[keep]
    mwhat1_trim <- mwhat1[keep]
    
    # calculate psi values for PLR
    psi_a_plr_array[keep, rep_idx, 1] <- -(w_trim - exhat_trim)^2
    psi_b_plr_array[keep, rep_idx, 1] <- (y_trim - mxhat_trim) * (w_trim - exhat_trim)
    
    # calculate psi values for AIPW
    psi_a_irm_array[keep, rep_idx, 1] <- -1
    psi_b_irm_array[keep, rep_idx, 1] <- mwhat1_trim - mwhat0_trim +
      w_trim * (y_trim - mwhat1_trim) / exhat_trim -
      (1 - w_trim) * (y_trim - mwhat0_trim) / (1 - exhat_trim)
  }
  
  psi_a_plr_array <- psi_a_plr_array[, 1:max_reps, , drop = FALSE]
  psi_b_plr_array <- psi_b_plr_array[, 1:max_reps, , drop = FALSE]
  psi_a_irm_array <- psi_a_irm_array[, 1:max_reps, , drop = FALSE]
  psi_b_irm_array <- psi_b_irm_array[, 1:max_reps, , drop = FALSE]
  
  psi_arrays_all_learners[[learner_name]] <- list(
    psi_a_plr_array = psi_a_plr_array,
    psi_b_plr_array = psi_b_plr_array,
    psi_a_irm_array = psi_a_irm_array,
    psi_b_irm_array = psi_b_irm_array
  )
  
}


# use psi arrays to compute aggregated estimates for each learner
all_estimates <- list()

for (learner_name in names(psi_arrays_all_learners)) {
  cat("\nComputing cumulative estimates for learner:", learner_name, "\n")
  
  psi <- psi_arrays_all_learners[[learner_name]]
  
  results_plr_df <- get_agg_estimates(psi$psi_a_plr_array, psi$psi_b_plr_array)
  results_irm_df <- get_agg_estimates(psi$psi_a_irm_array, psi$psi_b_irm_array)
  
  # annotate with learner and estimator identifiers
  results_plr_df$learner <- learner_name
  results_irm_df$learner <- learner_name
  results_plr_df$estimator <- "PLR"
  results_irm_df$estimator <- "IRM"
  
  # combine both estimators for this learner
  all_estimates[[learner_name]] <- rbind(results_plr_df, results_irm_df)
}

# combine results across all learners
final_estimates <- rbindlist(all_estimates)






# generate dynamic names depending on DGP
final_estimates_name <- paste0("final_estimates_", selected_dgp)
assign(final_estimates_name, final_estimates)

########## Plot ##########
# after getting final estimates for BOTH DGPs continue with combined Plotting 


plot_cumreps_by_dgp <- function(final_estimates, dgp_label, line_size = 1.1) {
  # compute CIs
  df_plot <- final_estimates %>%
    mutate(
      ci_low_mean = mean_theta - 1.96 * se_mean,
      ci_high_mean = mean_theta + 1.96 * se_mean,
      ci_low_median = median_theta - 1.96 * se_median,
      ci_high_median = median_theta + 1.96 * se_median,
      mean_theta_lower_seif = mean_theta - 1.96 * se_IF,
      mean_theta_upper_seif = mean_theta + 1.96 * se_IF,
    ) %>%
    mutate(
      learner = recode(learner,
                       "lasso" = "Lasso",
                       "random_forest" = "Random Forest",
                       "xgboost" = "Xgboost"
      ),
      estimator = recode(estimator,
                         "IRM" = "AIPW"
      ),
      learner = factor(learner, levels = c("Lasso", "Random Forest", "Xgboost")),
      estimator = factor(estimator, levels = c("PLR", "AIPW"))
    )
  
  # adapt limits
  if (dgp_label == "DGP_1" || grepl("1", dgp_label)) {
    xlim <- c(0, 150)
    ylim <- c(0.8, 0.2)      
    title_text <- "DGP 1"
  } else if (dgp_label == "DGP_2" || grepl("2", dgp_label)) {
    xlim <- c(0, 150)
    ylim <- c(0.435, 1.15)   
    title_text <- "DGP 2"
  } else {
    xlim <- c(0, 150)
    ylim <- range(df_plot$mean_theta, na.rm = TRUE)
    title_text <- dgp_label
  }
  
  ggp <- ggplot(df_plot, aes(x = cum_reps)) +
    geom_ribbon(
      aes(ymin = ci_low_mean, ymax = ci_high_mean, fill = "Mean-based"), 
      alpha = 0.18, inherit.aes = TRUE
    ) +
    geom_ribbon(
      aes(ymin = ci_low_median, ymax = ci_high_median, fill = "Median-based"), 
      alpha = 0.18, inherit.aes = TRUE
    ) +
    geom_ribbon(
      aes(ymin = mean_theta_lower_seif, ymax = mean_theta_upper_seif, fill = "IF-based"),
      alpha = 0.18, inherit.aes = TRUE
    ) +
    geom_line(aes(y = mean_theta, color = "Mean-based"), size = line_size) +
    geom_line(aes(y = ci_low_mean, color = "Mean-based"), linetype = "dotted", size = line_size) +
    geom_line(aes(y = ci_high_mean, color = "Mean-based"), linetype = "dotted", size = line_size) +
    geom_line(aes(y = median_theta, color = "Median-based"), size = line_size) +
    geom_line(aes(y = ci_low_median, color = "Median-based"), linetype = "dotted", size = line_size) +
    geom_line(aes(y = ci_high_median, color = "Median-based"), linetype = "dotted", size = line_size) +
    geom_line(aes(y = mean_theta_lower_seif, color = "IF-based"), linetype = "dotted", size = line_size) +
    geom_line(aes(y = mean_theta_upper_seif, color = "IF-based"), linetype = "dotted", size = line_size) +
    facet_grid(estimator ~ learner) +
    scale_color_manual(
      values = c(
        "Mean-based" = "#440154",
        "Median-based" = "#fde725",
        "IF-based" = "#21918c"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Mean-based" = "#440154",
        "Median-based" = "#fde725",
        "IF-based" = "#21918c"
      ),
      guide = "none"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = title_text,
      x = "Cross-fitting Repetition",
      y = "Estimate",
      color = "Aggregation Method:"
    ) +
    theme_minimal(base_size = 14)+
    theme(
      plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.title.position = "plot",
      axis.title.x    = element_text(size = 12),
      axis.title.y    = element_text(size = 12),
      axis.text.x     = element_text(size = 10),
      axis.text.y     = element_text(size = 10),
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 10),
      strip.text      = element_text(size = 12),
      panel.spacing   = unit(1.25, "lines")
    )
  
  
  return(ggp)
}




# solo plots
p1 <- plot_cumreps_by_dgp(final_estimates_DGP_1, "DGP_1")
p2 <- plot_cumreps_by_dgp(final_estimates_DGP_2, "DGP_2")

# combine plots
combined_plot <- p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(combined_plot)

