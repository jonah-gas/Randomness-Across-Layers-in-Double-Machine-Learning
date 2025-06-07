# check, install, and load necessary packages
packages <- c("data.table", "kableExtra")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

########## Data Preparation an Additional Analysis ##########


# load function
source("functions/get_estimates.R")

# specify input file and load results
dgp_choice <- "DGP_1"
# dgp_choice <- "DGP_2"
file_name <- paste0("results/sim2_results_", dgp_choice, ".rds")
results <- readRDS(file_name)

# define trimming
truncation_levels <- c(0, 0.01, 0.05, 0.1, 0.15, 0.2)

all_estimates <- list()

# loop over learners and compute aggregated estimates for different truncation levels
for (row_idx in 1:nrow(results)) {
  learner_name <- results$learner[row_idx]
  cat("\n==== Learner:", learner_name, " ====\n")
  
  # get nuisance estimates
  exhat_mat <- results$exhat[[row_idx]]
  mxhat_mat <- results$mxhat[[row_idx]]
  mwhat0_mat <- results$mwhat0[[row_idx]]
  mwhat1_mat <- results$mwhat1[[row_idx]]
  y <- results$y[[row_idx]]
  w <- results$w[[row_idx]]

  n_obs <- length(y)
  n_cv_rep <- ncol(exhat_mat)
  
  # loop over truncation levels
  for (trunc in truncation_levels) {
    psi_a_plr_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
    psi_b_plr_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
    psi_a_irm_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
    psi_b_irm_array <- array(NA, dim = c(n_obs, n_cv_rep, 1))
    
    for (rep_idx in 1:n_cv_rep) {
      exhat <- exhat_mat[, rep_idx]
      mxhat <- mxhat_mat[, rep_idx]
      mwhat0 <- mwhat0_mat[, rep_idx]
      mwhat1 <- mwhat1_mat[, rep_idx]
      
      keep <- (exhat > trunc) & (exhat < 1 - trunc)

      
      y_trim <- y[keep]
      w_trim <- w[keep]
      exhat_trim <- exhat[keep]
      mxhat_trim <- mxhat[keep]
      mwhat0_trim <- mwhat0[keep]
      mwhat1_trim <- mwhat1[keep]
      
      # compute psi
      psi_a_plr <- -(w_trim - exhat_trim)^2
      psi_b_plr <- (y_trim - mxhat_trim) * (w_trim - exhat_trim)
      psi_a_irm <- rep(-1, length(y_trim))
      psi_b_irm <- mwhat1_trim - mwhat0_trim +
        w_trim * (y_trim - mwhat1_trim) / exhat_trim -
        (1 - w_trim) * (y_trim - mwhat0_trim) / (1 - exhat_trim)
      
      psi_a_plr_array[keep, rep_idx, 1] <- psi_a_plr
      psi_b_plr_array[keep, rep_idx, 1] <- psi_b_plr
      psi_a_irm_array[keep, rep_idx, 1] <- psi_a_irm
      psi_b_irm_array[keep, rep_idx, 1] <- psi_b_irm
    }
    
    # compute estimates
    est_plr <- get_estimates(psi_a_plr_array, psi_b_plr_array)
    est_irm <- get_estimates(psi_a_irm_array, psi_b_irm_array)
    
    # compute average number of obs left
    mean_kept <- mean(colSums((exhat_mat > trunc) & (exhat_mat < 1 - trunc)))
    
    # save
    all_estimates[[length(all_estimates) + 1]] <- data.table(
      learner = learner_name,
      trunc_level = trunc,
      mean_obs_kept = mean_kept,
      mean_theta_plr = est_plr$mean_theta,
      se_IF_plr = est_plr$se_IF,
      mean_theta_irm = est_irm$mean_theta,
      se_IF_irm = est_irm$se_IF
    )
  }
}

# combine results
final_estimates_trunc_single_sim <- rbindlist(all_estimates)



########## Tables ##########

# nice names
nice_learner <- function(x) {
  switch(x,
         "random_forest" = "random forest",
         x)
}

# table function for AIPW results only
save_trunc_table_kable <- function(learner_name, dgp_choice, label_tag) {
  
  tab <- final_estimates_trunc_single_sim[learner == learner_name, .(
    `Truncation level` = trunc_level,
    `Mean observations kept` = mean_obs_kept,
    `\\(\\tilde{\\hat{\\theta}}^{\\text{mean}}\\)` = mean_theta_irm,
    `\\(\\text{se}_{\\text{IF}}\\)` = se_IF_irm
  )]
  
  file_out <- paste0("table_trunc_", label_tag, "_", dgp_choice, ".tex")
  Dgp_name <- ifelse(dgp_choice == "DGP_1", "DGP1", "DGP2")
  
  caption_title <- paste0(
    "Aggregated AIPW estimates for ", Dgp_name, 
    " under different propensity score truncation levels for ", 
    nice_learner(learner_name), "."
  )
  
  tab %>%
    kable("latex",
          booktabs = TRUE,
          linesep = "", 
          digits = 3,
          align = "c",
          caption = caption_title,
          label = paste0("tab:trunc_", label_tag, "_", dgp_choice),
          escape = FALSE,
    ) %>%
    kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 9) %>%
    writeLines(file_out)
}


# apply function
save_trunc_table_kable("xgboost",       dgp_choice, "xgboost")
save_trunc_table_kable("random_forest", dgp_choice, "rf")
save_trunc_table_kable("lasso",         dgp_choice, "lasso")
