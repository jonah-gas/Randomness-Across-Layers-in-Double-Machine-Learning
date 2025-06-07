# check, install, and load necessary packages
packages <- c("data.table", "kableExtra")
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
selected_dgp <- "DGP1_1"
# selected_dgp <- "DGP_2"
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


########## Tables ##########

# select repetitions to be described in tavle
selected_reps <- c(1, 5, 10, 20, 50, 100, 200, 300, 400, 500)
filtered <- as.data.table(final_estimates)[cum_reps %in% selected_reps]

learners <- unique(filtered$learner)

# nice names
learner_nice_names <- list(
  lasso = "lasso",
  random_forest = "random forest",
  xgboost = "xgboost"
)
estimator_nice_names <- list(
  PLR = "PLR",
  IRM = "AIPW"
)
format_dgp <- function(x) gsub("_", "", x)

# create tables and save as tex
for (l in learners) {
  for (est in c("PLR", "IRM")) {
    sub_table <- filtered[learner == l & estimator == est, .(
      C = cum_reps,
      `\\(\\tilde{\\hat{\\theta}}^{\\text{mean}}\\)` = mean_theta,
      `\\(\\tilde{\\hat{\\theta}}^{\\text{median}}\\)` = median_theta,
      `\\(\\text{se}_{\\text{mean}}\\)` = se_mean,
      `\\(\\text{se}_{\\text{median}}\\)` = se_median,
      `\\(\\text{se}_{\\text{IF}}\\)` = se_IF
    )][order(C)]
    
    filename <- paste0("table_", format_dgp(selected_dgp), "_", l, "_", est, ".tex")
    learner_nice <- learner_nice_names[[l]]
    estimator_nice <- estimator_nice_names[[est]]
    dgp_nice <- format_dgp(selected_dgp)
    
    tab_title <- paste0(
      "Detailed results of Simulation~2 for ",
      dgp_nice, ", ",
      learner_nice, ", ",
      estimator_nice, "."
    )
    
    kable_latex <- kable(
      sub_table,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      caption = tab_title,
      digits = 3,
      label = paste0("tab:sim2_", dgp_nice, "_", l, "_", est),
      col.names = colnames(sub_table),
      align = "c"
    ) %>%
      kable_styling(
        latex_options = c("hold_position"),
        font_size = 8,         
        position = "center"
      )
    
    cat(as.character(kable_latex), file = filename)
    cat("Saved LaTeX table to:", filename, "\n")
  }
}
