# check, install, and load necessary packages
packages <- c("data.table", "kableExtra")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


########## Ritzwoller and Romano Stopping Algorithm ##########

# specify input file and load results
selected_dgp <- "DGP_2"
input_file <- paste0("resutls/sim2_results_", selected_dgp, ".rds")
results_sim2 <- readRDS(input_file)

# config
xi <- 0.001
beta <- 0.05
trunc <- 0.01
z_val <- qnorm(1 - beta / 2)
cv_xi_beta <- (xi / (2 * z_val))^2
summary_table <- data.table()

# loop over alle learners
for (learner_name in unique(results_sim2$learner)) {
  cat("\n==== Processing learner:", learner_name, "====\n")
  row <- results_sim2[learner == learner_name]
  
  cat("  -- Simulation", 1, "\n")
  
  exhat_mat <- row$exhat[[1]]
  mxhat_mat <- row$mxhat[[1]]
  mwhat0_mat <- row$mwhat0[[1]]
  mwhat1_mat <- row$mwhat1[[1]]
  y <- row$y[[1]]
  w <- row$w[[1]]
  n_obs <- length(y)
  n_reps <- ncol(exhat_mat)
  
  for (estimator in c("PLR", "IRM")) {
    cat("    >> Estimator:", estimator, "\n")
    p_values <- numeric()
    theta_vals <- numeric()
    agg_p_vals <- numeric()
    se_vals <- numeric()
    g_stop <- NA
    
    for (g in 1:n_reps) {
      exhat <- exhat_mat[, g]
      keep <- (exhat > trunc) & (exhat < 1 - trunc)
      
      y_trim <- y[keep]
      w_trim <- w[keep]
      exhat_trim <- exhat[keep]
      
      if (estimator == "PLR") {
        mxhat_trim <- mxhat_mat[keep, g]
        psi_a <- -(w_trim - exhat_trim)^2
        psi_b <- (y_trim - mxhat_trim) * (w_trim - exhat_trim)
      } else if (estimator == "IRM") {
        m0 <- mwhat0_mat[, g]
        m1 <- mwhat1_mat[, g]
        m0_trim <- m0[keep]
        m1_trim <- m1[keep]
        psi_a <- rep(-1, length(y_trim))
        psi_b <- m1_trim - m0_trim +
          w_trim * (y_trim - m1_trim) / exhat_trim -
          (1 - w_trim) * (y_trim - m0_trim) / (1 - exhat_trim)
      }
      
      # compute theta and se for every new cross-fitting repetition
      theta_hat <- -sum(psi_b) / sum(psi_a)
      psi <- theta_hat * psi_a + psi_b
      sigma2 <- mean(psi^2) / mean(psi_a)^2
      se_hat <- sqrt(sigma2 / n_obs)
      
      # compute p-value
      p_value <- 1 - pnorm(theta_hat / se_hat)
      
      theta_vals[g] <- theta_hat
      p_values[g] <- p_value
      se_vals[g] <- se_hat
      
      # compute aggregate statistic
      if (g > 4) {
        agg_stat <- mean(p_values[1:g])
        print(agg_stat)
        agg_p_vals[g] <- agg_stat
        var_est <- (1 / g) * sum((p_values[1:g] - agg_stat)^2) / (g - 1)
        cat(sprintf("        g = %d | mean(p) = %.4f | Var(p) = %.6e\n", g, agg_stat, var_est))
        
        # stopping rule 
        if (var_est <= cv_xi_beta) {
          cat(sprintf("        ✅ Stop at g = %d (Variance threshold met)\n", g))
          g_stop <- g
          break
  
        }
      }
    }
     
    if (is.na(g_stop)) g_stop <- n_reps
    
    # save results
    summary_table <- rbind(summary_table, data.table(
      learner = learner_name,
      estimator = estimator,
      g_stop = g_stop,
      var_est = var_est,
      var_p = var(p_values[1:g_stop], na.rm = TRUE),
      mean_p = mean(p_values[1:g_stop], na.rm = TRUE),
      mean_theta = mean(theta_vals[1:g_stop], na.rm = TRUE)
    ), fill = TRUE)
  }
}


########## Table ##########

# nice names


summary_table_print <- as.data.frame(summary_table)
names(summary_table_print) <- c(
  "Learner",              # learner
  "Estimator",            # estimator
  "Stopping Repetition",  # g_stop
  "Variance (Threshold)", # var_est
  "Empirical Variance",   # var_p
  "Mean p-value",         # mean_p
  "Mean Estimate"         # mean_theta
)



summary_table_print$Estimator <- ifelse(
  summary_table_print$Estimator == "IRM", "AIPW", summary_table_print$Estimator
)


# change notation
summary_table_print$`Variance (Threshold)` <- formatC(summary_table_print$`Variance (Threshold)`, format = "e", digits = 2)
summary_table_print$`Empirical Variance`   <- formatC(summary_table_print$`Empirical Variance`, format = "e", digits = 2)
summary_table_print$`Mean p-value`         <- formatC(summary_table_print$`Mean p-value`, format = "e", digits = 2)
summary_table_print$Learner <- gsub("_", " ", summary_table_print$Learner)
summary_table_print$`Variance (Threshold)` <- NULL


# create and save table
kable(
  summary_table_print,
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Aggregation stopping points for different learners and estimators.",
  label = "tab:aggregation_pval_stopping",
  digits = 3,
  align = "c"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 8
  ) %>%
  cat(file = "aggregation_summary_pval.tex")


