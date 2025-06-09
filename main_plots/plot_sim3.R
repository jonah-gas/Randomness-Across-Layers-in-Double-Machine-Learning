# check, install, and load necessary packages
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "patchwork", "viridis", "scales", "tools")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


########## Boxplot ##########


# mapping for nice names
pretty_learner <- function(x) {
  recode(x,
                "lasso" = "Lasso",
                "random_forest" = "Random Forest",
                "xgboost" = "XGBoost"
  )
}

pretty_dgp <- function(x) {
  recode(x, "DGP_1" = "DGP 1", "DGP_2" = "DGP 2")
}

# plot function

plot_theta_by_dgp <- function(df_long, dgp_label, theta_true) {
  ggplot(df_long[dgp == dgp_label], aes(x = method_plot, y = theta_value, fill = theta_type)) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.85, outlier.shape = NA) +
    geom_hline(yintercept = theta_true, linetype = "dashed", color = "red") +
    scale_fill_viridis_d(name = "Aggregation:") +
    coord_cartesian(ylim = c(0, 2)) +
    labs(
      title = paste("DGP", substr(dgp_label, 5, 5)), 
      x = NULL,   
      y = "Point Estimate"
    ) +
    coord_cartesian(ylim = c(theta_true - 0.75, theta_true + 0.75)) +   # Y-Limit um true theta
    facet_wrap(~ learner, ncol = 3, labeller = as_labeller(pretty_learner)) +  # pretty_title entfernt
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(margin = margin(r = 10)),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10), 
      panel.spacing = unit(1.25, "lines"),
      strip.text = element_text(size = 12)
    )
}






# data preparation function
prep_dgp_data <- function(selected_dgp, theta_true) {
  input_filename <- paste0("results/sim3_results_", selected_dgp, ".rds")
  final_estimates <- readRDS(input_filename)
  df <- as.data.table(final_estimates)
  df[, dgp := selected_dgp]
  df_long <- melt(
    df[method %in% c("IRM", "PLR")],
    id.vars = c("learner", "method", "dgp"),
    measure.vars = c("mean_theta", "median_theta"),
    variable.name = "theta_type",
    value.name = "theta_value"
  )
  df_long[, theta_type := factor(theta_type,
                                 levels = c("mean_theta", "median_theta"),
                                 labels = c("Mean-based", "Median-based"))]
  df_long[, method_plot := fifelse(method == "IRM", "AIPW", method)]
  df_long[, theta_true := theta_true]
  return(df_long)
}

# load and combine DGPs
df_long_DGP1 <- prep_dgp_data("DGP_1", 0.5)
df_long_DGP2 <- prep_dgp_data("DGP_2", 1)
df_long_both <- rbind(df_long_DGP1, df_long_DGP2)

# solo plots
plot_dgp1 <- plot_theta_by_dgp(df_long_both, "DGP_1", 0.5)
plot_dgp2 <- plot_theta_by_dgp(df_long_both, "DGP_2", 1)

# combine plots 
combined_plot <- plot_dgp1 / plot_dgp2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)






########## MSE ##########

# preparation and plot function
plot_all_mse_faceted <- function(df, dgp_label = NULL, ncol = 3) {
  df_mse <- rbindlist(list(
    df[, .(
      bias2 = (mean(mean_theta) - unique(theta_true))^2,
      var = var(mean_theta),
      aggregation = "Mean-based"
    ), by = .(learner, method)],
    
    df[, .(
      bias2 = (mean(median_theta) - unique(theta_true))^2,
      var = var(median_theta),
      aggregation = "Median-based"
    ), by = .(learner, method)]
  ))
  
  
  df_mse_long <- melt(
    df_mse,
    id.vars = c("learner", "method", "aggregation"),
    measure.vars = c("bias2", "var"),  
    variable.name = "Component:",
    value.name = "MSE"
  )
  
  
  all_levels <- c("Mean\nAIPW", "Median\nAIPW", " ", "Mean\nPLR", "Median\nPLR")
  df_mse_long[, aggregation_clean := gsub("-based", "", aggregation)]
  df_mse_long[, method_plot := fifelse(method == "IRM", "AIPW", method)]
  df_mse_long[, agg_method := factor(
    paste(aggregation_clean, method_plot, sep = "\n"),
    levels = all_levels
  )]
  df_mse_long[, learner_facet := toTitleCase(gsub("_", " ", learner))]
  print(df_mse)
  print(df_mse_long)
  
  learners <- unique(df_mse_long$learner_facet)
  for (l in learners) {
    dummy <- data.table(
      learner = NA,
      method = NA,
      aggregation = NA,
      Component = NA,
      MSE = NA,
      aggregation_clean = NA,
      method_plot = NA,
      agg_method = factor(" ", levels = all_levels),
      learner_facet = l
    )
    df_mse_long <- rbind(df_mse_long, dummy, fill = TRUE)
  }
  
  global_y_max <- df_mse_long[!is.na(MSE), sum(MSE), by = .(learner, method_plot, aggregation_clean)][, max(V1)] * 1.1
  
  ggplot(df_mse_long, aes(
    x = agg_method,
    y = MSE,
    fill = factor(Component, levels = c("var", "bias2"))
  )) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.5,
      alpha = 0.9,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ learner_facet,
      ncol = ncol,
      strip.position = "top"
    ) +
    scale_fill_manual(
      name = "Component",
      values = c("var" = "#5ec962", "bias2" = "#3b528b"),
      labels = c("Variance", expression(Bias^2))
    ) +
    labs(
      title = dgp_label,
      x = NULL,
      y = "MSE"
    ) +
    coord_cartesian(ylim = c(0, global_y_max)) +
    scale_y_continuous(breaks = pretty_breaks(n=3))+
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size= 12),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.spacing = unit(1.25, "lines")
    )
}




# read in data
dt1 <- as.data.table(readRDS("results/sim3_results_DGP_1.rds"))
dt1[, theta_true := 0.5]
dt2 <- as.data.table(readRDS("results/sim3_results_DGP_2.rds"))
dt2[, theta_true := 1]

# solo plots 
p1 <- plot_all_mse_faceted(dt1, dgp_label = "DGP 1")
p2 <- plot_all_mse_faceted(dt2, dgp_label = "DGP 2")

#combine plots
combined_mse <- p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(combined_mse)




########## Coverage ##########

# preparation and plot function
plot_coverage_faceted <- function(
    df,
    dgp_label = NULL
) {
  df_coverage_all <- rbindlist(list(
    df[, .(
      learner,
      method,
      aggregation = "Mean",
      coverage = (theta_true >= (mean_theta - 1.96 * se_mean)) &
        (theta_true <= (mean_theta + 1.96 * se_mean))
    )],
    df[, .(
      learner,
      method,
      aggregation = "Median",
      coverage = (theta_true >= (median_theta - 1.96 * se_median)) &
        (theta_true <= (median_theta + 1.96 * se_median))
    )],
    df[, .(
      learner,
      method,
      aggregation = "IF",
      coverage = (theta_true >= (mean_theta - 1.96 * se_IF)) &
        (theta_true <= (mean_theta + 1.96 * se_IF))
    )]
  ))
  
  df_cov_summary <- df_coverage_all[, .(
    coverage = mean(coverage, na.rm = TRUE)
  ), by = .(learner, method, aggregation)]
  
  df_cov_summary[, learner_facet := toTitleCase(gsub("_", " ", learner))]
  df_cov_summary[, method_plot := fifelse(method == "IRM", "AIPW", method)]
  df_cov_summary[, aggregation := factor(aggregation, levels = c("Mean", "Median", "IF"))]
  
  ggplot(df_cov_summary, aes(
    x = aggregation,
    y = coverage,
    fill = aggregation
  )) +
    geom_hline(yintercept = c(0, 1), color = "grey90") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    geom_point(size = 5, shape = 21, stroke = 1, na.rm = TRUE) +
    facet_grid(cols = vars(learner_facet), rows = vars(method_plot)) +
    scale_fill_manual(
      values = c(
        "Mean" = "#440154",
        "Median" = "#fde725",
        "IF" = "#21918c"
      ),
      name = "Aggregation"
    ) +
    scale_y_continuous(limits = c(0.7, 1)) +
    labs(
      title = dgp_label,
      x = NULL,
      y = "Coverage"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.spacing = unit(1.25, "lines")
    )
}

# read in data
dt1 <- as.data.table(readRDS("results/sim3_results_DGP_1.rds"))
dt1[, theta_true := 0.5]
dt2 <- as.data.table(readRDS("results/sim3_results_DGP_2.rds"))
dt2[, theta_true := 1]

# solo plots
cov_block1 <- plot_coverage_faceted(dt1, dgp_label = "DGP 1")
cov_block2 <- plot_coverage_faceted(dt2, dgp_label = "DGP 2")

#combine plots
combined_cov <- cov_block1 / cov_block2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "none")
print(combined_cov)











