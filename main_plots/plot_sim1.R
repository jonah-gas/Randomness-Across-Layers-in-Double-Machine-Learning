# check, install, and load necessary packages
packages <- c("tidyverse", "patchwork")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))




# read in data
data_dgp1 <- read.csv("results/sim1_results_DGP1.csv")
data_dgp2 <- read.csv("results/sim1_results_DGP2.csv")

# prepare data
prep_fun <- function(df) {
  df %>%
    pivot_longer(
      cols = starts_with("share_variance"),
      names_to = "variance_component",
      values_to = "share"
    ) %>%
    mutate(
      component_value = share * variance_total,
      variance_component = factor(
        variance_component,
        levels = c("share_variance_dgp", "share_variance_cf_rep", "share_variance_learner_rep"),
        labels = c("Sample", "Cross-fitting", "Learner")
      )
    )
}

data_long1 <- prep_fun(data_dgp1)
data_long2 <- prep_fun(data_dgp2)

# plot function
plot_fun <- function(data_long, dgp_title, y_breaks = NULL) {
  if (is.null(y_breaks)) {
    y_min <- min(data_long$component_value, na.rm = TRUE)
    y_max <- max(data_long$component_value, na.rm = TRUE)
    y_breaks <- pretty(c(y_min, y_max), n = 3)
   
  }
  
  
  ggplot(data_long, aes(x = factor(n_obs), y = component_value, fill = variance_component)) +
    geom_bar(stat = "identity") +
    facet_grid(n_trees ~ n_folds,
               labeller = labeller(
                 n_trees = function(x) paste("Trees:", x),
                 n_folds = function(x) paste("Folds:", x)
               )) +
    labs(
      x = "Number of Observations",
      y = "Total Variance of Point Estimates",
      title = dgp_title,
      fill = "Variance Share\nComponent:"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size =8),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),     
      legend.title = element_text(size = 11),
      strip.text = element_text(size = 8),
      legend.text = element_text(size = 8),                        
      panel.spacing = unit(1.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 13)
    ) +
    scale_fill_manual(
      values = c(
        "Sample" = "#FDE725FF",
        "Cross-fitting" = "#22A884FF",
        "Learner" = "#7570b3FF"
      )
    ) +
    scale_y_continuous(breaks = y_breaks)
}


# build individual plots
plot1 <- plot_fun(data_long1, "DGP1", y_breaks = c(0, 0.04, 0.08))
plot2 <- plot_fun(data_long2, "DGP2")

# combine plots
combined_plot <- plot1 / plot2 +
  plot_layout(guides = "collect", heights = c(1, 1), ncol = 1) &
  theme(
    plot.title = element_text(size = 13) 
  )

print(combined_plot)
