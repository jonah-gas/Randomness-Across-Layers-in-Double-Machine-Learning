# check, install, and load necessary packages
packages <- c("data.table", "ggplot2", "gridExtra", "grid")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))



# load data
final_estimates_sim3 <- readRDS("results/sim3_results_DGP_2.rds")
final_estimates_add_ana <- readRDS("results/sim_appendix_c_results.rds")

# prepare data
df_1 <- as.data.table(final_estimates_sim3)
df_1 <- df_1[learner == "random_forest" & method == "IRM"]
summary_1 <- df_1[, .(
  `sd(mean)`   = sd(mean_theta),
  `sd(median)` = sd(median_theta),
  `IF`         = mean(se_IF)
)]
plot_dt_1 <- melt(summary_1, variable.name = "type", value.name = "value")

df_2_list <- lapply(final_estimates_add_ana, as.data.table)
df_2 <- rbindlist(df_2_list)
df_2 <- df_2[learner == "random_forest" & method == "IRM"]
summary_2 <- df_2[, .(
  `sd(mean)`   = sd(mean_theta),
  `sd(median)` = sd(median_theta),
  `IF`         = mean(se_IF)
)]
plot_dt_2 <- melt(summary_2, variable.name = "type", value.name = "value")




# nice names
pretty_breaks_n <- function(lim, n) {
  pretty(seq(lim[1], lim[2], length.out = 100), n = n)
}

# nice breaks
expand_range <- function(x, prop = 0.1) {
  r <- range(x)
  pad <- diff(r) * prop
  c(r[1] - pad, r[2] + pad)
}
n_breaks <- 6

# solo plots
p1 <- ggplot(plot_dt_1, aes(x = type, y = value)) +
  geom_point(size = 5, color = "black") +
  labs(x = NULL, y = "Standard Error") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = breaks1, limits = expand_range(plot_dt_1$value, 0.1)) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),          # <-- HIER Y-TICK-GRÖSSE
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "none"
  ) +
  ggtitle("n = 400")


p2 <- ggplot(plot_dt_2, aes(x = type, y = value)) +
  geom_point(size = 5, color = "black") +
  labs(x = NULL, y = "Standard Error") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = breaks2, limits = expand_range(plot_dt_2$value, 0.1)) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),          # <-- HIER Y-TICK-GRÖSSE
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "none"
  ) +
  ggtitle("n = 2000")



# combine plot
blank <- nullGrob()
grid.arrange(
  p1, blank, p2,
  ncol = 3,
  widths = c(1, 0.20, 1) 
)
