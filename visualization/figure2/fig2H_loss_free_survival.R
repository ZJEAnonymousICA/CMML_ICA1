# Figure 2H: Loss-free survival curves (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2H_loss_free_survival.R

find_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 1:10) {
    if (file.exists(file.path(cur, "visualization", "viz_helpers.R")) && dir.exists(file.path(cur, "data"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate project root (missing visualization/viz_helpers.R). Run from Improved_Model_BRs.")
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
start_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(find_root(start_dir))

source("visualization/viz_helpers.R")

surv_df <- read.csv("data/loss_survival_curves.csv") %>%
  dplyr::mutate(alpha_label = factor(sprintf("\u03b1=%.2f", alpha),
                                     levels = c("\u03b1=0.10", "\u03b1=0.45", "\u03b1=1.00")))
surv_tests <- read.csv("data/loss_survival_tests.csv")
p_mid_high <- surv_tests %>%
  dplyr::filter(abs(alpha_a - 0.45) < 1e-9, abs(alpha_b - 1.00) < 1e-9) %>%
  dplyr::pull(p_value)

p <- ggplot(surv_df, aes(x = time_days, y = survival_pct, colour = alpha_label)) +
  geom_step(linewidth = 0.8) +
  annotate("text", x = 3.15, y = 94,
           label = sprintf("log-rank:\n\u03b1=0.45 vs 1.00\np = %.3g", p_mid_high),
           family = "Arial", size = 3.0, hjust = 0) +
  scale_colour_manual(values = c("\u03b1=0.10" = "#A0CEC9",
                                 "\u03b1=0.45" = "#7E4C21",
                                 "\u03b1=1.00" = "#22215B"),
                      name = "") +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  scale_y_continuous(limits = c(40, 100), breaks = seq(40, 100, by = 10)) +
  labs(x = "Time (days)", y = "Loss-free runs (%)") +
  theme_abm() +
  theme(legend.position = "none")

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2H_loss_free_survival.png"), p, width = 5.4, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2H_loss_free_survival.png\n")
