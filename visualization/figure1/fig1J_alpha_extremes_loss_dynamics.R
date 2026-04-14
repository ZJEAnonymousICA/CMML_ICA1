# Figure 1J: BR5 alpha-extreme loss dynamics (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1J_alpha_extremes_loss_dynamics.R

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

alpha_extreme_dyn <- read.csv("data/main_figure1_alpha_extremes_regression_dynamics.csv")
alpha_extreme_dyn$alpha_label <- factor(alpha_extreme_dyn$alpha_label,
                                       levels = c("alpha = 0.00", "alpha = 1.00"))
alpha_extreme_dyn$branch <- factor(alpha_extreme_dyn$branch,
                                  levels = c("Proximal-first loss", "Distal-first loss"))

p <- ggplot(alpha_extreme_dyn,
            aes(x = time_days, y = loss_pct, colour = branch)) +
  geom_step(linewidth = 0.80) +
  facet_wrap(~ alpha_label, nrow = 1) +
  scale_colour_manual(values = c(
    "Proximal-first loss" = unname(pal_branch["proximal"]),
    "Distal-first loss" = unname(pal_branch["distal"])
  ), name = "") +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 60, by = 10)) +
  labs(x = "Time (days)", y = "Cumulative bifurcation loss (%)") +
  theme_abm() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 9)
  )

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1J_alpha_extremes_loss_dynamics.png"), p, width = 7.2, height = 3.0,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1J_alpha_extremes_loss_dynamics.png\n")
