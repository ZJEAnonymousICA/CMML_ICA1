# Figure 1G: BR5 joint-success sweep across alpha (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1G_br5_joint_success_sweep.R

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

sweep_summary <- read.csv("data/sweep_summary.csv")
plateau_summary <- read.csv("data/alpha_plateau_summary.csv")
peak_alpha <- plateau_summary$peak_alpha[1]
plateau_lo <- plateau_summary$plateau_alpha_lo[1]
plateau_hi <- plateau_summary$plateau_alpha_hi[1]

p <- ggplot(sweep_summary, aes(x = alpha, y = joint_success_pct)) +
  annotate("rect", xmin = plateau_lo, xmax = plateau_hi,
           ymin = -Inf, ymax = Inf, fill = "#2CA089", alpha = 0.07) +
  geom_ribbon(aes(ymin = joint_success_ci_lo, ymax = joint_success_ci_hi),
              fill = "#2CA089", alpha = 0.16) +
  geom_line(colour = "#2CA089", linewidth = 0.85) +
  geom_vline(xintercept = peak_alpha, colour = "#D36027",
             linetype = "dotted", linewidth = 0.45) +
  geom_point(data = sweep_summary %>% dplyr::filter(abs(alpha - peak_alpha) < 1e-9),
             colour = "#D36027", fill = "#D36027", size = 2.3) +
  geom_point(data = sweep_summary %>% dplyr::filter(abs(alpha - 0.45) < 1e-9),
             colour = "#1B4332", fill = "#1B4332", size = 2.3) +
  annotate("text", x = peak_alpha, y = 95,
           label = sprintf("alpha == %.2f*\" peak\"", peak_alpha),
           parse = TRUE, family = "Arial", size = 3.0) +
  annotate("text", x = (plateau_lo + plateau_hi) / 2, y = 8,
           label = sprintf("plateau: %.2f-%.2f", plateau_lo, plateau_hi),
           family = "Arial", size = 3.0) +
  annotate("text", x = 0.63, y = 72,
           label = expression(alpha == 0.45*" reference"),
           family = "Arial", size = 3.0) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = expression(alpha), y = "Joint success (%)") +
  theme_abm()

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1G_br5_joint_success_sweep.png"), p, width = 6.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1G_br5_joint_success_sweep.png\n")
