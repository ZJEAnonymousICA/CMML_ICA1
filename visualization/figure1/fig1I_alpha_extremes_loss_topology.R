# Figure 1I: BR5 alpha-extreme loss topology bar chart (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1I_alpha_extremes_loss_topology.R

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

alpha_extreme_bar <- read.csv("data/main_figure1_alpha_extremes_regression_bar.csv")
alpha_extreme_bar$alpha_label <- factor(alpha_extreme_bar$alpha_label,
                                       levels = c("alpha = 0.00", "alpha = 1.00"))
alpha_extreme_bar$branch <- factor(alpha_extreme_bar$branch,
                                  levels = c("Proximal-first loss", "Distal-first loss"))

p <- ggplot(alpha_extreme_bar, aes(x = alpha_label, y = loss_pct, fill = branch)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.62,
           colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c(
    "Proximal-first loss" = unname(pal_branch["proximal"]),
    "Distal-first loss" = unname(pal_branch["distal"])
  ), name = "") +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 60, by = 10)) +
  labs(x = NULL, y = "Runs with first loss in branch (%)") +
  theme_abm() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 9)
  )

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1I_alpha_extremes_loss_topology.png"), p, width = 6.0, height = 3.0,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1I_alpha_extremes_loss_topology.png\n")
