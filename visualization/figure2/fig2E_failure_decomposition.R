# Figure 2E: Failure decomposition across alpha (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2E_failure_decomposition.R

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

plateau_summary <- read.csv("data/alpha_plateau_summary.csv")
plateau_lo <- plateau_summary$plateau_alpha_lo[1]
plateau_hi <- plateau_summary$plateau_alpha_hi[1]

failure_df <- read.csv("data/alpha_failure_decomposition.csv")
regime_df <- dplyr::bind_rows(
  failure_df %>% dplyr::transmute(alpha = alpha, metric = "Topology preserved",
                                  value = topology_preserved_pct),
  failure_df %>% dplyr::transmute(alpha = alpha, metric = "Hierarchy among survivors",
                                  value = hierarchy_among_survivors_pct),
  failure_df %>% dplyr::transmute(alpha = alpha, metric = "Joint success",
                                  value = joint_success_pct)
)

p <- ggplot(regime_df, aes(x = alpha, y = value, colour = metric)) +
  annotate("rect", xmin = plateau_lo, xmax = plateau_hi,
           ymin = -Inf, ymax = Inf, fill = "#2CA089", alpha = 0.05) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Topology preserved" = "#445598",
                                 "Hierarchy among survivors" = "#D36027",
                                 "Joint success" = "#2CA089"),
                      name = "") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 20)) +
  labs(x = "\u03b1", y = "Day-5 proportion (%)") +
  theme_abm() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 10))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2E_failure_decomposition.png"), p, width = 6.6, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2E_failure_decomposition.png\n")
