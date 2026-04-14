# Figure 2B: Pressure-cell feedback (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2B_pressure_cell_feedback.R

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

pc_data <- read.csv("data/br5_stable_pressure_cells.csv")
p <- ggplot(pc_data, aes(x = time_days)) +
  geom_line(aes(y = dp_bifurc, colour = "Pressure drop"), linewidth = 1) +
  geom_line(aes(y = n_prox_mean * max(dp_bifurc) / max(n_prox_mean),
                colour = "Cell count (scaled)"), linewidth = 1) +
  scale_colour_manual(values = c("Pressure drop" = "#F18B89",
                                 "Cell count (scaled)" = "#465674"),
                      name = "") +
  scale_y_continuous(
    name = "Pressure drop (Pa)",
    sec.axis = sec_axis(~ . * max(pc_data$n_prox_mean) / max(pc_data$dp_bifurc),
                        name = "Mean cell count")
  ) +
  labs(x = "Time (days)") +
  theme_abm() +
  theme(axis.title.y.right = element_text(colour = "#465674"),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 9))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2B_pressure_cell_feedback.png"), p, width = 5.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2B_pressure_cell_feedback.png\n")
