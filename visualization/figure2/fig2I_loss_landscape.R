# Figure 2I: Loss landscape over alpha and time (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2I_loss_landscape.R

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

sweep_pct <- read.csv("data/sweep_loss_pct.csv")
p <- ggplot(sweep_pct, aes(x = alpha, y = time_days, fill = loss_pct)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlGn", direction = -1, name = "Loss (%)") +
  labs(x = "\u03b1", y = "Time (days)") +
  theme_abm() +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 9))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2I_loss_landscape.png"), p, width = 8.2, height = 2.6,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2I_loss_landscape.png\n")
