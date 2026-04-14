# Figure 2D: Stable vs unstable traces (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2D_stable_unstable_traces.R

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

br5 <- read.csv("data/br5_compare_multi_diam.csv")
br5_stable <- br5 %>% dplyr::filter(stable == 1)
br5_stab_summary <- br5_stable %>%
  dplyr::group_by(time_days) %>%
  dplyr::summarise(prox = mean(mean_diam_prox), dist = mean(mean_diam_dist),
                   .groups = "drop")
br5_unst <- read.csv("data/br5_unstable_mean_diam.csv")

p <- ggplot() +
  geom_line(data = br5_stab_summary,
            aes(x = time_days, y = prox, colour = "Stable proximal"),
            linewidth = 1) +
  geom_line(data = br5_stab_summary,
            aes(x = time_days, y = dist, colour = "Stable distal"),
            linewidth = 1) +
  geom_line(data = br5_unst,
            aes(x = time_days, y = mean_diam_prox, colour = "Unstable proximal"),
            linewidth = 1, linetype = "dashed") +
  geom_line(data = br5_unst,
            aes(x = time_days, y = mean_diam_dist, colour = "Unstable distal"),
            linewidth = 1, linetype = "dashed") +
  scale_colour_manual(values = c("Stable proximal" = "#B9BBBD",
                                 "Stable distal"   = "#DD6336",
                                 "Unstable proximal" = "#0088C1",
                                 "Unstable distal"   = "#805D76"),
                      name = "") +
  labs(x = "Time (days)", y = "Mean diameter (\u00b5m)") +
  theme_abm() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 9))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2D_stable_unstable_traces.png"), p, width = 6.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2D_stable_unstable_traces.png\n")
