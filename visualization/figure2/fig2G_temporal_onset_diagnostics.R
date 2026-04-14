# Figure 2G: Temporal onset diagnostics (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2G_temporal_onset_diagnostics.R

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

temporal_df <- dplyr::bind_rows(
  read.csv("data/time_to_imbalance_cumulative.csv"),
  read.csv("data/time_to_runaway_cumulative.csv")
) %>%
  dplyr::mutate(
    alpha_label = factor(sprintf("\u03b1=%.2f", alpha),
                         levels = c("\u03b1=0.10", "\u03b1=0.45", "\u03b1=1.00")),
    event = factor(event, levels = c("Hierarchy onset", "One-sided routing"))
  )

p <- ggplot(temporal_df, aes(x = time_days, y = cumulative_pct, colour = alpha_label)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ event, ncol = 2) +
  scale_colour_manual(values = c("\u03b1=0.10" = "#A0CEC9",
                                 "\u03b1=0.45" = "#7E4C21",
                                 "\u03b1=1.00" = "#22215B"),
                      name = "") +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Time (days)", y = "Runs reaching state (%)") +
  theme_abm() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2G_temporal_onset_diagnostics.png"), p, width = 7.4, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2G_temporal_onset_diagnostics.png\n")
