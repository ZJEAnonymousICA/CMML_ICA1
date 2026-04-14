# Figure 2F: First-loss topology bias across alpha (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2F_first_loss_topology_bias.R

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

bias_df <- read.csv("data/loss_topology_bias_summary.csv") %>%
  dplyr::transmute(
    alpha = alpha,
    `Proximal among losses` = proximal_share_of_losses_pct,
    `Distal among losses` = distal_share_of_losses_pct
  ) %>%
  tidyr::pivot_longer(cols = c(`Proximal among losses`, `Distal among losses`),
                      names_to = "metric", values_to = "share_pct")

p <- ggplot(bias_df, aes(x = alpha, y = share_pct, colour = metric)) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Proximal among losses" = "#4F95CE",
                                 "Distal among losses" = "#CE2627"),
                      name = "") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "\u03b1", y = "Share of failed runs (%)") +
  theme_abm() +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 10))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2F_first_loss_topology_bias.png"), p, width = 6.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2F_first_loss_topology_bias.png\n")
