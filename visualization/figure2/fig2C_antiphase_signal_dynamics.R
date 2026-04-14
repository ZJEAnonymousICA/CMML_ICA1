# Figure 2C: Anti-phase probability signal dynamics (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2C_antiphase_signal_dynamics.R

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

probs <- read.csv("data/br5_stable_probs.csv")
p <- ggplot(probs, aes(x = step)) +
  geom_line(aes(y = P_tau1, colour = "P_\u03c4"), linewidth = 1) +
  geom_line(aes(y = P_n1, colour = "P_n"), linetype = "dashed", linewidth = 1) +
  scale_colour_manual(values = c("P_\u03c4" = "#D3D478", "P_n" = "#87648F"),
                      name = "") +
  labs(x = "Time step", y = "Component value") +
  theme_abm() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 9))

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2C_antiphase_signal_dynamics.png"), p, width = 5.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2C_antiphase_signal_dynamics.png\n")
