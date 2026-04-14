# Figure 2A: Conductance amplification (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure2/fig2A_conductance_amplification.R

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

cond <- read.csv("data/br5_stable_conductance.csv") %>% dplyr::filter(Ncell > 0)
scale_const <- median(cond$G / (cond$Ncell^4))
ref_line <- tibble::tibble(Ncell = seq(min(cond$Ncell), max(cond$Ncell), length.out = 100)) %>%
  dplyr::mutate(G = scale_const * Ncell^4)

p <- ggplot(cond, aes(x = Ncell, y = G)) +
  geom_point(colour = "#64180A", alpha = 1, size = 2) +
  geom_line(data = ref_line, aes(x = Ncell, y = G),
            inherit.aes = FALSE, colour = "#CE2627", linewidth = 0.8) +
  annotate("text", x = 14.8, y = max(cond$G) / 4,
           label = expression(G %prop% n^4),
           family = "Arial", size = 3.3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Cell number, n", y = "Conductance, G") +
  theme_abm()

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig2A_conductance_amplification.png"), p, width = 4.2, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig2A_conductance_amplification.png\n")
