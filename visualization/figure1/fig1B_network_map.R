# Figure 1B: A-branch network map with segment-group labels (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1B_network_map.R

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

coords <- read.csv("data/segment_coords.csv")

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- ggplot(coords) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               linewidth = 1.2, lineend = "round", colour = "#CE2627") +
  coord_equal(expand = TRUE) +
  theme_abm(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  annotate("text", x = 5e7, y = 1.03e8, label = "A-branch network", fontface = "bold", size = 4.2) +
  annotate("text", x = 5e7, y = 9.2e7, label = "Distal (seg 20-39)", colour = "#CE2627", size = 3.4) +
  annotate("text", x = 3.0e7, y = 4.2e7, label = "Proximal (seg 5-14)", colour = "#CE2627", size = 3.4) +
  annotate("text", x = -4.0e6, y = 2.5e7, label = "Feeding\n(seg 0-4)", colour = "#4F95CE", size = 3.1, hjust = 0) +
  annotate("text", x = 4.0e7, y = 1.0e7, label = "Draining\n(seg 15-19)", colour = "#4F95CE", size = 3.1) +
  annotate("text", x = 0.0, y = 5.2e7, label = "Diverge\n(node 5)", size = 3.0, hjust = 1.1) +
  annotate("text", x = 1.0e8, y = 5.2e7, label = "Converge\n(node 15)", size = 3.0, hjust = -0.1)

ggsave(file.path(out_dir, "fig1B_network_map.png"), p, width = 6.6, height = 3.3,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1B_network_map.png\n")
