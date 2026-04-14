# Figure 1D: Network snapshots (cropped from the snapshot grid image)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1D_network_snapshots.R

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

library(png)
library(ggplot2)
library(grid)

trim_bbox <- function(img, thresh = 0.995) {
  # img: h x w x 3 (or 4)
  rgb <- img[, , 1:3, drop = FALSE]
  mask <- apply(rgb, c(1, 2), function(v) any(v < thresh))
  rows <- which(rowSums(mask) > 0)
  cols <- which(colSums(mask) > 0)
  if (!length(rows) || !length(cols)) {
    return(img)
  }
  img[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
}

img <- png::readPNG("figures/main_figure1_snapshots.png")
h <- dim(img)[1]
w <- dim(img)[2]

left <- img[, 1:floor(w * 0.50), , drop = FALSE]
left <- trim_bbox(left)

g <- grid::rasterGrob(left, interpolate = TRUE)
p <- ggplot() +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1D_network_snapshots.png"), p, width = 6.6, height = 3.3,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1D_network_snapshots.png\n")
