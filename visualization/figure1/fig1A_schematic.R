# Figure 1A: schematic workflow panel (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1A_schematic.R

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

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

box <- function(x0, y0, x1, y1, fill, label, cex = 1.0) {
  rect(x0, y0, x1, y1, col = fill, border = NA)
  rect(x0, y0, x1, y1, col = NA, border = "grey35", lwd = 1.0)
  text((x0 + x1) / 2, (y0 + y1) / 2, label, cex = cex, font = 2)
}

arrow_right <- function(x0, y0, x1, y1) {
  arrows(x0, y0, x1, y1, length = 0.10, angle = 20, lwd = 1.2, col = "grey25")
}

# Normalised device coordinates
x <- c(0.02, 0.26, 0.50, 0.74, 0.98)
y0 <- 0.18
y1 <- 0.86

png(file.path(out_dir, "fig1A_schematic.png"), width = 1900, height = 440, res = 200, type = "cairo")
op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
plot.new()
box(x[1], y0, x[2], y1, fill = "#EAF3FB", label = "Branching Rule (BR)\nBR1 / BR3 / BR5", cex = 1.00)
box(x[2], y0, x[3], y1, fill = "#EAF7EE", label = "Agent-Based Model (ABM)\nFlow + migration + polarity", cex = 1.00)
box(x[3], y0, x[4], y1, fill = "#EAF0FB", label = "Simulation\n5 days", cex = 1.00)
box(x[4], y0, x[5], y1, fill = "#FFF4E8", label = "Outcome\nTopology + hierarchy", cex = 1.00)
arrow_right(x[2] + 0.01, 0.52, x[2] + 0.08, 0.52)
arrow_right(x[3] + 0.01, 0.52, x[3] + 0.08, 0.52)
arrow_right(x[4] + 0.01, 0.52, x[4] + 0.08, 0.52)
par(op)
dev.off()

cat("Wrote figures/panels/fig1A_schematic.png\n")
