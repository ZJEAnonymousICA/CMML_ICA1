# Figure 1C: BR1/BR3/BR5 mean diameter traces (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1C_diameter_traces.R

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

br1 <- read.csv("data/br1_mean_diam.csv")
br3 <- read.csv("data/br3_compare_multi_diam.csv")
br5 <- read.csv("data/br5_compare_multi_diam.csv") %>% dplyr::filter(stable == 1)

p1 <- plot_diam_single(br1, title = "BR1") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

p3 <- plot_diam_mean_std(br3, title = "BR3") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

p5 <- plot_diam_mean_std(br5, title = "BR5") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

panel <- p1 + p3 + p5 + plot_layout(ncol = 3)

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "fig1C_diameter_traces.png"), panel,
       width = 9.6, height = 3.0, dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1C_diameter_traces.png\n")
