# Figure 1H: BR5 daily loss incidence vs alpha (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1H_daily_loss_incidence.R

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

daily_loss <- read.csv("data/main_figure1_daily_loss_by_alpha.csv")

day_cols <- c(
  "Day 1" = "#A6CEE2",
  "Day 2" = "#FAA41A",
  "Day 3" = "#CB79A7",
  "Day 4" = "#4F95CE",
  "Day 5" = "#8C8E90"
)
daily_loss$day_label <- factor(daily_loss$day_label,
                               levels = c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5"))

p <- ggplot(daily_loss, aes(x = alpha, y = loss_pct, colour = day_label)) +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = day_cols, name = "") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  labs(x = expression(alpha), y = "Bifurcation loss in day (%)") +
  theme_abm() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 9)
  )

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1H_daily_loss_incidence.png"), p, width = 6.0, height = 3.0,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1H_daily_loss_incidence.png\n")
