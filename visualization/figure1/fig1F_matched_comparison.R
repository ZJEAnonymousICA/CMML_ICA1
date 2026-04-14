# Figure 1F: Matched BR3 vs BR5 endpoint comparison (standalone)
#
# Run:
#   /usr/local/bin/Rscript --vanilla visualization/figure1/fig1F_matched_comparison.R

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

endpoint_panel <- read.csv("data/figure1d_endpoint_panel.csv")
endpoint_ann <- read.csv("data/figure1d_endpoint_annotations.csv")

endpoint_panel$metric <- factor(
  endpoint_panel$metric,
  levels = c("Topology preserved", "Hierarchy among survivors", "Joint success")
)
endpoint_ann$metric <- factor(endpoint_ann$metric, levels = levels(endpoint_panel$metric))

rule_cols <- c(BR3 = unname(pal_rules["BR3"]), BR5 = unname(pal_rules["BR5"]))
endpoint_ann_joint <- endpoint_ann %>% dplyr::filter(metric == "Joint success")

p <- ggplot(endpoint_panel, aes(x = metric, y = pct, colour = rule, group = rule)) +
  geom_linerange(aes(ymin = ci_lo, ymax = ci_hi),
                 position = position_dodge(width = 0.40),
                 linewidth = 0.75) +
  geom_point(position = position_dodge(width = 0.40), size = 4) +
  geom_text(data = endpoint_ann_joint,
            aes(x = metric, y = y, label = label),
            inherit.aes = FALSE, family = "Arial", size = 2.8) +
  scale_colour_manual(values = rule_cols, name = "") +
  scale_x_discrete(labels = c("Topology\npreserved",
                              "Hierarchy\namong survivors\n(stable only)",
                              "Joint\nsuccess")) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, by = 20)) +
  labs(x = NULL, y = "Seeds meeting criterion (%)") +
  theme_abm() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 9)
  )

out_dir <- file.path("figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig1F_matched_comparison.png"), p, width = 6.0, height = 3.2,
       dpi = 300, bg = "white")

cat("Wrote figures/panels/fig1F_matched_comparison.png\n")
