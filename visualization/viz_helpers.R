# Common plotting helpers for the Improved_Model_BRs panel scripts.
#
# All panel scripts are intended to be runnable standalone via:
#   /usr/local/bin/Rscript --vanilla visualization/figureX/<panel>.R
#
# This file centralises theme/palette helpers so each panel script can stay short
# while remaining reproducible.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ---- Custom theme ----
theme_abm <- function(base_size = 14) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text             = element_text(colour = "black"),
      axis.text        = element_text(colour = "black"),
      axis.ticks       = element_line(colour = "black"),
      panel.border     = element_rect(colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text       = element_text(colour = "black"),
      legend.text      = element_text(colour = "black"),
      legend.title     = element_text(colour = "black"),
      plot.title       = element_text(colour = "black", size = base_size + 1,
                                      face = "bold", hjust = 0)
    )
}

# ---- Colour palettes ----
pal_branch <- c(Proximal = "#4F95CE", Distal = "#CE2627",
                proximal = "#4F95CE", distal = "#CE2627")
pal_vessel <- c(feeding = "#A6CEE2", proximal = "#4F95CE",
                draining = "#FAA41A", distal = "#CB79A7")
pal_rules <- c(BR1 = "#4F95CE", BR2 = "#CE2627", BR3 = "#FAA41A",
               BR4 = "#CB79A7", BR5 = "#8C8E90")
pal_prob <- c(P_tau = "#445598", P_n = "#E6B481")
pal_prob2 <- c(P_tau = "#D8197E", P_n = "#39A2CE")

# ---- Diameter time-series helpers ----
plot_diam_single <- function(df, title = "", loss_time = NULL) {
  p <- ggplot(df, aes(x = time_days)) +
    geom_line(aes(y = mean_diam_prox, colour = "Proximal")) +
    geom_line(aes(y = mean_diam_dist, colour = "Distal")) +
    scale_colour_manual(values = pal_branch, name = "Branch") +
    labs(x = "Time (days)", y = "Mean diameter (µm)", title = title) +
    theme_abm() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )

  if (!is.null(loss_time) && loss_time > 0) {
    dt <- 10 / 3
    loss_day <- loss_time * dt / 24
    p <- p + geom_vline(xintercept = loss_day, linetype = "dashed",
                        colour = "grey40")
  }
  p
}

plot_diam_mean_std <- function(df, title = "") {
  summary_df <- df %>%
    group_by(step, time_days) %>%
    summarise(
      prox_mean = mean(mean_diam_prox),
      prox_sd   = sd(mean_diam_prox),
      dist_mean = mean(mean_diam_dist),
      dist_sd   = sd(mean_diam_dist),
      .groups = "drop"
    )
  ggplot(summary_df, aes(x = time_days)) +
    geom_ribbon(aes(ymin = prox_mean - prox_sd, ymax = prox_mean + prox_sd),
                fill = pal_branch["proximal"], alpha = 0.2) +
    geom_ribbon(aes(ymin = dist_mean - dist_sd, ymax = dist_mean + dist_sd),
                fill = pal_branch["distal"], alpha = 0.2) +
    geom_line(aes(y = prox_mean, colour = "Proximal")) +
    geom_line(aes(y = dist_mean, colour = "Distal")) +
    scale_colour_manual(values = pal_branch, name = "Branch") +
    labs(x = "Time (days)", y = "Mean diameter (µm)", title = title) +
    theme_abm()
}

