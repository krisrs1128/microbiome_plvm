#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Visualize a combined experiment across various NMF / fitting
## parameterizations.
##
## author: kriss1@stanford.edu

## ---- libraries-nmf-vis ----
## assumed running from NMF directory
library("data.table")
library("jsonlite")
library("ggplot2")
library("nmfSim")
library("plyr")
library("dplyr")
library("ggscaffold")

## ---- beta-reshape ----
## extract beta information from the fits
base_dir = "~/Desktop/microbiome_plvm"
#base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
nmf_dir <- file.path(base_dir, "src", "sim", "nmf")
fits_dir <- file.path(nmf_dir, "fits")
figure_dir <- file.path(base_dir, "doc", "figure")
dir.create(figure_dir, recursive = TRUE)

fits <- list.files(fits_dir, "fit-*", full.names = TRUE)
expers <- fromJSON(
  file.path(nmf_dir, "config.json"),
  simplifyVector = FALSE
)

beta_fits <- reshape_all_samples(
  fits,
  file.path(nmf_dir, "config.json"),
  "beta",
  c("j", "k")
)
beta_fits$method <- basename(as.character(beta_fits$method))
beta_fits$method <- beta_fits$method %>%
  revalue(
    c(
    "nmf_gamma_poisson.stan" = "GaP",
    "nmf_gamma_poisson_zero.stan" = "Z-GaP"
    )
  )

## ---- visualizebetas-prep -----
## Visualize the fitted betas, according to a few different simulation properties
plot_opts <- list(
  "x" = "value_1",
  "y" = "value_2",
  "fill" = "log(..level..)",
  "fill_type" = "gradient",
  "facet_terms" = c("N", "inference", "P"),
  "group" = "i",
  "alpha" = 0.05,
  "h" = 0.1,
  "mean_col" = "#e34a33",
  "x_lim" = c(0, 4),
  "y_lim" = c(0, 5),
  "text_size" = 4,
  "panel_border" = 0.1
)

## ---- visualize-zinf-betas-prep ----
zinf_data <- beta_fits
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
plot_opts$group <- "i"

## ---- zinf-betas-errors ----
zinf_data$a <- as.character(zinf_data$a)
zinf_data$zero_inf_prob <- as.character(zinf_data$zero_inf_prob)

perf <- zinf_data %>%
  mutate(
    a = substr(a, 1, 2),
    value_1 = sqrt(value_1),
    value_2 = sqrt(value_2),
    truth_1 = sqrt(truth_1),
    truth_2 = sqrt(truth_2)
  ) %>%
  group_by(j, inference, method, zero_inf_prob, a, N, P) %>%
  summarise(
    error = mean(sqrt((value_1 - truth_1) ^ 2 + (value_2 - truth_2) ^ 2)),
    error_bar = sd(value_1)
  )

theme_set(ggscaffold::min_theme(list(border_size = .7)))
method_cols <- c("#ae7664", "#64ae76", "#7664ae")
p <- ggplot(perf) +
  geom_abline(slope = 1, alpha = 0.6, size = 0.3) +
  geom_point(
    aes(x = error, y = error_bar, col = inference, shape = zero_inf_prob),
    size = 0.7, alpha = 0.6
  ) +
  facet_grid(P ~ N + a) +
  scale_color_manual(values = method_cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Error", y = "SD (k = 1)", col = "Inference")

ggsave(
  file.path(base_dir, "doc", "figure/beta_errors_nmf.pdf"),
  p,
  width = 5,
  height = 3
)

## ---- zinf-betas-contours ----
combined <- zinf_data %>%
  filter(P == 10)

plot_opts <- list(x = "sqrt(value_1)", y = "sqrt(value_2)",
                  group = "j", fill_type = "gradient", h = 0.1)

ggcontours(combined, plot_opts) +
  geom_text(
    data = combined %>%
      filter(iteration == 1),
    aes(
      x = sqrt(truth_1),
      y = sqrt(truth_2),
      label = j
    ),
    size = 2
  ) +
  geom_text(
    data = combined %>%
      group_by(j, N, a, b, zero_inf_prob, inference, method) %>%
      summarise(value_mean_1 = mean(value_1), value_mean_2 = mean(value_2)),
    aes(
      x = sqrt(value_mean_1),
      y = sqrt(value_mean_2),
      label = j
    ),
    size = 2, col = "#fc8d62") +
  facet_grid(method + inference + zero_inf_prob ~ N + a + b) +
  ylim(0, 3) +
  xlim(0, 3)
