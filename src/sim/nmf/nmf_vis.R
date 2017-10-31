#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Visualize a combined experiment across various NMF / fitting
## parameterizations.
##
## author: kriss1@stanford.edu

## ---- libraries-nmf-vis ----
## assumed running from NMF directory
library("jsonlite")
library("tidyverse")
library("nmfSim")
library("ggscaffold")
theme_set(min_theme(list(border_size = 0.7)))

plot_contours <- function(combined, plot_opts, ymax = 12, xmax = 12) {
  ggcontours(combined, plot_opts) +
    geom_segment(
      data = posterior_means,
      aes(
        x = sqrt(value_mean_1),
        y = sqrt(value_mean_2),
        xend = sqrt(truth_1),
        yend = sqrt(truth_2)
      ),
      size = 0.05,
      alpha = 0.5
    ) +
    geom_point(
      data = posterior_means,
      aes(
        x = sqrt(truth_1),
        y = sqrt(truth_2)
      ),
      size = 0.2,
      alpha = 0.2
    ) +
    geom_point(
      data = posterior_means,
      aes(
        x = sqrt(value_mean_1),
        y = sqrt(value_mean_2)
      ),
      size = 0.2,
      alpha = 0.5,
      col = "#fc8d62"
    ) +
    coord_fixed() +
    scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 12), expand = c(0, 0)) +
    facet_grid(method + inference ~ EN + zero_inf_prob) +
    labs(x = expression(sqrt(hat(beta)[1])), y = expression(sqrt(hat(beta)[2])))
}


## ---- beta-reshape ----
## extract beta information from the fits
base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
nmf_dir <- file.path(base_dir, "src", "sim", "nmf")
fits_dir <- file.path(nmf_dir, "fits")
figure_dir <- file.path(base_dir, "doc", "figure")
dir.create(figure_dir, recursive = TRUE)

fits <- list.files(fits_dir, "fit-*", full.names = TRUE)
expers <- fromJSON(
  file.path(nmf_dir, "config.json"),
  simplifyVector = FALSE
)

exper_ids <- paste0(sapply(expers, function(x) x$id), ".rda")
expers <- expers[exper_ids %in% basename(fits)]

beta_fits <- reshape_all_samples(
  fits,
  file.path(nmf_dir, "config.json"),
  "beta",
  c("j", "k")
) %>%
  mutate(
    inference = recode(inference, "gibbs" = "mcmc"),
    D = paste0("D = ", N),
    V = paste0("V = ", P)
  )

sort_levels <- function(x) {
  new_levels <- sort(as.numeric(gsub("[^0-9]+", "", unique(x))))
  prefix <- gsub("[0-9]+", "", x[1])
  factor(x, levels = paste0(prefix, new_levels))
}

beta_fits$D <- sort_levels(beta_fits$D)
beta_fits$V <- sort_levels(beta_fits$V)

beta_fits$method <- basename(as.character(beta_fits$method))
beta_fits$method <- beta_fits$method %>%
  recode(
    "nmf_gamma_poisson.stan" = "GaP",
    "nmf_gamma_poisson_zero.stan" = "Z-GaP"
  )
beta_fits$method <- factor(
  beta_fits$method,
  levels = c("bootstrap", "mcmc", "vb")
)

## ---- visualizebetas-prep -----
## Visualize the fitted betas, according to a few different simulation properties
plot_opts <- list(
  "x" = "value_1",
  "y" = "value_2",
  "fill" = "log(..level..)",
  "fill_type" = "gradient",
  "facet_terms" = c("D", "inference", "V"),
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
## zinf_data$a <- as.character(zinf_data$a)
zinf_data$zero_inf_prob <- as.character(zinf_data$zero_inf_prob)
zinf_data <- zinf_data %>%
  mutate(
    EN = paste0("E[N] = ", P * b * 20 / as.numeric(a))
  )

perf <- zinf_data %>%
  mutate(
    a = substr(a, 1, 2),
    value_1 = sqrt(value_1),
    value_2 = sqrt(value_2),
    truth_1 = sqrt(truth_1),
    truth_2 = sqrt(truth_2)
  ) %>%
  group_by(j, inference, method, zero_inf_prob, EN, a, b, D, V) %>%
  summarise(
    error = mean(sqrt((value_1 - truth_1) ^ 2 + (value_2 - truth_2) ^ 2)),
    error_bar = sd(value_1)
  )

method_cols <- c("#ae7664", "#64ae76", "#7664ae")
p <- ggplot(perf) +
  geom_abline(slope = 1, alpha = 0.6, size = 0.3) +
  geom_point(
    aes(x = error, y = error_bar, col = inference, shape = zero_inf_prob),
    size = 0.7, alpha = 0.6
  ) +
  facet_grid(method + V ~ D + EN) +
  scale_y_continuous(limits = c(0, 3)) +
  scale_x_continuous(limits = c(0, 30)) +
  scale_color_manual(values = method_cols) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 2)),
    shape = guide_legend(override.aes = list(alpha = 1, size = 2))
  ) +
  labs(x = "Root Mean Squared Error", y = "Standard Deviation (k = 1)", col = "Inference", shape = expression(p[0]))

ggsave(
  file.path(base_dir, "doc", "figure/beta_errors_nmf.png"),
  p,
  width = 5,
  height = 3
)

## ---- zinf-betas-contours ----
combined <- zinf_data %>%
  filter(P == 325)

posterior_means <- combined %>%
  group_by(j, D, a, b, zero_inf_prob, inference, method) %>%
  summarise(
    value_mean_1 = median(value_1),
    value_mean_2 = median(value_2),
    truth_1 = truth_1[1],
    truth_2 = truth_2[1]
  )

combined <- combined %>%
  filter(
    D == "D = 20",
    iteration > 400,
    sqrt(value_1) < 13,
    sqrt(value_2) < 13
  )
combined$D <- droplevels(combined$D)

plot_opts <- list(x = "sqrt(value_1)", y = "sqrt(value_2)",
                  group = "j", fill_type = "gradient", h = 1,
                  theme_opts = list(border_size = 0.7))

ggsave(
  file.path(base_dir, "doc", "figure/beta_contours_nmf_d20.png"),
  plot_contours(combined, plot_opts),
  width = 5.3,
  height = 6.5
)

combined <- zinf_data %>%
  filter(D == "D = 100")
combined <- combined %>%
  filter(
    iteration > 400,
    sqrt(value_1) < 13,
    sqrt(value_2) < 13
  )
combined$D <- droplevels(combined$D)

ggsave(
  file.path(base_dir, "doc", "figure/beta_contours_nmf_d100.png"),
  plot_contours(combined, plot_opts, 30, 30),
  width = 5.3,
  height = 6.5
)
