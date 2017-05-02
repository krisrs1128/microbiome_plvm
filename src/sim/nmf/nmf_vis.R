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

## ---- theta-reshape ----
## extract theta information from the fits
base_dir = "~/Desktop/microbiome_plvm"
#base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
nmf_dir <- file.path(base_dir, "src", "sim", "nmf")
fits_dir <- file.path(nmf_dir, "fits")
figure_dir <- file.path(base_dir, "doc", "figure")
dir.create(figure_dir, recursive = TRUE)

fits <- list.files(fits_dir, "p_10-*", full.names = TRUE)
expers <- fromJSON(
  file.path(nmf_dir, "config.json"),
  simplifyVector = FALSE
)


## Temporary realignment, to avoid having to rereun all the simulations
for (i in seq_along(fits)) {
  cur_samples <- get(load(fits[[i]]))
  if (dim(cur_samples$theta)[2] == 2) {
    cur_samples$theta <- aperm(cur_samples$theta, c(1, 3, 2))
    cur_samples$beta <- aperm(cur_samples$beta, c(1, 3, 2))
    save(cur_samples, file = fits[[i]])
  }
}

theta_fits <- reshape_all_samples(
  fits,
  file.path(nmf_dir, "config.json"),
  "theta",
  c("i", "k")
)
theta_fits$method <- basename(as.character(theta_fits$method))
theta_fits$method <- theta_fits$method %>%
  revalue(
    c(
    "nmf_gamma_poisson.stan" = "GaP",
    "nmf_gamma_poisson_zero.stan" = "Z-GaP"
    )
  )

## ---- visualizethetas-prep -----
## Visualize the fitted thetas, according to a few different simulation properties
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

## ---- visualize-zinf-thetas-prep ----
zinf_data <- theta_fits %>%
  filter(P == 10, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
plot_opts$group <- "i"

## ---- visualizezinfthetashist ----
zinf_data$a0 <- as.character(zinf_data$a0)
zinf_data$zero_inf_prob <- as.character(zinf_data$zero_inf_prob)

perf <- zinf_data %>%
  mutate(value_1 = sqrt(value_1), truth_1 = sqrt(truth_1), value_2 = sqrt(value_2), truth_2 = sqrt(truth_2)) %>%
  group_by(i, inference, method, zero_inf_prob, a0) %>%
  summarise(
    error = mean(sqrt((value_1 - truth_1) ^ 2 + (value_2 - truth_2) ^ 2)),
    error_bar = sqrt(det(cov(cbind(value_1, value_2))))
  )

ggplot(perf) +
  geom_point(aes(x = error, y = error_bar)) +
  facet_grid(inference ~ method + zero_inf_prob) +
  geom_abline(slope = 1) +
  coord_fixed()

x = zinf_data %>%
  filter(i == "6", inference == "gibbs", method == "GaP", zero_inf_prob == 0, a0 == "1") %>%
  select(value_1, value_2)
y = zinf_data %>%
  filter(i == "6", inference == "gibbs", method == "GaP", zero_inf_prob == 0, a0 == "0.2") %>%
  select(value_1, value_2)

ggplot(x) +
  geom_point(aes(x = value_1, y = value_2))
sqrt(det(cov(as.matrix(x))))
(sd(x$value_1) + sd(x$value_2) ) / 2
mean(sqrt(diag(cov(as.matrix(x)))))
sqrt(det(cov(as.matrix(x))))

head(zinf_data)
ggplot() +
  geom_text(
    data = zinf_data %>% sample_n(10000),
    aes(x = sqrt(value_1), y = sqrt(value_2), label = i), alpha = 0.2, size = 1) +
  geom_text(
    data = zinf_data %>%
      group_by(i, inference, method, zero_inf_prob, a) %>%
      summarise(truth_1 = truth_1[1], truth_2 = truth_2[1]),
    aes(x = sqrt(truth_1), y = sqrt(truth_2), label = i), alpha = 1, size = 3) +
  facet_grid(inference ~ zero_inf_prob + method + a)

ggplot() +
  geom_text(
    data = zinf_data %>%
      group_by(i, inference, method, zero_inf_prob, a) %>%
      summarise(value_1 = mean(value_1), value_2 = mean(value_2)),
    aes(x = sqrt(value_1), y = sqrt(value_2), label = i), alpha = 1, size = 3, col = "red") +
  geom_text(
    data = zinf_data %>%
      group_by(i, inference, method, zero_inf_prob, a) %>%
      summarise(truth_1 = truth_1[1], truth_2 = truth_2[1]),
    aes(x = sqrt(truth_1), y = sqrt(truth_2), label = i), alpha = 1, size = 3) +
  facet_grid(inference ~ zero_inf_prob + method + a)

mzinf_data <- melt_reshaped_samples(zinf_data)

p <- error_histograms(mzinf_data, plot_opts$facet_terms) +
  facet_grid(inference ~ zero_inf_prob + method)
ggsave(file.path(figure_dir, "visualizezinfthetashist-1.pdf"), p)
