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
dev.new(width = 4, height = 2)

## ---- theta-reshape ----
## extract theta information from the fits
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
  filter(P == 75, N == 100)
plot_opts$facet_terms <- c("zero_inf_prob", "inference", "method")
plot_opts$group <- "i"

## ---- visualizezinfthetashist ----
mzinf_data <- melt_reshaped_samples(zinf_data)
p <- error_histograms(mzinf_data, plot_opts$facet_terms) +
  facet_grid(inference ~ zero_inf_prob + method)
ggsave(file.path(figure_dir, "visualizezinfthetashist-1.pdf"), p)
