#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Script for visualizing output from LDA simulation and evaluation pipeline.
## Three main types of views: Boxplots of proportions estimates, across
## configurations, scatterplots of pairs of proportions estimates, and
## histograms of errors.
##
## author: kriss1@stanford.edu

## ---- libraries-boot-expers ----
library("feather")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("rstan")
library("ggplot2")
library("ggscaffold")
library("ldaSim")
base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
lda_dir <- file.path(base_dir, "src", "sim", "lda")

## ---- paths ----
output_path <- file.path(lda_dir, "fits")
metadata <- fread(file.path(output_path, "metadata.csv")) %>%
  unique()

## ---- beta-samples ----
beta <- get_truth_data(metadata, "beta") %>%
  rename(variable = v)
combined <- get_samples(metadata, "beta", c("iteration", "k", "variable")) %>%
  full_join(get_bootstraps(metadata, "beta")) %>%
  left_join(beta) %>%
  setDT()

## ---- beta-alignment ----
mcombined <- melt_reshaped_samples(combined)
mcombined <- rbind(
  align_posteriors(mcombined %>% filter(method %in% c("vb", "gibbs"))),
  align_bootstraps(mcombined %>% filter(method == "bootstrap"))
)

combined <- mcombined %>%
  gather(type, value, truth, estimate) %>%
  unite(temp, type, dimension) %>%
  spread(temp, value)

## ---- beta-contours-object ----
unique_V <- unique(mcombined$V)
p <- list()
for (i in seq_along(unique_V)) {
  p[[i]] <- experiment_contours(
    combined %>%
    filter_(sprintf("V == %s", unique_V[i]))
  )
}

## ---- betacontours1 ----
ggsave(file.path(base_dir, "doc", "figure/betacontours1-1.pdf", p[[1]])

## ---- betacontours2 ----
ggsave(file.path(base_dir, "doc", "figure/betacontours2-1.pdf", p[[2]])

## ---- betahistograms ----
p <- error_histograms(mcombined, c("V + method", "D + N"))
ggsave(file.path(base_dir, "doc", "figure/betahistograms-1.pdf", p)

## ---- beta-boxplots-object ----
p <- list()
for (i in seq_along(unique_V)) {
  p[[i]] <- experiment_boxplots(
    mcombined %>%
    filter_(sprintf("V == %s", unique_V[i]))
  )
}

## ---- betaboxplot1 ----
ggsave(file.path(base_dir, "doc", "figure/betaboxplot1-1.pdf", p[[1]])

## ---- betaboxplot2 ----
ggsave(file.path(base_dir, "doc", "figure/betaboxplot2-1.pdf", p[[2]])

## ---- theta-samples ----
theta <- get_truth_data(metadata, "theta", "i") %>%
  rename(variable = i)
combined <- get_samples(metadata, "theta", c("iteration", "variable", "k")) %>%
  full_join(get_bootstraps(metadata, "theta", "i")) %>%
  left_join(theta)

## ---- theta-alignment ----
mcombined <- melt_reshaped_samples(combined)
mcombined <- rbind(
  align_posteriors(mcombined %>% filter(method %in% c("vb", "gibbs"))),
  align_bootstraps(mcombined %>% filter(method == "bootstrap"))
)

combined <- mcombined %>%
  gather(type, value, truth, estimate) %>%
  unite(temp, type, dimension) %>%
  spread(temp, value)

## ---- theta-boxplot-object ----
unique_D <- unique(mcombined$D)
p <- list()
for (i in seq_along(unique_D)) {
  p[[i]] <- experiment_boxplots(
    mcombined %>%
    filter_(sprintf("D == %s", unique_D[i]), "dimension == 1")
  ) +
    facet_grid(
      N + V ~ variable,
      scale = "free_x"
    )
}

## ---- thetaboxplot1 ----
ggsave(file.path(base_dir, "doc", "figure/thetaboxplot1-1.pdf", p[[1]])

## ---- thetaboxplot2 ----
ggsave(file.path(base_dir, "doc", "figure/thetaboxplot2-1.pdf", p[[2]])

## ---- thetahistograms ----
p <- error_histograms(mcombined, c("D + method", "V + N"))
ggsave(file.path(base_dir, "doc", "figure/thetahistograms-1.pdf", p)
