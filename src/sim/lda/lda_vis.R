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
library("readr")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ldaSim")

## ---- paths ----
base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
lda_dir <- file.path(base_dir, "src", "sim", "lda")
output_path <- file.path(lda_dir, "fits")
metadata <- read_csv(file.path(output_path, "metadata.csv")) %>%
  unique() %>%
  mutate(file = file.path(lda_dir, "pipeline", file))

## ---- beta-samples ----
beta <- get_truth_data(metadata, "beta") %>%
  rename(variable = v)
combined <- get_samples(metadata, "beta", c("iteration", "k", "variable")) %>%
  full_join(get_bootstraps(metadata, "beta")) %>%
  left_join(beta)

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
p[[1]] <- p[[1]] +
  labs(x = expression(sqrt(hat(beta)[1])), y = expression(sqrt(hat(beta)[2]))) +
  facet_grid(N + D ~ method)

ggsave(
  file.path(base_dir, "doc", "figure/betacontours1-1.pdf"),
  p[[1]],
  width = 5,
  height = 3
)
