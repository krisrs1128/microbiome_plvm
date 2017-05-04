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
#base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
base_dir <- "~/Desktop/microbiome_plvm"
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
  facet_grid(method ~ D + N)

ggsave(
  file.path(base_dir, "doc", "figure/betacontours1-1.pdf"),
  p[[1]],
  width = 5,
  height = 3
)

## ---- relative-errors ----
perf <- combined %>%
  mutate(estimate_1 = sqrt(estimate_1), truth_1 = sqrt(truth_1), estimate_2 = sqrt(estimate_2), truth_2 = sqrt(truth_2)) %>%
  group_by(variable, method, D, V, N) %>%
  summarise(
    error = mean(sqrt((estimate_1 - truth_1) ^ 2 + (estimate_2 - truth_2) ^ 2)),
    error_bar = sd(estimate_1)
  )

theme_set(ggscaffold::min_theme(list(border_size = .7)))
method_cols <- c("#ae7664", "#64ae76", "#7664ae")
p <- ggplot(perf) +
  geom_abline(slope = 1, alpha = 0.6, size = 0.3) +
  geom_point(aes(x = error, y = error_bar, col = method), size = 0.7, alpha = 0.6) +
  scale_color_manual(values = method_cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Error", y = "SD (k = 1)", col = "Inference") +
  facet_grid(V ~ D + N)

ggsave(
  file.path(base_dir, "doc", "figure/beta_errors_lda.pdf"),
  p,
  width = 5,
  height = 3
)
