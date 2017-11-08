#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Script for visualizing output from LDA simulation and evaluation pipeline.
## Two main types of views: Kernel smoothed posterior views scatterplots of
## errors.
##
## author: kriss1@stanford.edu
## date: 10/23/2017

## ---- libraries-boot-expers ----
library("feather")
library("tidyverse")
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

combined <- combined %>%
  filter(iteration > 400)

## ---- beta-alignment ----
mcombined <- melt_reshaped_samples(combined)
aboot <- align_bootstraps(mcombined %>% filter(method == "bootstrap"))

mcombined <- rbind(
  align_posteriors(mcombined %>% filter(method %in% c("vb", "gibbs"))),
  aboot
)

combined <- mcombined %>%
  gather(type, value, truth, estimate) %>%
  unite(temp, type, dimension) %>%
  spread(temp, value) %>%
  ungroup() %>%
  mutate(
    method = recode(method, "gibbs" = "mcmc"),
    D = paste0("D = ", D),
    N = paste0("N = ", N),
    V = paste0("V = ", V)
  )

sort_levels <- function(x) {
  new_levels <- sort(as.numeric(gsub("[^0-9]+", "", unique(x))))
  prefix <- gsub("[0-9]+", "", x[1])
  factor(x, levels = paste0(prefix, new_levels))
}

combined$D <- sort_levels(combined$D)
combined$N <- sort_levels(combined$N)
combined$V <- sort_levels(combined$V)

## ---- beta-contours-object ----
unique_V <- unique(mcombined$V)
p <- list()
for (i in seq_along(unique_V)) {
  p[[i]] <- experiment_contours(
    combined %>%
    filter_(sprintf("V == 'V = %s'", unique_V[i]))
  )
}

## ---- betacontours1 ----
p[[1]] <- p[[1]] +
  scale_x_continuous(breaks = c(0, 0.05)) +
  labs(x = expression(sqrt(hat(beta)[1])), y = expression(sqrt(hat(beta)[2]))) +
  theme(panel.border = element_rect(fill = "transparent", size = 0.7)) +
  facet_grid(method ~ D + N)

p[[2]] <- p[[2]] +
  scale_x_continuous(breaks = c(0, 0.1)) +
  labs(x = expression(sqrt(hat(beta)[1])), y = expression(sqrt(hat(beta)[2]))) +
  theme(panel.border = element_rect(fill = "transparent", size = 0.7)) +
  facet_grid(method ~ D + N)

ggsave(
  file.path(base_dir, "doc", "figure/beta_contours_v325.png"),
  p[[1]],
  width = 5.2,
  height = 3.1,
  dpi = 450
)

ggsave(
  file.path(base_dir, "doc", "figure/beta_contours_v650.png"),
  p[[2]],
  width = 5.2,
  height = 3.1,
  dpi = 450
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
  scale_x_continuous(breaks = c(0, 0.01, 0.02)) +
  scale_color_manual(values = method_cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Root Mean Squared Error", y = "Standard Deviation (k = 1)", col = "Inference") +
  facet_grid(V ~ D + N)

ggsave(
  file.path(base_dir, "doc", "figure/beta_errors_lda.png"),
  p,
  width = 5.2,
  height = 3.1
)
