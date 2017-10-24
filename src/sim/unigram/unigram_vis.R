#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Visualization of output from the unigram simulation and evaluation pipeline.
## Generally comparable ot the LDA visualization script, but there are fewer
## parameters of interest.
##
## author: kriss1@stanford.edu
## date: 10/23/2017

###############################################################################
## Libraries and reading in data
###############################################################################
library("feather")
library("tidyverse")
library("ldaSim")
library("data.table")
library("abind")
theme_set(ggscaffold::min_theme(list(border_size = .7)))

base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
unigram_dir <- file.path(base_dir, "src", "sim", "unigram")
output_path <- file.path(unigram_dir, "unigram_fits")
metadata <- read_csv(
  file.path(output_path, "metadata.csv"),
  skip = 1,
  col_names = c("file", "D", "V", "N", "sigma0", "a0", "b0", "n_samples", "method", "iteration")
) %>%
  unique()

## get true underlying parameters
truth_paths <- metadata %>%
  filter(is.na(method), grepl("mu", file)) %>%
  select(file) %>%
  unlist()

mu <- feather_from_paths(truth_paths) %>%
  left_join(metadata)

ggplot(mu) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_point(
    aes(
      x = v, y = mu, group = i
    ),
    alpha = 0.2,
    size = 0.3
  ) +
  facet_grid(D ~ V, scales = "free")

## read in gibbs and vb samples
samples_paths <- metadata %>%
  filter(method %in% c("vb", "gibbs")) %>%
  select(file) %>%
  unlist()

lsamples <- list()

for (i in seq_along(samples_paths)) {
  cat(sprintf(
    "Processing sample %s [%s / %s] \n",
    samples_paths[i],
    i,
    length(samples_paths)
  ))
  fit <- get(load(samples_paths[i]))
  mu_i <- rstan::extract(fit)$mu
  lsamples[[i]] <- abind(
    apply(mu_i, c(2, 3), quantile),
    "mean" = array(colMeans(mu_i), c(1, dim(mu_i)[2:3])),
    along = 1
  )
}

names(lsamples) <- samples_paths
samples <- melt(lsamples)
colnames(samples) <- c("statistic", "i", "v", "mu", "file")
samples <- samples %>%
  left_join(metadata) %>%
  spread(statistic, mu)

combined <- mu %>%
  select(D, V, i, v, mu) %>%
  full_join(
    samples %>%
    select(method, N, D, V, i, v, `25%`, `50%`, `75%`)
  )

method_cols <- c("#ae7664", "#64ae76", "#7664ae")
ggplot(combined) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(
    slope = 1,
    intercept = 0,
    size = 0.6,
    alpha = 0.6
  ) +
  geom_pointrange(
    aes(
      x = mu,
      y = `50%`,
      ymin = `25%`,
      ymax = `75%`,
      col = method
    ),
    alpha = 0.01,
    size = 0.05,
    fatten = 0.01
  ) +
  scale_color_manual(values = method_cols) +
  coord_flip() +
  coord_fixed() +
  facet_grid(method + V ~ D + N) +
  ylim(-5, 20) +
  xlim(-5, 20)

## summary performance plot
perf <- combined %>%
  group_by(v, method, D, V, N) %>%
  summarise(
    error = mean(mu - `50%`),
    error_bar = sd(`50%`)
  )

ggplot(perf) +
  geom_abline(slope = 1, alpha = 0.6, size = 0.3) +
  geom_point(aes(x = error, y = error_bar, col = method), size = 0.7, alpha = 0.6) +
  scale_color_manual(values = method_cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Root Mean Squared Error", y = "Standard Deviation", col = "Inference") +
  facet_grid(V ~ D + N) +
  xlim(0, 25)

## study the bootstrap samples
bootstrap_paths <- metadata %>%
  filter(method == "bootstrap") %>%
  .[["file"]]

lbootstraps <- list()
for (i in seq_along(bootstrap_paths)) {
  cat(sprintf(
    "Processing bootstrap %s [%s / %s] \n",
    bootstrap_paths[i],
    i,
    length(bootstrap_paths)
  ))
  lbootstraps[[i]] <- read_feather(bootstrap_paths[i])
  lbootstraps[[i]]$file <- bootstrap_paths[i]
}

bootstraps <- bind_rows(lbootstraps) %>%
  left_join(metadata)
