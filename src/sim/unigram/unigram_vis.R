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
## Libraries and setup
###############################################################################
library("feather")
library("tidyverse")
library("ldaSim")
library("data.table")
theme_set(ggscaffold::min_theme(list(border_size = 0.7)))

base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
unigram_dir <- file.path(base_dir, "src", "sim", "unigram")
output_path <- file.path(unigram_dir, "unigram_fits")
metadata <- read_csv(
  file.path(output_path, "metadata.csv"),
  skip = 1,
  col_names = c("file", "D", "V", "N", "sigma0", "a0", "b0", "n_samples", "method", "iteration")
) %>%
  unique()

###############################################################################
## Read in all the data
###############################################################################
## get true underlying parameters
truth_paths <- metadata %>%
  filter(is.na(method), grepl("mu", file)) %>%
  select(file) %>%
  unlist()

mu <- feather_from_paths(truth_paths) %>%
  left_join(metadata)

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
  lsamples[[i]] <- apply(mu_i, c(2, 3), quantile, c(0.25, 0.5, 0.75)) %>%
    melt(
      varnames = c("quantile", "i", "v"),
      value.name = "mu"
    ) %>%
    as_data_frame()
  lsamples[[i]]$file <- samples_paths[i]
}

samples <- bind_rows(lsamples)
samples <- samples %>%
  left_join(metadata) %>%
  spread(quantile, mu)

## extract the bootstrap samples
bootstrap_paths <- metadata %>%
  filter(method == "bootstrap") %>%
  .[["file"]]

bootstrap_paths <- bootstrap_paths[1:5]
bootstraps <- feather_from_paths(bootstrap_paths)  %>%
  left_join(metadata) %>%
  group_by(i, v, D, V, N, sigma0, a0, b0, method) %>%
  do(
    data.frame(
      quantile = c("25%", "50%", "75%"),
      mu = quantile(.$mu, c(0.25, 0.5, 0.75))
    )
  ) %>%
  spread(quantile, mu)

combined_samples <- mu %>%
  select(D, V, i, v, mu) %>%
  full_join(
    samples %>%
    select(method, N, D, V, i, v, `25%`, `50%`, `75%`)
  )
combined_bootstrap <- mu %>%
  select(D, V, i, v, mu) %>%
  full_join(
    bootstraps %>%
    select(method, N, D, V, i, v, `25%`, `50%`, `75%`)
  )
combined <- bind_rows(
  combined_samples,
  combined_bootstrap
)

###############################################################################
## Visualize the performance of different methods
###############################################################################
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
  geom_linerange(
    aes(
      x = mu,
      ymin = `25%`,
      ymax = `75%`,
      col = method
    ),
    alpha = 0.2,
    size = 0.05,
  ) +
  scale_color_manual(values = method_cols) +
  guides(col = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  coord_flip() +
  coord_fixed() +
  facet_grid(method + V ~ D + N) +
  scale_x_continuous(limits = c(-3, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3, 20), expand = c(0, 0)) +
  labs(
    "x" = bquote(True~mu[tv]),
    "y" = bquote(Posterior~hat(mu[tv]))
  )

ggsave(
  file.path(base_dir, "doc", "figure/mu_intervals.png")
)

## summary performance plot
perf <- combined %>%
  group_by(v, method, D, V, N) %>%
  summarise(
    error = sqrt(mean((mu - `50%`) ^ 2)),
    error_bar = sd(`50%`)
  )

ggplot(perf) +
  geom_abline(slope = 1, alpha = 0.6, size = 0.3) +
  geom_point(aes(x = error, y = error_bar, col = method), size = 0.7, alpha = 0.6) +
  scale_color_manual(values = method_cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Root Mean Squared Error", y = "Standard Deviation", col = "Inference") +
  facet_grid(V ~ D + N) +
  scale_x_continuous(limits = c(0, 19), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 13), expand = c(0, 0))

ggsave(
  file.path(base_dir, "doc", "figure/mu_errors_unigram.png"),
  width = 5.2,
  height = 3.1
)
