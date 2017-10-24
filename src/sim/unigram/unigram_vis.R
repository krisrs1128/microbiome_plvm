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
  fit_params <- rstan::extract(fit)
  mu_i <- fit_params$mu
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

ggplot(samples) +
  geom_pointrange(
    aes(
      x = v, y = `50%`, ymin = `25%`, ymax = `75%`, col = N
    ),
    size = 0.5, alpha = 0.3, fatten = 0.1
  ) +
  facet_wrap(D ~ V, scale = "free_x") +
  ylim(-35, 35) +
  coord_fixed()

combined <- mu %>%
  select(D, V, i, v, mu) %>%
  full_join(
    samples %>%
    select(method, D, V, i, v, `25%`, `50%`, `75%`)
  )

ggplot(combined) +
  geom_point(
    aes(
      x = mu,
      y = `50%`
    ),
    alpha = 0.3,
    size = 0.5
  ) +
  facet_grid(D ~ V + method) +
  ylim(-20, 20)


## study the bootstrap samples
## bootstrap_paths <- metadata %>%
##   filter(method == "bootstrap") %>%
##   .[["file"]]
## rm(samples)

lbootstraps <- list()
for (i in seq_along(bootstrap_paths)) {
  cat(sprintf(
    "Processing sample %s [%s / %s] \n",
    samples_paths[i],
    i,
    length(samples_paths)
  ))

}

## boot <- feather_from_paths(bootstrap_paths)
curfit <- get(load("../unigram_fits/vb-042de8ce9ed1b5e9cfadafac738ec71f.RData"))

x <- read_feather("../unigram_fits/mu-05e91549a659b808c8bf3bae6b402c8c.feather")
