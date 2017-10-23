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

base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
unigram_dir <- file.path(base_dir, "src", "sim", "unigram")
output_path <- file.path(unigram_dir, "unigram_fits")
metadata <- read_csv(
  file.path(output_path, "metadata.csv"),
  skip = 1,
  col_names = c("file", "D", "V", "start_ix", "sigma0", "N", "?", "n_samples", "method")
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

lsamples <- quantiles
samples <- melt(
  lsamples,
  varnames = c("file", "statistic", "i", "v"),
  value.name = "mu"
)

## study the bootstrap samples
## bootstrap_paths <- metadata %>%
##   filter(method == "bootstrap") %>%
##   .[["file"]]
## rm(samples)

## boot <- feather_from_paths(bootstrap_paths)
