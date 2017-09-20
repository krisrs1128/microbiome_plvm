#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Simulate parametric bootstrap samples from a fitted dynamic unigram model.
##
## author: sankaran.kris@gmail.com
## date: 09/20/2017


###############################################################################
## Libraries and arguments setup
###############################################################################
library("rstan")
library("feather")
library("tidyverse")
library("reshape2")
source("../src/unigram_utils.R") # assumed running from pipeline dir
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(3141596)

args <- commandArgs(trailingOnly = TRUE)
argv <- list()
argv$output_dir <- args[[1]]
argv$start_iter <- as.integer(args[[2]])
argv$end_iter <- as.integer(args[[3]])
argv$fit_id <- args[[4]]
argv$input_path <- args[[5]]
argv$stan_path <- args[[6]]
argv$n_samples <- as.integer(args[[10]])
argv$N <- args[[7]]
argv$a0 <- as.numeric(args[[8]])
argv$b0 <- as.numeric(args[[9]])
dir.create(argv$output_dir)

###############################################################################
## Load fitted model and create replicates
###############################################################################
fit <- get(load(argv$input_path))
samples <- rstan::extract(fit)
mu_hat <- apply(samples$mu, c(2, 3), mean)

for (i in seq(argv$start_iter, argv$end_iter)) {
  cur_data <- unigram_data(argv$N, mu_hat)
  stan_data <- list(
    "V" = ncol(cur_data),
    "T" = nrow(cur_data),
    "N" = nrow(cur_data),
    "a0" = argv$a0,
    "b0" = argv$b0,
    "times" = seq_len(nrow(cur_data)),
    "times_mapping" = seq_len(nrow(cur_data)),
    "x" = cur_data
  )

  vb_fit <- vb(
    stan_model(argv$stan_path),
    stan_data,
    output_samples = argv$n_samples
  )
  mu_boot <- apply(rstan::extract(vb_fit)$mu, c(2, 3), mean) %>%
    melt(varnames = c("i", "v"), value.name = "mu")
  mu_path <- file.path(
    argv$output_dir,
    sprintf("mu-%s%s.feather", argv$fit_id, i)
  )
  write_feather(mu_boot, mu_path)

  metadata <- data.frame(
    "file" = mu_path,
    "D" = nrow(cur_data),
    "V" = ncol(cur_data),
    "N" = sum(cur_data[, 1]),
    "sigma0" = NA,
    "a0" = argv$a0,
    "b0" = argv$b0,
    "n_samples" = argv$n_samples,
    "method" = "bootstrap",
    "iteration" = i
  )

  metadata_path <- file.path(argv$output_dir, "..", "metadata.csv")
  write.table(
    metadata,
    file = metadata_path,
    append = TRUE,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(metadata_path)
  )
}
