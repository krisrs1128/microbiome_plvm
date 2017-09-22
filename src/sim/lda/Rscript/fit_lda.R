#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Fit an LDA Model using STAN.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[[1]]
gen_id <- args[[2]]
data_path <- args[[3]]
fit_method <- args[[4]]
n_samples <- as.integer(args[[5]])
K <- as.integer(args[[6]])
alpha <- as.numeric(args[[7]])
gamma <- as.numeric(args[[8]])

## ---- libraries ----
library("feather")
library("dplyr")
library("rstan")
library("reshape2")
library("ldaSim")

## ---- get-data ----
n <- read_feather(data_path) %>%
    dcast(i ~ v) %>%
    select(-i) %>%
    as.matrix()

stan_data <- list(
    "K" = K,
    "V" = ncol(n),
    "D" = nrow(n),
    "n" = n,
    "alpha" = rep(alpha, K),
    "gamma" = rep(gamma, ncol(n))
)

## ---- fit-model ----
stan_path <- file.path(.libPaths()[1], "ldaSim", "extdata", "lda.stan")
if (tolower(fit_method) == "vb") {
    fit <- vb(
      stan_model(stan_path),
      data = stan_data,
      iter = 1000,
      output_samples = n_samples
    )
} else if (tolower(fit_method) == "gibbs") {
    fit <- stan(
      stan_path,
      data = stan_data,
      chains = 1,
      warmup = 1000,
      iter = 1000 + n_samples
    )
} else {
    stop("fit_method must be either 'gibbs' or 'vb'")
}

## ---- save ----
output_path <- file.path(
  output_dir,
  paste0(fit_method, "-", gen_id, ".RData")
)
save(fit, file = output_path)

## ---- update-metadata ----
metadata <- data.frame(
  "file" = output_path,
  "D" = nrow(n),
  "V" = ncol(n),
  "N" = sum(n) / nrow(n),
  "K" = K,
  "alpha0" = NA,
  "gamma0" = NA,
  "alpha_fit" = alpha,
  "gamma_fit" = gamma,
  "n_samples" = n_samples,
  "method" = fit_method,
  "iteration" = NA
)

metadata_path <- file.path(output_dir, "metadata.csv")
write.table(
  metadata,
  file = metadata_path,
  append = TRUE,
  sep = ",",
  row.names = FALSE,
  col.names = !file.exists(metadata_path)
)
