#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulate parameters according to an LDA model, using command line arguments
args <- commandArgs(trailingOnly=TRUE)
output_dir <- args[[1]]
output_id <- args[[2]]
D <- as.integer(args[[3]])
V <- as.integer(args[[4]])
K <- as.integer(args[[5]])
alpha0 <- as.numeric(args[[6]])
gamma0 <- as.numeric(args[[7]])

## ---- libraries ----
library("feather")
library("reshape2")
library("ldaSim")
set.seed(3141596)

## ---- simulate ----
params <- generate_params(D, K, V, alpha0, gamma0)
dir.create(output_dir)
beta_path <- file.path(output_dir, paste0("beta-", output_id, ".feather"))
write_feather(
  melt(
    params$beta,
    varnames = c("v", "k"),
    value.name = "beta"
  ),
  beta_path
)

theta_path <- file.path(output_dir, paste0("theta-", output_id, ".feather"))
write_feather(
  melt(
    params$theta,
    varnames = c("i", "k"),
    value.name = "theta"
  ),
  theta_path
)

## ---- update-metadata ----
metadata <- data.frame(
  "file" = c(beta_path, theta_path),
  "D" = D,
  "V" = V,
  "N" = NA,
  "K" = K,
  "alpha0" = alpha0,
  "gamma0" = gamma0,
  "alpha_fit" = NA,
  "gamma_fit" = NA,
  "n_samples" = NA,
  "method" = NA,
  "iteration" = NA
)

write.table(
  metadata,
  file = file.path(output_dir, "metadata.csv"),
  append = TRUE,
  sep = ",",
  row.names = FALSE,
  col.names = !file.exists(file.path(output_dir, "metadata.csv"))
)
