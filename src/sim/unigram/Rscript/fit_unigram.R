#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Fit a dynamic unigram model using STAN.
##
## author: sankaran.kris@gmail.com
## date: 09/19/2017

###############################################################################
## arguments and libraries
###############################################################################
library("feather")
library("tidyverse")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[[1]]
fit_method <- args[[2]]
stan_path <- args[[3]]
fit_id <- args[[4]]
data_path <- args[[5]]
n_samples <- as.integer(args[[6]])
a0 <- as.numeric(args[[7]])
b0 <- as.numeric(args[[8]])

###############################################################################
## prepare model input
###############################################################################
X <- read_feather(data_path) %>%
  spread(v, x) %>%
  select(-i) %>%
  as.matrix()

stan_data <- list(
  "V" = ncol(X),
  "T" = nrow(X),
  "N" = nrow(X),
  "a0" = a0,
  "b0" = b0,
  "times" = seq_len(nrow(X)),
  "times_mapping" = seq_len(nrow(X)),
  "x" = X
)

###############################################################################
## fit the model
###############################################################################
if (tolower(fit_method) == "vb") {
  fit <- vb(
    stan_model(stan_path),
    data = stan_data,
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

###############################################################################
## save results
###############################################################################
output_path <- file.path(
  output_dir,
  sprintf("%s-%s.RData", fit_method, fit_id)
)
save(fit, file = output_path)

###############################################################################
##  update the metdata
###############################################################################
metadata <- data.frame(
  "file" = output_path,
  "D" = nrow(X),
  "V" = ncol(X),
  "N" = sum(X[, 1]),
  "sigma0" = NA,
  "a0" = a0,
  "b0" = b0,
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
