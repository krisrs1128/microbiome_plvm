#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulate parameters according to an LDA model, using command line arguments

args <- commandArgs(trailingOnly=TRUE)
output_dir <- args[[1]]
output_id <- args[[2]]
N <- as.integer(args[[3]])
beta_path <- args[[4]]
theta_path <- args[[5]]

## ---- libraries ----
library("feather")
library("tidyverse")
library("ldaSim")
set.seed(3141596)

## ---- simulate ----
beta <- read_feather(beta_path) %>%
  spread(k, beta) %>%
  select(-v) %>%
  as.matrix()

theta <- read_feather(theta_path) %>%
  spread(k, theta) %>%
  select(-i) %>%
  as.matrix()

n <- generate_data(N, theta, beta) %>%
  as_data_frame() %>%
  rownames_to_column("i") %>%
  gather(v, n, -i) %>%
  mutate(
    i = as.integer(i),
    v = as.integer(gsub("V", "", v))
  )

output_path <- file.path(output_dir, paste0("n-", output_id, ".feather"))
write_feather(
  as_data_frame(n),
  output_path
)

## ---- update-metadata ----
metadata <- data.frame(
  "file" = output_path,
  "D" = max(n$i),
  "V" = max(n$v),
  "N" = N,
  "K" = ncol(beta),
  "alpha0" = NA,
  "gamma0" = NA,
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
