#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Simulate parameters according to the dynamic unigrams model, using command
## line arguments.
##
## author: sankaran.kris@gmail.com
## date: 09/19/2017

###############################################################################
## setup arguments + functions
###############################################################################
library("reshape2")
library("feather")
args <- commandArgs(trailingOnly=TRUE)
argv <- list()
argv$output_dir <- args[[1]]
argv$gen_id <- args[[2]]
argv$D <- as.integer(args[[3]])
argv$V <- as.integer(args[[4]])
argv$sigma0 <- as.numeric(args[[5]])

#' Simulate unigram parameters
unigram_params <- function(D, V, sigma0) {
  mu <- matrix(nrow = D, ncol = V)
  mu[1, ] <- rnorm(V, 0, sigma0)
  for (i in seq_len(D - 1)) {
    mu[i + 1, ] <- rnorm(V, mu[i, ], sigma0)
  }

  mu
}

###############################################################################
## simulate and write results to file
###############################################################################
mu <- unigram_params(argv$D, argv$V, argv$sigma0)
dir.create(argv$output_dir)
mu_path <- file.path(argv$output_dir, sprintf("mu-%s.feather", argv$gen_id))

write_feather(
  melt(
    mu,
    varnames = c("i", "v"),
    value.name = "mu"
  ),
  mu_path
)

###############################################################################
## update the metaata
###############################################################################
metadata <- data.frame(
  "file" = mu_path,
  "D" = argv$D,
  "V" = argv$V,
  "N" = NA,
  "sigma0" = argv$sigma0,
  "sigma_fit" = NA,
  "n_samples" = NA,
  "method" = NA,
  "iteration" = NA
)

write.table(
  metadata,
  file = file.path(argv$output_dir, "metadata.csv"),
  append = TRUE,
  sep = ",",
  row.names = FALSE,
  col.names = !file.exists(file.path(argv$output_dir, "metadata.csv"))
)
