#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Utilities used throughout the unigram simulation pipeline.
##
## author: sankaran.kris@gmail.com
## date: 09/19/2017

softmax <- function(x) {
  exp(x) / log(sum(exp(x)))
}

#' Simulate from a unigram model
unigram_data <- function(N, mu_mat) {
  X <- matrix(nrow = nrow(mu_mat), ncol = ncol(mu_mat))

  for (i in seq_len(nrow(mu_mat))) {
    X[i, ] <- rmultinom(1, N, prob = softmax(mu_mat[i, ]))[, 1]
  }

  X
}
