#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Simulate data according to the dynamic unigrams model, given simulated
## parameters $\mu_{t}$.
##
## author: sankaran.kris@gmail.com
## date: 09/19/2017

###############################################################################
## setup arguments + functions
###############################################################################
library("feather")
library("tidyverse")
library("reshape2")
source("../Rscript/unigram_utils.R") # assumed running from pipeline dir

args <- commandArgs(trailingOnly = TRUE)
argv <- list()
argv$output_dir <- args[[1]]
argv$gen_id <- args[[2]]
argv$N <- as.numeric(args[[3]])
argv$mu_path <- args[[4]]

###############################################################################
## Simulate and write results to file
###############################################################################
mu_mat <- read_feather(argv$mu_path) %>%
  spread(v, mu) %>%
  select(-i)
X <- unigram_data(argv$N, mu_mat)

dir.create(argv$output_dir)
x_path <- file.path(argv$output_dir, sprintf("x-%s.feather", argv$gen_id))

write_feather(
  melt(
    X,
    varnames = c("i", "v"),
    value.name = "x"
  ),
  x_path
)
