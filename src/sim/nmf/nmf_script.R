#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## We don't want to run all the experiment configurations in parallel or series,
## we need something intermediate. This script runs a subset of the experiments
## and can be called by a SLURM submitter.
##
## author: kriss1@stanford.edu

## ---- libraries ----
library("jsonlite")
library("nmfSim")

## ---- parse-args ----
args <- commandArgs(trailingOnly = TRUE)
expers <- fromJSON(args[[1]], simplifyVector = TRUE, simplifyDataFrame = FALSE)
subset_ix <- as.integer(args[[2]])

for (i in seq_along(expers)) {
  if (expers[[i]]$batch != subset_ix) next
  set.seed(01112017)
  output_path <- file.path(expers[[i]]$output_dir, paste0(expers[[i]]$id, ".rda"))
  if (file.exists(output_path)) {
    warning(paste("File", output_path, "exists, not overwriting."))
    next
  }

  cur_data <- nmf_sim(expers[[i]]$sim_opts)
  cur_fit <- fit_model(cur_data$y, expers[[i]]$model_opts)
  save(cur_fit, file = output_path)
}
