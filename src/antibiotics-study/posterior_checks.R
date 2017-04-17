#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## While posterior checks for each method individually are contained in those
## methods' R scripts, to visualize the posterior check values of the different
## methods side-by-side we read the figure data separately and visualize here.

library("plyr")
library("dplyr")
library("ggplot2")
library("feather")
library("tidyr")
source("./posterior_check_funs.R")
theme_set(ggscaffold::min_theme(list(text_size = 10, subtitle_size = 11)))
set.seed(11242016)

## ---- read-reshape ----
input_paths <- list.files("../../data/figure-input/", full.names = TRUE)
input_data <- lapply(
  input_paths,
  read_feather
)

input_types <- data_frame(
  basename = gsub("\\.feather", "", basename(input_paths))
) %>%
  separate(basename, c("method", "data"), "-")

data_types <- unique(input_types$data)
merged_data <- list()
for (i in seq_along(data_types)) {
  cur_ix <- which(input_types$data == data_types[i])
  cur_data <- list()
  for (j in seq_along(cur_ix)) {
    cur_data[[j]] <- input_data[[cur_ix[j]]]
    cur_data[[j]]$method <- input_types$method[cur_ix[j]]
  }
  merged_data[[data_types[i]]] <- do.call(rbind, cur_data)
}

## ---- plot ----
p <- posterior_checks_plots(merged_data)

output_dir <- "../../doc/figure/"
dir.create(output_dir, recursive = TRUE)

ggsave(sprintf("%s/posterior_check_hists.png", output_dir), p[["hists"]], width = 6, height = 3.5)
ggsave(sprintf("%s/posterior_check_quantiles.png", output_dir), p[["quantiles"]], width = 3.5, height = 2)
ggsave(sprintf("%s/posterior_check_margins.png", output_dir), p[["margins"]], width = 6, height = 3.5)
ggsave(sprintf("%s/posterior_check_ts.png", output_dir), p[["ts"]], width = 6, height = 4)
ggsave(sprintf("%s/posterior_check_scores.png", output_dir), p[["scores"]], width = 7, height = 1.5)
ggsave(sprintf("%s/posterior_check_loadings.png", output_dir), p[["loadings"]], width = 6, height = 3.5)
ggsave(sprintf("%s/posterior_check_evals.png", output_dir), p[["evals"]], width = 6, height = 3.5)
