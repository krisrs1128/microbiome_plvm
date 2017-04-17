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
theme_set(ggscaffold::min_theme())
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
p <- posterior_checks_plots(
  merged_data,
  "../../doc/figure/",
  width = 5, height = 3
)
