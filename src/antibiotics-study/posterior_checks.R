#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## While posterior checks for each method individually are contained in those
## methods' R scripts, to visualize the posterior check values of the different
## methods side-by-side we read the figure data separately and visualize here.

library("argparser")
parser <- arg_parser("Perform dynamic unigrams on the antibiotics dataset")
parser <- add_argument(parser, "--subject", help = "Subject on which to perform analysis", default = "F")
argv <- parse_args(parser)

## ---- setup ----
library("plyr")
library("dplyr")
library("ggplot2")
library("scales")
library("feather")
library("tidyr")
source("./posterior_check_funs.R")
theme_set(ggscaffold::min_theme(list(text_size = 10, subtitle_size = 11)))
set.seed(11242016)

## ---- read-reshape ----
input_paths <- list.files(
  "../../data/figure-input",
  sprintf("*-%s-*", argv$subject),
  full.names = TRUE
)
input_data <- lapply(
  input_paths,
  read_feather
)

input_types <- data_frame(
  basename = gsub("\\.feather", "", basename(input_paths))
) %>%
  separate(basename, c("method", "subject", "data"), "-")

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

ggsave(sprintf("%s/posterior_check_quantiles-%s.png", output_dir, argv$subject), p[["quantiles"]], width = 4.5, height = 2.2)
ggsave(sprintf("%s/posterior_check_margins-%s.png", output_dir, argv$subject), p[["margins"]], width = 6, height = 3.5)
ggsave(sprintf("%s/posterior_check_ts-%s.png", output_dir, argv$subject), p[["ts"]], width = 6, height = 4)
ggsave(
  sprintf("%s/posterior_check_scores-loadings-%s.png", output_dir, argv$subject),
  grid.arrange(p[["scores"]], p[["loadings"]], ncol = 2),
  width = 7.3, height = 6.4
)
ggsave(sprintf("%s/posterior_check_evals-%s.png", output_dir, argv$subject), p[["evals"]], width = 6, height = 2.5)
