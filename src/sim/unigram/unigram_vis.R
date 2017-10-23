#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Visualization of output from the unigram simulation and evaluation pipeline.
## Generally comparable ot the LDA visualization script, but there are fewer
## parameters of interest.
##
## author: kriss1@stanford.edu
## date: 10/23/2017

###############################################################################
## Libraries and reading in data
###############################################################################
library("feather")
library("tidyverse")
library("ldaSim")
library("data.table")

base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
unigram_dir <- file.path(base_dir, "src", "sim", "unigram")
output_path <- file.path(unigram_dir, "unigram_fits")
metadata <- read_csv(file.path(output_path, "metadata.csv")) %>%
  unique() %>%
  mutate(file = file.path(unigram_dir, "pipeline", file))


## get true underlying parameters
truth_paths <- metadata %>%
  filter(is.na(method), grepl("mu", file)) %>%
  select(file) %>%
  unlist()

mu <- feather_from_paths(truth_paths) %>%
  dcast(file + i ~ v, value = "mu")

## read in gibbs and vb samples
samples_paths <- metadata %>%
  filter(method %in% c("vb", "gibbs")) %>%
  select(file) %>%
  unlist()
