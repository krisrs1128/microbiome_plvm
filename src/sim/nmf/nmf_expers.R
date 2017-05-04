#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiments with NMF, under different parameterizations and fitting
## techniques. Here, we consider the parameterizations in the parameterizations
## table and then save the resulting fits in a fits/ directory. These are
## visualized separately.
##
## author: kriss1@stanford.edu

## ---- libraries-nmf-expers ----
library("jsonlite")
library("nmfSim")

## ---- configuration ----
## create the configuration JSON file
#base_dir <- Sys.getenv("MICROBIOME_PLVM_DIR")
base_dir <- "~/Desktop/microbiome_plvm/"
nmf_dir <- file.path(base_dir, "src", "sim", "nmf")
config_path <- file.path(nmf_dir, "config.json")
stan_path <- file.path(.libPaths()[1], "nmfSim", "extdata")
fits_dir = file.path(nmf_dir, "fits")
dir.create(fits_dir, recursive = TRUE)

sim_factors <- list(
  "N" = c(20, 100),
  "P" = 10,
  "b" = c(1, 4/7, 1/3),
  "a" = c(1, 10/7, 5/3),
  "zero_inf_prob" = c(0, 0.2)
)

sim_factors_high <- sim_factors
sim_factors_high$P <- 50
sim_factors_high$a <- c(5/3, 4/3, 1)
sim_factors_high$b <- c(1/3, 2/3, 1)

model_factors <- list(
  "inference" = c("gibbs", "vb", "bootstrap"),
  "method" = c(
    file.path(stan_path, "nmf_gamma_poisson.stan"),
    file.path(stan_path, "nmf_gamma_poisson_zero.stan")
  )
)

config_df <- rbind(
  expand.grid(c(sim_factors, model_factors)),
  expand.grid(c(sim_factors_high, model_factors))
)

write_configs(
  config_df,
  n_batches = 3,
  config_path = config_path,
  output_dir = fits_dir
)

## ---- submit-jobs ----
## loop over unique values in the "batch" field of the json file
configs <- fromJSON(config_path, simplifyVector = FALSE)
batches <- sapply(configs, function(x) { x$batch })

rscript_file <- file.path(nmf_dir, "nmf_script.R")
for (i in seq_along(unique(batches))) {
  rscript_cmd <- paste("Rscript", rscript_file, config_path, i)
  system(paste(rscript_cmd, "&"))
}

sims_complete <- FALSE
while(!sims_complete) {
  Sys.sleep(60)
  completed <- list.files(fits_dir)
  cat(sprintf("Completed fits %s \n", paste0(completed, collapse = "\t")))
  if (length(completed) == length(configs)) {
    sims_complete <- TRUE
  }
}
