#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# RStan code to run LDA on antibiotics dataset.
# Based on
# https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf
# page 157

library("argparser")
parser <- arg_parser("Perform LDA on the antibiotics dataset")
parser <- add_argument(parser, "--subject", help = "Subject on which to perform analysis", default = "F")
parser <- add_argument(parser, "--K", help = "Number of LDA topics", default = 2)
argv <- parse_args(parser)

## ---- setup ----
library("rstan")
library("reshape2")
library("stringr")
library("tidyr")
library("dplyr")
library("tibble")
library("phyloseq")
library("feather")
theme_set(theme_bw(base_size = 7))
source("./posterior_check_funs.R")
dir.create("../../data/fits/", recursive = TRUE)
dir.create("../../data/figure-input/", recursive = TRUE)
dir.create("../../doc/figure/", recursive = TRUE)
set.seed(11242016)
options(mc.cores = parallel::detectCores())

## ---- get-data ----
abt <- get(load("../../data/antibiotics-study/abt.rda"))
abt <- abt %>%
  subset_samples(ind == argv$subject)

releveled_sample_data <- abt %>%
  sample_data %>%
  mutate(
    condition = recode(
      condition,
      "Pre Cp" = "Pre",
      "1st Cp" = "1st Course",
      "1st WPC" = "1st Course",
      "2nd Cp" = "2nd Course",
      "2nd WPC" = "2nd Course",
      "Post Cp" = "Post"
    )
  )

rownames(releveled_sample_data) <- abt %>%
  sample_data %>%
  rownames
sample_data(abt) <- releveled_sample_data

## ---- histograms ---
transformed_counts <- data_frame(
  count = c(get_taxa(abt), asinh(get_taxa(abt))),
  transformation = c(
    rep("original", ntaxa(abt) * nsamples(abt)),
    rep("asinh", ntaxa(abt) * nsamples(abt))
  )
)

p <- ggplot(transformed_counts) +
  geom_histogram(aes(x = count)) +
  facet_grid(. ~ transformation, scale = "free_x")
ggsave("../../doc/figure/histograms-1.png", p)

## ---- heatmaps ----
x_order <- names(sort(taxa_sums(abt)))
y_order <- names(sort(sample_sums(abt)))

## ---- lda ----
x <- t(get_taxa(abt))
out_ix <- c(160, 343, 891, 1036, 1045, 1086, 1116, 1131, 1417, 1463, 1659, 1895)

dimnames(x) <- NULL
stan_data <- list(
  K = argv$K,
  V = ncol(x[, -out_ix]),
  D = nrow(x),
  n = x[, -out_ix],
  alpha = rep(1, argv$K),
  gamma = rep(0.3, ncol(x[, -out_ix]))
)

start_fit <- Sys.time()
f <- stan_model(file = sprintf("../stan/lda_counts.stan", argv$K))
stan_fit <- vb(
  f,
  data = stan_data,
  output_samples = 1000,
  eta = 0.1,
  adapt_engaged = FALSE
)
cat(sprintf(
  "Finished in %f minutes\n",
  Sys.time() - start_fit, 4)
)

save(
  stan_fit,
  file = sprintf(
      "../../data/fits/lda-%s-%s-%s.rda",
      argv$subject,
      argv$K,
      gsub("[:|| ||-]", "", Sys.time())
  )
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)

## ---- simulate-holdout ----
thetai <- samples$theta[1,,]
x_out <- x[, out_ix]

B <- nrow(samples$theta)
lx_hat <- matrix(nrow=B, ncol=nrow(x_out))
for (i in seq_len(B)) {
  lx_hat[i, ] <- predict(lm(asinh(x_out[, 2]) ~ samples$theta[i,,]))
}

plot(asinh(x_out[, 7]), type="l")
for (i in seq_len(B)) {
  points(lx_hat[i, ], col=rgb(0, 0, 0, 0.1))
}

plot(asinh(x_out[, 7]))

## ---- posterior-checks ----
checks_data <- posterior_checks_input(
  x,
  samples$x_sim,
  sprintf("../../data/figure-input/lda-%s-%s", argv$subject, argv$K)
)
