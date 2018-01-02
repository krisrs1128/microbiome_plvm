#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Predictive likelihood on Train / Test split
##
## The idea of this approach is to hold out some of the species when fitting the
## LDA model. We can then try to estimate the predictive likelihood on the held
## out data set, p(x^{*} | x) by drawing samples from the approximate posterior
## and averaging the log likelihood of the observed samples given those drawn
## parameter values. One nuisance is that the posterior for local variables
## (theta's) on the test data is the same as the prior for those variables,
## since they are unrelated to the training data. For those, we just have to
## draw from the prior. More sophisticated approaches are documented in Wallach
## et. al. "Evaluation Methods for Topic Models", but none of these seem to
## apply to the case where we've marginalized individual topic assignments.
##
## author: sankaran.kris@gmail.com
## date: 01/02/2017

library("argparser")
parser <- arg_parser("Evaluate LDA on the antibiotics dataset")
parser <- add_argument(parser, "--subject", help = "Subject on which to perform analysis", default = "F")
argv <- parse_args(parser)

###############################################################################
## Setup / libraries
###############################################################################
library("rstan")
library("tidyverse")
library("feather")
library("phyloseq")
library("MCMCpack")

abt <- get(load("../../data/antibiotics-study/abt.rda"))
abt <- abt %>%
  subset_samples(ind == argv$subject)

###############################################################################
## Fit to just the training data
###############################################################################

train_ix <- sample(1:nsamples(abt), 0.85 * nsamples(abt))
x <- t(get_taxa(abt))[train_ix, ]
dimnames(x) <- NULL
stan_data <- list(
  K = 4,
  V = ncol(x),
  D = nrow(x),
  n = x,
  alpha = rep(1, 4),
  gamma = rep(0.5, ncol(x))
)

## fit to just the training data
f <- stan_model(file = "../stan/lda_counts.stan")
stan_fit <- vb(
  f,
  data = stan_data,
  output_samples = 1000,
  eta = 1,
  adapt_engaged = FALSE
)

###############################################################################
## Evaluate log probabilities on the test data
###############################################################################
test_data <- stan_data
test_data$n <- t(get_taxa(abt))[-train_ix, ]
test_data$D <- nrow(test_data$n)
dummy_fit <- stan(fit = stan_fit, data = test_data, iter = 1)

post_params <- rstan::extract(stan_fit)
beta_star <- post_params$beta
theta_star <- post_params$theta

## sometimes constraints aren't exactly respected -- normalize just in case
for (i in 1:dim(beta_star)[1]) {
  for (j in 1:dim(beta_star)[2]) {
    beta_star[i, j, ] <- beta_star[i, j, ] / sum(beta_star[i, j, ])
  }

  for (j in 1:dim(theta_star)[2]) {
    theta_star[i, j, ] <- theta_star[i, j, ] / sum(theta_star[i, j, ])
  }
}

## Loop over posterior samples to approximate predictive likelihood
log_probs <- list(
  train = vector(length = nrow(beta_star)),
  test = vector(length = nrow(beta_star))
)

for (i in seq_len(nrow(beta_star))) {
  if (i %% 20 == 0) {
    message("Processing posterior sample ", i)
  }

  ## likelihood on test data
  theta_sim <- rdirichlet(test_data$D, test_data$alpha)
  pars <- unconstrain_pars(dummy_fit, list(theta = theta_sim, beta = beta_star[i,, ]))
  log_probs$test[i] <- log_prob(dummy_fit, upars = pars) / test_data$D

  ## likelihood on train data
  pars <- unconstrain_pars(stan_fit, list(theta = theta_star[i,, ], beta = beta_star[i,, ]))
  log_probs$train[i] <- log_prob(stan_fit, upars = pars) / stan_data$D
}

## visualize -- as expected train is much higher likelihood. Also, test is much
## larger variance, presumably because sampling theta from prior
log_prob_df <- bind_rows(log_probs) %>%
  gather()

ggplot(log_prob_df) +
  geom_histogram(
    aes(x = value, fill = key),
    binwidth = 100
  )
