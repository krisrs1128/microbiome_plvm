/*
 * Dynamic unigram model
 *
 * This models documents whose word proportions evolve according to a
 * softmax on an ordinary random walk. This is the one-topic version
 * of the kalman-filter based dynamic topic model.
 */

data {
  int<lower=0> N; // number of samples
  int<lower=0> V; // number of words
  int<lower=0> T; // number of unique times
  real<lower=0> a0; // prior on sigma
  real<lower=0> b0; // prior on sigma

  real times[T]; // unique times
  int<lower=0> times_mapping[N]; // times associated to each sample
  int<lower=0> x[N, V]; // word counts
}

parameters {
  vector[V] mu[T];
  real<lower=0> sigma2;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(sigma2);
}

model {
  sigma2 ~ inv_gamma(a0, b0);

  for (i in 1:(T - 1)) {
    mu[i + 1] ~ normal(mu[i], sqrt(times[i + 1] - times[i]) * sigma);
  }

  for (i in 1:N) {
    x[i] ~ multinomial(softmax(mu[times_mapping[i]]));
  }
}

generated quantities {
  int<lower=0> x_sim[N, V]; // simulated word counts, for posterior checking
  for (i in 1:N) {
    x_sim[i] = multinomial_rng(softmax(mu[times_mapping[i]]), sum(x[i]));
  }
}
