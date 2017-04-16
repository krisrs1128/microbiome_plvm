/*
 * Dynamic Topic Model
 *
 * This models documents using a sequence of cluster and topic
 * parameters that evolve over time. See Blei & Lafferty 2006.
 */

data {
  int<lower=0> N; // number of samples
  int<lower=0> V; // number of words
  int<lower=0> T; // number of unique times
  int<lower=0> K; // number of topics
  real<lower=0> sigma_hyper[2]; // hyperparameters for topic evolution on simplex
  real<lower=0> delta_hyper[2]; // hyperparamters for membership evolution on simplex

  real times[T]; // unique times
  int<lower=0> times_mapping[N]; // times associated to each sample
  int<lower=0> x[N, V]; // word counts
}

parameters {
  vector[V] mu[T, K]; // logitted vocabulary probabilities
  vector[K] alpha[T]; // logitted cluster probabilities
  real<lower=0> sigma2; // topic evolution rate
  real<lower=0> delta2; // membership evolution rate
}

transformed parameters {
  vector[V] beta[T, K];
  vector[K] theta[T];
  real<lower=0> sigma;
  real<lower=0> delta;

  for (i in 1:T) {
    for (k in 1:K) {
      beta[i, k] = softmax(mu[i, k]);
    }
    theta[i] = softmax(alpha[i]);
  }

  sigma = sqrt(sigma2);
  delta = sqrt(delta2);
}

model {
  sigma2 ~ inv_gamma(sigma_hyper[1], sigma_hyper[2]);
  delta2 ~ inv_gamma(delta_hyper[1], delta_hyper[2]);

  for (i in 1:(T - 1)) {
    for (k in 1:K) {
      mu[i + 1, k] ~ normal(mu[i, k], sqrt(times[i + 1] - times[i]) * sigma);
    }
    alpha[i + 1] ~ normal(alpha[i], sqrt(times[i + 1] - times[i]) * delta);
  }

  for (i in 1:N) {
    vector[V] gamma;
    gamma = beta[times_mapping[i], 1] * theta[times_mapping[i]][1];
    for (k in 2:K) {
      gamma = gamma + beta[times_mapping[i], k] * theta[times_mapping[i]][k];
    }

    x[i] ~ multinomial(gamma);
  }
}

generated quantities {
  int<lower=0> x_sim[N, V]; // simulated word counts, for posterior checking
  for (i in 1:N) {
    vector[V] gamma;
    gamma = beta[times_mapping[i], 1] * theta[times_mapping[i]][1];
    for (k in 2:K) {
      gamma = gamma + beta[times_mapping[i], k] * theta[times_mapping[i]][k];
    }

    x_sim[i] = multinomial_rng(gamma, sum(x[i]));
  }
}
