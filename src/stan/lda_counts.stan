data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma;
}

parameters {
  simplex[K] theta[D]; // topic mixtures
  simplex[V] beta[K]; // word dist for k^th topic
}

model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);
  }

  for (k in 1:K) {
    beta[k] ~ dirichlet(gamma);
  }

  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    n[d] ~ multinomial(eta);
  }

}

generated quantities {
  int<lower=0> x_sim[D, V]; // simulated word counts, for posterior checking

  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    x_sim[d] = multinomial_rng(eta, sum(n[d]));
  }
}
