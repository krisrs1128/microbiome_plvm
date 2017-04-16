/*
* Nonnegative Matrix Factorization usinga Gamma-Poisson Model
*
* reference:https://arxiv.org/pdf/1506.03431.pdf
*/

data{
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  int<lower=0> y[N, P];
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0, upper=1> zero_inf_prob;
}

parameters{
  vector<lower=0>[K] theta[N]; // scores
  vector<lower=0>[K] beta[P]; // latent factors
}

model{
  for(i in 1:N) {
    theta[i] ~ gamma(a,b); // componentwise gamma
  }

  for (j in 1:P) {
    beta[j] ~ gamma(c,d);//componentwise gamma
  }

  for (i in 1:N) {
    for (j in 1:P) {
      if (y[i, j] == 0) {
        target += log_sum_exp(bernoulli_lpmf(1 | zero_inf_prob),
                              bernoulli_lpmf(0 | zero_inf_prob) + poisson_lpmf(y[i, j] | theta[i]' * beta[j]));

      } else {
        target += bernoulli_lpmf(0 | zero_inf_prob) + poisson_lpmf(y[i, j] | theta[i]' * beta[j]);
      }
    }
  }
}
