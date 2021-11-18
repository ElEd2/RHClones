data {
int N;
vector<lower = 0>[N] age;
int<lower = 0> counts[N];
int<lower = 0> num_crypts[N];
}

parameters {
vector<lower = 0>[N] a_i;
//vector[N] b_i;
// real<upper = max(age)> x_intercept;
real x_intercept;
real<lower = 0> pop_mu;
real<lower = 0> pop_sd;
}

transformed parameters {
  // vector<lower = 0, upper = 1>[N]   p_i; // These are the 
  vector[N]   p_i; // These are the 
  p_i = a_i .* (age-x_intercept);
  for(i in 1:N){
    if(p_i[i]>1) p_i[i] = 1-1e-16;
    if(p_i[i]<0) p_i[i] = 1e-16;
  }
}

model {
// Likelihood
x_intercept  ~ normal(0, 10);
pop_mu       ~ gamma(1e-2, 1e-2);
pop_sd       ~ gamma(1e-2, 1e-2);

// Hierarchical model of slopes
a_i ~ normal(pop_mu, pop_sd);
counts ~ binomial(num_crypts, p_i);
}

generated quantities {
  int y_pred[N];
  vector[N] log_lik;
  // For the log likelihood
  for (n in 1:N) {
    // preferred Stan syntax as of version 2.10.0
    log_lik[n] = binomial_lpmf(counts[n] | num_crypts[n], p_i[n]);
  }

  for (n in 1:N) {
    // preferred Stan syntax as of version 2.10.0
    y_pred[n] = binomial_rng(num_crypts[n], p_i[n]);
  }

}







