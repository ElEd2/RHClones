data {
int N;
vector<lower = 0>[N] age;
int<lower = 0> counts[N];
int<lower = 0> num_crypts[N];
}

parameters {
vector<lower = 0, upper = 1>[N] intercept_i;
real<lower = 0> pop_mu;
real<lower = 0> pop_sd;
}

model {
// Likelihood
pop_mu ~ gamma(2, 0.1);
pop_sd ~ gamma(2, 0.1);

// Hierarchical model of slopes
intercept_i ~ normal(pop_mu, pop_sd);
counts      ~ binomial(num_crypts, intercept_i);
}




