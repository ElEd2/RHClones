data {
  int N; // Number of patients
  int max_patch; // maximum patch size
  vector<lower = 0>[N] age;
  vector<lower = 0>[N] num_crypts;
  int<lower = 0> counts[N, max_patch];
  // row_vector counts[M];
}

parameters {
  // real<lower = 0> rho;
  vector<lower = 0>[N]      rho;
  real<lower = 0>        pop_mu;
  real<lower = 0>        pop_sd;
  real<lower = 0>        pop_nu;
}

model {
  vector[max_patch]  prob_i;
  int                n_vals;
  real               p2_rnd;
  // Likelihood
  pop_mu ~ normal(0, 1e-2);
  pop_sd ~ normal(0, 1e-2); //gamma(2, 0.1);
  pop_nu ~  gamma(2,  0.1);
  
  
  // Hierarchical model of parameter rho
  // rho ~ gamma(2, 0.1);
  // rho ~ normal(pop_mu, pop_sd);
  
  rho ~ student_t(pop_nu, pop_mu, pop_sd);
  
  // Likelihood
  for (j in 1:N){
    for (k in 1:max_patch){
      prob_i[k] = (1-exp(-rho[j]*age[j]))^k/(rho[j]*age[j]*k);
    } 
    // Correction for random double counts ==========
    n_vals = counts[j][1] + 2*counts[j][2];
    p2_rnd = 6.0*(n_vals-1)/(2*num_crypts[j]);
    if(p2_rnd<0) p2_rnd = 0;
    prob_i[2] = prob_i[2] + p2_rnd;
    prob_i[1] = 1-(sum(prob_i)-prob_i[1]); // Adjust prob of patch = 1 ; 1 - sum(prob[2:max_patch])
    // ==============================================
    prob_i = prob_i/sum(prob_i);
    counts[j] ~ multinomial(prob_i); 
  }
}


