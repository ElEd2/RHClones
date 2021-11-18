#' Fit model to patch sizes
#'
#' @export
#' @param pat_data data frame containing all clone info.
#' @param num_chains Number of MCMC chains to run default 4.
#' @param multicore Use multicore boolean, defaults to TRUE.
#' @return A list with two MCMC tibbles one with per patient info and one with population level info
#'
clonal_mark_fission_outpat = function(pat_data, num_chains = 4, multicore = T)
{
  # Check whether to run single or multi core
  if(multicore){
    options(mc.cores = parallel::detectCores())
    num_cores_use = getOption("mc.cores", 1L)
  }else{
    num_cores_use = 1
  }
  # Check that there are clones to fit!
  pat_filt      = pat_data %>% dplyr::filter(monoclonal>0)
  if(nrow(pat_filt)==0){
    print("No clones!")
    return(NA)
  }
  patch_mat     = pat_filt %>% dplyr::select(dplyr::contains("full")) %>% as.matrix(.)
  stan_data     = list(N = nrow(patch_mat), age = pat_filt$Age, counts = patch_mat, max_patch = ncol(patch_mat), num_crypts = pat_filt$num_crypts)
  starting_vals = list(rho = rep(1e-3, stan_data$N), pop_mu = 1e-3, pop_sd = 1e-3*0.1)

  # run model
  fission_fit = rstan::sampling(stanmodels$crypt_fission_inference, data = stan_data, iter = 4e3, chains = num_chains,
                         verbose = T,  thin = 5, init = rep(list(starting_vals),num_chains),
                         control = list(max_treedepth=70))
  ## Check Rhat
  # if(any(summary(fission_fit)$summary[,10] > 1.1)) message("CONVERGENCE PROBLEMS!!!!!")

  # Pop parameters
  pop_pars_fission = c("pop_nu", "pop_mu", "pop_sd")
  tbl_fission = dplyr::bind_cols(rstan::extract(fission_fit, pop_pars_fission))

  # Per pat parameter
  tbl_fission_indv = dplyr::as_tibble(rstan::extract(fission_fit, "rho")[[1]])
  per_pat_fission = tbl_fission_indv %>%
    tidyr::gather(pat, val) %>%
    dplyr::group_by(pat) %>%
    dplyr::summarise(rho_mid = stats::quantile(val, 0.5),
                     rho_min = stats::quantile(val, 0.975),
                     rho_max = stats::quantile(val, 0.025))

  pat_filt = pat_filt %>% dplyr::mutate(pat = paste0("V", 1:dplyr::n())) %>% dplyr::left_join(per_pat_fission)
  list(pop = tbl_fission, indv = pat_filt)
}


#' Calculate probability of chosen patch size or smaller
#'
#' @export
#' @param rho fission rate
#' @param age patient Age
#' @param size_measured patch size of interest
#' @return Probability of patch size or less
#'
calc_prob_below = function(rho, age, size_measured){
  max_patch = 30
  prob_i = rep(0, max_patch)
  for (k in 1:max_patch){
    prob_i[k] = (1-exp(-rho*age))^k/(rho*age*k);
  }
  cumsum(prob_i)[size_measured]
}

#' Calculate probability of chosen patch size
#'
#' @export
#' @param rho fission rate
#' @param age patient Age
#' @param size_measured patch size of interest
#' @return Probability of patch size
#'
calc_prob_at = function(rho, age, size_measured){
  max_patch = 30
  prob_i = rep(0, max_patch)
  for (k in 1:max_patch){
    prob_i[k] = (1-exp(-rho*age))^k/(rho*age*k);
  }
  prob_i[size_measured]
}



