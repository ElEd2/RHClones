# library(gridExtra)
# library(tidyverse)
# library(stringr)
# library(scales)
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

fancy_scientific = function(l_all){
  exponent = -4
  out_val = NULL
  for(l in l_all)
  {
    l = format(l, scientific = TRUE)
    get_base    = stringr::str_split(string = l, pattern = "e")[[1]]
    current_exp = as.numeric(get_base[2])
    diff_exp    = current_exp-exponent
    current_non_exp = as.numeric(get_base[1])
    non_exp = format(current_non_exp*10^(diff_exp), scientific = F)
    wanted_num = paste0(non_exp, "%*%10^", exponent)  # return this as an expression
    out_val = c(out_val, parse(text=wanted_num))
  }
  out_val
}

getFit = function(time_vals, mcmc_mark, interval = 0.9){
  num_tp = length(time_vals)
  result = matrix(0, num_tp, 3)
  for(ii in 1:num_tp){
    time_i       = time_vals[ii]
    y_vals       = mcmc_mark$mean_slope*(time_i-mcmc_mark$x_intercept)
    result[ii,]  = stats::quantile(y_vals, c(0.025, 0.5, 0.975))
  }
  # bind_rows(result)
  dplyr::tibble(ymin = result[,1], ymid = result[,2], ymax = result[,3], time = time_vals)
}


#' Plot monoclonals and partials and fit
#'
#' @export
#' @param pat_data tibble with patient data.
#' @param col_mono Colour for monoclonal
#' @param col_part Colour for partials
#' @param mcmc_fit Include fit
#' @param with_partials Include partials
#' @return grid plot
#'
plot_clone_data = function(pat_data, col_mono = "steelblue4", col_part = "steelblue2", mcmc_fit = NULL, with_partials = F)
{
  pp_reg = pat_data %>% dplyr::mutate(freq_part = partial/num_crypts, freq_mono = monoclonal/num_crypts) %>%
    ggplot2::ggplot(aes(x = Age)) + ggplot2::geom_point( aes(size = num_crypts,  y = freq_mono), col = col_mono, alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(0.7,3), name = "Number of Crypts", label = scales::comma) + ggplot2::theme(legend.position = c(0.15,0.68)) +
    ggplot2::xlab("Patient Age") + ggplot2::ylab("Clone frequency") + ggplot2::xlim(0, 100)+ ggplot2::scale_y_continuous(labels=fancy_scientific)
  if(!is.null(mcmc_fit)){
    slopes_CI = mcmc_fit %>% getFit(25:100, .)
    pp_reg =  pp_reg + ggplot2::geom_line(data = slopes_CI, aes(y = ymid,  x = time), col = "red") +
      ggplot2::geom_ribbon(data = slopes_CI, aes(ymin = ymin,  ymax = ymax, x = time), alpha = 0.2)
  }

  if(with_partials){
    pp_reg =  pp_reg + ggplot2::geom_point(aes(y = freq_part), col = col_part, pch = 15, alpha = 0.7)
  }
  return(pp_reg)
}

#' Plot partials
#'
#' @export
#' @param pat_data tibble with patient data.
#' @param col_part Colour for partials
#' @return plot
#'
plot_clone_only_part = function(pat_data, col_part = "steelblue2")
{
  pp_reg = pat_data %>% dplyr::mutate(freq_part = partial/num_crypts, freq_mono = monoclonal/num_crypts) %>%
    ggplot2::ggplot(ggplot2::aes(x = Age)) + ggplot2::geom_point( ggplot2::aes(size = num_crypts,  y = freq_part), col = col_part,  pch = 15, alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(0.7,3), name = "Number of Crypts", label = scales::comma) + ggplot2::theme(legend.position = c(0.15,0.78)) +
    ggplot2::xlab("Patient Age") + ggplot2::ylab("Clone frequency") + ggplot2::xlim(0, 100)+ ggplot2::scale_y_continuous(labels=fancy_scientific)
  pp_reg
}

#' Plot monoclonals and partials separately
#'
#' @export
#' @param pat_data tibble with patient data.
#' @param col_mono Colour for monoclonal
#' @param col_part Colour for partials
#' @return grid plot
#'
plot_full_part_sep = function(pat_data, col_mono = "steelblue4", col_part = "steelblue2")
{
  pp1 = plot_clone_data(pat_data, col_mono, col_part, with_partials = T)
  pp2 = plot_clone_only_part(pat_data, col_part) + ggplot2::theme(legend.position="none")
  gridExtra::arrangeGrob(pp1, pp2, heights = c(0.6, 0.4))
}

#' Fit model to monoclonals and partials
#'
#' @export
#' @param count_mono Vector of monoclonal counts.
#' @param count_part Vector of partial clone counts.
#' @param crypts_counted Vector of crypt counts.
#' @param Age Vector of Ages.
#' @param num_chains Number of MCMC chains to run default 4.
#' @param multicore Use multicore boolean, defaults to TRUE.
#' @return A tibble with the MCMC output of monoclonal and partial clone model fit
#'
clonal_mark_regression = function(count_mono, count_part, crypts_counted, Age, num_chains = 4, multicore = T)
{
  # Check whether to run single or multi core
  if(multicore){
    options(mc.cores = parallel::detectCores())
    num_cores_use = getOption("mc.cores", 1L)
  }else{
    num_cores_use = 1
  }

  pat_data      = dplyr::tibble(Age, count_mono = as.integer(count_mono),
                                count_part = as.integer(count_part), crypts_counted = as.integer(crypts_counted),
                                freq_mono = count_mono/crypts_counted, freq_part = count_part/crypts_counted)

  coef_vals     = abs(stats::coef(stats::lm(data = pat_data, freq_mono~Age)))
  mean_part     = mean(pat_data$count_part/pat_data$crypts_counted) + 1e-10
  starting_vals = list(a_i = rep(coef_vals[[2]], nrow(pat_data)), x_intercept = 0)
  starting_vals_part = list(intercept_i = rep(mean_part, nrow(pat_data)), pop_mu =mean_part, pop_sd = mean_part*0.1)

  stan_data      = list(N = nrow(pat_data), age = pat_data$Age, num_crypts = pat_data$crypts_counted,  counts = pat_data$count_mono)
  stan_data_part = list(N = nrow(pat_data), age = pat_data$Age, num_crypts = pat_data$crypts_counted,  counts = pat_data$count_part)

  # run model
  mono_fit = rstan::sampling(stanmodels$pop_model_normal, data = stan_data, iter = 1e4, chains = num_chains,
                      verbose = F,  thin = 5, init = rep(list(starting_vals),num_chains), cores = num_cores_use)

  part_fit = rstan::sampling(stanmodels$partials_normal, data = stan_data_part, iter = 1e4, chains = num_chains, verbose = F,  thin = 5,
                      init = rep(list(starting_vals_part),num_chains), control = list(adapt_delta = 0.9999, stepsize = 0.0001), cores = num_cores_use)

  pars_mono = c("x_intercept", "pop_mu", "pop_sd")
  pars_part = c("pop_mu", "pop_sd")

  tbl_reg  = dplyr::bind_cols(rstan::extract(mono_fit, pars_mono), rstan::extract(part_fit, pars_part))
  colnames(tbl_reg) = c("x_intercept", "mean_slope", "sd_slope", "partials_mean", "partials_sd")
  tbl_reg
}

