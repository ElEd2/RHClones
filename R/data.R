#' MCMC tibble with regression parameters for HDAC6, KDM6A, MAOA, mPAS, NONO, P53, PTEN, STAG2
#'
#' @format A data frame with 6 variables: \code{x_intercept}, \code{mean_slope},
#'   \code{sd_slope}, \code{partials_mean}, \code{partials_sd},
#'   \code{clonal_mark}
"all_slope_part"

#' MCMC tibble with fission population estimates for NONO, P53, mPAS, MAOA, Stag2
#'
#' @format A data frame with 3 variables: \code{pop_mu}, \code{pop_sd},
#'   \code{clonal_mark}
#' \describe{
#' \item{pop_mu}{Mean fission of population}
#' \item{pop_sd}{Standard deviation of fission in population}
#' \item{clonal_mark}{One of several clonal marks}
#' }
"fission_marks"
