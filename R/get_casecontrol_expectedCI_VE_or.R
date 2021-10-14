#' Expected lower and expected upper limit of the confidence intervals of strain and vaccine-specific efficacy based on odds ratio.
#' @description
#' The function `get_cohort_expectedCI_VE_or` simulates confindence intervals for strain and vaccine-specific efficacy (VE) for a given sample size.
#' The efficacy is defined as `VE = 1 - odds ratio`, and the function returns expected lower and expected upper confidence interval limit
#' for both absolute and relative VE.
#'
#' @param anticipated_VE_for_each_brand_and_strain a matrix of vaccine efficacy of each vaccine (row) against each strain (column). Each value must be a real number between 0 and 1.
#' @param brand_proportions_in_vaccinated a vector denoting the proportion in which vaccines are given in the vaccinated subjects of the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_vaccine_coverage the proportion of the study cohort that will be vaccinated. It should be a real number between 0 and 1.
#' @param proportion_strains_in_unvaccinated_cases a vector of the proportions in which each strain is expected to be present in the unvaccinated and infected subjects in the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param controls_per_case the number of controls to be sampled per case in the study cohort. It should be a real number between 0 and Inf.
#' @param calculate_relative_VE a logical indicating if calculations should also be done for relative vaccine efficacy (default `TRUE`).
#' @param alpha controls the width `[100*alpha/2, 100*(1 - alpha/2)]%` of the confidence interval. It is a numeric that must take a value between 0 and 1.
#' @param confounder_adjustment_Rsquared we use this parameter to adjust the calculations for potential confounders using the methodology proposed by Hsieh and Lavori (2000). It represents the amount of variance (R^2) explained in a regression model where vaccination status is the outcome and confounders of interest are predictors. It is a numeric that must take a value between 0 (no adjustment for confounders) and 1.
#' @param prob_missing_data to adjust the calculations for non-informative and random subject loss to follow-up/dropout. it should take a numeric value between 0 and 1.
#' @param total_cases a vector of the number of cases in our cohort for which calculations should be done.
#' @param nsims total number of Monte Carlo simulations conducted.
#'
#' In this function efficacy is defined as `VE = 1 - odds ratio`, where 'odds ratio' is
#' the ratio of odds of being vaccinated with a particular vaccine versus being unvaccinated among cases and the same odds in controls.
#' When the groups being compared are a particular vaccine versus placebo then we call the VE
#' as the absolute VE of the vaccine. For `M` vaccines there are `M` absolute VE, one each for the `M` vaccines.
#' When the groups being compared are a particular vaccine versus another vaccine then we call the VE
#' as the relative VE of the vaccines, for a particular strain. For `M` vaccines and `I` strains there are `I x 2 x utils::combn(M, 2)`
#' permutations of relative VE of two vaccines against the same strain.
#'
#' We first transform the user inputs for `I` strains and `M` vaccines into two conditional tables, one for cases, and another for controls.
#' First an `I x (M + 1)` cross table of the probability of being unvaccinated or vaccinated with a vaccine, given that the subject is a case. The sum of the cells of this cases table is equal to 1.
#' Second a `1 x (M + 1)` vector of the probability of being unvaccinated or vaccinated with a vaccine, given that the subject is a controls. The sum of this controls vector is equal to 1.
#' The first column in both of these corresponds to subjects who are unvaccinated and the remaining `M` columns correspond to subjects who are vaccinated with a particular vaccine.
#'
#' The next step is to simulate the data. To speed up our computations we sample an `I x (M + 1)`
#' cross table of cases data from a multinomial distribution with probabilities taken from our cases table.
#' The total cases sampled in the cross table are are `total_cases * (1 - prob_missing_data)`.
#' We also sample a vector of `1 x (M + 1)` controls data from a multinomial distribution with probabilities taken from our controls vector.
#' The total controls sampled in the cross table are are `controls_per_case * total_cases * (1 - prob_missing_data)`.
#' We then estimate the absolute and relative VE of each vaccine using odds ration based on the sampled data.
#' The confidence intervals with widths `[100*alpha/2, 100*(1 - alpha/2)]%` are obtained using normal approximation
#' to the distribution of log of odds ratio (Morris and Gardner, 1988).
#' To adjust for confounders, the standard-error used in the confidence interval is rescaled to `SE/(1 - confounder_adjustment_Rsquared)` (Hsieh and Lavori, 2000)
#' We repeat this procedure `nsims` times, and in each such simulation we obtain `nsims` confidence intervals.
#'
#' To conduct simulations faster all the calculations are done without using for loops.
#' Instead we use a three-dimensional R arrays, with one-dimension for `nsims`,
#' another for the sample size vector `total_cases`,
#' and another for a vector containing the flattened cross table of simulated data on cases and controls.
#'
#' @return A data frame consisting of the input parameters, absolute and relative VE combinations,
#' and the expected lower and expected upper width of the confidence intervals for each absolute and relative VE combination.
#'
#' @examples As an example we recommend running the function without passing any parameter to it.
#' The default scenario is for three vaccines and three pathogen strains.
#'
#' @references
#' 1. Hsieh, F. Y., & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled clinical trials, 21(6), 552-560.
#' 2. Morris, J. A., & Gardner, M. J. (1988). Statistics in medicine: Calculating confidence intervals for relative risks (odds ratios) and standardised ratios and rates. British medical journal (Clinical research ed.), 296(6632), 1313.
#' @importFrom utils combn
#' @export
get_casecontrol_expectedCI_VE_or = function(anticipated_VE_for_each_brand_and_strain=
                                              matrix(data=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), nrow = 3, ncol = 3, byrow = T,
                                                     dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                            brand_proportions_in_vaccinated =
                                              c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                            overall_vaccine_coverage=0.3,
                                            proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                            controls_per_case=2,
                                            calculate_relative_VE = T,
                                            alpha=0.05,
                                            confounder_adjustment_Rsquared = 0,
                                            prob_missing_data = 0.1,
                                            total_cases=seq(50,500, 10),
                                            nsims = 500){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_vaccine_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  total_total_case_settings = length(total_cases)
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  missing_data_adjusted_total_controls = missing_data_adjusted_total_cases * controls_per_case

  case_control_full_tables = get_case_control_full_tables(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                                          proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  flat_full_table_cases = c(case_control_full_tables$full_table_cases)

  ##c(full_table_cases) converts the table into a vector, going column by column
  cell_counts_case_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_cases, rmultinom, n=nsims,
                                                               prob=flat_full_table_cases)),
                                dim = c((total_vaccine_brands+1)*total_case_strains, nsims, total_total_case_settings))
  #Haldane's correction
  cell_counts_case_sims[cell_counts_case_sims==0] = 0.5

  #cell_count_control_sims has dimension total_vaccine_brands x nsims
  cell_counts_control_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_controls, rmultinom, n=nsims,
                                                                  prob=case_control_full_tables$full_vector_controls)),
                                   dim = c(total_vaccine_brands+1, nsims, total_total_case_settings))
  #Haldane's correction
  cell_counts_control_sims[cell_counts_control_sims==0] = 0.5

  UNVACCINATED_ROW = 1
  brand1_control_indices = relative_VE_combn[BRAND1,] + UNVACCINATED_ROW
  brand1_case_indices = brand1_control_indices + (relative_VE_combn[STRAIN,]-1)*nrow(case_control_full_tables$full_table_cases)
  brand2_control_indices = relative_VE_combn[BRAND2,] + UNVACCINATED_ROW
  brand2_case_indices =  brand2_control_indices + (relative_VE_combn[STRAIN,]-1)*nrow(case_control_full_tables$full_table_cases)

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  cases_vaccine1 = cell_counts_case_sims[brand1_case_indices,,,drop=F]
  cases_vaccine2 = cell_counts_case_sims[brand2_case_indices,,,drop=F]
  controls_vaccine1 = cell_counts_control_sims[brand1_control_indices,,,drop=F]
  controls_vaccine2 = cell_counts_control_sims[brand2_control_indices,,,drop=F]

  estimate = (cases_vaccine1/cases_vaccine2) / (controls_vaccine1/controls_vaccine2)
  standard_error = sqrt(1/cases_vaccine1 + 1/cases_vaccine2 + 1/controls_vaccine1 + 1/controls_vaccine2)
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  odds_vaccinated_vaccine1 = flat_full_table_cases[brand1_case_indices]/case_control_full_tables$full_vector_controls[brand1_control_indices]
  odds_vaccinated_vaccine2 = flat_full_table_cases[brand2_case_indices]/case_control_full_tables$full_vector_controls[brand2_control_indices]
  anticipated_VE = 1 - odds_vaccinated_vaccine1/odds_vaccinated_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[BRAND1,]),total_total_case_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                          no = paste("Brand", relative_VE_combn[BRAND2,]),
                                          yes = "Unvaccinated"), total_total_case_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN,]), total_total_case_settings),
                   total_cases = rep(total_cases, each=ncol(relative_VE_combn)),
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit),
                   alpha = alpha,
                   overall_vaccine_coverage = overall_vaccine_coverage,
                   controls_per_case = controls_per_case,
                   confounder_adjustment_Rsquared = confounder_adjustment_Rsquared,
                   prob_missing_data = prob_missing_data
  )

  ret = cbind(ret,
              t(data.frame(brand_proportions_in_vaccinated,
                           row.names = paste0('brand_prop_', 1:total_vaccine_brands))),
              t(data.frame(proportion_strains_in_unvaccinated_cases,
                           row.names = paste0('strain_prop_', 1:total_case_strains))))

  return(ret)
}
