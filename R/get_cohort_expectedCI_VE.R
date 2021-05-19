#' @importFrom utils combn
#' @export
get_cohort_expectedCI_VE = function(anticipated_VE_for_each_brand_and_strain=
                                      matrix(data=c(0.8, 0.5, 0.3, 0.3, 0.5, 0.8, 0.9, 0.5, 1), nrow = 3, ncol = 3, byrow = F,
                                             dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                    brand_proportions_in_vaccinated =
                                      c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                    overall_vaccine_coverage=0.3,
                                    proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                    overall_attack_rate_in_unvaccinated = 0.1,
                                    calculate_relative_VE = T,
                                    alpha=0.05,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1,
                                    total_subjects=seq(1000, 10000, 25),
                                    nsims = 500){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  full_table = get_cohort_full_table(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage, overall_attack_rate_in_unvaccinated,
                                     proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  STRAIN_ROW = 3
  CONTROL_ROW=total_case_strains + 1

  #c(full_table) converts the table into a vector, going column by column
  cell_counts_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_subjects, rmultinom, n=nsims, prob=c(full_table))),
                           dim = c((total_vaccine_brands+1)*(total_case_strains+1), nsims, total_total_subject_settings))

  #cell_counts_sims = rmultinom(n=nsims, size=missing_data_adjusted_total_subjects, prob = c(full_table))
  #Haldane's correction
  cell_counts_sims[cell_counts_sims==0] = 0.5

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  cases_vaccine1 = cell_counts_sims[(relative_VE_combn[1,]-1)*(total_case_strains +1) + relative_VE_combn[STRAIN_ROW,],,,drop=F]
  cases_vaccine2 = cell_counts_sims[(relative_VE_combn[2,]-1)*(total_case_strains +1) + relative_VE_combn[STRAIN_ROW,],,,drop=F]
  controls_vaccine1 = cell_counts_sims[(relative_VE_combn[1,]-1)*(total_case_strains +1) + CONTROL_ROW,,,drop=F]
  controls_vaccine2 = cell_counts_sims[(relative_VE_combn[2,]-1)*(total_case_strains +1) + CONTROL_ROW,,,drop=F]

  estimate = (cases_vaccine1/(cases_vaccine1 + controls_vaccine1))/(cases_vaccine2/(cases_vaccine2 + controls_vaccine2))
  standard_error = sqrt(1/cases_vaccine1 + 1/cases_vaccine2 - 1/(cases_vaccine1 + controls_vaccine1) - 1/(cases_vaccine2 + controls_vaccine2))
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  relative_risk_for_each_brand_and_strain = cbind(1 - t(anticipated_VE_for_each_brand_and_strain), 1)
  relative_risk_vaccine1 = relative_risk_for_each_brand_and_strain[(relative_VE_combn[1,]-1)*total_case_strains + relative_VE_combn[STRAIN_ROW,]]
  relative_risk_vaccine2 = relative_risk_for_each_brand_and_strain[(relative_VE_combn[2,]-1)*total_case_strains + relative_VE_combn[STRAIN_ROW,]]
  anticipated_VE = 1 - relative_risk_vaccine1/relative_risk_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[1,]),total_total_subject_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[2,]==total_vaccine_brands+1,
                                          no = paste("Brand", relative_VE_combn[2,]),
                                          yes = "Unvaccinated"), total_total_subject_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN_ROW,]), total_total_subject_settings),
                   total_subjects = rep(total_subjects, each=ncol(relative_VE_combn)),
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit))

  return(ret)
}
