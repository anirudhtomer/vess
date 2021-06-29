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

  relative_VE_combn = get_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

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
