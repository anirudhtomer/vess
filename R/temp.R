#' @importFrom utils combn
#' @export
get_cohort_expectedCI_VE_irr = function(anticipated_VE_for_each_brand_and_strain=
                                          matrix(data=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1), nrow = 3, ncol = 3, byrow = F,
                                                 dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                        brand_proportions_in_vaccinated =
                                          c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                        overall_vaccine_coverage=0.3,
                                        proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                        overall_attack_rate_in_unvaccinated = 0.1,
                                        study_period=1,
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

  all_tables = get_cohort_full_table_irr(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                         overall_attack_rate_in_unvaccinated, study_period,
                                         proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  #c(full_table) converts the table into a vector, going column by column
  flat_full_table_prob = c(all_tables$full_table_prob)
  flat_event_rates_table = c(all_tables$conditional_table_incidence_rate)

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  cell_counts_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_subjects, rmultinom, n=nsims, prob=flat_full_table_prob)),
                           dim = c((1+total_vaccine_brands)*(1+total_case_strains), nsims, total_total_subject_settings))

  person_time_sims = array(dim = dim(cell_counts_sims))

  CONTROL_ROW = 1
  CONTROL_INDICES = CONTROL_ROW + (0:total_vaccine_brands)*(total_case_strains+1)
  person_time_sims[CONTROL_INDICES,,] = cell_counts_sims[CONTROL_INDICES,,,drop=F] * study_period

  trunc_exp_means = trunc_exp_mean(event_rate = flat_event_rates_table, trunc_limit = study_period)
  #this happens when VE = 100%, that is event rate is zero. consequently the person time becomes NaN
  #i want to set that NaN to 0 because 100% VE means no events and 0 person time contributed
  trunc_exp_means[is.nan(trunc_exp_means)] = 0
  cases_person_time_sum_mean = cell_counts_sims[-CONTROL_INDICES,,,drop=F] * trunc_exp_means
  cases_person_time_sum_sd = sqrt(cell_counts_sims[-CONTROL_INDICES,,,drop=F] * study_period^2/12)
  person_time_sims[-CONTROL_INDICES,,] = array(
    data = pmax(
      rnorm(
        n = prod(dim(cases_person_time_sum_mean)),
        mean = c(cases_person_time_sum_mean),
        sd = c(cases_person_time_sum_sd)
      ),
      0
    ),
    dim = dim(cases_person_time_sum_mean)
  )

  browser()
  #summarize person time by summing it up per brand. that becomes the denominator
  # person_time_sims = apply(person_time_sims, c(2,3), FUN = function(x){
  #   cumsum_overall = c(0, cumsum(x))
  #   a = seq(1, length(x), by = (1 + total_case_strains))
  #   b = seq(1+total_case_strains, length(x), by = (1 + total_case_strains))
  #   return(cumsum_overall[b + 1] - cumsum_overall[a])
  # })
  person_time_sims = array(
    c(
      apply(person_time_sims, c(3), function(x){
        cumsum_overall = rbind(0,colCumsums(x))
        a = seq(1, nrow(x), by = (1 + total_case_strains))
        b = seq(1+total_case_strains, nrow(x), by = (1 + total_case_strains))
        return(cumsum_overall[b + 1,] - cumsum_overall[a,])
      })
    ),
    dim = c(1+total_vaccine_brands, nsims, total_total_subject_settings)
  )
  #think about event rate zero: person time should be study period
  cell_counts_sims[cell_counts_sims==0] = 0.5

  UNVACCINATED_EVENTS_ROW = 1

  brand1_total_events_indices = relative_VE_combn[BRAND1,]*(1 + total_case_strains) + relative_VE_combn[STRAIN,] + CONTROL_ROW
  brand1_person_time_indices = relative_VE_combn[BRAND1,] + UNVACCINATED_EVENTS_ROW
  brand2_total_events_indices = relative_VE_combn[BRAND2,]*(1 + total_case_strains) + relative_VE_combn[STRAIN,] + CONTROL_ROW
  brand2_person_time_indices = relative_VE_combn[BRAND2,] + UNVACCINATED_EVENTS_ROW

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  brand1_total_events = cell_counts_sims[brand1_total_events_indices,,,drop=F]
  brand1_person_time = person_time_sims[brand1_person_time_indices,,,drop=F]
  brand2_total_events = cell_counts_sims[brand2_total_events_indices,,,drop=F]
  brand2_person_time = person_time_sims[brand2_person_time_indices,,,drop=F]

  estimate = (brand1_total_events/brand1_person_time)/(brand2_total_events/brand2_person_time)
  standard_error = sqrt(1/brand1_total_events + 1/brand2_total_events)
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  rate_table = rbind(NA, all_tables$conditional_table_incidence_rate)
  rate_vaccine1 = rate_table[brand1_total_events_indices]
  rate_vaccine2 = rate_table[brand2_total_events_indices]
  anticipated_VE = 1 - rate_vaccine1/rate_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[BRAND1,]), total_total_subject_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                          no = paste("Brand", relative_VE_combn[BRAND2,]),
                                          yes = "Unvaccinated"), total_total_subject_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN,]), total_total_subject_settings),
                   total_subjects = rep(total_subjects, each=ncol(relative_VE_combn)),
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit),
                   alpha = alpha,
                   overall_vaccine_coverage = overall_vaccine_coverage,
                   overall_attack_rate_in_unvaccinated = overall_attack_rate_in_unvaccinated,
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
