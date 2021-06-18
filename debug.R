get_cohort_full_table_irr = function(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                     overall_attack_rate_in_unvaccinated, study_period,
                                     proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated){
  #We are going to create a table with dimensions (row x columns) = (total_case_strains + 1) x (total_vaccine_brands + 1)
  #The extra 1 column is for unvaccinated
  #The extra 1 row is for controls

  #the last column of unvaccinated is given by
  prob_unvaccinated = 1 - overall_vaccine_coverage
  prob_vaccinated_with_brands = brand_proportions_in_vaccinated * overall_vaccine_coverage

  incidence_rate_unvaccinated =  - log(1 - overall_attack_rate_in_unvaccinated) / study_period
  incidence_rate_vaccinated  = (1 - anticipated_VE_for_each_brand_and_strain) * incidence_rate_unvaccinated

  full_table = cbind('unvaccinated'=c('p'=prob_unvaccinated, 'strain1'=incidence_rate_unvaccinated),
                     rbind('p'=prob_vaccinated_with_brands, t(incidence_rate_vaccinated)))

  full_table = rbind(full_table, 1 - exp(-full_table[2,] * study_period))

  return(full_table)
}

get_cohort_expectedCI_VE_irr = function(anticipated_VE_for_each_brand_and_strain=
                                          matrix(data=c(0.7, 0.4, 0.1),
                                                 nrow = 3, ncol = 1,
                                                 byrow = T,
                                                 dimnames = list(paste0('brand', 1:3), paste0('strain', 1))),
                                        brand_proportions_in_vaccinated =
                                          c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                        overall_vaccine_coverage=0.3,
                                        proportion_strains_in_unvaccinated_cases = c('strain1'=1),
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

  full_table = get_cohort_full_table_irr(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                         overall_attack_rate_in_unvaccinated, study_period,
                                         proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  #c(full_table) converts the table into a vector, going column by column
  flat_full_table = c(full_table)

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  event_probs = full_table[1,] * full_table[3,]
  names(event_probs) = paste0(names(event_probs), '_events')
  non_event_probs = full_table[1,] * (1-full_table[3,])
  names(non_event_probs) = paste0(names(non_event_probs), '_non_events')
  all_probs = c(event_probs, non_event_probs)

  cell_counts_sims = array(
    data = do.call(
      'cbind',
      lapply(missing_data_adjusted_total_subjects,
             FUN = function(n_subjects){

               N = rmultinom(n=nsims, size = n_subjects, prob=all_probs)

               total_event_matrix = N[1:(1+total_vaccine_brands),]
               total_non_event_matrix = N[(1+total_vaccine_brands + 1):nrow(N),]

               flat_total_events =c(total_event_matrix)

               #see simulations for the reason of this choice
               flat_events_person_time = pmax(
                 rnorm(
                   n = length(flat_total_events),
                   mean = flat_total_events*study_period/2,
                   sd = sqrt(flat_total_events * study_period^2/12)
                 ),
                 0
               )

               events_person_time_matrix = matrix(data = flat_events_person_time,
                                                  nrow = nrow(total_event_matrix),
                                                  byrow = F)

               total_person_time_matrix = events_person_time_matrix + total_non_event_matrix * study_period
               rownames(total_person_time_matrix) = paste0(names(event_probs), '_person_time')

               return(rbind(total_event_matrix, total_non_event_matrix, total_person_time_matrix))
             }
      )
    ), dim = c((1+total_vaccine_brands)*3, nsims, total_total_subject_settings)
  )

  UNVACCINATED_EVENTS_ROW = 1
  UNVACCINATED_PERSON_TIME_ROW = (1+total_vaccine_brands) * 2 + UNVACCINATED_EVENTS_ROW

  brand1_total_events_indices = relative_VE_combn[BRAND1,] + UNVACCINATED_EVENTS_ROW
  brand1_person_time_indices = relative_VE_combn[BRAND1,] + UNVACCINATED_PERSON_TIME_ROW
  brand2_total_events_indices = relative_VE_combn[BRAND2,] + UNVACCINATED_EVENTS_ROW
  brand2_person_time_indices = relative_VE_combn[BRAND2,] + UNVACCINATED_PERSON_TIME_ROW

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  brand1_total_events = cell_counts_sims[brand1_total_events_indices,,,drop=F]
  brand1_person_time = cell_counts_sims[brand1_person_time_indices,,,drop=F]
  brand2_total_events = cell_counts_sims[brand2_total_events_indices,,,drop=F]
  brand2_person_time = cell_counts_sims[brand2_person_time_indices,,,drop=F]

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

  rate_vaccine1 = full_table[2, brand1_total_events_indices]
  rate_vaccine2 = full_table[2, brand2_total_events_indices]
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
