source("sample_size/cohort_util.R")

#all parameters are scalars
#coverage: coverage in overall population
#nsims: larger value give more accurate results
cohort_precision_overall_VE = function(anticipated_VE=0.8, alpha=0.05, 
                                       overall_vaccine_coverage=0.3, 
                                       attack_rate_unvaccinated=0.5,
                                       confounder_adjustment_Rsquared = 0,
                                       prob_missing_data = 0.1, 
                                       prob_getting_swabbed_given_ili_sari = 0.5,
                                       ili_sari_cumulative_prob=NA, 
                                       ili_sari_incidence_rate=NA, 
                                       study_period_length = NA, 
                                       total_subjects=500,
                                       nsims = 500){
  
  #I do it first because it checks if default NA parameters have been given required value
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  CI_limits = fill_cohort_VE_expected_CI_limits(anticipated_VE, alpha, 
                                                overall_vaccine_coverage, 
                                                attack_rate_unvaccinated,
                                                confounder_adjustment_Rsquared, 
                                                prob_missing_data, total_subjects, nsims)
  
  catchment = total_subjects / (prob_getting_swabbed_given_ili_sari * ili_sari_cumulative_prob)
  
  return(c(CI_limits, catchment=catchment))
}

#all parameters are scalars
#coverage: coverage in overall population
#nsims: larger value give more accurate results
cohort_precision_brand_VE = function(anticipated_VE=0.8, alpha=0.05, 
                                     brand_proportion=0.8,
                                     overall_vaccine_coverage=0.3, 
                                     attack_rate_unvaccinated=0.5,
                                     confounder_adjustment_Rsquared = 0,
                                     prob_missing_data = 0.1, 
                                     prob_getting_swabbed_given_ili_sari = 0.5,
                                     ili_sari_cumulative_prob=NA, 
                                     ili_sari_incidence_rate=NA, 
                                     study_period_length = NA, 
                                     total_subjects=500,
                                     nsims = 500){
  
  #I do it first because it checks if default NA parameters have been given required value
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  brand_coverage = brand_proportion *overall_vaccine_coverage
  brand_and_unvaccinated_coverage = brand_coverage + 1-overall_vaccine_coverage
  CI_limits = fill_cohort_VE_expected_CI_limits(anticipated_VE, alpha, 
                                                overall_vaccine_coverage=brand_coverage / brand_and_unvaccinated_coverage, 
                                                attack_rate_unvaccinated,
                                                confounder_adjustment_Rsquared, 
                                                prob_missing_data, 
                                                total_subjects = total_subjects * brand_and_unvaccinated_coverage, nsims)
  
  catchment = total_subjects / (prob_getting_swabbed_given_ili_sari * ili_sari_cumulative_prob)
  
  return(c(CI_limits, catchment=catchment))
}
