source("sample_size/util.R")

#all parameters are scalars
#coverage: coverage in overall population
#nsims: larger value give more accurate results
cohort_precision_VE = function(anticipated_VE=0.8, alpha=0.05, 
                               brand_proportion=0.8,
                               overall_vaccine_coverage=0.3, 
                               attack_rate_unvaccinated=0.5,
                               confounder_adjustment_Rsquared = 0,
                               prob_missing_data = 0.1, 
                               prob_getting_swabbed_given_ili_sari = 0.5,
                               ili_sari_symptom_prob=NA, 
                               ili_sari_symptom_incidence_rate=NA, 
                               study_period_length = NA, 
                               total_subjects=500,
                               nsims = 500){
  
  #I do it first because it checks if default NA parameters have been given required value
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)
  
  brand_coverage = brand_proportion *overall_vaccine_coverage
  brand_and_unvaccinated_coverage = brand_coverage + 1-overall_vaccine_coverage
  
  total_subjects = total_subjects * brand_and_unvaccinated_coverage
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))
  
  overall_vaccine_coverage=brand_coverage / brand_and_unvaccinated_coverage
  
  relative_risk = 1 - anticipated_VE
  p_a_vaccinated_case = attack_rate_unvaccinated * overall_vaccine_coverage * relative_risk
  p_b_vaccinated_control = overall_vaccine_coverage * (1- attack_rate_unvaccinated*relative_risk)
  p_c_unvaccinated_case = attack_rate_unvaccinated * (1-overall_vaccine_coverage)
  p_d_unvaccinated_control = (1-attack_rate_unvaccinated) * (1-overall_vaccine_coverage)
  
  upp = rep(NA, nsims)
  low = rep(NA, nsims)
  
  for(j in 1:nsims){
    temp = c(rmultinom(n=1, size=missing_data_adjusted_total_subjects, prob = c(p_a_vaccinated_case, p_b_vaccinated_control, p_c_unvaccinated_case, p_d_unvaccinated_control)))
    
    a_vaccinated_case = temp[1]
    b_vaccinated_control = temp[2]
    c_unvaccinated_case = temp[3]
    d_unvaccinated_control = temp[4]
    
    estimate = (a_vaccinated_case/(a_vaccinated_case + b_vaccinated_control)) / (c_unvaccinated_case/(c_unvaccinated_case + d_unvaccinated_control))
    se = sqrt(1/a_vaccinated_case + 1/c_unvaccinated_case - 1/(a_vaccinated_case + b_vaccinated_control) - 1/(c_unvaccinated_case + d_unvaccinated_control))
    confounder_adjusted_se = se / (1-confounder_adjustment_Rsquared)
    low = exp(log(estimate) + qnorm(alpha/2) * confounder_adjusted_se)
    upp = exp(log(estimate) + qnorm(1-alpha/2) * confounder_adjusted_se)
  }
  
  catchment = total_subjects / (prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob)
  
  return(c(avg_lower_limit_VE = mean(1-upp, na.rm = T), 
           avg_upper_limit_VE = mean(1-low, na.rm = T),
           catchment=catchment))
}
