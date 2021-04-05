source("sample_size/util.R")

#all parameters are scalars
#coverage: coverage in overall population
#nsims: larger value give more accurate results
cohort_precision_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                               overall_brand_proportions = c(0.3, 0.5, 0.2),
                               overall_vaccine_coverage=0.3, 
                               attack_rate_unvaccinated = 0.1,
                               alpha=0.05, 
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
  
  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }
  
  total_vaccines = length(anticipated_brand_VEs)
  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))
  
  #first n-1 numbers are 
  relative_risks = c(1 - anticipated_brand_VEs)
  prob_vaccinated_case = attack_rate_unvaccinated * brand_vaccine_coverages * relative_risks
  #made a dummy variable for ease of reading the code and the math
  relative_risk_placebo = 1
  prob_unvaccinated_case = attack_rate_unvaccinated * (1-overall_vaccine_coverage) * relative_risk_placebo
  prob_vaccinated_control = brand_vaccine_coverages * (1- attack_rate_unvaccinated*relative_risks)
  prob_unvaccinated_control = (1-attack_rate_unvaccinated) * (1-overall_vaccine_coverage)
  
  upp = matrix(ncol=total_vaccines, nrow = nsims, data = NA)
  low = matrix(ncol=total_vaccines, nrow = nsims, data = NA)
  
  CASE=1
  CONTROL=2
  for(j in 1:nsims){
    cell_counts = c(rmultinom(n=1, size=missing_data_adjusted_total_subjects, 
                              prob = c(prob_vaccinated_case, prob_unvaccinated_case, 
                                       prob_vaccinated_control, prob_unvaccinated_control)))
    #2 x n table, row 1 for cases, row 2 for controls. nth column for unvaccinated and other columns for different vaccines
    cell_counts = matrix(cell_counts, nrow=2, byrow = T)
    
    estimates = (cell_counts[CASE, 1:total_vaccines]/colSums(cell_counts[, 1:total_vaccines])) / (cell_counts[CASE, total_vaccines+1]/sum(cell_counts[, total_vaccines+1]))
    standard_errors = sqrt(
      cell_counts[CASE, 1:total_vaccines]^-1 + 
        cell_counts[CASE, total_vaccines+1]^-1 - 
        (cell_counts[CASE, 1:total_vaccines] + cell_counts[CONTROL, 1:total_vaccines])^-1 - 
        (cell_counts[CASE, total_vaccines+1] + cell_counts[CONTROL, total_vaccines+1])^-1
    )
    
    confounder_adjusted_standard_errors = standard_errors / (1-confounder_adjustment_Rsquared)
    low[j,] = exp(log(estimates) + qnorm(alpha/2) * confounder_adjusted_standard_errors)
    upp[j,] = exp(log(estimates) + qnorm(1-alpha/2) * confounder_adjusted_standard_errors)
  }
  
  avg_lower_limit_VE = apply(1-upp, MARGIN = 2, FUN = mean, na.rm=T)
  avg_upper_limit_VE = apply(1-low, MARGIN = 2, FUN = mean, na.rm=T)
  
  catchment = total_subjects / (prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob)
  
  return(list(avg_lower_limit_VE=avg_lower_limit_VE, avg_upper_limit_VE=avg_upper_limit_VE, catchment=catchment))
}
