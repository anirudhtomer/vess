source("sample_size/util.R")

get_casecontrol_catchment = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                     overall_brand_proportions = c(0.3, 0.5, 0.2),
                                     overall_vaccine_coverage=0.3, 
                                     attack_rate_unvaccinated = 0.1,
                                     prob_getting_swabbed_given_ili_sari = 0.5,
                                     ili_sari_symptom_prob=NA, 
                                     ili_sari_symptom_incidence_rate=NA,
                                     study_period_length = NA,
                                     total_cases = 500){
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)
  
  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  prob_unvaccinated_case = attack_rate_unvaccinated * (1 - overall_vaccine_coverage)
  prob_vaccinated_case = brand_vaccine_coverages * (1 + (1-attack_rate_unvaccinated)/((1-anticipated_brand_VEs) * attack_rate_unvaccinated))^-1
  prob_case = sum(prob_unvaccinated_case, prob_vaccinated_case)
  
  return(total_cases / (prob_case * prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob))
}