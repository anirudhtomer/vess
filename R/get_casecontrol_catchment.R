#' @export
get_casecontrol_catchment = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                     overall_brand_proportions = c(0.3, 0.5, 0.2),
                                     overall_vaccine_coverage=0.3,
                                     overall_attack_rate_in_unvaccinated = 0.1,
                                     prob_getting_swabbed_given_ili_sari = 0.5,
                                     ili_sari_symptom_prob=NA,
                                     ili_sari_symptom_incidence_rate=NA,
                                     study_period_length = NA,
                                     total_cases=seq(50,500, 10)){
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)

  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  prob_unvaccinated_case = overall_attack_rate_in_unvaccinated * (1 - overall_vaccine_coverage)
  prob_vaccinated_case = brand_vaccine_coverages * (1 + (1-overall_attack_rate_in_unvaccinated)/((1-anticipated_brand_VEs) * overall_attack_rate_in_unvaccinated))^-1
  prob_case = sum(prob_unvaccinated_case, prob_vaccinated_case)

  catchment = total_cases / (prob_case * prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob)
  return(data.frame(total_cases=total_cases,
                    catchment = catchment,
                    overall_attack_rate_in_unvaccinated = overall_attack_rate_in_unvaccinated,
                    prob_getting_swabbed_given_ili_sari = prob_getting_swabbed_given_ili_sari,
                    ili_sari_symptom_incidence_rate = ili_sari_symptom_incidence_rate,
                    study_period_length = study_period_length,
                    ili_sari_symptom_prob = ili_sari_symptom_prob)
  )
}
