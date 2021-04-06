#' @export
get_cohort_catchment = function(prob_getting_swabbed_given_ili_sari = 0.5,
                                ili_sari_symptom_prob=NA,
                                ili_sari_symptom_incidence_rate=NA,
                                study_period_length = NA,
                                total_subjects = 500){
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)

  return(total_subjects / (prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob))
}
