source("sample_size/cohort_util.R")

#all parameters are vectors
#coverage: coverage in overall population
cohort_mindet_overall_VE = function(power=0.8, alpha=0.05, 
                                    overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                    attack_rate_unvaccinated=0.5,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1, 
                                    prob_getting_swabbed_given_ili_sari = 0.5,
                                    ili_sari_cumulative_prob=NA, 
                                    ili_sari_incidence_rate=NA, 
                                    study_period_length = NA, 
                                    total_subjects=seq(500,5000,500)){
  
  #I do it first because it checks if default NA parameters have been given required value
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  ret_mindet_VE = expand.grid(power=power, alpha=alpha, 
                              overall_vaccine_coverage=overall_vaccine_coverage, 
                              attack_rate_unvaccinated=attack_rate_unvaccinated, 
                              confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                              prob_missing_data=prob_missing_data, 
                              total_subjects=total_subjects)
  ret_mindet_VE$mindet_overall_VE = apply(ret_mindet_VE, MARGIN = 1, FUN = get_cohort_mindetVE)
  
  ret_catchment = expand.grid(prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                              ili_sari_cumulative_prob=ili_sari_cumulative_prob, 
                              total_subjects=total_subjects)
  ret_catchment$catchment = ret_catchment$total_subjects / (ret_catchment$prob_getting_swabbed_given_ili_sari * ret_catchment$ili_sari_cumulative_prob)
  
  return(list(mindet_VE=ret_mindet_VE, catchment=ret_catchment))
}

#all parameters are vectors
#coverage: coverage in overall population
#total_subjects: total subjects in the total population, not just for brand
cohort_mindet_brand_VE = function(power=0.8, alpha=0.05, 
                                  brand_proportion = c(0.3, 0.5, 1),
                                  overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                  attack_rate_unvaccinated=0.5,
                                  confounder_adjustment_Rsquared = 0,
                                  prob_missing_data = 0.1, 
                                  prob_getting_swabbed_given_ili_sari = 0.5,
                                  ili_sari_cumulative_prob=NA, 
                                  ili_sari_incidence_rate=NA, 
                                  study_period_length = NA, 
                                  total_subjects=seq(500,5000,500)){
  
  #I do it first because it checks if default NA parameters have been given required value
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  ret_mindet_VE = expand.grid(power=power, alpha=alpha, 
                              brand_proportion = brand_proportion,
                              overall_vaccine_coverage=overall_vaccine_coverage, 
                              attack_rate_unvaccinated=attack_rate_unvaccinated, 
                              confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                              prob_missing_data=prob_missing_data, 
                              total_subjects=total_subjects)
  
  ret_mindet_VE$mindet_overall_VE = apply(ret_mindet_VE, MARGIN = 1, FUN = function(df_row){
    brand_coverage = df_row['brand_proportion']*df_row['overall_vaccine_coverage']
    brand_and_unvaccinated_coverage = brand_coverage + 1-df_row['overall_vaccine_coverage']
    df_row['overall_vaccine_coverage'] = brand_coverage / brand_and_unvaccinated_coverage
    df_row['total_subjects'] = df_row['total_subjects'] * brand_and_unvaccinated_coverage
    get_cohort_mindetVE(df_row)
  })
  
  ret_catchment = expand.grid(prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                              ili_sari_cumulative_prob=ili_sari_cumulative_prob, 
                              total_subjects=total_subjects)
  
  ret_catchment$catchment = ret_catchment$total_subjects / (ret_catchment$prob_getting_swabbed_given_ili_sari * ret_catchment$ili_sari_cumulative_prob)
  
  return(list(mindet_VE=ret_mindet_VE, catchment=ret_catchment))
}