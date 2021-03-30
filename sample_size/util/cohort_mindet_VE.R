get_ili_sari_cumulative_prob=function(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length){
  if(is.na(ili_sari_cumulative_prob)){
    if(!is.na(ili_sari_incidence_rate) & !is.na(study_period_length)){
      rate = expand.grid(ili_sari_incidence_rate, study_period_length)
      return(1 - exp(-rate[,1] * rate[,2]))
    }else{
      stop("Non NA values required for either 'ili_sari_cumulative_prob' or for both 'ili_sari_incidence_rate' and 'study_period_length'")
    }
  }else{
    return(ili_sari_cumulative_prob)
  }
}

VE_calculator = function(df_row){
  confounder_loss_follow_adjustment = (1-df_row['prob_loss_follow_up']) * (1-df_row['confounder_adjustment_Rsquared'])
  #the 'n' parameter need not be integer for this API.
  epi.sscohortc(irexp1 = NA, irexp0 = df_row['attack_rate_unvaccinated'], 
                n = df_row['total_subjects'] * confounder_loss_follow_adjustment,
                power = df_row['power'], 
                r = df_row['overall_vaccine_coverage']/(1-df_row['overall_vaccine_coverage']), 
                design = 1, sided.test = 2, conf.level = 1-df_row['alpha'])$irr[1]  
}

#all parameters are vectors
#coverage: coverage in overall population
cohort_mindet_overall_VE = function(power=0.8, alpha=0.05, 
                                    overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                    attack_rate_unvaccinated=0.5,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_loss_follow_up = 0.1, 
                                    prob_getting_swabbed_given_ili_sari = 0.5,
                                    ili_sari_cumulative_prob=NA, 
                                    ili_sari_incidence_rate=NA, 
                                    study_period_length = NA, 
                                    total_subjects=seq(500,5000,500)){
  
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  ret = expand.grid(power=power, alpha=alpha, 
                    overall_vaccine_coverage=overall_vaccine_coverage, 
                    attack_rate_unvaccinated=attack_rate_unvaccinated, 
                    confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                    prob_loss_follow_up=prob_loss_follow_up, 
                    prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                    ili_sari_cumulative_prob=ili_sari_cumulative_prob, 
                    total_subjects=total_subjects)
  
  ret$mindet_overall_VE = apply(ret, MARGIN = 1, FUN = VE_calculator)
  ret$catchment = ret$total_subjects / (ret$prob_getting_swabbed_given_ili_sari * ret$ili_sari_cumulative_prob)
  return(ret)
}

#all parameters are vectors
#coverage: coverage in overall population
#total_subjects: total subjects in the total population, not just for brand
cohort_mindet_brand_VE = function(power=0.8, alpha=0.05, 
                                  brand_proportion = c(0.3, 0.5, 1),
                                  overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                  attack_rate_unvaccinated=0.5,
                                  confounder_adjustment_Rsquared = 0,
                                  prob_loss_follow_up = 0.1, 
                                  prob_getting_swabbed_given_ili_sari = 0.5,
                                  ili_sari_cumulative_prob=NA, 
                                  ili_sari_incidence_rate=NA, 
                                  study_period_length = NA, 
                                  total_subjects=seq(500,5000,500)){
  ili_sari_cumulative_prob = get_ili_sari_cumulative_prob(ili_sari_cumulative_prob, ili_sari_incidence_rate, study_period_length)
  
  ret = expand.grid(power=power, alpha=alpha, 
                    brand_proportion = brand_proportion,
                    overall_vaccine_coverage=overall_vaccine_coverage, 
                    attack_rate_unvaccinated=attack_rate_unvaccinated, 
                    confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                    prob_loss_follow_up=prob_loss_follow_up, 
                    prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                    ili_sari_cumulative_prob=ili_sari_cumulative_prob, 
                    total_subjects=total_subjects)
  
  ret$mindet_overall_VE = apply(ret, MARGIN = 1, FUN = function(df_row){
    brand_coverage = df_row['brand_proportion']*df_row['overall_vaccine_coverage']
    brand_and_unvaccinated_coverage = brand_coverage + 1-df_row['overall_vaccine_coverage']
    df_row['overall_vaccine_coverage'] = brand_coverage / brand_and_unvaccinated_coverage
    df_row['total_subjects'] = df_row['total_subjects'] * brand_and_unvaccinated_coverage
    VE_calculator(df_row)
  })
  ret$catchment = ret$total_subjects / (ret$prob_getting_swabbed_given_ili_sari * ret$ili_sari_cumulative_prob)
  return(ret)
}