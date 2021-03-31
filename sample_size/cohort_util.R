source('sample_size/util.R')

get_cohort_mindetVE = function(df_row){
  confounder_loss_follow_adjustment = (1-df_row['prob_missing_data']) * (1-df_row['confounder_adjustment_Rsquared'])
  #the 'n' parameter need not be integer for this API.
  epi.sscohortc(irexp1 = NA, irexp0 = df_row['attack_rate_unvaccinated'], 
                n = df_row['total_subjects'] * confounder_loss_follow_adjustment,
                power = df_row['power'], 
                r = df_row['overall_vaccine_coverage']/(1-df_row['overall_vaccine_coverage']), 
                design = 1, sided.test = 2, conf.level = 1-df_row['alpha'])$irr[1]  
}

fill_cohort_VE_expected_CI_limits = function(anticipated_VE, alpha, 
                                             overall_vaccine_coverage, 
                                             attack_rate_unvaccinated,
                                             confounder_adjustment_Rsquared, 
                                             prob_missing_data, total_subjects, nsims){
  
  loss_follow_up_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))
  
  overall_vaccine_coverage = overall_vaccine_coverage
  
  relative_risk = 1 - anticipated_VE
  p_a_vaccinated_case = attack_rate_unvaccinated * overall_vaccine_coverage * relative_risk
  p_b_vaccinated_control = overall_vaccine_coverage * (1- attack_rate_unvaccinated*relative_risk)
  p_c_unvaccinated_case = attack_rate_unvaccinated * (1-overall_vaccine_coverage)
  p_d_unvaccinated_control = (1-attack_rate_unvaccinated) * (1-overall_vaccine_coverage)
  
  upp = rep(NA, nsims)
  low = rep(NA, nsims)
  
  for(j in 1:nsims){
    temp = c(rmultinom(n=1, size=loss_follow_up_adjusted_total_subjects, prob = c(p_a_vaccinated_case, p_b_vaccinated_control, p_c_unvaccinated_case, p_d_unvaccinated_control)))
    
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
  
  return(c(avg_lower_limit_VE = mean(1-upp, na.rm = T),
           avg_upper_limit_VE = mean(1-low, na.rm = T)))
}