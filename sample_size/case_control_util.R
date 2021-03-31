source('sample_size/util.R')

get_casecontrol_mindetVE = function(df_row){
  confounder_loss_follow_adjustment = (1-df_row['prob_missing_data']) * (1-df_row['confounder_adjustment_Rsquared'])
  #the 'n' parameter need not be integer for this API.
  total_subjects = df_row['total_cases'] * (df_row['controls_per_case'] + 1) * confounder_loss_follow_adjustment
  
  epi.sscc(OR = NA, p0 = df_row['overall_vaccine_coverage'], n = total_subjects,
           power = df_row['power'], r = df_row['controls_per_case'],
           sided.test = 2, conf.level = 1-df_row['alpha'], 
           method = "unmatched", fleiss = FALSE)$OR[1]
}

fill_casecontrol_VE_expected_CI_limits = function(df_row, nsims){
  
  upp = rep(NA, nsims)
  low = rep(NA, nsims)
  
  total_controls = df_row$total_cases * df_row$controls_per_case
  for(j in 1:nsims){
    
    statusvaccine_control = rbinom(total_controls, 1, df$overall_vaccine_coverage)
    
    b_yesvaccine_yescontrol = sum(statusvaccine_control)
    d_novaccine_yescontrol = total_controls-b_yesvaccine_yescontrol
    
    a_divided_c = (1-df$anticipated_VE) * (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
    
    a_yesvaccine_yescase = sum(rbinom(df_row$total_cases, 1, 1 / (1 + a_divided_c^-1)))
    c_novaccine_yescase = df_row$total_cases - a_yesvaccine_yescase
    
    vaccine = c(rep(1, a_yesvaccine_yescase), rep(0, c_novaccine_yescase),
                rep(1, b_yesvaccine_yescontrol), rep(0,d_novaccine_yescontrol))
    hospitalized = c(rep(1, df_row$total_cases), rep(0, total_controls))
    
    estimate = (a_yesvaccine_yescase/c_novaccine_yescase) / (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
    se = sqrt(1/a_yesvaccine_yescase + 1/c_novaccine_yescase + 1/b_yesvaccine_yescontrol + 1/d_novaccine_yescontrol)
    confounder_adjusted_se = se / (1-df_row$confounder_adjustment_Rsquared)
    low[j] = exp(log(estimate[j]) + qnorm(alpha/2) * se)
    upp[j] = exp(log(estimate[j]) + qnorm(1 - alpha/2) * se)
  }
  
  df_row$avg_lower_limit_VE = mean(1-upp, na.rm = T)
  df_row$avg_upper_limit_VE = mean(1-low, na.rm = T)
  return(df_row)
}