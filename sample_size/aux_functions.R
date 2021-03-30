##########################################
# Auxiliary function for the functions   #
# Casecontrol_overall_precision          #
# and                                    #
# Casecontrol_brandspecific_precision    #
##########################################

smoothAndRound = function(dat, xname){
  fit = smooth.spline(x = dat[,xname], y = dat$VE_lower, df = SPLINE_DF)
  return(round(pmax(0, stats:::predict.smooth.spline(object = fit, x = dat[,xname])$y), 2))
}

findpar=function(a,b,c,d,R){
  x1 = (1/2)* sqrt(( (a^2)* R^2 + 2*a*b*R - 2*a*c*R^2 + 2*a*c*R + b^2 + 2*b*c*R - 2*b*c + (c^2)* R^2 - 2*(c^2)*R + c^2)/((R - 1)^2)) + (a*R + b + c*R - c)/(2*(R - 1))
  x2 =a-x1;
  x3 =c-x1;
  x4 =d-x2;
  return(list(x1,x2,x3,x4))
}

getVELowerCohortPrecision = function(alpha, df){
  
  df$VE = df$VE/100
  df$AR = df$AR/100
  df$Coverage = df$Coverage/100
  
  df$RR = 1 - df$VE
  df$p_a_vaccinated_case = df$AR * df$Coverage * df$RR
  df$p_b_vaccinated_control = df$Coverage * (1- df$AR*df$RR)
  df$p_c_unvaccinated_case = df$AR * (1-df$Coverage)
  df$p_d_unvaccinated_control = (1-df$AR) * (1-df$Coverage)
  
  cl = makeCluster(MAX_CORES)
  registerDoParallel(cl)
  cachePresence_and_VE_lower = foreach(i=1:nrow(df), .combine='rbind', .export = 'cohort_cache') %dopar% {
    
    matches = cohort_cache$Coverage == df$Coverage[i] & cohort_cache$SS == df$SS[i] & cohort_cache$AR == df$AR[i] & cohort_cache$VE == df$VE[i]
    if(any(matches)){
      return(c(T, cohort_cache$VE_lower[matches]))
    }
    
    set.seed(2020 + i)
    n_studies_back_intime = 2000
    estimate = rep(NA, n_studies_back_intime)
    upp = rep(NA, n_studies_back_intime)
    low = rep(NA, n_studies_back_intime)
    
    for(j in 1:n_studies_back_intime){
      temp = c(rmultinom(n=1, size=df$SS[i], c(df$p_a_vaccinated_case[i], df$p_b_vaccinated_control[i], df$p_c_unvaccinated_case[i], df$p_d_unvaccinated_control[i])))
      
      a_vaccinated_case = temp[1]
      b_vaccinated_control = temp[2]
      c_unvaccinated_case = temp[3]
      d_unvaccinated_control = temp[4]
      
      estimate[j] = (a_vaccinated_case/(a_vaccinated_case + b_vaccinated_control)) / (c_unvaccinated_case/(c_unvaccinated_case + d_unvaccinated_control))
      se = sqrt(1/a_vaccinated_case + 1/c_unvaccinated_case - 1/(a_vaccinated_case + b_vaccinated_control) - 1/(c_unvaccinated_case + d_unvaccinated_control))
      low[j] = exp(log(estimate[j]) + qnorm(alpha/2) * se)
      upp[j] = exp(log(estimate[j]) + qnorm(1-alpha/2) * se)
    }
    
    return(c(F, 100 * max(0, 1-mean(upp[!is.infinite(upp)], na.rm = T))))
  }
  stopCluster(cl)
  df$VE_lower = cachePresence_and_VE_lower[,2]
  cohort_cache <<- rbind(cohort_cache, df[!cachePresence_and_VE_lower[,1],!colnames(df) %in% c("prop_brand")])
  return(df$VE_lower)
}

getVELowerCaseControlPrecision = function(alpha, df){
  cl = makeCluster(MAX_CORES)
  registerDoParallel(cl)
  cachePresence_and_VE_lower = foreach(i=1:nrow(df), .combine='rbind', .export = 'casecontrol_cache') %dopar% {
    
    matches = casecontrol_cache$Coverage == df$Coverage[i] & casecontrol_cache$nr_cases == df$nr_cases[i] & casecontrol_cache$r == df$r[i] & casecontrol_cache$VE == df$VE[i]
    if(any(matches)){
      return(c(T, casecontrol_cache$VE_lower[matches]))
    }
    
    set.seed(2020 + i)
    control_case_ratio = df$r[i]
    total_cases = df$nr_cases[i]
    total_controls = total_cases * control_case_ratio
    
    n_studies_back_intime = 2000
    estimate = rep(NA, n_studies_back_intime)
    upp = rep(NA, n_studies_back_intime)
    low = rep(NA, n_studies_back_intime)
    
    for(j in 1:n_studies_back_intime){
      
      statusvaccine_control = rbinom(total_controls, 1, df$Coverage[i]/100)
      
      b_yesvaccine_yescontrol = sum(statusvaccine_control)
      d_novaccine_yescontrol = total_controls-b_yesvaccine_yescontrol
      
      a_divided_c = (1-df$VE[i]/100) * (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
      
      a_yesvaccine_yescase = sum(rbinom(total_cases, 1, 1 / (1 + a_divided_c^-1)))
      c_novaccine_yescase = total_cases - a_yesvaccine_yescase
      
      vaccine = c(rep(1, a_yesvaccine_yescase), rep(0, c_novaccine_yescase),
                  rep(1, b_yesvaccine_yescontrol), rep(0,d_novaccine_yescontrol))
      hospitalized = c(rep(1, total_cases), rep(0, total_controls))
      
      #tt = fisher.test(table(vaccine, hospitalized), or = 1)
      estimate[j] = (a_yesvaccine_yescase/c_novaccine_yescase) / (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
      se = sqrt(1/a_yesvaccine_yescase + 1/c_novaccine_yescase + 1/b_yesvaccine_yescontrol + 1/d_novaccine_yescontrol)
      low[j] = exp(log(estimate[j]) + qnorm(alpha/2) * se)
      upp[j] = exp(log(estimate[j]) + qnorm(1 - alpha/2) * se)
    }
    #return(list(mean(estimate, na.rm = T), mean(low, na.rm = T), mean(upp, na.rm = T)))
    return(c(F, 100 * max(0, 1-mean(upp[!is.infinite(upp)], na.rm = T))))
  }
  stopCluster(cl)
  df$VE_lower = cachePresence_and_VE_lower[,2]
  casecontrol_cache <<- rbind(casecontrol_cache, df[!cachePresence_and_VE_lower[,1],!colnames(df) %in% c("prop_brand")])
  return(df$VE_lower)
}