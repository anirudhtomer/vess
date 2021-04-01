source("sample_size/util.R")

get_casecontrol_mindetVE = function(df_row){
  confounder_and_missing_data_adjustment = (1-df_row['prob_missing_data']) * (1-df_row['confounder_adjustment_Rsquared'])
  #the 'n' parameter need not be integer for this API.
  total_subjects = df_row['total_cases'] * (df_row['controls_per_case'] + 1) * confounder_and_missing_data_adjustment
  
  epi.sscc(OR = NA, p0 = df_row['overall_vaccine_coverage'], n = total_subjects,
           power = df_row['power'], r = df_row['controls_per_case'],
           sided.test = 2, conf.level = 1-df_row['alpha'], 
           method = "unmatched", fleiss = FALSE)$OR[1]
}


#first something important regarding vaccine coverage in controls
#P(Vaccination) = P(Vaccination | Case) P(Case) + P(Vaccination | Control) (1-P(Case))
#In general P(Case) is really low for diseases like influenza or even dengue. hardly 1-3% in the population
#Hence overall_vaccine_coverage is almost equal to vaccine coverage in controls

#second regarding catchment
#Unlike cohort scenario an additional filter for subjects is they need to be cases
#controls are chosen after cases are chosen. hence in calculation of catchment
#i need to know the probability of being a case
#given you have ILI symptoms and that you are swabbed

#all parameters are vectors
#coverage: coverage in overall population
casecontrol_mindet_overall_VE = function(power=0.8, alpha=0.05, 
                                         overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                         attack_rate_unvaccinated = c(0.01),
                                         controls_per_case=2,
                                         confounder_adjustment_Rsquared = 0,
                                         prob_missing_data = 0.1, 
                                         prob_getting_swabbed_given_ili_sari = 0.5,
                                         ili_sari_symptom_prob=NA, 
                                         ili_sari_symptom_incidence_rate=NA, 
                                         study_period_length = NA, 
                                         total_cases=seq(500,5000,500)){
  
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)
  
  ret = expand.grid(power=power, alpha=alpha, 
                    overall_vaccine_coverage=overall_vaccine_coverage, 
                    controls_per_case=controls_per_case, 
                    confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                    prob_missing_data=prob_missing_data, 
                    prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                    ili_sari_symptom_prob=ili_sari_symptom_prob, 
                    attack_rate_unvaccinated=attack_rate_unvaccinated, 
                    total_cases=total_cases)
  
  ret$mindet_overall_VE = apply(ret, MARGIN = 1, FUN = get_casecontrol_mindetVE)
  
  ret$prob_unvaccinated_case = ret$attack_rate_unvaccinated * (1 - ret$overall_vaccine_coverage)
  ret$prob_vaccinated_case = ret$overall_vaccine_coverage * (1 + (1-ret$attack_rate_unvaccinated)/((1-ret$mindet_overall_VE) * ret$attack_rate_unvaccinated))^-1
  
  ret$catchment = ret$total_cases / ((ret$prob_unvaccinated_case + ret$prob_vaccinated_case)*ret$prob_getting_swabbed_given_ili_sari * ret$ili_sari_symptom_prob)
  return(ret)
}

#all parameters are vectors
#coverage: coverage in overall population
#total_cases: total subjects in the total population, not just for brand
casecontrol_mindet_brand_VE = function(power=0.8, alpha=0.05, 
                                       brand_proportion = c(0.3, 0.5, 1),
                                       overall_vaccine_coverage=c(0.3, 0.5, 0.7), 
                                       attack_rate_unvaccinated = c(0.01),
                                       controls_per_case=2,
                                       confounder_adjustment_Rsquared = 0,
                                       prob_missing_data = 0.1, 
                                       prob_getting_swabbed_given_ili_sari = 0.5,
                                       ili_sari_symptom_prob=NA, 
                                       ili_sari_symptom_incidence_rate=NA, 
                                       study_period_length = NA, 
                                       total_cases=seq(500,5000,500)){
  
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)
  
  ret = expand.grid(power=power, alpha=alpha, 
                    brand_proportion = brand_proportion,
                    overall_vaccine_coverage=overall_vaccine_coverage, 
                    controls_per_case=controls_per_case, 
                    confounder_adjustment_Rsquared=confounder_adjustment_Rsquared, 
                    prob_missing_data=prob_missing_data, 
                    prob_getting_swabbed_given_ili_sari=prob_getting_swabbed_given_ili_sari, 
                    ili_sari_symptom_prob=ili_sari_symptom_prob, 
                    total_cases=total_cases)
  
  ret$mindet_overall_VE = apply(ret, MARGIN = 1, FUN = function(df_row){
    brand_coverage = df_row['brand_proportion']*df_row['overall_vaccine_coverage']
    brand_and_unvaccinated_coverage = brand_coverage + 1-df_row['overall_vaccine_coverage']
    df_row['overall_vaccine_coverage'] = brand_coverage / brand_and_unvaccinated_coverage
    
    #I am multiplying with total cases here, but it is actually the total subjects I should multiply to
    #But the total_subjects = total_cases * (controls_per_case + 1), in the get_casecontrol_mindetVE function
    #hence this multiplcation here makes sense
    df_row['total_cases'] = df_row['total_cases'] * brand_and_unvaccinated_coverage
    get_casecontrol_mindetVE(df_row)
  })
  
  #This is wrong but I have got to fix it
  ret$catchment = ret$total_cases / (ret$prob_getting_swabbed_given_ili_sari * ret$ili_sari_symptom_prob)
  return(ret)
}