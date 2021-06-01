# smoothAndRound = function(dat, xname){
#   fit = stats::smooth.spline(x = dat[,xname], y = dat$VE_lower, df = SPLINE_DF)
#   return(round(pmax(0, stats:::predict.smooth.spline(object = fit, x = dat[,xname])$y), 2))
# }

BRAND1 = 'Brand 1'
BRAND2 = 'Brand 2'
STRAIN = 'Strain'
CONTROLS = 'controls'
UNVACCINATED = 'unvaccinated'

get_ili_sari_symptom_prob=function(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length){
  if(is.na(ili_sari_symptom_prob)){
    if(!is.na(ili_sari_symptom_incidence_rate) & !is.na(study_period_length)){
      rate = expand.grid(ili_sari_symptom_incidence_rate, study_period_length)
      return(1 - exp(-rate[,1] * rate[,2]))
    }else{
      stop("Non NA values required for either 'ili_sari_symptom_prob' or for both 'ili_sari_symptom_incidence_rate' and 'study_period_length'")
    }
  }else{
    return(ili_sari_symptom_prob)
  }
}

check_input = function(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases){
  if(!sum(brand_proportions_in_vaccinated, na.rm = T)==1){
    stop("Sum of the values of the parameter 'brand_proportions_in_vaccinated' should be equal to 1")
  }

  if(!sum(proportion_strains_in_unvaccinated_cases, na.rm = T)==1){
    stop("Sum of the values of the parameter 'proportion_strains_in_unvaccinated_cases' should be equal to 1")
  }

  if(nrow(anticipated_VE_for_each_brand_and_strain)!=length(brand_proportions_in_vaccinated)){
    stop("Total number of rows of the parameter 'anticipated_VE_for_each_brand_and_strain' should be equal to length of the parameter 'brand_proportions_in_vaccinated'")
  }

  if(ncol(anticipated_VE_for_each_brand_and_strain)!=length(proportion_strains_in_unvaccinated_cases)){
    stop("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_strain' should be equal to length of the parameter 'proportion_strains_in_unvaccinated_cases'")
  }
}

get_comparison_combinations = function(total_vaccine_brands, total_case_strains, calculate_relative_VE){
  relative_VE_combn = matrix(c(1:total_vaccine_brands, rep(0, total_vaccine_brands)), byrow = T, nrow = 2)
  if(calculate_relative_VE & total_vaccine_brands>1){
    relative_VE_combn = cbind(combn(total_vaccine_brands, 2), relative_VE_combn)
  }
  relative_VE_combn = rbind(relative_VE_combn[,rep(1:ncol(relative_VE_combn), each=total_case_strains), drop=F],
                            rep(1:total_case_strains, ncol(relative_VE_combn)))

  rownames(relative_VE_combn) = c(BRAND1, BRAND2, STRAIN)
  colnames(relative_VE_combn) = paste('Combination' , 1:ncol(relative_VE_combn))
  return(relative_VE_combn)
}

get_cohort_full_table = function(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage, overall_attack_rate_in_unvaccinated,
                                 proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated){
  #We are going to create a table with dimensions (row x columns) = (total_case_strains + 1) x (total_vaccine_brands + 1)
  #The extra 1 column is for unvaccinated
  #The extra 1 row is for controls

  #the last column of unvaccinated is given by
  prob_unvaccinated = 1 - overall_vaccine_coverage
  prob_case_given_unvaccinated = overall_attack_rate_in_unvaccinated

  #the cell number [total_case_strains + 1, total_vaccine_brands + 1]
  prob_control_given_unvaccinated = 1 - prob_case_given_unvaccinated
  prob_control_and_unvaccinated = prob_control_given_unvaccinated * prob_unvaccinated

  #the cell numbers [1:total_case_strains, total_vaccine_brands + 1]
  #P(strain & case | unvaccinated) = p(strain | case, unvaccinated) * p(case|unvaccinated)
  prob_case_each_strain_given_unvaccinated = proportion_strains_in_unvaccinated_cases * prob_case_given_unvaccinated
  prob_case_each_strain_and_unvaccinated = prob_case_each_strain_given_unvaccinated * prob_unvaccinated

  #So now we are done with the last column of our table. the remaining columns are for vaccines
  #sum of columns 1:total_vaccine_brands
  prob_vaccinated_with_brands = brand_proportions_in_vaccinated * overall_vaccine_coverage

  relative_risk_for_each_brand_and_strain = 1 - t(anticipated_VE_for_each_brand_and_strain)
  risk_for_each_brand_and_strain_given_unvaccinated = prob_case_each_strain_and_unvaccinated/(prob_case_each_strain_and_unvaccinated + prob_control_and_unvaccinated)
  risk_for_each_brand_and_strain_given_vaccinated = relative_risk_for_each_brand_and_strain * risk_for_each_brand_and_strain_given_unvaccinated
  odds_for_each_brand_and_strain_given_vaccinated = (risk_for_each_brand_and_strain_given_vaccinated^-1 - 1)^-1
  prob_control_each_brand_given_vaccinated = prob_vaccinated_with_brands / (1 + colSums(odds_for_each_brand_and_strain_given_vaccinated))
  prob_case_each_strain_and_each_brand_given_vaccinated = t(odds_for_each_brand_and_strain_given_vaccinated) * prob_control_each_brand_given_vaccinated

  full_table = cbind('unvaccinated'=c('controls'=prob_control_and_unvaccinated, prob_case_each_strain_and_unvaccinated),
                     rbind('controls'=prob_control_each_brand_given_vaccinated, t(prob_case_each_strain_and_each_brand_given_vaccinated)))

  return(full_table)
}


get_case_control_full_tables = function(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                       proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated){
  prob_vaccinated_each_brand = brand_proportions_in_vaccinated * overall_vaccine_coverage
  #P(Vaccination) = P(Vaccination | Case) P(Case) + P(Vaccination | Control) (1-P(Case))
  #In general P(Case) is really low for diseases like influenza or even dengue. hardly 1-3% in the population
  #Hence overall_vaccine_coverage is almost equal to vaccine coverage in controls
  prob_vaccinated_each_brand_given_control = prob_vaccinated_each_brand
  prob_unvaccinated_given_control = 1 - overall_vaccine_coverage

  #the idea in a case-control study is that we generate two sets of multinomial data
  #1. for controls
  #2. for cases
  #the counts remain constant

  odds_vaccinated_for_each_brand_and_strain = 1 - anticipated_VE_for_each_brand_and_strain
  prob_case_and_unvaccinated = sum(proportion_strains_in_unvaccinated_cases * (1 + colSums(odds_vaccinated_for_each_brand_and_strain * prob_vaccinated_each_brand_given_control)/prob_unvaccinated_given_control))^-1
  prob_unvaccinated_case_each_strain = prob_case_and_unvaccinated * proportion_strains_in_unvaccinated_cases
  prob_vaccinated_case_each_strain_and_brand = t(t(odds_vaccinated_for_each_brand_and_strain * prob_vaccinated_each_brand_given_control / prob_unvaccinated_given_control) * prob_unvaccinated_case_each_strain)

  return(list(full_table_cases=rbind('unvaccinated'=prob_unvaccinated_case_each_strain, prob_vaccinated_case_each_strain_and_brand),
              full_vector_controls = c('unvaccinated'=prob_unvaccinated_given_control, prob_vaccinated_each_brand_given_control)))
}

