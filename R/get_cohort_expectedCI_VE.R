#' @importFrom utils combn
#' @export
get_cohort_expectedCI_VE = function(anticipated_VE_for_each_brand_and_strain=
                                      list('brand1'=c('strain1'=0.8, 'strain2'=0.3, 'strain3'=0.9),
                                           'brand2'=c('strain1'=0.5, 'strain2'=0.5, 'strain3'=0.5),
                                           'brand3'=c('strain1'=0.3, 'strain2'=0.8, 'strain3'=0.1)),
                                    brand_proportions_in_vaccinated =
                                      c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                    overall_vaccine_coverage=0.3,
                                    proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                    overall_attack_rate_in_unvaccinated = 0.1,
                                    calculate_relative_VE = T,
                                    alpha=0.05,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1,
                                    total_subjects=seq(1000, 10000, 25),
                                    nsims = 500){

  if(!sum(brand_proportions_in_vaccinated, na.rm = T)==1){
    stop("Sum of the values of the parameter 'brand_proportions_in_vaccinated' should be equal to 1")
  }

  if(!sum(proportion_strains_in_unvaccinated_cases, na.rm = T)==1){
    stop("Sum of the values of the parameter 'proportion_strains_in_unvaccinated_cases' should be equal to 1")
  }

  if(length(anticipated_VE_for_each_brand_and_strain)!=length(brand_proportions_in_vaccinated)){
    stop("Length of the parameter 'anticipated_VE_for_each_brand_and_strain' should be equal to length of the parameter 'brand_proportions_in_vaccinated'")
  }

  total_total_subject_settings = length(total_subjects)

  #Assuming that the all subjects have the same chance of dropout due to one or another reason
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  #For each vaccine we know in what proportion are they given in the general population
  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  #For each strain we know in what proportion are they affect the total number of cases
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = matrix(c(1:total_vaccine_brands, rep(total_vaccine_brands+1, total_vaccine_brands)), byrow = T, nrow = 2)
  if(calculate_relative_VE & total_vaccine_brands>1){
    relative_VE_combn = cbind(combn(total_vaccine_brands, 2), relative_VE_combn)
  }
  relative_VE_combn = rbind(relative_VE_combn[,rep(1:ncol(relative_VE_combn), each=total_case_strains)], rep(1:total_case_strains, total_case_strains))

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
  #P(strain & case | unvaccinated) p(strain | case, unvaccinated) * p(case|unvaccinated)
  prob_case_each_strain_given_unvaccinated = proportion_strains_in_unvaccinated_cases * prob_case_given_unvaccinated
  prob_case_each_strain_and_unvaccinated = prob_case_each_strain_given_unvaccinated * prob_unvaccinated

  #So now we are done with the last column of our table. the remaining columns are for vaccines
  #sum of columns 1:total_vaccine_brands
  prob_vaccinated_with_brands = brand_proportions_in_vaccinated * overall_vaccine_coverage

  relative_risk_for_each_brand_and_strain = 1 - do.call('cbind', anticipated_VE_for_each_brand_and_strain)
  risk_for_each_brand_and_strain_given_unvaccinated = prob_case_each_strain_and_unvaccinated/(prob_case_each_strain_and_unvaccinated + prob_control_and_unvaccinated)
  risk_for_each_brand_and_strain_given_vaccinated = relative_risk_for_each_brand_and_strain * risk_for_each_brand_and_strain_given_unvaccinated
  odds_for_each_brand_and_strain_given_vaccinated = (risk_for_each_brand_and_strain_given_vaccinated^-1 - 1)^-1
  prob_control_each_brand_given_vaccinated = prob_vaccinated_with_brands / (1 + colSums(odds_for_each_brand_and_strain_given_vaccinated))
  prob_case_each_strain_and_each_brand_given_vaccinated = t(odds_for_each_brand_and_strain_given_vaccinated) * prob_control_each_brand_given_vaccinated

  full_table = cbind(rbind(t(prob_case_each_strain_and_each_brand_given_vaccinated), 'controls'=prob_control_each_brand_given_vaccinated),
                     'unvaccinated'=c(prob_case_each_strain_and_unvaccinated, prob_control_and_unvaccinated))

  #first row is for vaccine 1, second row is for vaccine 2, third is for strain
  STRAIN_ROW = 3
  CONTROL_ROW=ncol(full_table)

  #c(full_table) converts the table into a vector, going column by column
  cell_counts_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_subjects, rmultinom, n=nsims, prob=c(full_table))),
                           dim = c((total_vaccine_brands+1)*(total_case_strains+1), nsims, total_total_subject_settings))

  #cell_counts_sims = rmultinom(n=nsims, size=missing_data_adjusted_total_subjects, prob = c(full_table))
  #Haldane's correction
  cell_counts_sims[cell_counts_sims==0] = 0.5

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  cases_vaccine1 = cell_counts_sims[(relative_VE_combn[1,]-1)*(total_case_strains +1) + relative_VE_combn[STRAIN_ROW,],,,drop=F]
  cases_vaccine2 = cell_counts_sims[(relative_VE_combn[2,]-1)*(total_case_strains +1) + relative_VE_combn[STRAIN_ROW,],,,drop=F]
  controls_vaccine1 = cell_counts_sims[(relative_VE_combn[1,]-1)*(total_case_strains +1) + CONTROL_ROW,,,drop=F]
  controls_vaccine2 = cell_counts_sims[(relative_VE_combn[2,]-1)*(total_case_strains +1) + CONTROL_ROW,,,drop=F]

  estimate = (cases_vaccine1/(cases_vaccine1 + controls_vaccine1))/(cases_vaccine2/(cases_vaccine2 + controls_vaccine2))
  standard_error = sqrt(1/cases_vaccine1 + 1/cases_vaccine2 - 1/(cases_vaccine1 + controls_vaccine1) - 1/(cases_vaccine2 + controls_vaccine2))
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  relative_risk_for_each_brand_and_strain = cbind(relative_risk_for_each_brand_and_strain,1)
  relative_risk_vaccine1 = relative_risk_for_each_brand_and_strain[(relative_VE_combn[1,]-1)*total_case_strains + relative_VE_combn[STRAIN_ROW,]]
  relative_risk_vaccine2 = relative_risk_for_each_brand_and_strain[(relative_VE_combn[2,]-1)*total_case_strains + relative_VE_combn[STRAIN_ROW,]]
  anticipated_VE = 1 - relative_risk_vaccine1/relative_risk_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Vaccine", relative_VE_combn[1,]),total_total_subject_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[2,]==total_vaccine_brands+1,
                                          no = paste("Vaccine", relative_VE_combn[2,]),
                                          yes = "Unvaccinated"), total_total_subject_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN_ROW,]), total_total_subject_settings),
                   total_subjects = rep(total_subjects, each=ncol(relative_VE_combn)),
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit))

  return(ret)
}
