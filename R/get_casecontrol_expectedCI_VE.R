#' @importFrom utils combn
#' @export
get_casecontrol_expectedCI_VE = function(anticipated_VE_for_each_brand_and_strain=
                                           matrix(data=c(0.8, 0.5, 0.3, 0.3, 0.5, 0.8, 0.9, 0.5, 1), nrow = 3, ncol = 3, byrow = F,
                                                  dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                         brand_proportions_in_vaccinated =
                                           c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                         overall_vaccine_coverage=0.3,
                                         proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                         controls_per_case=2,
                                         calculate_relative_VE = T,
                                         alpha=0.05,
                                         confounder_adjustment_Rsquared = 0,
                                         prob_missing_data = 0.1,
                                         total_cases=seq(50,500, 10),
                                         nsims = 500){

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

  total_total_case_settings = length(total_cases)

  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  missing_data_adjusted_total_controls = missing_data_adjusted_total_cases * controls_per_case

  #For each vaccine we know in what proportion are they given in the general population
  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  #For each strain we know in what proportion are they affect the total number of cases
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = matrix(c(1:total_vaccine_brands, rep(total_vaccine_brands+1, total_vaccine_brands)), byrow = T, nrow = 2)
  if(calculate_relative_VE & total_vaccine_brands>1){
    relative_VE_combn = cbind(combn(total_vaccine_brands, 2), relative_VE_combn)
  }
  relative_VE_combn = rbind(relative_VE_combn[,rep(1:ncol(relative_VE_combn), each=total_case_strains), drop=F],
                            rep(1:total_case_strains, ncol(relative_VE_combn)))

  #We are going to create a table with dimensions (row x columns) = (total_vaccine_brands + 1) x (total_case_strains + 1)
  #Note that this table is the transpose of what we did in the cohort study.
  #The extra 1 column is for unvaccinated
  #The extra 1 row is for controls

  prob_vaccinated_each_brand = brand_proportions_in_vaccinated * overall_vaccine_coverage
  #P(Vaccination) = P(Vaccination | Case) P(Case) + P(Vaccination | Control) (1-P(Case))
  #In general P(Case) is really low for diseases like influenza or even dengue. hardly 1-3% in the population
  #Hence overall_vaccine_coverage is almost equal to vaccine coverage in controls
  prob_vaccinated_each_brand_given_control = prob_vaccinated_each_brand
  prob_unvaccinated_given_control = 1 - overall_vaccine_coverage

  STRAIN_ROW = 3
  UNVACCINATED_ROW = total_vaccine_brands+1
  #the idea in a case-control study is that we generate two sets of multinomial data
  #1. for controls
  #2. for cases
  #the counts remain constant

  odds_vaccinated_for_each_brand_and_strain = 1 - anticipated_VE_for_each_brand_and_strain
  prob_case_and_unvaccinated = sum(proportion_strains_in_unvaccinated_cases * (1 + colSums(odds_vaccinated_for_each_brand_and_strain * prob_vaccinated_each_brand_given_control)/prob_unvaccinated_given_control))^-1
  prob_unvaccinated_case_each_strain = prob_case_and_unvaccinated * proportion_strains_in_unvaccinated_cases
  prob_vaccinated_case_each_strain_and_brand = t(t(odds_vaccinated_for_each_brand_and_strain * prob_vaccinated_each_brand_given_control / prob_unvaccinated_given_control) * prob_unvaccinated_case_each_strain)

  full_table_cases = rbind(prob_vaccinated_case_each_strain_and_brand, 'unvaccinated'=prob_unvaccinated_case_each_strain)

  ##c(full_table_cases) converts the table into a vector, going column by column
  cell_counts_case_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_cases, rmultinom, n=nsims, prob=c(full_table_cases))),
                                dim = c((total_vaccine_brands+1)*total_case_strains, nsims, total_total_case_settings))
  #Haldane's correction
  cell_counts_case_sims[cell_counts_case_sims==0] = 0.5

  #cell_count_control_sims has dimension total_vaccine_brands x nsims
  cell_counts_control_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_controls, rmultinom, n=nsims, prob=c(prob_vaccinated_each_brand_given_control, prob_unvaccinated_given_control))),
                                   dim = c(total_vaccine_brands+1, nsims, total_total_case_settings))
  #Haldane's correction
  cell_counts_control_sims[cell_counts_control_sims==0] = 0.5

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  cases_vaccine1 = cell_counts_case_sims[(relative_VE_combn[STRAIN_ROW,]-1)*(total_vaccine_brands+1) + relative_VE_combn[1,],,,drop=F]
  cases_vaccine2 = cell_counts_case_sims[(relative_VE_combn[STRAIN_ROW,]-1)*(total_vaccine_brands+1) + relative_VE_combn[2,],,,drop=F]
  controls_vaccine1 = cell_counts_control_sims[relative_VE_combn[1,],,,drop=F]
  controls_vaccine2 = cell_counts_control_sims[relative_VE_combn[2,],,,drop=F]

  estimate = (cases_vaccine1/cases_vaccine2) / (controls_vaccine1/controls_vaccine2)
  standard_error = sqrt(1/cases_vaccine1 + 1/cases_vaccine2 + 1/controls_vaccine1 + 1/controls_vaccine2)
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  odds_vaccinated_for_each_brand_and_strain = rbind(odds_vaccinated_for_each_brand_and_strain,1)
  odds_vaccinated_vaccine1 = odds_vaccinated_for_each_brand_and_strain[(relative_VE_combn[STRAIN_ROW,]-1)*(total_vaccine_brands+1) + relative_VE_combn[1,]]
  odds_vaccinated_vaccine2 = odds_vaccinated_for_each_brand_and_strain[(relative_VE_combn[STRAIN_ROW,]-1)*(total_vaccine_brands+1) + relative_VE_combn[2,]]
  anticipated_VE = 1 - odds_vaccinated_vaccine1/odds_vaccinated_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Vaccine", relative_VE_combn[1,]),total_total_case_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[2,]==total_vaccine_brands+1,
                                          no = paste("Vaccine", relative_VE_combn[2,]),
                                          yes = "Unvaccinated"), total_total_case_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN_ROW,]), total_total_case_settings),
                   total_cases = rep(total_cases, each=ncol(relative_VE_combn)),
                   controls_per_case = controls_per_case,
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit))

  return(ret)
}
