source("sample_size/util.R")
#first something important regarding vaccine coverage in controls
#P(Vaccination) = P(Vaccination | Case) P(Case) + P(Vaccination | Control) (1-P(Case))
#In general P(Case) is really low for diseases like influenza or even dengue. hardly 1-3% in the population
#Hence overall_vaccine_coverage is almost equal to vaccine coverage in controls

#all parameters are scalars except anticipated_brand_VEs and overall_brand_proportions
#nsims: larger value give more accurate results
#Important
#we have a n x 2 table for this situation, with n-1 rows for n-1 vaccines and 1 row for unvaccinated people
#the 2 columns are for cases and controls. our aim is to obtain a cell count for each cell
#there is a controls:case ratio which is selected on a population level and not on vaccine level

casecontrol_precision_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                    overall_brand_proportions = c(0.3, 0.5, 0.2),
                                    overall_vaccine_coverage=0.3, 
                                    alpha=0.05, 
                                    controls_per_case=2,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1, 
                                    total_cases=500,
                                    nsims = 500){
  
  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }
  
  total_vaccines = length(anticipated_brand_VEs)
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  total_controls = missing_data_adjusted_total_cases * controls_per_case
  
  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  #assuming that coverage in controls is same as overall coverage (see top of this file)
  cell_prob_control = c(brand_vaccine_coverages, 1-overall_vaccine_coverage)
  
  upp = matrix(ncol=total_vaccines, nrow = nsims, data = NA)
  low = matrix(ncol=total_vaccines, nrow = nsims, data = NA)
  
  for(j in 1:nsims){
    
    #first n-1 cells are of different vaccines, the n-th cell is of unvaccinated people
    cell_count_controls = c(rmultinom(n = 1, size = total_controls, prob = cell_prob_control))
    
    cell_prob_unvaccinated_case = (1 + sum(cell_count_controls[1:total_vaccines] * (1-anticipated_brand_VEs))/cell_count_controls[total_vaccines+1])^-1  
    cell_prob_vaccinated_cases = (1-anticipated_brand_VEs) * cell_count_controls[1:total_vaccines] * cell_prob_unvaccinated_case / cell_count_controls[total_vaccines+1]
    
    cell_count_cases = c(rmultinom(n = 1, size = missing_data_adjusted_total_cases, prob = c(cell_prob_vaccinated_cases, cell_prob_unvaccinated_case)))
    
    estimates = (cell_count_cases[1:total_vaccines]/cell_count_cases[total_vaccines+1]) / (cell_count_controls[1:total_vaccines]/cell_count_controls[total_vaccines+1])
    standard_errors = sqrt(cell_count_cases[1:total_vaccines]^-1 + cell_count_cases[total_vaccines+1]^-1 + cell_count_controls[1:total_vaccines]^-1 + cell_count_controls[total_vaccines+1]^-1)
    confounder_adjusted_standard_errors = standard_errors / (1-confounder_adjustment_Rsquared)
    low[j,] = exp(log(estimates) + qnorm(alpha/2) * confounder_adjusted_standard_errors)
    upp[j,] = exp(log(estimates) + qnorm(1 - alpha/2) * confounder_adjusted_standard_errors)
  }
  
  avg_lower_limit_VE = apply(1-upp, MARGIN = 2, FUN = mean, na.rm=T)
  avg_upper_limit_VE = apply(1-low, MARGIN = 2, FUN = mean, na.rm=T)
  
  return(list(avg_lower_limit_VE=avg_lower_limit_VE, avg_upper_limit_VE=avg_upper_limit_VE))
}
