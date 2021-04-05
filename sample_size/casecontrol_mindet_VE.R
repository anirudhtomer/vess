source("sample_size/util.R")

#all parameters are vectors
#coverage: coverage in overall population
#total_cases: total subjects in the total population, not just for brand
casecontrol_mindet_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                 overall_brand_proportions = c(0.3, 0.5, 0.2),
                                 overall_vaccine_coverage=0.3, power=0.8, alpha=0.05, 
                                 attack_rate_unvaccinated = c(0.01),
                                 controls_per_case=2,
                                 confounder_adjustment_Rsquared = 0,
                                 prob_missing_data = 0.1, 
                                 total_cases=500){
  
  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }
  
  total_vaccines = length(anticipated_brand_VEs)
  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  total_controls = missing_data_adjusted_total_cases * controls_per_case
  
  cell_prob_control = c(brand_vaccine_coverages, 1-overall_vaccine_coverage)
  cell_count_controls = total_controls * cell_prob_control
  cell_prob_unvaccinated_case = (1 + sum(cell_count_controls[1:total_vaccines] * (1-anticipated_brand_VEs))/cell_count_controls[total_vaccines+1])^-1  
  cell_prob_vaccinated_cases = (1-anticipated_brand_VEs) * cell_count_controls[1:total_vaccines] * cell_prob_unvaccinated_case / cell_count_controls[total_vaccines+1]
  cell_count_cases = total_cases * c(cell_prob_vaccinated_cases, cell_prob_unvaccinated_case)
  
  mindet_VE = sapply(1:total_vaccines, FUN = function(vaccine_nr){
    cases = cell_count_cases[vaccine_nr] + cell_count_cases[total_vaccines+1]
    controls = cell_count_controls[vaccine_nr] + cell_count_controls[total_vaccines+1]
    
    epi.sscc(OR = NA, p0 = brand_vaccine_coverages[vaccine_nr], 
             n =  controls + cases,
             power = power, r = controls/cases,
             sided.test = 2, conf.level = 1-alpha, 
             method = "unmatched", fleiss = FALSE)$OR[1]
  })
  
  return(mindet_VE)
}