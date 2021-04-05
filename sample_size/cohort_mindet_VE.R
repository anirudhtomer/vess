source("sample_size/util.R")

#all parameters are scalars except overall_brand_proportions for different vaccines
#coverage: coverage in overall population
#nsims: larger value give more accurate results
cohort_mindet_VE = function(overall_brand_proportions = c(0.3, 0.5, 0.2),
                            overall_vaccine_coverage=0.3, 
                            attack_rate_unvaccinated = 0.1,
                            poweer=0.8,
                            alpha=0.05, 
                            confounder_adjustment_Rsquared = 0,
                            prob_missing_data = 0.1, 
                            total_subjects=500){
  
  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }
  
  total_vaccines = length(anticipated_brand_VEs)
  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))
  
  mindet_VE = sapply(brand_vaccine_coverages, FUN = function(brand_coverage){
    brand_and_unvaccinated_coverage = brand_coverage + 1-overall_vaccine_coverage
    brand_subpopulation_coverage = brand_coverage / brand_and_unvaccinated_coverage
    total_subpopulation_subjects = missing_data_adjusted_total_subjects * brand_and_unvaccinated_coverage
    
    #the 'n' parameter need not be integer for this API.
    epi.sscohortc(irexp1 = NA, irexp0 = attack_rate_unvaccinated, 
                  n = total_subpopulation_subjects * (1-confounder_adjustment_Rsquared),
                  power = power(), 
                  r = brand_subpopulation_coverage/(1-brand_subpopulation_coverage), 
                  design = 1, sided.test = 2, conf.level = 1-alpha)$irr[1]  
  })
  
  return(mindet_VE)
}