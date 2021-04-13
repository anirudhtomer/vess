get_cohort_mindet_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                overall_brand_proportions = c(0.3, 0.5, 0.2),
                                overall_vaccine_coverage=0.3,
                                attack_rate_unvaccinated = 0.1,
                                calculate_for_relative_VE = T,
                                power=0.8,
                                alpha=0.05,
                                confounder_adjustment_Rsquared = 0,
                                prob_missing_data = 0.1,
                                total_subjects=500){

  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }else if(length(anticipated_brand_VEs)!=length(overall_brand_proportions)){
    stop("Length of anticipated brand VE should be equal to length of overall brand proportions")
  }

  total_vaccines = length(anticipated_brand_VEs)

  relative_VE_combn = matrix(c(1:total_vaccines, rep(total_vaccines+1, total_vaccines)), byrow = T, nrow = 2)
  if(calculate_for_relative_VE & total_vaccines>1){
    relative_VE_combn = cbind(combn(total_vaccines, 2), relative_VE_combn)
  }

  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  #the last index is for placebo (unvaccinated)
  coverages = c(overall_brand_proportions * overall_vaccine_coverage, 1-overall_vaccine_coverage)
  attack_rates = attack_rate_unvaccinated * (1 - c(anticipated_brand_VEs, 0))

  mindet_VE = apply(relative_VE_combn, MARGIN = 2, function(vaccines){
    group_coverage = sum(coverages[vaccines])
    subpopulation_coverage = coverages[vaccines[1]] / group_coverage

    #the 'n' parameter need not be integer for this API.
    1 - epi.sscohortc(irexp1 = NA, irexp0 = attack_rates[vaccines[2]],
                      n = missing_data_adjusted_total_subjects * group_coverage * (1-confounder_adjustment_Rsquared),
                      power = power,
                      r = subpopulation_coverage/(1-subpopulation_coverage),
                      design = 1, sided.test = 2, conf.level = 1-alpha)$irr[1]
  })

  anticipated_brand_VEs = c(anticipated_brand_VEs, 0)

  ret = data.frame(Vaccine1=paste("Vaccine", relative_VE_combn[1,]),
                   Vaccine2=ifelse(relative_VE_combn[2,]==total_vaccines+1, no = paste("Vaccine", relative_VE_combn[2,]), yes = "Unvaccinated"),
                   anticipated_VE = apply(relative_VE_combn, 2, FUN = function(x){
                     1 - (1-anticipated_brand_VEs[x[1]])/(1-anticipated_brand_VEs[x[2]])
                   }),
                   mindet_VE = mindet_VE)
  return(ret)
}
