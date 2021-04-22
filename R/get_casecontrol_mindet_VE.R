#' @importFrom epiR epi.sscc
#' @importFrom utils combn
#' @export
get_casecontrol_mindet_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                     overall_brand_proportions = c(0.3, 0.5, 0.2),
                                     overall_vaccine_coverage=0.1,
                                     controls_per_case=2,
                                     calculate_for_relative_VE=T,
                                     power=0.8,
                                     alpha=0.05,
                                     confounder_adjustment_Rsquared = 0,
                                     prob_missing_data = 0.1,
                                     total_cases=500){
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

  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  total_controls = missing_data_adjusted_total_cases * controls_per_case

  cell_prob_control = c(brand_vaccine_coverages, 1-overall_vaccine_coverage)
  cell_count_controls = total_controls * cell_prob_control
  cell_prob_unvaccinated_case = (1 + sum(cell_count_controls[1:total_vaccines] * (1-anticipated_brand_VEs))/cell_count_controls[total_vaccines+1])^-1
  cell_prob_vaccinated_cases = (1-anticipated_brand_VEs) * cell_count_controls[1:total_vaccines] * cell_prob_unvaccinated_case / cell_count_controls[total_vaccines+1]
  cell_count_cases = missing_data_adjusted_total_cases * c(cell_prob_vaccinated_cases, cell_prob_unvaccinated_case)

  mindet_VE = apply(relative_VE_combn, MARGIN = 2, FUN = function(vaccines){
    cases = sum(cell_count_cases[vaccines]) * (1-confounder_adjustment_Rsquared)
    controls = sum(cell_count_controls[vaccines]) * (1-confounder_adjustment_Rsquared)

    ret = try(1 - epi.sscc(OR = NA,
                           p0 = cell_prob_control[vaccines[1]] / (cell_prob_control[vaccines[2]] + cell_prob_control[vaccines[1]]) ,
                           n =  controls + cases,
                           power = power, r = controls/cases,
                           sided.test = 2, conf.level = 1-alpha,
                           method = "unmatched", fleiss = FALSE)$OR[1], silent = T)
    if(inherits(ret, "try-error")){
      return(NA)
    }else{
      return(ret)
    }
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

