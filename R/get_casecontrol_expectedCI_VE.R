#' @importFrom utils combn
#' @export
get_casecontrol_expectedCI_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                         overall_brand_proportions = c(0.3, 0.5, 0.2),
                                         overall_vaccine_coverage=0.3,
                                         controls_per_case=2,
                                         calculate_for_relative_VE = T,
                                         alpha=0.05,
                                         confounder_adjustment_Rsquared = 0,
                                         prob_missing_data = 0.1,
                                         total_cases=500,
                                         nsims = 500){
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

  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  total_controls = missing_data_adjusted_total_cases * controls_per_case

  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  #assuming that coverage in controls is same as overall coverage (see top of this file)
  cell_prob_control = c(brand_vaccine_coverages, 1-overall_vaccine_coverage)

  upp = matrix(ncol=ncol(relative_VE_combn), nrow = nsims, data = NA)
  low = matrix(ncol=ncol(relative_VE_combn), nrow = nsims, data = NA)

  for(j in 1:nsims){

    #first n-1 cells are of different vaccines, the n-th cell is of unvaccinated people
    cell_count_controls = c(rmultinom(n = 1, size = total_controls, prob = cell_prob_control))

    cell_prob_unvaccinated_case = (1 + sum(cell_count_controls[1:total_vaccines] * (1-anticipated_brand_VEs))/cell_count_controls[total_vaccines+1])^-1
    cell_prob_vaccinated_cases = (1-anticipated_brand_VEs) * cell_count_controls[1:total_vaccines] * cell_prob_unvaccinated_case / cell_count_controls[total_vaccines+1]
    cell_count_cases = c(rmultinom(n = 1, size = missing_data_adjusted_total_cases, prob = c(cell_prob_vaccinated_cases, cell_prob_unvaccinated_case)))

    for(k in 1:ncol(relative_VE_combn)){
      v1index = relative_VE_combn[1, k]
      v2index = relative_VE_combn[2, k]

      estimate = (cell_count_cases[v1index]/cell_count_cases[v2index]) / (cell_count_controls[v1index]/cell_count_controls[v2index])
      standard_error = sqrt(cell_count_cases[v1index]^-1 + cell_count_cases[v2index]^-1 + cell_count_controls[v1index]^-1 + cell_count_controls[v2index]^-1)
      confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)
      low[j,k] = exp(log(estimate) + qnorm(alpha/2) * confounder_adjusted_standard_error)
      upp[j,k] = exp(log(estimate) + qnorm(1 - alpha/2) * confounder_adjusted_standard_error)
    }
  }

  anticipated_brand_VEs = c(anticipated_brand_VEs, 0)

  ret = data.frame(Vaccine1=paste("Vaccine", relative_VE_combn[1,]),
                   Vaccine2=ifelse(relative_VE_combn[2,]==total_vaccines+1, no = paste("Vaccine", relative_VE_combn[2,]), yes = "Unvaccinated"),
                   anticipated_VE = apply(relative_VE_combn, 2, FUN = function(x){
                     1 - (1-anticipated_brand_VEs[x[1]])/(1-anticipated_brand_VEs[x[2]])
                   }),
                   avg_lower_limit = apply(1-upp, MARGIN = 2, FUN = mean, na.rm=T),
                   avg_upper_limit = apply(1-low, MARGIN = 2, FUN = mean, na.rm=T))

  return(ret)
}
