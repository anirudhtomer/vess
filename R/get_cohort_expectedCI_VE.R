#' @importFrom utils combn
#' @export
get_cohort_expectedCI_VE = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                    overall_brand_proportions = c(0.3, 0.5, 0.2),
                                    overall_vaccine_coverage=0.3,
                                    attack_rate_unvaccinated = 0.1,
                                    calculate_for_relative_VE = T,
                                    alpha=0.05,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1,
                                    total_subjects=500,
                                    nsims = 500){

  if(!sum(overall_brand_proportions, na.rm = T)==1){
    stop("Sum of brand proportions should be equal to 1")
  }

  total_vaccines = length(anticipated_brand_VEs)

  relative_VE_combn = matrix(c(1:total_vaccines, rep(total_vaccines+1, total_vaccines)), byrow = T, nrow = 2)
  if(calculate_for_relative_VE){
    relative_VE_combn = cbind(combn(total_vaccines, 2), relative_VE_combn)
  }

  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  #first n-1 numbers are
  relative_risks = c(1 - anticipated_brand_VEs)
  prob_vaccinated_case = attack_rate_unvaccinated * brand_vaccine_coverages * relative_risks
  #made a dummy variable for ease of reading the code and the math
  relative_risk_placebo = 1
  prob_unvaccinated_case = attack_rate_unvaccinated * (1-overall_vaccine_coverage) * relative_risk_placebo
  prob_vaccinated_control = brand_vaccine_coverages * (1- attack_rate_unvaccinated*relative_risks)
  prob_unvaccinated_control = (1-attack_rate_unvaccinated) * (1-overall_vaccine_coverage)

  upp = matrix(ncol=ncol(relative_VE_combn), nrow = nsims, data = NA)
  low = matrix(ncol=ncol(relative_VE_combn), nrow = nsims, data = NA)

  CASE=1
  CONTROL=2
  for(j in 1:nsims){
    cell_counts = c(rmultinom(n=1, size=missing_data_adjusted_total_subjects,
                              prob = c(prob_vaccinated_case, prob_unvaccinated_case,
                                       prob_vaccinated_control, prob_unvaccinated_control)))
    #2 x n table, row 1 for cases, row 2 for controls. nth column for unvaccinated and other columns for different vaccines
    cell_counts = matrix(cell_counts, nrow=2, byrow = T)

    for(k in 1:ncol(relative_VE_combn)){
      v1index = relative_VE_combn[1, k]
      v2index = relative_VE_combn[2, k]

      estimate = (cell_counts[CASE, v1index]/sum(cell_counts[, v1index])) / (cell_counts[CASE, v2index]/sum(cell_counts[, v2index]))
      standard_error = sqrt(
        cell_counts[CASE, v1index]^-1 +
          cell_counts[CASE, v2index]^-1 -
          (cell_counts[CASE, v1index] + cell_counts[CONTROL, v1index])^-1 -
          (cell_counts[CASE, v2index] + cell_counts[CONTROL, v2index])^-1
      )
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
