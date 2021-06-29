#' @importFrom epiR epi.sscc
#' @importFrom utils combn
#' @export
get_casecontrol_mindet_VE_or = function(anticipated_VE_for_each_brand_and_strain=
                                          matrix(data=c(0.8, 0.5, 0.3, 0.3, 0.5, 0.8, 0.9, 0.5, 1), nrow = 3, ncol = 3, byrow = F,
                                                 dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                        brand_proportions_in_vaccinated =
                                          c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                        overall_vaccine_coverage=0.3,
                                        proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                        controls_per_case=2,
                                        calculate_relative_VE=T,
                                        power=0.8,
                                        alpha=0.05,
                                        confounder_adjustment_Rsquared = 0,
                                        prob_missing_data = 0.1,
                                        total_cases=seq(50,500, 10)){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  total_total_case_settings = length(total_cases)
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  missing_data_adjusted_total_controls = missing_data_adjusted_total_cases * controls_per_case

  CONTROL_COL = 1
  UNVACCINATED_ROW = 1

  case_control_full_tables = get_case_control_full_tables(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                                          proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  full_table = cbind(case_control_full_tables$full_table_cases, case_control_full_tables$full_vector_controls)

  mindet_VE = t(apply(relative_VE_combn, MARGIN = 2, FUN = function(comparison_set){
    strain_index = comparison_set[STRAIN] + CONTROL_COL
    vaccine1_index = comparison_set[BRAND1] + UNVACCINATED_ROW
    vaccine2_index = comparison_set[BRAND2] + UNVACCINATED_ROW

    control_probs = full_table[c(vaccine1_index, vaccine2_index), CONTROL_COL]
    case_probs = full_table[c(vaccine1_index, vaccine2_index), strain_index]

    coverage_subpopulation = control_probs[1] / sum(control_probs)

    case_counts = missing_data_adjusted_total_cases * sum(case_probs) * (1-confounder_adjustment_Rsquared)
    control_counts = missing_data_adjusted_total_controls * sum(control_probs) * (1-confounder_adjustment_Rsquared)
    sapply(1:total_total_case_settings, FUN = function(i){
      ret = try(1 - epi.sscc(OR = NA,
                             p0 = coverage_subpopulation,
                             n =  control_counts[i] + case_counts[i],
                             power = power, r = control_counts[i]/case_counts[i],
                             sided.test = 2, conf.level = 1-alpha,
                             method = "unmatched", fleiss = FALSE)$OR[1], silent = T)
      if(inherits(ret, "try-error")){
        return(NA)
      }else{
        return(ret)
      }
    })
  }))

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[BRAND1,]),total_total_case_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                          no = paste("Brand", relative_VE_combn[BRAND2,]),
                                          yes = "Unvaccinated"), total_total_case_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN,]), total_total_case_settings),
                   total_cases = rep(total_cases, each=ncol(relative_VE_combn)),
                   mindet_VE = c(mindet_VE),
                   alpha = alpha,
                   power = power,
                   overall_vaccine_coverage = overall_vaccine_coverage,
                   controls_per_case = controls_per_case,
                   confounder_adjustment_Rsquared = confounder_adjustment_Rsquared,
                   prob_missing_data = prob_missing_data
  )

  ret = cbind(ret,
              t(data.frame(brand_proportions_in_vaccinated,
                           row.names = paste0('brand_prop_', 1:total_vaccine_brands))),
              t(data.frame(proportion_strains_in_unvaccinated_cases,
                           row.names = paste0('strain_prop_', 1:total_case_strains))))

  return(ret)
}

