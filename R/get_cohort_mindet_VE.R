#' @importFrom epiR epi.sscc
#' @importFrom utils combn
#' @export
get_cohort_mindet_VE = function(anticipated_VE_for_each_brand_and_strain=
                                  matrix(data=c(0.8, 0.5, 0.3, 0.3, 0.5, 0.8, 0.9, 0.5, 1), nrow = 3, ncol = 3, byrow = F,
                                         dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                brand_proportions_in_vaccinated =
                                  c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                overall_vaccine_coverage=0.3,
                                proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                overall_attack_rate_in_unvaccinated = 0.1,
                                calculate_relative_VE = T,
                                power=0.8,
                                alpha=0.05,
                                confounder_adjustment_Rsquared = 0,
                                prob_missing_data = 0.1,
                                total_subjects=seq(1000, 10000, 25)){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  full_table = get_cohort_full_table(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage, overall_attack_rate_in_unvaccinated,
                                     proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  CONTROL_ROW = 1
  UNVACCINATED_COL = 1

  mindet_VE = t(apply(relative_VE_combn, MARGIN = 2, function(comparison_set){
    sub_table = full_table[c(0, comparison_set[STRAIN]) + CONTROL_ROW, c(comparison_set[BRAND1], comparison_set[BRAND2]) + UNVACCINATED_COL]

    vaccine1_coverage = sum(sub_table[,1])
    group_coverage = sum(sub_table)
    subpopulation_coverage = vaccine1_coverage / group_coverage

    #the 'n' parameter need not be integer for this API.
    sapply(missing_data_adjusted_total_subjects * group_coverage * (1-confounder_adjustment_Rsquared),
           function(n){
             ret = try(1 - epi.sscohortc(irexp1 = NA, irexp0 = sub_table[2, 2]/sum(sub_table[,2]),
                                         n = n,
                                         power = power,
                                         r = subpopulation_coverage/(1-subpopulation_coverage),
                                         design = 1, sided.test = 2, conf.level = 1-alpha)$irr[1], silent = T)
             if(inherits(ret, "try-error")){
               return(NA)
             }else{
               return(ret)
             }
           })
  }))

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[BRAND1,]),total_total_subject_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                          no = paste("Brand", relative_VE_combn[BRAND2,]),
                                          yes = "Unvaccinated"), total_total_subject_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN,]), total_total_subject_settings),
                   total_subjects = rep(total_subjects, each=ncol(relative_VE_combn)),
                   mindet_VE = c(mindet_VE),
                   alpha = alpha,
                   power = power,
                   overall_vaccine_coverage = overall_vaccine_coverage,
                   vaccine_1_prop = rep(brand_proportions_in_vaccinated[relative_VE_combn[BRAND1,]], total_total_subject_settings),
                   vaccine_2_prop = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                               no = relative_VE_combn[BRAND2,],
                                               yes = NA), total_total_subject_settings),
                   strain_prop = rep(proportion_strains_in_unvaccinated_cases[relative_VE_combn[STRAIN,]], total_total_subject_settings),
                   overall_attack_rate_in_unvaccinated = overall_attack_rate_in_unvaccinated,
                   confounder_adjustment_Rsquared = confounder_adjustment_Rsquared,
                   prob_missing_data = prob_missing_data
  )
  return(ret)
}
