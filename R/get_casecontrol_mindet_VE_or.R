#' Minimum detectable variant and vaccine-specific efficacy based on odds ratio.
#' @description
#' The function `get_cohort_mindet_VE_or` calculates the minimum detectable variant and vaccine-specific efficacy (VE) for a given sample size, power, and type-I error rate.
#' The efficacy is defined as `VE = 1 - odds ratio`, and the function returns
#' both absolute and relative minimum detectable VE. The NULL hypothesis is that `VE = 0` and
#' the alternate hypothesis in that `VE is not equal to 0`.
#'
#' @param anticipated_VE_for_each_brand_and_variant a matrix of vaccine efficacy of each vaccine (row) against each variant (column). Each value must be a real number between 0 and 1.
#' @param brand_proportions_in_vaccinated a vector denoting the proportion in which vaccines are given in the vaccinated subjects of the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_vaccine_coverage the proportion of the study cohort that will be vaccinated. It should be a real number between 0 and 1.
#' @param proportion_variants_in_unvaccinated_cases a vector of the proportions in which each variant is expected to be present in the unvaccinated and infected subjects in the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param controls_per_case the number of controls to be sampled per case in the study cohort. It should be a real number between 0 and Inf.
#' @param calculate_relative_VE a logical indicating if calculations should also be done for relative vaccine efficacy (default `TRUE`).
#' @param power the power to detect the VE. It is equal to 1 - Type-II error. It is a numeric that must take a value between 0 and 1.
#' @param alpha Type-I error probability. It is a numeric that must take a value between 0 and 1.
#' @param confounder_adjustment_Rsquared we use this parameter to adjust the calculations for potential confounders using the methodology proposed by Hsieh and Lavori (2000). It represents the amount of variance (R^2) explained in a regression model where vaccination status is the outcome and confounders of interest are predictors. It is a numeric that must take a value between 0 (no adjustment for confounders) and 1.
#' @param prob_missing_data to adjust the calculations for non-informative and random subject loss to follow-up/dropout. it should take a numeric value between 0 and 1.
#' @param total_cases a vector of the number of cases in our cohort for which calculations should be done.
#'
#' @details
#' In this function efficacy is defined as `VE = 1 - odds ratio`, where 'odds ratio' is
#' the ratio of odds of being vaccinated with a particular vaccine versus being unvaccinated among cases and the same odds in controls.
#' When the groups being compared are a particular vaccine versus placebo then we call the VE
#' as the absolute VE of the vaccine. For `M` vaccines there are `M` absolute VE, one each for the `M` vaccines.
#' When the groups being compared are a particular vaccine versus another vaccine then we call the VE
#' as the relative VE of the vaccines, for a particular variant. For `M` vaccines and `I` variants there are `I x 2 x utils::combn(M, 2)`
#' permutations of relative VE of two vaccines against the same variant.
#'
#' We first transform the user inputs for `I` variants and `M` vaccines into two conditional tables, one for cases, and another for controls.
#' First an `I x (M + 1)` cross table of the probability of being unvaccinated or vaccinated with a vaccine, given that the subject is a case. The sum of the cells of this cases table is equal to 1.
#' Second a `1 x (M + 1)` vector of the probability of being unvaccinated or vaccinated with a vaccine, given that the subject is a controls. The sum of this controls vector is equal to 1.
#' The first column in both of these corresponds to subjects who are unvaccinated and the remaining `M` columns correspond to subjects who are vaccinated with a particular vaccine.
#'
#' In general, while calculating minimum detectable VE it is not necessary to ask for `anticipated_VE_for_each_brand_and_variant`.
#' The reason we need it in the multiple variant and multiple vaccine scenario is explained next.
#'
#' After we obtain the cases table and controls vector our next step is to calculate the minimum detectable VE for each absolute and relative VE combination.
#' For absolute VE we extract a `1 x 2` vector of probabilities of being a control larger `1 x (M + 1)` controls vector. The two columns are for the vaccine of interest and placebo.
#' We also extract a `1 x 2` vector of probabilities of being a case from the larger `I x (M + 1)` controls vector.
#' The row is the row corresponding to the variant of interest and the two columns are for the vaccine of interest and placebo.
#' The effective vaccine coverage is then rescaled probability of being vaccinated in the selected `1 x 2` controls vector.
#' The effective number of cases are `total_cases * sub_cases_vector_probability_sum * (1 - prob_missing_data) * (1 - confounder_adjustment_Rsquared)`.
#' Here `sub_cases_vector_probability_sum` is the sum of the cells of the `1 x 2` sub vector of cases.
#' The effective number of controls are `total_cases * controls_per_case * sub_cases_vector_probability_sum * (1 - prob_missing_data) * (1 - confounder_adjustment_Rsquared)`.
#' Here `sub_cases_vector_probability_sum` is the sum of the cells of the `1 x 2` sub vector of controls.
#' We then simply pass these adjusted parameters to the function `epi.sscc` of the R package `epiR` to obtain the minimum detectable absolute VE.
#' For relative VE the process is similar to absolute VE, except that instead of placebo we select a `1 x 2` sub-vector where
#' the columns are two vaccines of interest.
#'
#' @return A data frame consisting of the input parameters, absolute and relative VE combinations,
#' and the minimum detectable VE for each absolute and relative VE combination.
#'
#' @examples As an example we recommend running the function without passing any parameter to it.
#' The default scenario is for three vaccines and three pathogen variants.
#'
#' @references
#' 1. Hsieh, F. Y., & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled clinical trials, 21(6), 552-560.
#' 2. Morris, J. A., & Gardner, M. J. (1988). Statistics in medicine: Calculating confidence intervals for relative risks (odds ratios) and standardised ratios and rates. British medical journal (Clinical research ed.), 296(6632), 1313.
#' @importFrom epiR epi.sscc
#' @importFrom utils combn
#' @export
get_casecontrol_mindet_VE_or = function(anticipated_VE_for_each_brand_and_variant=
                                          matrix(data=c(0.8, 0.5, 0.3, 0.3, 0.5, 0.8, 0.9, 0.5, 1), nrow = 3, ncol = 3, byrow = F,
                                                 dimnames = list(paste0('brand', 1:3), paste0('variant', 1:3))),
                                        brand_proportions_in_vaccinated =
                                          c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                        overall_vaccine_coverage=0.3,
                                        proportion_variants_in_unvaccinated_cases = c('variant1'=0.6, 'variant2'=0.3, 'variant3'=0.1),
                                        controls_per_case=2,
                                        calculate_relative_VE=T,
                                        power=0.8,
                                        alpha=0.05,
                                        confounder_adjustment_Rsquared = 0,
                                        prob_missing_data = 0.1,
                                        total_cases=seq(50,500, 10)){

  check_input(anticipated_VE_for_each_brand_and_variant, brand_proportions_in_vaccinated, proportion_variants_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_variants = length(proportion_variants_in_unvaccinated_cases)

  relative_VE_combn = get_vaccine_comparison_combinations(total_vaccine_brands, total_case_variants, calculate_relative_VE)

  total_total_case_settings = length(total_cases)
  missing_data_adjusted_total_cases = round(total_cases * (1-prob_missing_data))
  missing_data_adjusted_total_controls = missing_data_adjusted_total_cases * controls_per_case

  CONTROL_COL = 1
  UNVACCINATED_ROW = 1

  case_control_full_tables = get_case_control_full_tables(anticipated_VE_for_each_brand_and_variant, overall_vaccine_coverage,
                                                          proportion_variants_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  full_table = cbind(case_control_full_tables$conditional_vector_controls, case_control_full_tables$full_table_cases)

  mindet_VE = t(apply(relative_VE_combn, MARGIN = 2, FUN = function(comparison_set){
    browser()
    variant_index = comparison_set[VARIANT] + CONTROL_COL
    vaccine1_index = comparison_set[BRAND1] + UNVACCINATED_ROW
    vaccine2_index = comparison_set[BRAND2] + UNVACCINATED_ROW

    control_probs = full_table[c(vaccine1_index, vaccine2_index), CONTROL_COL]
    case_probs = full_table[c(vaccine1_index, vaccine2_index), variant_index]

    coverage_subpopulation = control_probs[1] / sum(control_probs)

    if(comparison_set[BRAND2] == 0){
      or_index = 1
    }else{
      anticipated_VEs = anticipated_VE_for_each_brand_and_variant[comparison_set[c(BRAND1, BRAND2)],
                                                                 comparison_set[VARIANT]]
      #irr index 1 is lower limit of incidence rate ratio and to be chosen when VE is 0 to 100%
      #irr index 2 is upper limit of incidence rate ratio and to be chosen when VE is between -100% and 0%
      or_index = ifelse(anticipated_VEs[2] > anticipated_VEs[1], 2, 1)
    }

    case_counts = missing_data_adjusted_total_cases * sum(case_probs) * (1-confounder_adjustment_Rsquared)
    control_counts = missing_data_adjusted_total_controls * sum(control_probs) * (1-confounder_adjustment_Rsquared)
    sapply(1:total_total_case_settings, FUN = function(i){
      ret = try(1 - epi.sscc(OR = NA,
                             p0 = coverage_subpopulation,
                             n =  control_counts[i] + case_counts[i],
                             power = power, r = control_counts[i]/case_counts[i],
                             sided.test = 2, conf.level = 1-alpha,
                             method = "unmatched", fleiss = FALSE)$OR[or_index], silent = T)
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
                   variant = rep(paste("Variant", relative_VE_combn[VARIANT,]), total_total_case_settings),
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
              t(data.frame(proportion_variants_in_unvaccinated_cases,
                           row.names = paste0('variant_prop_', 1:total_case_variants))))

  return(ret)
}

