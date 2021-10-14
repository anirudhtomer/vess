#' Minimum detectable strain and vaccine-specific efficacy based on incidence rate ratio.
#' @description
#' The function `get_cohort_mindet_VE_irr` calculates the minimum detectable strain and vaccine-specific efficacy (VE) for a given sample size, power, and type-I error rate.
#' The efficacy is defined as `VE = 1 - incidence rate ratio`, and the function returns
#' both absolute and relative minimum detectable VE. The NULL hypothesis is that `VE = 0` and
#' the alternate hypothesis in that `VE is not equal to 0`.
#'
#' @param anticipated_VE_for_each_brand_and_strain a matrix of vaccine efficacy of each vaccine (row)  against each strain (column). Each value must be a real number between 0 and 1.
#' @param brand_proportions_in_vaccinated a vector denoting the proportion in which vaccines are given in the vaccinated subjects of the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_vaccine_coverage the proportion of the study cohort that will be vaccinated. It should be a real number between 0 and 1.
#' @param proportion_strains_in_unvaccinated_cases a vector of the proportions in which each strain is expected to be present in the unvaccinated and infected subjects in the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_attack_rate_in_unvaccinated the proportion of the study cohort that is expected to infected over the study period. It should be a real number between 0 and 1.
#' @param study_period the study period (epochs) should be numeric value greater than 0.
#' @param calculate_relative_VE a logical indicating if calculations should also be done for relative vaccine efficacy (default `TRUE`).
#' @param power the power to detect the VE. It is equal to 1 - Type-II error. It is a numeric that must take a value between 0 and 1.
#' @param alpha controls the width `[100*alpha/2, 100*(1 - alpha/2)]%` of the confidence interval. It is a numeric that must take a value between 0 and 1.
#' @param confounder_adjustment_Rsquared we use this parameter to adjust the calculations for potential confounders using the methodology proposed by Hsieh and Lavori (2000). It represents the amount of variance (R^2) explained in a regression model where vaccination status is the outcome and confounders of interest are predictors. It is a numeric that must take a value between 0 (no adjustment for confounders) and 1.
#' @param prob_missing_data to adjust the calculations for non-informative and random subject loss to follow-up/dropout. it should take a numeric value between 0 and 1.
#' @param total_subjects a vector of study cohort size for which calculations should be done.
#'
#' @details
#' In this function efficacy is defined as `VE = 1 - incidence rate ratio`, where 'incidence rate ratio' is
#' the ratio of incidence rates of being a case of a particular strain/variant among the groups being compared.
#' When the groups being compared are a particular vaccine versus placebo then we call the VE
#' as the absolute VE of the vaccine. For `M` vaccines there are `M` absolute VE, one each for the `M` vaccines.
#' When the groups being compared are a particular vaccine versus another vaccine then we call the VE
#' as the relative VE of the vaccines, for a particular strain. For `M` vaccines and `I` strains there are `I x 2 x utils::combn(M, 2)`
#' permutations of relative VE of two vaccines against the same strain.
#'
#' We first transform the user inputs for `I` strains and `M` vaccines into a `(I + 1) x (M + 1)` cross table of
#' cumulative-incidences of being a case or a control over the study period. The overall sum of all cumulative-incidences,
#' i.e., all cells, of this table is 1. The first row of our cumulative-incidence table contain cumulative-incidence of being a control.
#' The first column corresponds to subjects who are unvaccinated.
#' Thus, the cell `{1,1}` contains the probabitity (cumulative-incidence) that over the study period a subject will be a control and unvaccinated.
#' The remaining `Ì` rows correspond to subjects who are cases of a particular strain/variant of the pathogen,
#' and the remaining `M` columns correspond to subjects who are vaccinated with a particular vaccine. In the
#' person-time based estimator of incidence rates the number of cases only constitute the numerator.
#' The denominator is person-time contributed by subjects vaccinated with a particular vaccine (or placebo).
#' Hence, in addition to the table of cumulative-incidences, using the user input we also
#' obtain a `I x (M + 1)` table of incidence rates.
#'
#' In general, while calculating minimum detectable VE it is not necessary to ask for `anticipated_VE_for_each_brand_and_strain`.
#' The reason we need it in the multiple variant and mulitple vaccine scenario is explained next.
#'
#' After we obtain the cumulative incidence table our next step is to calculate the minimum detectable VE for each absolute and relative VE combination.
#' For absolute VE we extract a `2 x 2` sub-table of cumulative-incidences from the larger `(I + 1) x (M + 1)` cross table of
#' cumulative-incidences. The two rows of this table are for the strain of interest and controls.
#' The two columns are for the vaccine of interest and placebo. The effective sample size for this sub-table and the absolute VE is then
#' `total_subjects * sub_table_probability_sum * (1 - prob_missing_data) * (1 - confounder_adjustment_Rsquared)`.
#' Here `sub_table_probability_sum` is the sum of the cells of the sub_table.
#' This sum corresponds to the percentage of the `total_subjects` that can be used for a specific absolute VE calculation.
#' The effective overall coverage in this sub-table is the rescaled probability of being vaccinated in the `2 x 2` sub-table.
#' The effective incidence rate of the strain in this comparison is the strain-specific incidence rate of the strain taken
#' from the table of incidence rates.
#' We then simply pass these adjusted parameters to the function `epi.sscohortt` of the R package `epiR` to obtain the minimum detectable absolute VE.
#' For relative VE the process is similar to absolute VE, except that instead of placebo we select a `2 x 2` sub-table where
#' the columns are two vaccines of interest.
#'
#' @return A data frame consisting of the input parameters, absolute and relative VE combinations,
#' and the minimum detectable VE for each absolute and relative VE combination.
#'
#' @examples As an example we recommend running the function without passing any parameter to it.
#' The default scenario is for three vaccines and three pathogen strains.
#'
#' @references
#' 1. Hsieh, F. Y., & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled clinical trials, 21(6), 552-560.
#' 2. Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
#' @importFrom epiR epi.sscohortt
#' @importFrom utils combn
#' @export
get_cohort_mindet_VE_irr = function(anticipated_VE_for_each_brand_and_strain=
                                      matrix(data=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1), nrow = 3, ncol = 3, byrow = F,
                                             dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                    brand_proportions_in_vaccinated =
                                      c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                    overall_vaccine_coverage=0.3,
                                    proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                    overall_attack_rate_in_unvaccinated = 0.1,
                                    study_period=1,
                                    calculate_relative_VE = T,
                                    power=0.8,
                                    alpha=0.05,
                                    confounder_adjustment_Rsquared = 0,
                                    prob_missing_data = 0.1,
                                    total_subjects=seq(1000, 10000, 25)){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_vaccine_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  all_tables = get_cohort_full_table_irr(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                         overall_attack_rate_in_unvaccinated, study_period,
                                         proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  full_table_prob = all_tables$full_table_prob
  event_rates_table = all_tables$conditional_table_incidence_rate

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  CONTROL_ROW = 1
  UNVACCINATED_COL = 1

  mindet_VE = t(apply(relative_VE_combn, MARGIN = 2, function(comparison_set){
    sub_table = full_table_prob[c(0, comparison_set[STRAIN]) + CONTROL_ROW, c(comparison_set[BRAND1], comparison_set[BRAND2]) + UNVACCINATED_COL]
    event_rate_comparison_group = event_rates_table[comparison_set[STRAIN], comparison_set[BRAND2] + UNVACCINATED_COL]

    vaccine1_coverage = sum(sub_table[,1])
    group_coverage = sum(sub_table)
    subpopulation_coverage = vaccine1_coverage / group_coverage

    if(comparison_set[BRAND2] == 0){
      irr_index = 1
    }else{
      anticipated_VEs = anticipated_VE_for_each_brand_and_strain[comparison_set[c(BRAND1, BRAND2)],
                                                                 comparison_set[STRAIN]]
      #irr index 1 is lower limit of incidence rate ratio and to be chosen when VE is 0 to 100%
      #irr index 2 is upper limit of incidence rate ratio and to be chosen when VE is between -100% and 0%
      irr_index = ifelse(anticipated_VEs[2] > anticipated_VEs[1], 2, 1)
    }

    #the 'n' parameter need not be integer for this API.
    sapply(missing_data_adjusted_total_subjects * group_coverage * (1-confounder_adjustment_Rsquared),
           function(n){
             ret = try(1 - epi.sscohortt(irexp1 = NA, irexp0 = event_rate_comparison_group,
                                         FT = study_period, n = n, power = power,
                                         r = subpopulation_coverage/(1-subpopulation_coverage),
                                         design = 1, sided.test = 2, conf.level = 1-alpha
             )$irr[irr_index], silent = T)

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
                   overall_attack_rate_in_unvaccinated = overall_attack_rate_in_unvaccinated,
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
