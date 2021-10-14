#' Expected lower and expected upper limit of the confidence intervals of strain and vaccine-specific efficacy based on incidence rate ratio.
#' @description
#' The function `get_cohort_expectedCI_VE_irr` simulates confindence intervals for strain and vaccine-specific efficacy (VE) for a given sample size.
#' The efficacy is defined as `VE = 1 - incidence rate ratio`, and the function returns expected lower and expected upper confidence interval limit for both absolute and relative VE.
#'
#' @param anticipated_VE_for_each_brand_and_strain a matrix of vaccine efficacy of each vaccine (row)  against each strain (column). Each value must be a real number between 0 and 1.
#' @param brand_proportions_in_vaccinated a vector denoting the proportion in which vaccines are given in the vaccinated subjects of the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_vaccine_coverage the proportion of the study cohort that will be vaccinated. It should be a real number between 0 and 1.
#' @param proportion_strains_in_unvaccinated_cases a vector of the proportions in which each strain is expected to be present in the unvaccinated and infected subjects in the study cohort. Each value of this vector must be a real number between 0 and 1 and the sum of the values of this vector must be equal to 1.
#' @param overall_attack_rate_in_unvaccinated the proportion of the study cohort that is expected to infected over the study period. It should be a real number between 0 and 1.
#' @param study_period the study period (epochs) should be numeric value greater than 0.
#' @param calculate_relative_VE a logical indicating if calculations should also be done for relative vaccine efficacy (default `TRUE`).
#' @param alpha controls the width `[100*alpha/2, 100*(1 - alpha/2)]%` of the confidence interval. It is a numeric that must take a value between 0 and 1.
#' @param confounder_adjustment_Rsquared we use this parameter to adjust the calculations for potential confounders using the methodology proposed by Hsieh and Lavori (2000). It represents the amount of variance (R^2) explained in a regression model where vaccination status is the outcome and confounders of interest are predictors. It is a numeric that must take a value between 0 (no adjustment for confounders) and 1.
#' @param prob_missing_data to adjust the calculations for non-informative and random subject loss to follow-up/dropout. it should take a numeric value between 0 and 1.
#' @param total_subjects a vector of study cohort size for which calculations should be done.
#' @param nsims total number of Monte Carlo simulations conducted.
#'
#' @details
#' #' In this function efficacy is defined as `VE = 1 - incidence rate ratio`, where 'incidence rate ratio' is
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
#' and the remaining `M` columns correspond to subjects who are vaccinated with a particular vaccine.
#' The next step is to simulate the data. To speed up our computations we sample an `(I + 1) x (M + 1)`
#' cross table of data from a multinomial distribution with probabilities taken from our cumulative-incidence table.
#' The total subjects sampled in the cross table are are `total_subjects * (1 - prob_missing_data)`.
#'
#' In the person-time based estimator of incidence rates the number of cases only constitute the numerator.
#' The denominator is person-time contributed by subjects vaccinated with a particular vaccine (or placebo).
#' Hence, in addition to the table of cumulative-incidences, using the user input we also obtain a `I x (M + 1)` table of incidence rates.
#' We use this table to calculate the person-time contribution of subjects over `study_period`.
#' In this table each cell contains the incidence rate of being a case of a particular strain given their vaccination status.
#' The naive method to calculate the overall person-time contribution of subjects vaccinated with a particular vaccine is to
#' simulate event-time of each subject using an exponential distribution (R function rexp). Subsequently, one can add them up over the subjects.
#' But this is computationally very slow. Hence we instead exploit the central limit theorem to directly sample the
#' 'sum of the person-time of all subjects with a certain vaccination status'.
#' For this, assume that the `study_period=t` and the sum of the event rates
#' for `I` strains is `p = p1 + p2 + ... + pI` per unit time, where `p1, p2, ..., pI` are the individual
#' event rates of the `I` strains. Then `exp(-pt)%` subjects will not be infected with the strain
#' and contribute `t` units time each. The remaining `100 - exp(-pt)%` subjects will
#' obtain the event and individually contribute an event time whose sum can be directly sampled from a
#' normal distribution (central limit theorem). This is because the sum of event times of `K`
#' eventful subjects follows a normal distributon with mean equal to
#' `K x mean of a truncated exponential distribution in the interval` `[0, t]`
#' and variance equal to the `K x variance of a truncated exponential distribution in the interval` `[0, t]`.
#'
#' We then estimate the absolute and relative VE of each vaccine using the counts and person-times based on the sampled data.
#' For this purpose we also reduce the sample size by a factor `prob_missing_data * confounder_adjustment_Rsquared`.
#' The confidence intervals with widths `[100*alpha/2, 100*(1 - alpha/2)]%` are obtained using normal approximation
#' to the distribution of log of incidence rate ratio (Rothman KJ, 2012).
#' To adjust for confounders, the standard-error used in the confidence interval is rescaled to `SE/(1 - confounder_adjustment_Rsquared)` (Hsieh and Lavori, 2000)
#' We repeat this procedure `nsims` times, and in each such simulation we obtain `nsims` confidence intervals.
#'
#' To conduct simulations faster all the calculations are done without using for loops.
#' Instead we use a three-dimensional R arrays, with one-dimension for `nsims`,
#' another for the sample size vector `total_subjects`,
#' and another for a vector containing the flattened cross table of simulated data on cases and controls.
#'
#' @return A data frame consisting of the input parameters, absolute and relative VE combinations, and the expected lower and expected upper width of the confidence intervals for each absolute and relative VE.
#' @references
#' 1. Hsieh, F. Y., & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled clinical trials, 21(6), 552-560.
#' 2. Rothman KJ (2012) Epidemiology: An Introduction. 2nd Ed., Oxford University Press, Oxford.
#' @importFrom matrixStats colCumsums
#' @importFrom utils combn
#' @export
get_cohort_expectedCI_VE_irr = function(anticipated_VE_for_each_brand_and_strain=
                                              matrix(data=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1), nrow = 3, ncol = 3, byrow = F,
                                                     dimnames = list(paste0('brand', 1:3), paste0('strain', 1:3))),
                                            brand_proportions_in_vaccinated =
                                              c('brand1'=0.3, 'brand2'=0.5, 'brand3'=0.2),
                                            overall_vaccine_coverage=0.3,
                                            proportion_strains_in_unvaccinated_cases = c('strain1'=0.6, 'strain2'=0.3, 'strain3'=0.1),
                                            overall_attack_rate_in_unvaccinated = 0.1,
                                            study_period=1,
                                            calculate_relative_VE = T,
                                            alpha=0.05,
                                            confounder_adjustment_Rsquared = 0,
                                            prob_missing_data = 0.1,
                                            total_subjects=seq(1000, 10000, 25),
                                            nsims = 500){

  check_input(anticipated_VE_for_each_brand_and_strain, brand_proportions_in_vaccinated, proportion_strains_in_unvaccinated_cases)

  total_vaccine_brands = length(brand_proportions_in_vaccinated)
  total_case_strains = length(proportion_strains_in_unvaccinated_cases)

  relative_VE_combn = get_vaccine_comparison_combinations(total_vaccine_brands, total_case_strains, calculate_relative_VE)

  all_tables = get_cohort_full_table_irr(anticipated_VE_for_each_brand_and_strain, overall_vaccine_coverage,
                                         overall_attack_rate_in_unvaccinated, study_period,
                                         proportion_strains_in_unvaccinated_cases, brand_proportions_in_vaccinated)

  #c(full_table) converts the table into a vector, going column by column
  flat_full_table_prob = c(all_tables$full_table_prob)
  flat_overall_event_rates = colSums(all_tables$conditional_table_incidence_rate)

  total_total_subject_settings = length(total_subjects)
  missing_data_adjusted_total_subjects = round(total_subjects * (1-prob_missing_data))

  cell_counts_sims = array(data = do.call('cbind', lapply(missing_data_adjusted_total_subjects, rmultinom, n=nsims, prob=flat_full_table_prob)),
                           dim = c((1+total_vaccine_brands)*(1+total_case_strains), nsims, total_total_subject_settings))

  trunc_exp_means = trunc_exp_mean(event_rate = flat_overall_event_rates, trunc_limit = study_period)
  trunc_exp_vars = trunc_exp_var(event_rate = flat_overall_event_rates, trunc_limit = study_period)
  #this happens when VE = 100%, that is event rate is zero. consequently the person time becomes NaN
  #i want to set that NaN to 0 because 100% VE means no events and 0 person time contributed
  trunc_exp_means[is.nan(trunc_exp_means)] = 0
  trunc_exp_vars[is.nan(trunc_exp_vars)] = 0

  CONTROL_ROW = 1
  CONTROL_INDICES = CONTROL_ROW + (0:total_vaccine_brands)*(total_case_strains+1)
  #cell_counts_sims[CONTROL_INDICES,,,drop=F] is the counts per vaccine among controls
  person_time_sims = cell_counts_sims[CONTROL_INDICES,,,drop=F] * study_period

  #summarizing the counts per vaccine among cases. so want to add over all strains per vaccine
  cases_counts_sims = array(
    c(
      apply(cell_counts_sims[-CONTROL_INDICES,,,drop=F], c(3), function(x){
        cumsum_overall = rbind(0,colCumsums(x))
        a = seq(1, nrow(x), by = total_case_strains)
        b = seq(total_case_strains, nrow(x), by = total_case_strains)
        return(cumsum_overall[b + 1,] - cumsum_overall[a,])
      })
    ),
    dim = c(1+total_vaccine_brands, nsims, total_total_subject_settings)
  )

  cases_person_time_sum_mean = cases_counts_sims * trunc_exp_means
  cases_person_time_sum_sd = sqrt(cases_counts_sims * trunc_exp_vars)

  #the pmax is because rnorm can also sample negative person-time
  #think about event rate zero against a strain: person time should be study period
  person_time_sims = person_time_sims + array(
    data = pmax(
      rnorm(
        n = prod(dim(cases_person_time_sum_mean)),
        mean = c(cases_person_time_sum_mean),
        sd = c(cases_person_time_sum_sd)
      ),
      0
    ),
    dim = dim(person_time_sims)
  )

  cell_counts_sims[cell_counts_sims==0] = 0.5

  UNVACCINATED_EVENTS_ROW = 1

  brand1_total_events_indices = relative_VE_combn[BRAND1,]*(1 + total_case_strains) + relative_VE_combn[STRAIN,] + CONTROL_ROW
  brand1_person_time_indices = relative_VE_combn[BRAND1,] + UNVACCINATED_EVENTS_ROW
  brand2_total_events_indices = relative_VE_combn[BRAND2,]*(1 + total_case_strains) + relative_VE_combn[STRAIN,] + CONTROL_ROW
  brand2_person_time_indices = relative_VE_combn[BRAND2,] + UNVACCINATED_EVENTS_ROW

  #dimensions of the following are ncol(relative_VE_combn) x nsims
  brand1_total_events = cell_counts_sims[brand1_total_events_indices,,,drop=F]
  brand1_person_time = person_time_sims[brand1_person_time_indices,,,drop=F]
  brand2_total_events = cell_counts_sims[brand2_total_events_indices,,,drop=F]
  brand2_person_time = person_time_sims[brand2_person_time_indices,,,drop=F]

  estimate = (brand1_total_events/brand1_person_time)/(brand2_total_events/brand2_person_time)
  standard_error = sqrt(1/brand1_total_events + 1/brand2_total_events)
  confounder_adjusted_standard_error = standard_error / (1-confounder_adjustment_Rsquared)

  log_estimate = log(estimate)
  quantile_times_adjusted_standard_error = qnorm(1-alpha/2) * confounder_adjusted_standard_error
  low = exp(log_estimate - quantile_times_adjusted_standard_error)
  upp = exp(log_estimate + quantile_times_adjusted_standard_error)

  #expected_VE = rowMeans(1-estimate, na.rm=T)
  expected_VE = apply(X = 1-estimate, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_lower_limit = apply(X = 1-upp, MARGIN = 3, FUN = rowMeans, na.rm=T)
  avg_upper_limit = apply(X = 1-low, MARGIN = 3, FUN = rowMeans, na.rm=T)

  rate_table = rbind(NA, all_tables$conditional_table_incidence_rate)
  rate_vaccine1 = rate_table[brand1_total_events_indices]
  rate_vaccine2 = rate_table[brand2_total_events_indices]
  anticipated_VE = 1 - rate_vaccine1/rate_vaccine2

  ret = data.frame(vaccine_1 = rep(paste("Brand", relative_VE_combn[BRAND1,]), total_total_subject_settings),
                   vaccine_2 = rep(ifelse(relative_VE_combn[BRAND2,]==0,
                                          no = paste("Brand", relative_VE_combn[BRAND2,]),
                                          yes = "Unvaccinated"), total_total_subject_settings),
                   strain = rep(paste("Strain", relative_VE_combn[STRAIN,]), total_total_subject_settings),
                   total_subjects = rep(total_subjects, each=ncol(relative_VE_combn)),
                   anticipated_VE = anticipated_VE,
                   expected_VE = c(expected_VE),
                   bias_VE = c(expected_VE - anticipated_VE),
                   avg_lower_limit = c(avg_lower_limit),
                   avg_upper_limit = c(avg_upper_limit),
                   alpha = alpha,
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
