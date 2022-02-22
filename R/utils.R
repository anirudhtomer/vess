# smoothAndRound = function(dat, xname){
#   fit = stats::smooth.spline(x = dat[,xname], y = dat$VE_lower, df = SPLINE_DF)
#   return(round(pmax(0, stats:::predict.smooth.spline(object = fit, x = dat[,xname])$y), 2))
# }

BRAND1 = 'Brand 1'
BRAND2 = 'Brand 2'
VARIANT = 'Variant'
BRAND = 'Brand'
VARIANT1 = 'Variant 1'
VARIANT2 = 'Variant 2'
CONTROLS = 'controls'
UNVACCINATED = 'unvaccinated'

# itexp <- function(u, rate, t) { -log(1-u*(1-exp(-t*rate)))/rate }
# rtexp <- function(n, rate, t) { itexp(runif(n), rate, t) }
trunc_exp_mean = function(event_rate, trunc_limit){
  lambda = 1/event_rate
  k=trunc_limit/lambda

  return(lambda * (1-(k+1)*exp(-k))/(1-exp(-k)))
}

trunc_exp_var = function(event_rate, trunc_limit){
  lambda = 1/event_rate
  k=trunc_limit/lambda

  ret = (2 * lambda^2)/(1 - exp(-k))
  ret = ret * (1 - exp(-k)*(k^2 + 2*k + 2)*0.5)

  return(ret - trunc_exp_mean(event_rate, trunc_limit)^2)
}

get_ili_sari_symptom_prob=function(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length){
  if(is.na(ili_sari_symptom_prob)){
    if(!is.na(ili_sari_symptom_incidence_rate) & !is.na(study_period_length)){
      rate = expand.grid(ili_sari_symptom_incidence_rate, study_period_length)
      return(1 - exp(-rate[,1] * rate[,2]))
    }else{
      stop("Non NA values required for either 'ili_sari_symptom_prob' or for both 'ili_sari_symptom_incidence_rate' and 'study_period_length'")
    }
  }else{
    return(ili_sari_symptom_prob)
  }
}

check_input = function(anticipated_VE_for_each_brand_and_variant, brand_proportions_in_vaccinated, proportion_variants_in_unvaccinated_cases){
  if(!sum(brand_proportions_in_vaccinated, na.rm = T)==1){
    stop("Sum of the values of the parameter 'brand_proportions_in_vaccinated' should be equal to 1")
  }

  if(!sum(proportion_variants_in_unvaccinated_cases, na.rm = T)==1){
    stop("Sum of the values of the parameter 'proportion_variants_in_unvaccinated_cases' should be equal to 1")
  }

  if(nrow(anticipated_VE_for_each_brand_and_variant)!=length(brand_proportions_in_vaccinated)){
    stop("Total number of rows of the parameter 'anticipated_VE_for_each_brand_and_variant' should be equal to length of the parameter 'brand_proportions_in_vaccinated'")
  }

  if(ncol(anticipated_VE_for_each_brand_and_variant)!=length(proportion_variants_in_unvaccinated_cases)){
    stop("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_variant' should be equal to length of the parameter 'proportion_variants_in_unvaccinated_cases'")
  }
}

get_vaccine_comparison_combinations = function(total_vaccine_brands, total_case_variants, calculate_relative_VE=F){
  relative_VE_combn = matrix(c(1:total_vaccine_brands, rep(0, total_vaccine_brands)), byrow = T, nrow = 2)
  if(calculate_relative_VE & total_vaccine_brands>1){
    combinations = combn(total_vaccine_brands, 2)
    relative_VE_combn = cbind(relative_VE_combn, combinations, combinations[c(2,1),])
  }
  relative_VE_combn = rbind(relative_VE_combn[,rep(1:ncol(relative_VE_combn), each=total_case_variants), drop=F],
                            rep(1:total_case_variants, ncol(relative_VE_combn)))

  rownames(relative_VE_combn) = c(BRAND1, BRAND2, VARIANT)
  colnames(relative_VE_combn) = paste('Combination' , 1:ncol(relative_VE_combn))
  return(relative_VE_combn)
}

get_variant_comparison_combinations = function(total_vaccine_brands, total_case_variants){
  if(total_case_variants>1){
    relative_VE_combn = combn(total_case_variants, 2)
    relative_VE_combn = rbind(relative_VE_combn[,rep(1:ncol(relative_VE_combn), each=total_vaccine_brands), drop=F],
                              rep(1:total_vaccine_brands, ncol(relative_VE_combn)))

    rownames(relative_VE_combn) = c(VARIANT1, VARIANT2, BRAND)
    relative_VE_combn = cbind(relative_VE_combn, relative_VE_combn[c(2,1,3),])
    colnames(relative_VE_combn) = paste('Combination' , 1:ncol(relative_VE_combn))
  }else{
    relative_VE_combn = NULL
  }
  return(relative_VE_combn)
}

get_cohort_full_table_crr = function(anticipated_VE_for_each_brand_and_variant, overall_vaccine_coverage,
                                     overall_attack_rate_in_unvaccinated,
                                     proportion_variants_in_unvaccinated_cases, brand_proportions_in_vaccinated){
  #We are going to create a table with dimensions (row x columns) = (total_case_variants + 1) x (total_vaccine_brands + 1)
  #The extra 1 column is for unvaccinated
  #The extra 1 row is for controls

  p_v_0 = 1 - overall_vaccine_coverage
  p_c_any_given_v_0 = overall_attack_rate_in_unvaccinated
  p_c_0_given_v_0 = 1 - p_c_any_given_v_0
  p_c_i_given_v_0 = proportion_variants_in_unvaccinated_cases * p_c_any_given_v_0
  p_v_m = brand_proportions_in_vaccinated * overall_vaccine_coverage

  p_c_i_given_v_m = (1 - t(anticipated_VE_for_each_brand_and_variant)) * p_c_i_given_v_0
  p_c_0_given_v_m = 1 -colSums(p_c_i_given_v_m)

  conditional_table = rbind('controls'=p_c_0_given_v_m, p_c_i_given_v_m)
  conditional_table = cbind('unvaccinated'=c('controls'=p_c_0_given_v_0, p_c_i_given_v_0), conditional_table)

  full_table = t(t(conditional_table) * c(p_v_0, p_v_m))

  return(list(full_table=full_table, conditional_table = conditional_table))
}

#irr is incidence rate ratio
get_cohort_full_table_irr = function(anticipated_VE_for_each_brand_and_variant, overall_vaccine_coverage,
                                     overall_attack_rate_in_unvaccinated, study_period,
                                     proportion_variants_in_unvaccinated_cases, brand_proportions_in_vaccinated){

  #read variant as cause, and hazard as incidence rate
  prob_unvaccinated = 1 - overall_vaccine_coverage
  prob_case_given_unvaccinated = overall_attack_rate_in_unvaccinated
  prob_control_given_unvaccinated = 1 - prob_case_given_unvaccinated

  overall_incidence_rate_given_unvaccinated = (- log(1 - prob_case_given_unvaccinated) / study_period)
  variant_specific_incidence_rate_given_unvaccinated = overall_incidence_rate_given_unvaccinated * proportion_variants_in_unvaccinated_cases
  prob_case_each_variant_given_unvaccinated = prob_case_given_unvaccinated * proportion_variants_in_unvaccinated_cases

  relative_rate_ratio_for_each_brand_and_variant = 1 - t(anticipated_VE_for_each_brand_and_variant)

  variant_specific_incidence_rate_given_vaccinated_with_brand = relative_rate_ratio_for_each_brand_and_variant * variant_specific_incidence_rate_given_unvaccinated
  overall_incidence_rate_given_vaccinated_with_brand = colSums(variant_specific_incidence_rate_given_vaccinated_with_brand)
  prob_case_given_vaccinated_with_brand = 1 - exp(-overall_incidence_rate_given_vaccinated_with_brand * study_period)
  prob_control_given_vaccinated_with_brand = 1 - prob_case_given_vaccinated_with_brand
  prob_case_each_variant_given_vaccinated_with_brand = (t(variant_specific_incidence_rate_given_vaccinated_with_brand) / overall_incidence_rate_given_vaccinated_with_brand) * prob_case_given_vaccinated_with_brand

  conditional_prob_table = rbind('unvaccinated'=c('controls'=prob_control_given_unvaccinated, prob_case_each_variant_given_unvaccinated),
                                 cbind('controls'=prob_control_given_vaccinated_with_brand, prob_case_each_variant_given_vaccinated_with_brand))
  prob_vaccinated_with_brands = brand_proportions_in_vaccinated * overall_vaccine_coverage

  full_table_prob = t(conditional_prob_table * c(prob_unvaccinated, prob_vaccinated_with_brands))
  conditional_table_incidence_rate = cbind('unvaccinated'=variant_specific_incidence_rate_given_unvaccinated,
                                           variant_specific_incidence_rate_given_vaccinated_with_brand)

  return(list(full_table_prob=full_table_prob,
              conditional_table_incidence_rate=conditional_table_incidence_rate))
}

get_case_control_full_tables = function(anticipated_VE_for_each_brand_and_variant, overall_vaccine_coverage,
                                        proportion_variants_in_unvaccinated_cases, brand_proportions_in_vaccinated){

  prob_vaccinated_each_brand = brand_proportions_in_vaccinated * overall_vaccine_coverage
  #P(Vaccination) = P(Vaccination | Case) P(Case) + P(Vaccination | Control) (1-P(Case))
  #In general P(Case) is really low for diseases like influenza or even dengue. hardly 1-3% in the population
  #Hence overall_vaccine_coverage is almost equal to vaccine coverage in controls
  prob_vaccinated_each_brand_given_control = prob_vaccinated_each_brand
  prob_unvaccinated_given_control = 1 - overall_vaccine_coverage

  #the idea in a case-control study is that we generate two sets of multinomial data
  #1. for controls
  #2. for cases
  #the counts remain constant

  odds_vaccinated_for_each_brand_and_variant = 1 - anticipated_VE_for_each_brand_and_variant
  prob_case_and_unvaccinated = sum(proportion_variants_in_unvaccinated_cases * (1 + colSums(odds_vaccinated_for_each_brand_and_variant * prob_vaccinated_each_brand_given_control)/prob_unvaccinated_given_control))^-1
  prob_unvaccinated_case_each_variant = prob_case_and_unvaccinated * proportion_variants_in_unvaccinated_cases
  prob_vaccinated_case_each_variant_and_brand = t(t(odds_vaccinated_for_each_brand_and_variant * prob_vaccinated_each_brand_given_control / prob_unvaccinated_given_control) * prob_unvaccinated_case_each_variant)

  full_table_cases=rbind('unvaccinated'=prob_unvaccinated_case_each_variant, prob_vaccinated_case_each_variant_and_brand)
  conditional_vector_controls = c('unvaccinated'=prob_unvaccinated_given_control, prob_vaccinated_each_brand_given_control)

  return(
    list(
      full_table_cases = full_table_cases,
      conditional_table_cases = t(t(full_table_cases) / colSums(full_table_cases)),
      conditional_vector_controls = conditional_vector_controls
    )
  )
}

get_allornone_ve_estimand_components = function(anticipated_VE_for_each_brand_and_variant, incidence_rate_unvaccinated, study_period){

  lambda_i = incidence_rate_unvaccinated
  Lambda = sum(lambda_i)
  Theta_m = rowSums(anticipated_VE_for_each_brand_and_variant)

  pr_c_0_v_0 = exp(-Lambda*study_period)
  overall_attack_rate_unvaccinated = 1 - pr_c_0_v_0
  pr_c_i_v_0 = overall_attack_rate_unvaccinated * lambda_i/Lambda
  e_y_v_0 = overall_attack_rate_unvaccinated / Lambda

  total_brands = nrow(anticipated_VE_for_each_brand_and_variant)
  total_variants = length(incidence_rate_unvaccinated)
  combination_list = lapply(X = 1:total_variants, FUN = combn, x=total_variants)

  pr_c_0_v_m = pr_c_0_v_0 * (1 - Theta_m) +
    sapply(
      X = 1:total_brands, FUN = function(brand){
        thetas = anticipated_VE_for_each_brand_and_variant[brand,]
        counter = 1
        summation = 0
        for(i in 1:length(combination_list)){
          for(j in 1:ncol(combination_list[[i]])){
            combination = combination_list[[i]][,j]
            summation = summation + exp(-(Lambda - sum(lambda_i[combination])) * study_period) * thetas[counter]
            counter = counter + 1
          }
        }
        return(summation)
      }
    )

  pr_c_i_v_m = lambda_i %*% t(overall_attack_rate_unvaccinated * (1 - Theta_m) / Lambda)
  rownames(pr_c_i_v_m) = paste0('variant', 1:total_variants)

  for(variant in 1:nrow(pr_c_i_v_m)){
    for(brand in 1:ncol(pr_c_i_v_m)){
      thetas = anticipated_VE_for_each_brand_and_variant[brand,]
      combination_list = lapply(X = 1:total_variants, FUN = combn, x=total_variants)
      counter = 1
      for(i in 1:length(combination_list)){
        for(j in 1:ncol(combination_list[[i]])){
          combination = combination_list[[i]][,j]
          if(!variant %in% combination){
            pr_c_i_v_m[variant, brand] = pr_c_i_v_m[variant, brand] +
              if(Lambda - sum(lambda_i[combination]) == 0){
                lambda_i[variant] * study_period * thetas[counter]
              }else{
                lambda_i[variant]*(1 - exp(-(Lambda - sum(lambda_i[combination])) * study_period)) * thetas[counter] / (Lambda - sum(lambda_i[combination]))
              }

          }
          counter = counter + 1
        }
      }
    }
  }

  e_y_v_m = overall_attack_rate_unvaccinated * (1 - Theta_m) / Lambda +
    sapply(
      X = 1:total_brands, FUN = function(brand){
        thetas = anticipated_VE_for_each_brand_and_variant[brand,]
        combination_list = lapply(X = 1:total_variants, FUN = combn, x=total_variants)
        counter = 1
        summation = 0
        for(i in 1:length(combination_list)){
          for(j in 1:ncol(combination_list[[i]])){
            combination = combination_list[[i]][,j]
            summation = summation +
              if((Lambda - sum(lambda_i[combination])) == 0){
                study_period * thetas[counter]
              }else{
                (1 - exp(-(Lambda - sum(lambda_i[combination])) * study_period)) * thetas[counter] / (Lambda - sum(lambda_i[combination]))
              }

            counter = counter + 1
          }
        }
        return(summation)
      }
    )


  absolute_ve_true = matrix(data=0, nrow = total_variants, ncol = total_brands)
  for(variant in 1:nrow(absolute_ve_true)){
    for(brand in 1:ncol(absolute_ve_true)){
      thetas = anticipated_VE_for_each_brand_and_variant[brand,]
      combination_list = lapply(X = 1:total_variants, FUN = combn, x=total_variants)
      counter = 1
      for(i in 1:length(combination_list)){
        for(j in 1:ncol(combination_list[[i]])){
          combination = combination_list[[i]][,j]
          if(variant %in% combination){
            absolute_ve_true[variant, brand] = absolute_ve_true[variant, brand] + thetas[counter]
          }
          counter = counter + 1
        }
      }
    }
  }
  rownames(absolute_ve_true) = paste0("variant",1:total_variants)
  colnames(absolute_ve_true) = paste0("brand",1:total_brands)

  return(list(pr_c_0_v_0 = pr_c_0_v_0, pr_c_i_v_0 = pr_c_i_v_0, e_y_v_0 = e_y_v_0,
              pr_c_0_v_m = pr_c_0_v_m, pr_c_i_v_m = pr_c_i_v_m, e_y_v_m = e_y_v_m,
              theta_im = 1 - absolute_ve_true))
}

get_leaky_ve_estimand_components = function(anticipated_VE_for_each_brand_and_variant,
                                            incidence_rate_unvaccinated,
                                            study_period){
  lambda_i = incidence_rate_unvaccinated
  Lambda = sum(lambda_i)
  pr_c_0_v_0 = exp(-Lambda*study_period)
  overall_attack_rate_unvaccinated = 1 - pr_c_0_v_0
  pr_c_i_v_0 = overall_attack_rate_unvaccinated * lambda_i/Lambda
  e_y_v_0 = overall_attack_rate_unvaccinated / Lambda

  theta_im = t(1 - anticipated_VE_for_each_brand_and_variant)
  Theta_m = colSums(theta_im * lambda_i) / Lambda
  Theta_mLambda = Theta_m * Lambda
  pr_c_0_v_m = exp(-Theta_mLambda * study_period)
  overall_attack_rate_vaccinated_m = 1 - pr_c_0_v_m

  pr_c_i_v_m = t(t(theta_im * lambda_i) * (overall_attack_rate_vaccinated_m / Theta_mLambda))
  e_y_v_m = overall_attack_rate_vaccinated_m / Theta_mLambda

  return(list(pr_c_0_v_0 = pr_c_0_v_0, pr_c_i_v_0 = pr_c_i_v_0, e_y_v_0 = e_y_v_0,
              pr_c_0_v_m = pr_c_0_v_m, pr_c_i_v_m = pr_c_i_v_m, e_y_v_m = e_y_v_m,
              theta_im = theta_im))
}

get_ve_from_components = function(components){
  pr_c_0_v_0 = components$pr_c_0_v_0
  pr_c_i_v_0 = components$pr_c_i_v_0
  e_y_v_0 = components$e_y_v_0
  pr_c_0_v_m = components$pr_c_0_v_m
  pr_c_i_v_m = components$pr_c_i_v_m
  e_y_v_m = components$e_y_v_m
  theta_im = components$theta_im

  total_variants = nrow(theta_im)
  total_brands = ncol(theta_im)

  absolute_ve_true = 1 - theta_im
  absolute_ve_irr = 1 - (pr_c_i_v_m / pr_c_i_v_0) * (e_y_v_0 / e_y_v_m)
  absolute_ve_crr = 1 - (pr_c_i_v_m / pr_c_i_v_0)
  absolute_ve_or = 1 - (pr_c_i_v_m / pr_c_i_v_0) * (pr_c_0_v_0/pr_c_0_v_m)

  variant_comparisons = combn(total_variants, 2)
  variant_comparisons = cbind(variant_comparisons, variant_comparisons[c(2,1),])
  variant_comparisons = rbind(variant_comparisons[, rep(1:ncol(variant_comparisons), each=total_brands), drop=F],
                              rep(1:total_brands, ncol(variant_comparisons)), NA, NA, NA, NA)
  rownames(variant_comparisons) = c("variant1", "variant2", "vaccine", "relative_ve_true", "relative_ve_irr", "relative_ve_crr", "relative_ve_or")
  variant_comparisons[c("relative_ve_true","relative_ve_irr", "relative_ve_crr", "relative_ve_or"),] = apply(X = variant_comparisons, MARGIN = 2, FUN = function(comparison){
    variant1 = comparison['variant1']
    variant2 = comparison['variant2']
    vaccine = comparison['vaccine']

    thetas = theta_im[c(variant1, variant2), vaccine]

    c(
      1 - thetas[1]/thetas[2],
      rep(1 - (pr_c_i_v_m[variant1, vaccine]/pr_c_i_v_m[variant2, vaccine])/ (pr_c_i_v_0[variant1]/pr_c_i_v_0[variant2]), 3)
    )
  })

  vaccine_comparisons = combn(total_brands, 2)
  vaccine_comparisons = cbind(vaccine_comparisons, vaccine_comparisons[c(2,1),])
  vaccine_comparisons = rbind(vaccine_comparisons[, rep(1:ncol(vaccine_comparisons), each=total_variants), drop=F],
                              rep(1:total_variants, ncol(vaccine_comparisons)), NA, NA, NA, NA)
  rownames(vaccine_comparisons) = c("vaccine1", "vaccine2", "variant", "relative_ve_true", "relative_ve_irr", "relative_ve_crr", "relative_ve_or")

  vaccine_comparisons[c("relative_ve_true","relative_ve_irr", "relative_ve_crr", "relative_ve_or"),] = apply(X = vaccine_comparisons, MARGIN = 2, FUN = function(comparison){
    vaccine1 = comparison['vaccine1']
    vaccine2 = comparison['vaccine2']
    variant = comparison['variant']

    thetas = theta_im[variant, c(vaccine1, vaccine2)]
    theta_ratio = thetas[1]/thetas[2]

    c(1 - theta_ratio,
      1 - (pr_c_i_v_m[variant, vaccine1]/pr_c_i_v_m[variant, vaccine2]) / (e_y_v_m[vaccine2]/e_y_v_m[vaccine1]),
      1 - (pr_c_i_v_m[variant, vaccine1]/pr_c_i_v_m[variant, vaccine2]),
      1 - (pr_c_i_v_m[variant, vaccine1]/pr_c_i_v_m[variant, vaccine2]) / (pr_c_0_v_m[vaccine2]/pr_c_0_v_m[vaccine1])
    )
  })

  return(list(absolute_ve = list(absolute_ve_true=absolute_ve_true,
                                 absolute_ve_irr=absolute_ve_irr,
                                 absolute_ve_crr=absolute_ve_crr,
                                 absolute_ve_or=absolute_ve_or),
              relative_ve_across_variant_given_vaccine = variant_comparisons,
              relative_ve_across_vaccines_given_variant = vaccine_comparisons))
}

get_casecontrol_catchment_or = function(anticipated_brand_VEs=c(0.8, 0.5, 0.3),
                                        overall_brand_proportions = c(0.3, 0.5, 0.2),
                                        overall_vaccine_coverage=0.3,
                                        overall_attack_rate_in_unvaccinated = 0.1,
                                        prob_getting_swabbed_given_ili_sari = 0.5,
                                        ili_sari_symptom_prob=NA,
                                        ili_sari_symptom_incidence_rate=NA,
                                        study_period_length = NA,
                                        total_cases=seq(50,500, 10)){
  ili_sari_symptom_prob = get_ili_sari_symptom_prob(ili_sari_symptom_prob, ili_sari_symptom_incidence_rate, study_period_length)

  brand_vaccine_coverages = overall_brand_proportions * overall_vaccine_coverage
  prob_unvaccinated_case = overall_attack_rate_in_unvaccinated * (1 - overall_vaccine_coverage)
  prob_vaccinated_case = brand_vaccine_coverages * (1 + (1-overall_attack_rate_in_unvaccinated)/((1-anticipated_brand_VEs) * overall_attack_rate_in_unvaccinated))^-1
  prob_case = sum(prob_unvaccinated_case, prob_vaccinated_case)

  catchment = total_cases / (prob_case * prob_getting_swabbed_given_ili_sari * ili_sari_symptom_prob)
  return(data.frame(total_cases=total_cases,
                    catchment = catchment,
                    overall_attack_rate_in_unvaccinated = overall_attack_rate_in_unvaccinated,
                    prob_getting_swabbed_given_ili_sari = prob_getting_swabbed_given_ili_sari,
                    ili_sari_symptom_incidence_rate = ili_sari_symptom_incidence_rate,
                    study_period_length = study_period_length,
                    ili_sari_symptom_prob = ili_sari_symptom_prob)
  )
}
