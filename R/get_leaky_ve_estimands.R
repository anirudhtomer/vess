#' Time-dependency of variant and vaccine-specific efficacy (VE) estimands based on incidence rate ratio, cumulative-incidence ratio, and odds ratio
#' @description
#' The function `get_leaky_ve_estimands` provides the VE estimands based on incidence rate ratio, cumulative-incidence ratio, and odds ratio
#' for a leaky vaccine given the actual efficacy of the vaccine. The function returns both absolute and relative VE estimands.
#'
#' @param anticipated_VE_for_each_brand_and_variant a matrix of vaccine efficacy of each vaccine (row) against each variant (column). Each value must be a real number between 0 and 1.
#' @param incidence_rate_unvaccinated a vector denoting the incidence rates of each variant in the unvaccinated subjects.
#' @param study_period the study period (epochs) should be numeric value greater than 0.
#'
#' @details
#' A leaky vaccine is characterized by its ability to reduce the force of infection `lambda` in unvaccinated subjects to factor
#' `theta * lambda` in vaccinated subjects, where `0 <= 1 - theta <=1 `is the true VE of the leaky vaccine. This true VE
#'  is to be entered in the matrix anticipated_VE_for_each_brand_and_variant.
#'
#' Smith et al. (1984) have shown that when a vaccine has a leaky action mechanism (Halloran et al., 2010, page 132),
#' then VE calculated as `VE = 1 - cumulative-incidence ratio` depends on the length of the study period `t`.
#' Specifically, when `t -> Infinity` then `1 - cumulative-incidence ratio -> 0` and thus it does not reflect the
#' biological effect of the vaccine for a given exposure to the pathogen. The aim of this function is to
#' extend upon the findings of Smith et al. (1984) for multiple variant/variants and vaccines. And to help
#' practitioners understand why certains VE measures may not reflect the biological VE.
#' More details can be found in our related paper.

#' @return A list with three elements, namely `absolute_ve`, `relative_ve_across_variant_given_vaccine`,
#' and `relative_ve_across_vaccines_given_variant`.

#' The list element `absolute_ve` pertains to the absolute efficacy of the vaccines.
#' It contains four elements, namely `absolute_ve_true`, `absolute_ve_irr`,
#' `absolute_ve_cir`, and `absolute_ve_or`. Here `absolute_ve_true` is same as the input parameter
#' `anticipated_VE_for_each_brand_and_variant`, and indicates the true biological efficacy of the vaccine.
#' Whereas, the remaining three are the VE estimands based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_variant_given_vaccine` pertains to the relative VE of a vaccine against two different
#' variants/variants. It is a matrix with 7 rows and multiple columns. Each column being one of the comparison sets of variants and vaccines.
#' Rows 1 and 2 indicate the two variants of interest for comparison and row 3 has the vaccine.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_vaccines_given_variant` pertains to the relative VE of two vaccines against the
#' same variants/variants. It is a matrix with 7 rows and multiple columns.
#' Each column being one of the comparison sets of variants and vaccines.
#' Rows 1 and 2 indicate the two vaccines of interest for comparison and row 3 has the variant.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.

#'
#' @examples As an example we recommend running the function without passing any parameter to it.
#' The default scenario is for three vaccines and three pathogen variants.
#'
#' @references
#' 1. Halloran, M. E., Longini, I. M., Struchiner, C. J., and Longini, I. M. (2010).Design and analysis of vaccine studies, volume 18. Springer.
#' 2. Smith, P., Rodrigues, L., and Fine, P. (1984). Assessment of the protective efficacy of vaccines against common diseases using case-control and cohort studies.International journal of epidemiology, 13(1):87–93.
#' @importFrom utils combn
#' @export
get_leaky_ve_estimands = function(anticipated_VE_for_each_brand_and_variant=
                                    matrix(data=c(0.1, 0.3, 0.4, 0.6, 0.7, 0.9), nrow = 2, ncol = 3, byrow = F,
                                           dimnames = list(paste0('brand', 1:2), paste0('variant', 1:3))),
                                  incidence_rate_unvaccinated = c(0.05, 0.1, 0.2),
                                  study_period = 1){

  if(any(is.na(incidence_rate_unvaccinated)) | !is.numeric(incidence_rate_unvaccinated) | any(incidence_rate_unvaccinated < 0)){
    stop("The values of the parameter 'incidence_rate_unvaccinated' should be more than 0")
  }

  if(ncol(anticipated_VE_for_each_brand_and_variant)!=length(incidence_rate_unvaccinated)){
    stop("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_variant' should be equal to length of the parameter 'incidence_rate_unvaccinated'")
  }

  components = get_leaky_ve_estimand_components(anticipated_VE_for_each_brand_and_variant, incidence_rate_unvaccinated, study_period)
  return(get_ve_from_components(components))
}
