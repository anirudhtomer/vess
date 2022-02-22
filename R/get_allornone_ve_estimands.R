#' Time-dependency of variant and vaccine-specific efficacy (VE) estimands based on incidence rate ratio, cumulative-risk ratio, and odds ratio
#' @description
#' The function `get_allornone_ve_estimands` provides the VE estimands based on incidence rate ratio, cumulative-risk ratio, and odds ratio
#' for an all-or-none vaccine given the actual efficacy of the vaccine. The function returns both absolute and relative VE estimands.
#'
#' @param anticipated_VE_for_each_brand_and_variant a matrix of proportion of subjects immune to certain combination (column) of variants given a certain vaccine (row). Each value must be a real number between 0 and 1 and sum of each row should not exceed 1.
#' @param incidence_rate_unvaccinated a vector denoting the incidence rates of each variant in the unvaccinated subjects.
#' @param study_period the study period (epochs) should be numeric value greater than 0.
#'
#' @details
#' To understand the parameter `anticipated_VE_for_each_brand_and_variant` and the action of all-or-none vaccines,
#' assume that there are only three variants `1`,`2`, and `3` circulating. The subjects vaccinated with the vaccine `1`
#' are represented in row 1 of the matrix `anticipated_VE_for_each_brand_and_variant[1, ]`. Specifically,
#' `anticipated_VE_for_each_brand_and_variant[1, 1]` is the proportion of subjects who become immune to the
#' variant `1` but not to the other two variants. Among these subjects infections occur with the combined incidence rates
#' of the variants `2` and `3`. `anticipated_VE_for_each_brand_and_variant[1, 2]` and `anticipated_VE_for_each_brand_and_variant[1, 3]`
#' hold similar interpretation. Next, let `anticipated_VE_for_each_brand_and_variant[1, 4]`, be the proportion of subjects who become
#' immune to both variants `1` and `2` but not to variant `3`. Among these subjects infection with variant `3` happens with
#' the incidence rate of variant `3`. We can define `anticipated_VE_for_each_brand_and_variant[1, 5]` and
#' `anticipated_VE_for_each_brand_and_variant[1, 6]`, similarly. Lastly, `anticipated_VE_for_each_brand_and_variant[1, 7]` are
#' the proportion of subjects who become immune to all three variants.
#' The remaining proportion of subjects do not become immune to any variant despite vaccination.
#' For them infections happens as in as in the placebo subjects. This proportion is not entered in the matrix
#' `anticipated_VE_for_each_brand_and_variant` because this proportion can be obtained as
#' `1 - rowSums(anticipated_VE_for_each_brand_and_variant)`.
#'
#' Smith et al. (1984) have shown that when a vaccine has an all-or-none action mechanism (Halloran et al., 2010, page 132),
#' then VE calculated as `VE = 1 - incidence rate ratio` depends on the length of the study period `t`.
#' Specifically, when `t -> Infinity` then `1 - incidence rate ratio -> 1` and thus it does not reflect the
#' biological effect of the vaccine for a given exposure to the pathogen. The aim of this function is to
#' extend upon the findings of Smith et al. (1984) for multiple variant/variants and vaccines. And to help
#' practitioners understand why certain VE measures may not reflect the biological VE.
#' More details can be found in our related paper.

#' @return A list with three elements, namely `absolute_ve`, `relative_ve_across_variant_given_vaccine`,
#' and `relative_ve_across_vaccines_given_variant`.

#' The list element `absolute_ve` pertains to the absolute efficacy of the vaccines.
#' It contains four elements, namely `absolute_ve_true`, `absolute_ve_irr`,
#' `absolute_ve_crr`, and `absolute_ve_or`. Here `absolute_ve_true` is same as the input parameter
#' `anticipated_VE_for_each_brand_and_variant`, and indicates the true biological efficacy of the vaccine.
#' Whereas, the remaining three are the VE estimands based on incidence rate ratio, cumulative-risk ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_variant_given_vaccine` pertains to the relative VE of a vaccine against two different
#' variants/variants. It is a matrix with 7 rows and multiple columns. Each column being one of the comparison sets of variants and vaccines.
#' Rows 1 and 2 indicate the two variants of interest for comparison and row 3 has the vaccine.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-risk ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_vaccines_given_variant` pertains to the relative VE of two vaccines against the
#' same variants/variants. It is a matrix with 7 rows and multiple columns.
#' Each column being one of the comparison sets of variants and vaccines.
#' Rows 1 and 2 indicate the two vaccines of interest for comparison and row 3 has the variant.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-risk ratio,
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
get_allornone_ve_estimands = function(anticipated_VE_for_each_brand_and_variant=
                                        matrix(data=c(0.04, 0.1, 0.03, 0.05, 0.03, 0.07, 0.4, 0.05, 0.18, 0.01, 0.08, 0.04, 0.07, 0.2, 0.08, 0.05, 0.07, 0.01, 0.02, 0.13, 0.1),
                                               nrow = 3, ncol = 7, byrow = T,
                                               dimnames = list(paste0('brand', 1:3), paste0('variant', c('1','2','3','12','13','23','123')))),
                                      incidence_rate_unvaccinated = c(0.05, 0.1, 0.2),
                                      study_period = 1){

  if(any(is.na(incidence_rate_unvaccinated)) | !is.numeric(incidence_rate_unvaccinated) | any(incidence_rate_unvaccinated < 0)){
    stop("The values of the parameter 'incidence_rate_unvaccinated' should be more than 0")
  }

  total_comp_expected = sum(sapply(lapply(1:length(incidence_rate_unvaccinated), combn, x=length(incidence_rate_unvaccinated)), ncol))

  if(ncol(anticipated_VE_for_each_brand_and_variant)!=total_comp_expected){
    stop(paste0("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_variant' should be equal to", total_comp_expected))
  }

  components = get_allornone_ve_estimand_components(anticipated_VE_for_each_brand_and_variant, incidence_rate_unvaccinated, study_period)
  return(get_ve_from_components(components))
}
