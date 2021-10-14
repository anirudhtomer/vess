#' Time-dependency of strain and vaccine-specific efficacy (VE) estimands based on incidence rate ratio, cumulative-incidence ratio, and odds ratio
#' @description
#' The function `get_allornone_ve_estimands` provides the VE estimands based on incidence rate ratio, cumulative-incidence ratio, and odds ratio
#' for an all-or-none vaccine given the actual efficacy of the vaccine. The function returns both absolute and relative VE estimands.
#'
#' @param anticipated_VE_for_each_brand_and_strain a matrix of proportion of subjects immune to certain combination (column) of strains given a certain vaccine (row). Each value must be a real number between 0 and 1 and sum of each row should not exceed 1.
#' @param incidence_rate_unvaccinated a vector denoting the incidence rates of each strain in the unvaccinated subjects.
#' @param study_period the study period (epochs) should be numeric value greater than 0.
#'
#' @details
#' To understand the parameter `anticipated_VE_for_each_brand_and_strain` and the action of all-or-none vaccines,
#' assume that there are only three strains `1`,`2`, and `3` circulating. The subjects vaccinated with the vaccine `1`
#' are represented in row 1 of the matrix `anticipated_VE_for_each_brand_and_strain[1, ]`. Specifically,
#' `anticipated_VE_for_each_brand_and_strain[1, 1]` is the proportion of subjects who become immune to the
#' strain `1` but not to the other two strains. Among these subjects infections occur with the combined incidence rates
#' of the strains `2` and `3`. `anticipated_VE_for_each_brand_and_strain[1, 2]` and `anticipated_VE_for_each_brand_and_strain[1, 3]`
#' hold similar interpretation. Next, let `anticipated_VE_for_each_brand_and_strain[1, 4]`, be the proportion of subjects who become
#' immune to both strains `1` and `2` but not to strain `3`. Among these subjects infection with variant `3` happens with
#' the incidence rate of strain `3`. We can define `anticipated_VE_for_each_brand_and_strain[1, 5]` and
#' `anticipated_VE_for_each_brand_and_strain[1, 6], similarly. Lastly, `anticipated_VE_for_each_brand_and_strain[1, 7]` are
#' the proportion of subjects who become immune to all three strains.
#' The remaining proportion of subjects do not become immune to any variant despite vaccination.
#' For them infections happens as in as in the placebo subjects. This proportion is not entered in the matrix
#' `anticipated_VE_for_each_brand_and_strain` because this proportion can be obtained as
#' `1 - rowSums(anticipated_VE_for_each_brand_and_strain)`.
#'
#' Smith et al. (1984) have shown that when a vaccine has an all-or-none action mechanism (Halloran et al., 2010, page 132),
#' then VE calculated as `VE = 1 - incidence rate ratio` depends on the length of the study period `t`.
#' Specifically, when `t -> Infinity` then `1 - incidence rate ratio -> 1` and thus it does not reflect the
#' biological effect of the vaccine for a given exposure to the pathogen. The aim of this function is to
#' extend upon the findings of Smith et al. (1984) for multiple strain/variants and vaccines. And to help
#' practitioners understand why certains VE measures may not reflect the biological VE.
#' More details can be found in our related paper.

#' @return A list with three elements, namely `absolute_ve`, `relative_ve_across_strain_given_vaccine`,
#' and `relative_ve_across_vaccines_given_strain`.

#' The list element `absolute_ve` pertains to the absolute efficacy of the vaccines.
#' It contains four elements, namely `absolute_ve_true`, `absolute_ve_irr`,
#' `absolute_ve_cir`, and `absolute_ve_or`. Here `absolute_ve_true` is same as the input parameter
#' `anticipated_VE_for_each_brand_and_strain`, and indicates the true biological efficacy of the vaccine.
#' Whereas, the remaining three are the VE estimands based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_strain_given_vaccine` pertains to the relative VE of a vaccine against two different
#' strains/variants. It is a matrix with 7 rows and multiple columns. Each column being one of the comparison sets of strains and vaccines.
#' Rows 1 and 2 indicate the two strains of interest for comparison and row 3 has the vaccine.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.
#'
#' The list element `relative_ve_across_vaccines_given_strain` pertains to the relative VE of two vaccines against the
#' same strains/variants. It is a matrix with 7 rows and multiple columns.
#' Each column being one of the comparison sets of strains and vaccines.
#' Rows 1 and 2 indicate the two vaccines of interest for comparison and row 3 has the strain.
#' The remaining rows show the true relative VE, and the relative VE based on incidence rate ratio, cumulative-incidence ratio,
#' and odds ratio, respectively.

#'
#' @examples As an example we recommend running the function without passing any parameter to it.
#' The default scenario is for three vaccines and three pathogen strains.
#'
#' @references
#' 1. Halloran, M. E., Longini, I. M., Struchiner, C. J., and Longini, I. M. (2010).Design and analysis of vaccine studies, volume 18. Springer.
#' 2. Smith, P., Rodrigues, L., and Fine, P. (1984). Assessment of the protective efficacy of vaccines against common diseases using case-control and cohort studies.International journal of epidemiology, 13(1):87–93.
#' @importFrom utils combn
#' @export
get_allornone_ve_estimands = function(anticipated_VE_for_each_brand_and_strain=
                                        matrix(data=c(0.04, 0.1, 0.03, 0.05, 0.03, 0.07, 0.4, 0.05, 0.18, 0.01, 0.08, 0.04, 0.07, 0.2, 0.08, 0.05, 0.07, 0.01, 0.02, 0.13, 0.1),
                                               nrow = 3, ncol = 7, byrow = T,
                                               dimnames = list(paste0('brand', 1:3), paste0('strain', c('1','2','3','12','13','23','123')))),
                                      incidence_rate_unvaccinated = c(0.05, 0.1, 0.2),
                                      study_period = 1){

  if(any(is.na(incidence_rate_unvaccinated)) | !is.numeric(incidence_rate_unvaccinated) | any(incidence_rate_unvaccinated < 0)){
    stop("The values of the parameter 'incidence_rate_unvaccinated' should be more than 0")
  }

  total_comp_expected = sum(sapply(lapply(1:length(incidence_rate_unvaccinated), combn, x=length(incidence_rate_unvaccinated)), ncol))

  if(ncol(anticipated_VE_for_each_brand_and_strain)!=total_comp_expected){
    stop(paste0("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_strain' should be equal to", total_comp_expected))
  }

  components = get_allornone_ve_estimand_components(anticipated_VE_for_each_brand_and_strain, incidence_rate_unvaccinated, study_period)
  return(get_ve_from_components(components))
}
