#' @importFrom utils combn
#' @export
get_leaky_ve_estimands = function(anticipated_VE_for_each_brand_and_strain=
                                    matrix(data=c(0.1, 0.3, 0.4, 0.6, 0.7, 0.9), nrow = 2, ncol = 3, byrow = F,
                                           dimnames = list(paste0('brand', 1:2), paste0('strain', 1:3))),
                                  incidence_rate_unvaccinated = c(0.05, 0.1, 0.2),
                                  study_period = 1){

  if(any(is.na(incidence_rate_unvaccinated)) | !is.numeric(incidence_rate_unvaccinated) | any(incidence_rate_unvaccinated < 0)){
    stop("The values of the parameter 'incidence_rate_unvaccinated' should be more than 0")
  }

  if(ncol(anticipated_VE_for_each_brand_and_strain)!=length(incidence_rate_unvaccinated)){
    stop("Total number of columns of the parameter 'anticipated_VE_for_each_brand_and_strain' should be equal to length of the parameter 'incidence_rate_unvaccinated'")
  }

  components = get_leaky_ve_estimand_components(anticipated_VE_for_each_brand_and_strain, incidence_rate_unvaccinated, study_period)
  return(get_ve_from_components(components))
}
