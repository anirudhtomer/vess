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
