# smoothAndRound = function(dat, xname){
#   fit = stats::smooth.spline(x = dat[,xname], y = dat$VE_lower, df = SPLINE_DF)
#   return(round(pmax(0, stats:::predict.smooth.spline(object = fit, x = dat[,xname])$y), 2))
# }

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
