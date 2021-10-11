# library(epiR)
# library(ggplot2)
#
# source("R/utils.R")
# source("R/get_cohort_mindet_VE_cir.R")
# source("R/get_cohort_mindet_VE_irr.R")
#
# vemat = matrix(data=c(0.4, 0.2), nrow = 2, ncol = 1, byrow = F,
#                dimnames = list(paste0('brand', 1:2), paste0('strain', 1)))
# brand_prop = c('brand1'=0.8, 'brand2'=0.2)
# total_sub = seq(10000, 50000, by = 500)
# coverage = 0.5
# attack_rate = 0.1
#
# browser()
# tt_cir = get_cohort_mindet_VE_cir(anticipated_VE_for_each_brand_and_strain = vemat,
#                                 brand_proportions_in_vaccinated = brand_prop,
#                                 overall_vaccine_coverage = coverage,
#                                 proportion_strains_in_unvaccinated_cases = 1,
#                                 overall_attack_rate_in_unvaccinated = attack_rate,
#                                 calculate_relative_VE = T,
#                                 power = 0.8, alpha = 0.05, total_subjects = total_sub)
#
# tt_irr = get_cohort_mindet_VE_irr(anticipated_VE_for_each_brand_and_strain = vemat,
#                                   brand_proportions_in_vaccinated = brand_prop,
#                                   overall_vaccine_coverage = coverage,
#                                   proportion_strains_in_unvaccinated_cases = 1,
#                                   overall_attack_rate_in_unvaccinated = attack_rate,
#                                   calculate_relative_VE = T,
#                                   power = 0.8, alpha = 0.05, total_subjects = total_sub)
#
# tt_cir$type = "Incidence risk"
# tt_irr$type = "Incidence rate"
#
# tt = rbind(tt_cir, tt_irr)
# tt$group = paste(tt$vaccine_1, ' vs. ',tt$vaccine_2, ' for ', tt$strain)
#
# library(ggplot2)
# ggplot(data=tt) + geom_line(aes(x = total_subjects, y=mindet_VE, group = type, color=type)) +
#   facet_grid(.~group) + ylim(0,1) + theme_bw() + theme(legend.position = 'bottom')
#
#
# plot_df = expand.grid(subpopulation_coverage = seq(0.1, 0.9, by = 0.2),
#                       irexp0 = seq(0.1, 0.9, 0.2), total_sub = total_sub,
#                       type = c("Incidence rate", "Incidence risk"), mindet_VE = NA)
#
# plot_df$r = plot_df$subpopulation_coverage/(1-plot_df$subpopulation_coverage)
#
# plot_df$mindet_VE[plot_df$type=="Incidence rate"] = c(sapply(X = 1:nrow(plot_df), FUN = function(i){
#   ret = try(1 - epi.sscohortt(
#     irexp1 = NA, irexp0 = plot_df$irexp0[i],
#     FT = 1, n = plot_df$total_sub[i], power = 0.8,
#     r = plot_df$r[i], design = 1, sided.test = 2, conf.level = 1-0.05,
#   )$irr[1], silent = T)
#   if(inherits(ret, "try-error")){
#     return(NA)
#   }else{
#     return(ret)
#   }
# }))
#
# plot_df$mindet_VE[plot_df$type=="Incidence risk"] = c(sapply(X = 1:nrow(plot_df), FUN = function(i){
#   ret = try(1 - epi.sscohortc(
#     irexp1 = NA, irexp0 = plot_df$irexp0[i],
#     n = plot_df$total_sub[i], power = 0.8,
#     r = plot_df$r[i], design = 1, sided.test = 2, conf.level = 1-0.05
#   )$irr[1], silent = T)
#   if(inherits(ret, "try-error")){
#     return(NA)
#   }else{
#     return(ret)
#   }
# }))
#
# ggplot(data=plot_df) + geom_line(aes(x = total_sub, y=mindet_VE, group = type, color=type)) +
#   facet_grid(paste("Cov", subpopulation_coverage, "(", round(r, 2), ")")~paste("irexp0", irexp0)) +
#   ylim(0,1) + theme_bw() + theme(legend.position = 'bottom')
#
