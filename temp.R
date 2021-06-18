library(ggplot2)
library(tidyverse)

FONT = 15

itexp <- function(u, rate, t) { -log(1-u*(1-exp(-t*rate)))/rate }
rtexp <- function(n, rate, t) { itexp(runif(n), rate, t) }

N = 10000
settings_df = expand.grid(
                          rate = c(0.0001, 0.001, 0.01, 0.1, 1),
                          study_period = c(1, 3, 5, 7)
)

plots = lapply(1:nrow(settings_df), FUN = function(i){

  rate = settings_df$rate[i]
  study_period = settings_df$study_period[i]

  rtexp_data = rtexp(n = N, rate = rate, study_period)
  rexp_data =rexp(n = N / (1-exp(-rate*study_period)), rate = rate)
  rexp_data = rexp_data[rexp_data<=study_period]
  runif_data = runif(n = N, min = 0, max = study_period)

  plot_df = data.frame(rate=rate, study_period = study_period,
                       data=c(rtexp_data, rexp_data, runif_data),
                       dist = c(rep('rtexp', N), rep('rexp', length(rexp_data)), rep('runif', N)))

  ggplot(plot_df) + geom_density(aes(data, color=dist)) +
    theme_bw() + ggtitle(paste0('study_period = ',study_period, ', rate = ', rate)) + ylim(0, 1)
})

dev.off()
ggpubr::ggarrange(plotlist = plots, ncol = 5, nrow = 4, common.legend = T,
                  legend = 'bottom')


#reasonable assumption of unif dist, and rtexp works fine

#next check sampling dist of sum
settings_df = expand.grid(study_period = c(1, 3, 5, 7),
                          n = c(2, 5, 10, 25, 50),
                          rate = c(0.01, 0.1)
)

plots = lapply(1:nrow(settings_df), FUN = function(i){

  rate = settings_df$rate[i]
  study_period = settings_df$study_period[i]
  n = settings_df$n[i]

  rtexp_data = sapply(1:N, function(i){
    #sum(rtexp(n = n, rate = rate, study_period))
    mean(rtexp(n = n, rate = rate, study_period))
  })

  runif_data = sapply(1:N, function(i){
    #sum(runif(n = n, min = 0, max = study_period))
    mean(runif(n = n, min = 0, max = study_period))
  })

  #normal_dist using var = n * var(unif)
  #rnorm_data = rnorm(n = N, mean = n*study_period/2, sd = sqrt(n * study_period^2/12))
  rnorm_data = rnorm(n = N, mean = study_period/2, sd = sqrt(study_period^2/(12*n)))

  plot_df = data.frame(rate=rate, study_period = study_period,n =n,
                       data=c(rtexp_data, runif_data, rnorm_data),
                       dist = c(rep('rtexp', N), rep('runif', N), rep('rnorm', N)))

  ggplot(plot_df) + geom_density(aes(data, color=dist)) +
    theme_bw() +
    ggtitle(paste0('n = ',n, ', study_period = ',study_period, ', rate = ', rate))
})
dev.off()
ggpubr::ggarrange(plotlist = plots, ncol = 4, nrow = 5, common.legend = T,
                  legend = 'bottom')
