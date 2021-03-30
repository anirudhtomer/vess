library(epitools)
library(doParallel)

rm(list = ls())

df = expand.grid(coverage = seq(0.05, 0.99, 0.005),
                 vaccine_effectiveness = seq(0.05, 0.99, 0.01),
                 control_case_ratio = c(1:10,20),
                 total_cases = seq(50, 1000, 25))

df = expand.grid(coverage = c(0.05, 0.1, 0.15, 0.2),
                 vaccine_effectiveness = 0.8,
                 control_case_ratio = c(1),
                 total_cases = seq(25, 1000, 25))

t1 = Sys.time()
cl = makeCluster(detectCores()-1)
registerDoParallel(cl)
temp = foreach(i=1:nrow(df)) %dopar% {
  
  set.seed(2020 + i)
  control_case_ratio = df$control_case_ratio[i]
  total_cases = df$total_cases[i]
  total_controls = total_cases * control_case_ratio
  
  n_studies_back_intime = 1000
  estimate = rep(NA, n_studies_back_intime)
  upp = rep(NA, n_studies_back_intime)
  low = rep(NA, n_studies_back_intime)
  
  for(j in 1:n_studies_back_intime){
    
    statusvaccine_control = rbinom(total_controls, 1, df$coverage[i])
    
    b_yesvaccine_yescontrol = sum(statusvaccine_control)
    d_novaccine_yescontrol = total_controls-b_yesvaccine_yescontrol
    
    a_divided_c = (1-df$vaccine_effectiveness[i]) * (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
    
    a_yesvaccine_yescase = sum(rbinom(total_cases, 1, 1 / (1 + a_divided_c^-1)))
    c_novaccine_yescase = total_cases - a_yesvaccine_yescase
    
    vaccine = c(rep(1, a_yesvaccine_yescase), rep(0, c_novaccine_yescase), 
                rep(1, b_yesvaccine_yescontrol), rep(0,d_novaccine_yescontrol))
    hospitalized = c(rep(1, total_cases), rep(0, total_controls))
    
    #tt = fisher.test(table(vaccine, hospitalized), or = 1)
    estimate[j] = (a_yesvaccine_yescase/c_novaccine_yescase) / (b_yesvaccine_yescontrol/d_novaccine_yescontrol)
    se = sqrt(1/a_yesvaccine_yescase + 1/c_novaccine_yescase + 1/b_yesvaccine_yescontrol + 1/d_novaccine_yescontrol)
    low[j] = exp(log(estimate[j]) + qnorm(0.025) * se)
    upp[j] = exp(log(estimate[j]) + qnorm(0.975) * se)
  }
  return(list(mean(estimate, na.rm = T), mean(low, na.rm = T), mean(upp, na.rm = T)))
}
stopCluster(cl)
t2=Sys.time()

df$mean_estimate = 1-sapply(temp, "[[", 1)
#df$maria_low = sapply(temp, "[[", 4)
df$mean_low = pmax(1-sapply(temp, "[[", 3), 0)
df$mean_upp = 1-sapply(temp, "[[", 2)


ggplot(data=df) +
  geom_line(aes(x=total_cases, y=mean_low, group=factor(control_case_ratio), color=factor(control_case_ratio))) +
  geom_hline(aes(yintercept = vaccine_effectiveness)) + 
  geom_label(aes(x=250, y=vaccine_effectiveness, label=vaccine_effectiveness)) +
  facet_grid(paste0("Cov=",coverage)~paste0("VE=",vaccine_effectiveness)) +
  theme_bw() + theme(text=element_text(size=14), legend.position = 'bottom') +
  ylim(0,1) +
  xlab("Total cases") + ylab("Lower 95% CI") + labs(color='Controls per case')


ggplot(data=df) +
  geom_line(aes(x=total_cases, y=maria_low, group=factor(control_case_ratio), color=factor(control_case_ratio))) +
  geom_hline(aes(yintercept = vaccine_effectiveness)) + 
  geom_label(aes(x=250, y=vaccine_effectiveness, label=vaccine_effectiveness)) +
  facet_grid(paste0("Cov=",coverage)~paste0("VE=",vaccine_effectiveness)) +
  ylim(0,1) +
  theme_bw() + theme(text=element_text(size=14), legend.position = 'bottom') +
  xlab("Total cases") + ylab("Lower 95% CI") + labs(color='Controls per case')

ggplot(data=df) +
  geom_line(aes(x=total_cases, y=(mean_low-maria_low), group=factor(control_case_ratio), color=factor(control_case_ratio))) +
  geom_hline(aes(yintercept = vaccine_effectiveness)) + 
  geom_label(aes(x=250, y=vaccine_effectiveness, label=vaccine_effectiveness)) +
  facet_grid(paste0("Cov=",coverage)~paste0("VE=",vaccine_effectiveness)) +
  theme_bw() + theme(text=element_text(size=14), legend.position = 'bottom') +
  xlab("Total cases") + ylab("Absolute difference in lower 95% CI") + labs(color='Controls per case')



