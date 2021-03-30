
###########################################
# Case control-overall-min det VE         #
###########################################

Casecontrol_overall_mindetVE<-function(power=0.8,                      #Power (Input=numeric)
                                       alpha=0.05,                     #Significance level (Input=numeric)
                                       design=1,                       #Design-effect - inflates the variance of the parameter estimates to allow for correlations among clustered observations (Input=numeric)
                                       COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                       r=1,                            #Controls per case (Input=numeric)
                                       ncases=seq(100, 4000, 100))     #Total number of cases in the study (Input=vector)
{
  clevel=1-0.05
  #--- pre-allocate output matrix
  n_cov <- length(COV)
  n_ncases <- length(ncases)
  nn <- n_cov * n_ncases
  output <- array(numeric(4*nn), dim=c(nn,4))
  #--- perform the sample size calculation; what is the minimal detectable odds ratio given power, alpha, coverage, allocation ratio and number of cases
  cnt <- 0
  for (i in 1:n_cov){
    for (ii in 1:n_ncases) {
      cnt <- cnt + 1
      cov <- COV[i]/100
      N <- (1 + r) * ncases[ii]                   # total sample size: cases + controls
      if (ncases[ii]<200 & cov==0.05){
        output[cnt,4]<-NA
      }
      else{
        t <- epi.sscc(OR = NA, p0 = cov, n = N, power = power, r = r, rho=0, design = design, sided.test = 2, conf.level = clevel, method = "unmatched", fleiss = FALSE)
        output[cnt,4] <- pmin((1 - t$OR[1]),1) * 100
      }
      output[cnt,1] <- cov * 100
      output[cnt,2] <- ncases[ii]
      output[cnt,3] <- r
    }
  }
  #--- output table
  output2 <- data.frame(output)
  output2[,1] <- as.factor(output2[,1])  
  colnames(output2) <- c("Coverage", "nr_cases", "r", "minVE")
  return(output2)
}
