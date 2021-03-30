###########################################
# Cohort-brand specific-min det VE        #
###########################################
Cohort_brandspecific_mindetVE<-function(power=0.8,                #Power (Input=numeric)
                                  alpha=0.05,                     #Significance level (Input=numeric)
                                  design=1,                       #Design-effect - inflates the variance of the parameter estimates to allow for correlations among clustered observations (Input=numeric)
                                  COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                  AR=5,                           #attack rate in the unvaccinated (Input=numeric)
                                  N=seq(1000,50000,200),          #Total number of subjects in the study (Input=vector)
                                  propbrand=c(10,20,40,60,80,90)) #Proportion of the brand (Input=vector/numeric)    
{
  clevel=1-alpha
  pbrand=propbrand/100
  #--- pre-allocate output matrix
  n_cov <- length(COV)
  n_n <- length(N)
  n_pbrand <- length(propbrand) 
  nn <- n_cov * n_n * n_pbrand
  output <- array(numeric(5*nn), dim=c(nn,5))
  #--- perform the sample size calculation; what is the minimal detectable risk ratio given power, alpha, coverage, baseline risk and sample size
  cnt <- 0
  for (i in 1:n_cov){
    for (ii in 1:n_n) {
      for (iii in 1:n_pbrand) {
        cnt <- cnt + 1
        pbrandi <- pbrand[iii]
        covi <- COV[i]/100
        covi2 <- (pbrandi * covi)/(1 - covi + pbrandi * covi)   
        r <- covi2/(1-covi2)                                    
        n <- N[ii]                                              
        Ni <- (1 - covi + pbrandi * covi) * n                   
        t <- epi.sscohortc(NA, AR/100, n = Ni, power = power, r = r, design = design, sided.test = 2, conf.level = clevel)
        output[cnt,1] <- AR
        output[cnt,2] <- covi*100
        output[cnt,3] <- n
        output[cnt,4] <- pmin(1 - t$irr[1],1) * 100
        output[cnt,5] <- pbrandi*100
      }#iii
    }#ii
  }#i
  #--- output table
  output2 <- data.frame(output)
  output2[,2] <- as.factor(output2[,2])  
  colnames(output2) <- c("AR","Coverage", "SS","minVE","prop_brand")
  return(output2)
}