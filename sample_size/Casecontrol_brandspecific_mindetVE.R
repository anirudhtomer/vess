Casecontrol_brandspecific_mindetVE<-function(power=0.8,                #Power (Input=numeric)
                                       alpha=0.05,                     #Significance level (Input=numeric)
                                       design=1,                       #Design-effect - inflates the variance of the parameter estimates to allow for correlations among clustered observations (Input=numeric)
                                       COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                       r=1,                            #Controls per case (Input=numeric)
                                       ncases=seq(100, 4000, 100),     #Total number of cases in the study (Input=vector)
                                       propbrand=c(10,20,40,60,80,90)) #Proportion of the brand (Input=vector/numeric) 
{
  clevel=1-alpha
  pbrand=propbrand/100
  #--- pre-allocate output matrix
  n_cov <- length(COV)
  n_ncases <- length(ncases)
  n_pbrand <- length(pbrand) 
  nn <- n_cov * n_ncases * n_pbrand
  output <- array(numeric(5*nn), dim=c(nn,5))
  #--- perform the sample size calculation; what is the minimal detectable odds ratio given power, alpha, coverage, allocation ratio and number of cases
  cnt <- 0
  for (i in 1:n_cov){
    for (ii in 1:n_ncases) {
      for (iii in 1:n_pbrand) {
        cnt <- cnt + 1
        covi <- COV[i]/100
        ncasesi <- ncases[ii]  
        pbrandi <- pbrand[iii]
        covi2 <- (pbrandi * covi)/(1 - covi + pbrandi * covi)    # brand-specific coverage, recalculated after excluding other brands
        Ni <- (1 - covi + pbrandi * covi)*(1+r)*ncasesi           # Total sample size (cases + controls), after excluding other brands
        t <- try(epi.sscc(OR = NA, p0 = covi2, n = Ni, power = power, r = r, rho=0, design = design, sided.test = 2, conf.level = clevel, method = "unmatched", fleiss = FALSE),silent=T)
        output[cnt,1] <- covi*100
        output[cnt,2] <- pbrandi * 100
        output[cnt,3] <- ncasesi
        output[cnt,4] <- r
        if(class(t) != "try-error")
          output[cnt,5] <- pmin((1 - t$OR[1]),1) * 100
        if(output[cnt,5]==0)
          output[cnt,5]=NA
      }
    }
  }
  #--- output table
  output2 <- data.frame(output)
  output2[,1] <- as.factor(output2[,1])  
  colnames(output2) <- c("Coverage", "prop_brand", "nr_cases", "r", "minVE")
  return(output2)
}