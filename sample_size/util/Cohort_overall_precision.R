
###########################################
# Cohort-overall-precision                #
###########################################

Cohort_overall_precision<-function(alpha=0.05,                #Significance level (Input=numeric)
                                   VE=c(20,50,70),            #True vaccine effectiveness(Input=vector/numeric)
                                   COV=c(5, 20, 50, 70),      #Coverage (Input=vector/numeric)
                                   AR=5,                      #Attack rate in the unvaccinated (Input=numeric)
                                   N=seq(1000,50000,200))     #Total number of subjects in the study (Input=vector)
{
  
  
  df = expand.grid(Coverage = COV, SS = N, AR = AR, VE=VE)
  df$VE_lower = getVELowerCohortPrecision(alpha, df)
  
  df = df[order(df$Coverage, df$SS),]
  df$VE_lower = unlist(by(data = df, INDICES = df$Coverage, FUN = smoothAndRound, xname="SS"))
  
  return(df)
  
  # #--- pre-allocate output matrix
  # n_VE <- length(VE)
  # n_cov <- length(COV)
  # n_n <- length(N)
  # nn <- n_VE * n_cov * n_n
  # output <- array(numeric(5*nn), dim=c(nn,5))
  # #--- perform the sample size calculation; what is the precision (lower bound of the CI) for a given number of cases
  # cnt <- 0
  # for (i in 1:n_VE){
  #   for (ii in 1:n_cov){
  #     for (iii in 1:n_n) {
  #       cnt <- cnt + 1
  #       # given parameters
  #       VEi <- VE[i]/100
  #       covi <- COV[ii]/100
  #       ni <- N[iii]
  #       ARi  <-AR/100
  #       casesunexposedi    <- ni*(1-covi)*ARi
  #       casesexposedi      <- ni*covi*ARi*(1-VEi)
  #       unexposed <-ni*(1-covi)
  #       exposed <-ni*covi
  #       RRi <-1-VEi
  #       selogRRi <-sqrt(1/casesunexposedi +1/casesexposedi - (1/exposed +1/unexposed))
  #       UCLRRi <- exp(log(RRi)+qnorm(1-alpha/2)*selogRRi)
  #       LCLVEi <- pmax(1-UCLRRi,0)
  #       # output
  #       output[cnt,1] <- covi*100
  #       output[cnt,2] <- ni
  #       output[cnt,3] <- AR
  #       output[cnt,4] <- VEi * 100
  #       output[cnt,5] <- LCLVEi * 100
  #     }#iii
  #   }#ii
  # }#i
  # #--- output table
  # output2 <- data.frame(output)
  # colnames(output2) <- c("Coverage", "SS", "AR", "VE", "VE_lower")
  # return(output2)
}
