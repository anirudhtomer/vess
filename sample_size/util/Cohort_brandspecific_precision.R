
###########################################
# Cohort-brand specific-precision         #
###########################################

Cohort_brandspecific_precision<-function(alpha=0.05,               #Significance level (Input=numeric)
                                   VE=c(20,50,70),                 #True vaccine effectiveness(Input=vector/numeric)
                                   COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                   AR=5,                           #Attack rate in the unvaccinated (Input=numeric)
                                   N=seq(1000,50000,200),          #Total number of subjects in the study (Input=vector)
                                   propbrand=c(10,20,40,60,80,90)) #Proportion of the brand (Input=vector/numeric) 
{
  
  df = expand.grid(Coverage = COV, SS = N, AR = AR, prop_brand=propbrand, VE=VE)
  full_cov = df$Coverage

  df$Coverage = 100 * df$prop_brand * df$Coverage / (100* (100 - df$Coverage) + df$prop_brand * df$Coverage)
  df$VE_lower = getVELowerCohortPrecision(alpha, df)

  df$Coverage = full_cov
  
  df = df[order(df$prop_brand, df$SS),]
  df$VE_lower = unlist(by(data = df, INDICES = df$prop_brand, FUN = smoothAndRound, xname="SS"))
  
  return(df)
  
  # pbrand=propbrand/100
  # ARi=AR/100
  # #--- pre-allocate output matrix
  # n_VE <- length(VE)
  # n_cov <- length(COV)
  # n_n <- length(N)
  # n_pbrand <- length(pbrand)
  # nn <- n_VE * n_cov * n_n *n_pbrand
  # output <- array(numeric(6*nn), dim=c(nn,6))
  # #--- perform the sample size calculation; what is the precision (lower bound of the CI) for a given number of cases
  # cnt <- 0
  # for (i in 1:n_VE){
  #   for (ii in 1:n_cov){
  #     for (iii in 1:n_n) {
  #       for (iv in 1:n_pbrand) {
  #         cnt <- cnt + 1
  #         # given parameters
  #         VEi <- VE[i]/100
  #         covi <- COV[ii]/100
  #         ni <- N[iii]
  #         pbrandi <-pbrand[iv]
  #         covi2 <- (pbrandi * covi)/(1 - covi + pbrandi * covi)   # brand-specific coverage, recalculated after excluding other brands
  #         ni2 <- (1 - covi + pbrandi * covi) * ni                 # Total sample size (cases + controls), after excluding other brands
  #         casesunexposedi    <- ni2*(1-covi2)*ARi
  #         casesexposedi      <- ni2*   covi2 *ARi*(1-VEi)
  #         unexposedi <-ni2*(1-covi2)
  #         exposedi <-ni2*covi2
  #         RRi <-1-VEi
  #         selogRRi <-sqrt(1/casesunexposedi +1/casesexposedi-(1/unexposedi + 1/exposedi))
  #         UCLRRi <- exp(log(RRi)+qnorm(1-alpha/2)*selogRRi)
  #         LCLVEi <- pmax(1-UCLRRi,0)
  #         # output
  #         output[cnt,1] <- covi*100
  #         output[cnt,2] <- ni
  #         output[cnt,3] <- AR
  #         output[cnt,4] <- pbrandi*100
  #         output[cnt,5] <- VEi * 100
  #         output[cnt,6] <- LCLVEi * 100
  #       }#iv
  #     }#iii
  #   }#ii
  # }#i
  # #--- output table
  # output2 <- data.frame(output)
  # colnames(output2) <- c("Coverage", "SS", "AR", "prop_brand","VE", "VE_lower")
  # return(output2)
}