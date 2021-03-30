
###########################################
# Case control-brand specific-precision   #
###########################################

Casecontrol_brandspecific_precision<-function(alpha=0.05,                     #Significance level (Input=numeric)
                                              VE=c(20,50,70),                 #True vaccine effectiveness(Input=vector/numeric)
                                              COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                              r=1,                            #Controls per case (Input=numeric)
                                              ncases=seq(100, 4000, 100),     #Total number of cases in the study (Input=vector) 
                                              propbrand=c(10,20,40,60,80,90)) #Proportion of the brand (Input=vector/numeric) 
{
  
  df = expand.grid(Coverage = COV, nr_cases = ncases, r = r, prop_brand=propbrand, VE=VE)
  full_cov = df$Coverage
  df$Coverage = 100 * df$prop_brand * df$Coverage / (100* (100 - df$Coverage) + df$prop_brand * df$Coverage)

  df$VE_lower = getVELowerCaseControlPrecision(alpha, df)
  df$Coverage = full_cov
  
  df = df[order(df$prop_brand, df$nr_cases),]
  df$VE_lower = unlist(by(data = df, INDICES = df$prop_brand, FUN = smoothAndRound, xname="nr_cases"))

  return(df)
  
  # pbrand=propbrand/100
  # #--- pre-allocate output matrix
  # n_VE <- length(VE)
  # n_cov <- length(COV)
  # n_ncases <- length(ncases)
  # n_r <- length(r)
  # n_pbrand <- length(pbrand)
  # nn <- n_VE * n_cov * n_ncases* n_r *n_pbrand
  # output <- array(numeric(6*nn), dim=c(nn,6))
  # #--- perform the sample size calculation; what is the precision (lower bound of the CI) for a given number of cases
  # cnt <- 0
  # for (i in 1:n_VE){
  #   for (ii in 1:n_cov){
  #     for (iii in 1:n_ncases) {
  #       for (v in 1:n_pbrand) {
  #         cnt <- cnt + 1
  #         # given parameters
  #         VEi <- VE[i]/100
  #         covi <- COV[ii]/100
  #         ncasesi <- ncases[iii]
  #         ri <- r
  #         pbrandi <-pbrand[v]
  #         ORi <- 1 - VEi
  #         ncontrolsi <- ncasesi * ri
  #         ni <- ncasesi + ncontrolsi
  #         vaci <-ni*covi
  #         nonvaci <-ni*(1-covi) #It does not change with brand specific precision
  #         # calculate original 2x2 table's parameters
  #         a<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[1])
  #         b<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[2])
  #         c<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[3])
  #         d<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[4])
  #         a2<-pbrandi*a
  #         b2<-pbrandi*b
  #         vaci2 <-ni*covi*pbrandi
  #         ni2 <- nonvaci+vaci2
  #         ncasesi2 <-a2+c
  #         ncontrolsi2 <- ni2-ncasesi2
  #         # calculate lower CI of VE
  #         SE_logOR <- sqrt(1/a2 + 1/b2 + 1/c + 1/d)
  #         OR_UL <- exp(log(ORi) + qnorm(1-alpha/2) * SE_logOR)
  #         VE_LL <- pmax(1 - OR_UL,0)
  #         # output
  #         output[cnt,1] <- covi*100
  #         output[cnt,2] <- ncasesi
  #         output[cnt,3] <- ri
  #         output[cnt,4] <- pbrandi*100
  #         output[cnt,5] <- VEi * 100
  #         output[cnt,6] <- VE_LL * 100
  #       }
  #     }
  #   }
  # }
  # #--- output table
  # output2 <- data.frame(output)
  # colnames(output2) <- c("Coverage", "nr_cases", "r","prop_brand","VE", "VE_lower")
  # return(output2)
}