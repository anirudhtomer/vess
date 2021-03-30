###########################################
# Case control-overall-min det VE         #
###########################################
Casecontrol_overall_precision<-function(alpha=0.05,                     #Significance level (Input=numeric)
                                        VE=c(20,50,70),                 #True vaccine effectiveness(Input=vector/numeric)
                                        COV=c(5, 20, 50, 70),           #Coverage (Input=vector/numeric)
                                        r=1,                            #Controls per case (Input=numeric)
                                        ncases=seq(100, 4000, 100))     #Total number of cases in the study (Input=vector)
{
  df = expand.grid(Coverage = COV, nr_cases = ncases, r = r, VE=VE)
  df$VE_lower = getVELowerCaseControlPrecision(alpha, df)
  
  df = df[order(df$Coverage, df$nr_cases),]
  df$VE_lower = unlist(by(data = df, INDICES = df$Coverage, FUN = smoothAndRound, xname="nr_cases"))
  
  return(df)
  
  # #--- pre-allocate output matrix
  # n_VE <- length(VE)
  # n_cov <- length(COV)
  # n_ncases <- length(ncases)
  # n_r <- length(r)
  # nn <- n_VE * n_cov * n_ncases * n_r
  # output <- array(numeric(5*nn), dim=c(nn,5))
  # #--- perform the sample size calculation; what is the precision (lower bound of the CI) for a given number of cases
  # cnt <- 0
  # for (i in 1:n_VE){
  #   for (ii in 1:n_cov) {
  #     for (iii in 1:n_ncases) {
  #       cnt <- cnt + 1
  #       # given parameters
  #       VEi <- VE[i]/100
  #       covi <- COV[ii]/100
  #       ncasesi <- ncases[iii]
  #       ri <- r
  #       ncontrolsi <- ncasesi * ri
  #       ni <- ncasesi + ncontrolsi
  #       ORi <- 1 - VEi
  #       vaci <-ni*covi
  #       nonvaci <-ni*(1-covi)
  #       # calculate 2x2 table's parameters
  #       a<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[1])
  #       b<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[2])
  #       c<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[3])
  #       d<-as.numeric(findpar(a=vaci,b=nonvaci,c=ncasesi,d=ncontrolsi,R=ORi)[4])
  #       # calculate lower CI of VE
  #       SE_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)
  #       OR_UL <- exp(log(ORi) + qnorm(1-alpha/2) * SE_logOR)
  #       VE_LL <- pmax(1 - OR_UL,0)
  #       # output
  #       output[cnt,1] <- covi*100
  #       output[cnt,2] <- ncasesi
  #       output[cnt,3] <- ri
  #       output[cnt,4] <- VEi * 100
  #       output[cnt,5] <- VE_LL * 100
  #     }#iii
  #   }#ii
  # }#i
  # #--- output table
  # output2 <- data.frame(output)
  # colnames(output2) <- c("Coverage", "nr_cases",  "r","VE", "VE_lower")
  # return(output2)
}
