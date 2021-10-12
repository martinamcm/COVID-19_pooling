####################################
# Outbreak simulation for HK RCHEs #
####################################


library(deSolve)
library(EpiModel)
library(gtable)
library(grid)


#' Renewal Equation
#' with seasonal forcing and
#' vital dynamics 
#' @param N subpartition of time unit. N=1 <=> SIR discrete, N=large <=> SIR continuous.
#' @param R0 Basic reproduction number
#' @param pop_size Population size (=1)
#' @param I.init Initial infectious individuals
#' @param mu Immgration and death rate (balanced)
#' @param prd Period of the seasonal forcing
#' @param amp Amplitude of the seasonal forcing
#' @param horizon Time horizon of the simulation
#' 
RE.simul.N <- function(N,
                       pop_size, 
                       I.init,
                       R0,
                       tau,
                       rel,
                       gamma.sh,
                       gamma.sc,
                       mu, 
                       horizon){
  inc <- vector()
  S <- vector()
  C <- vector()
  inc[1] <- I.init
  S[1] <- pop_size - I.init
  C[1] <- 0
  
  for(t in 2:(1 + horizon * N))
  {
    z <- 0
    for(j in 1:(t-1))
    {
      GI_j <- tran_curve_H((j-1),tau,rel,gamma.sh,gamma.sc,R0)
      z <- z + GI_j * inc[t-j] * exp(-mu/N * j) 
    } 
    R0_force <- R0 
    
    inc.tmp <- S[t-1] / pop_size * R0_force * z 
    inc[t] <- min(inc.tmp, S[t-1])
    
    S[t-1] <- S[t-1] + mu/N * (1 - S[t-1])
    S[t] <- max(0, S[t-1] - inc[t])
    
    C[t] <- C[t-1] + inc[t]
  }
  tt <- seq(0, 1+ horizon, by = 1/N )
  return(data.frame(time = tt[1:length(S)], S=S, inc=inc))
}


RE.simul.N.notest <- function(N,
                              pop_size, 
                              I.init,
                              R0,
                              tau,
                              gamma.sh,
                              gamma.sc,
                              mu, 
                              horizon){
  inc <- vector()
  S <- vector()
  C <- vector()
  inc[1] <- I.init
  S[1] <- pop_size - I.init
  C[1] <- 0
  
  for(t in 2:(1 + horizon * N))
  {
    z <- 0
    for(j in 1:(t-1))
    {
      GI_j <- tran_curve_notest((j-1),tau,gamma.sh,gamma.sc,R0)
      z <- z + GI_j * inc[t-j] * exp(-mu/N * j) 
    } 
    R0_force <- R0 
    inc.tmp <- S[t-1] / pop_size * R0_force * z  
    inc[t] <- min(inc.tmp, S[t-1])
    
    S[t-1] <- S[t-1] + mu/N * (1 - S[t-1])
    S[t] <- max(0, S[t-1] - inc[t])
    
    C[t] <- C[t-1] + inc[t]
  }
  tt <- seq(0, 1+ horizon, by = 1/N )
  return(data.frame(time = tt[1:length(S)], S=S, inc=inc))
}


## Assume 50 staff in RCHE with 100 beds tested once every 4 weeks ind testing 
disc <- 1
RCHEpop <- 150
no.inf <- 1
timeframe <- 30
BDrate <- 0
sh <- sh_Park <- 1.960923
sc <- sc_Park <- 4.334694
#sh <- sh_Li = 4.8685
#sc <- sc_Li = 1.821917
 shWu = 2.419753
 scWu = 2.892857

det_symp <- function(tau,se,n,inc){
  time <- min((1/se-0.5)*tau,(1/0.64)*inc)+1
  return(time)
}

sens_HCWp1 <- sens_HCW[sens_HCW$Prevalence==2e-04,]

################################################################################
# Cumulative number of cases under different R0 and different testing resources
################################################################################

###Scenario 1
R0vect <- c(1.8,2.5,3.6)

tau_ind <- 7

### Define number of cases across time for different R0 values when inc=7.76, testing frequency=28days

notestarray <- sapply(R0vect,function(x){RE.simul.N.notest(disc,RCHEpop,no.inf,x,tau_ind,sh,sc,BDrate,timeframe)})
indarray <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind,1,sh,sc,BDrate,timeframe)})
npool2array <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],sh,sc,BDrate,timeframe)})
npool5array <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],sh,sc,BDrate,timeframe)})
npool10array <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],sh,sc,BDrate,timeframe)})

n.array <- c(notestarray,indarray,npool2array,npool5array,npool10array)

CRl=RCHEpop-c(notestarray[,1]$S,indarray[,1]$S,npool2array[,1]$S,npool5array[,1]$S,npool10array[,1]$S)
CR=RCHEpop-c(notestarray[,2]$S,indarray[,2]$S,npool2array[,2]$S,npool5array[,2]$S,npool10array[,2]$S)
CRh=RCHEpop-c(notestarray[,3]$S,indarray[,3]$S,npool2array[,3]$S,npool5array[,3]$S,npool10array[,3]$S)

### Scenario 2

tau_ind2 <- 14

notestarray2 <- sapply(R0vect,function(x){RE.simul.N.notest(disc,RCHEpop,no.inf,x,tau_ind,sh,sc,BDrate,timeframe)})
indarray2 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind,1,sh,sc,BDrate,timeframe)})
npool2array2 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],sh,sc,BDrate,timeframe)})
npool5array2 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],sh,sc,BDrate,timeframe)})
npool10array2 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],sh,sc,BDrate,timeframe)})

CRl.2=RCHEpop-c(notestarray2[,1]$S,indarray2[,1]$S,npool2array2[,1]$S,npool5array2[,1]$S,npool10array2[,1]$S)
CR.2=RCHEpop-c(notestarray2[,2]$S,indarray2[,2]$S,npool2array2[,2]$S,npool5array2[,2]$S,npool10array2[,2]$S)
CRh.2=RCHEpop-c(notestarray2[,3]$S,indarray2[,3]$S,npool2array2[,3]$S,npool5array2[,3]$S,npool10array2[,3]$S)

### Scenario 3

tau_ind3 <- 21

notestarray3 <- sapply(R0vect,function(x){RE.simul.N.notest(disc,RCHEpop,no.inf,x,tau_ind,sh,sc,BDrate,timeframe)})
indarray3 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind,1,sh,sc,BDrate,timeframe)})
npool2array3 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==2],sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],sh,sc,BDrate,timeframe)})
npool5array3 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==5],sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],sh,sc,BDrate,timeframe)})
npool10array3 <- sapply(R0vect,function(x){RE.simul.N(disc,RCHEpop,no.inf,x,tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==10],sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],sh,sc,BDrate,timeframe)})

CRl.3=RCHEpop-c(notestarray3[,1]$S,indarray3[,1]$S,npool2array3[,1]$S,npool5array3[,1]$S,npool10array3[,1]$S)
CR.3=RCHEpop-c(notestarray3[,2]$S,indarray3[,2]$S,npool2array3[,2]$S,npool5array3[,2]$S,npool10array3[,2]$S)
CRh.3=RCHEpop-c(notestarray3[,3]$S,indarray3[,3]$S,npool2array3[,3]$S,npool5array3[,3]$S,npool10array3[,3]$S)

tran_curve_H(a,tau,rel,gamma.sh,gamma.sc,R0)

testing <- c(rep("notest",dim(notest)[1]),rep("ind",dim(ind)[1]),rep("npool2",dim(npool2)[1]),rep("npool5",dim(npool4)[1]),rep("npool10",dim(npool8)[1]))
time=rep(ind$time,5)





############################################################################################
# Expected time of detection under different incubation periods and percentage asymptomatic
############################################################################################

#Stochastic incubation & asymptomatic evaluation?

## Tau=28
incvect <- c(5.2,6.7,8.2)

Enotest <- sapply(incvect,function(x){(1/0.64)*x})
Eind <- sapply(incvect,function(x){median(det_symp(tau_ind,sensind,n,x))})
Enpool2 <- sapply(incvect,function(x){median(det_symp(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],n,x))})
Enpool4 <- sapply(incvect,function(x){median(det_symp(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==4],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==4],n,x))})
Enpool8 <- sapply(incvect,function(x){median(det_symp(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==8],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==8],n,x))})
Enpool14 <- sapply(incvect,function(x){median(det_symp(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==14],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==14],n,x))})

Exp.Rl <- c(Enotest[1],Eind[1],Enpool2[1],Enpool4[1],Enpool8[1],Enpool14[1])
Exp.R <- c(Enotest[2],Eind[2],Enpool2[2],Enpool4[2],Enpool8[2],Enpool14[2])
Exp.Rh <- c(Enotest[3],Eind[3],Enpool2[3],Enpool4[3],Enpool8[3],Enpool14[3])

## Tau=14
Enotest.2 <- sapply(incvect,function(x){(1/0.64)*x})
Eind.2 <- sapply(incvect,function(x){median(det_symp(tau_ind2,sensind,n,x))})
Enpool2.2 <- sapply(incvect,function(x){median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],n,x))})
Enpool4.2 <- sapply(incvect,function(x){median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==4],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==4],n,x))})
Enpool8.2 <- sapply(incvect,function(x){median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==8],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==8],n,x))})
Enpool14.2 <- sapply(incvect,function(x){median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==14],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==14],n,x))})

Exp2.Rl <- c(Enotest.2[1],Eind.2[1],Enpool2.2[1],Enpool4.2[1],Enpool8.2[1],Enpool14.2[1])
Exp2.R <- c(Enotest.2[2],Eind.2[2],Enpool2.2[2],Enpool4.2[2],Enpool8.2[2],Enpool14.2[2])
Exp2.Rh <- c(Enotest.2[3],Eind.2[3],Enpool2.2[3],Enpool4.2[3],Enpool8.2[3],Enpool14.2[3])

##Tau=7
Enotest.3 <- sapply(incvect,function(x){(1/0.64)*x})
Eind.3 <- sapply(incvect,function(x){median(det_symp(tau_ind3,sensind,n,x))})
Enpool2.3 <- sapply(incvect,function(x){median(det_symp(tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==2],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],n,x))})
Enpool4.3 <- sapply(incvect,function(x){median(det_symp(tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==4],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==4],n,x))})
Enpool8.3 <- sapply(incvect,function(x){median(det_symp(tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==8],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==8],n,x))})
Enpool14.3 <- sapply(incvect,function(x){median(det_symp(tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==14],sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==14],n,x))})

Exp3.Rl <- c(Enotest.3[1],Eind.3[1],Enpool2.3[1],Enpool4.3[1],Enpool8.3[1],Enpool14.3[1])
Exp3.R <- c(Enotest.3[2],Eind.3[2],Enpool2.3[2],Enpool4.3[2],Enpool8.3[2],Enpool14.3[2])
Exp3.Rh <- c(Enotest.3[3],Eind.3[3],Enpool2.3[3],Enpool4.3[3],Enpool8.3[3],Enpool14.3[3])


### Cases for Tau=28
times <- sequence(0,30,1)

Cases1.Rla <- sapply(Exp.Rl,function(x){CRl[ceiling(x)+1]})
Cases1.Ra <- sapply(Exp.Rl,function(x){CR[ceiling(x)+1]})
Cases1.Rha <- sapply(Exp.Rl,function(x){CRh[ceiling(x)+1]})

Cases1.Rlb <- sapply(Exp.R,function(x){CRl[ceiling(x)+1]})
Cases1.Rb <- sapply(Exp.R,function(x){CR[ceiling(x)+1]})
Cases1.Rhb <- sapply(Exp.R,function(x){CRh[ceiling(x)+1]})

Cases1.Rlc <- sapply(Exp.Rh,function(x){CRl[ceiling(x)+1]})
Cases1.Rc <- sapply(Exp.Rh,function(x){CR[ceiling(x)+1]})
Cases1.Rhc <- sapply(Exp.Rh,function(x){CRh[ceiling(x)+1]})

#### Cases for Tau=14

Cases2.Rla <- sapply(Exp2.Rl,function(x){CRl.2[ceiling(x)+1]})
Cases2.Ra <- sapply(Exp2.Rl,function(x){CR.2[ceiling(x)+1]})
Cases2.Rha <- sapply(Exp2.Rl,function(x){CRh.2[ceiling(x)+1]})

Cases2.Rlb <- sapply(Exp2.R,function(x){CRl.2[ceiling(x)+1]})
Cases2.Rb <- sapply(Exp2.R,function(x){CR.2[ceiling(x)+1]})
Cases2.Rhb <- sapply(Exp2.R,function(x){CRh.2[ceiling(x)+1]})

Cases2.Rlc <- sapply(Exp2.Rh,function(x){CRl.2[ceiling(x)+1]})
Cases2.Rc <- sapply(Exp2.Rh,function(x){CR.2[ceiling(x)+1]})
Cases2.Rhc <- sapply(Exp2.Rh,function(x){CRh.2[ceiling(x)+1]})

### Cases for Tau=7

Cases3.Rla <- sapply(Exp3.Rl,function(x){CRl.3[ceiling(x)+1]})
Cases3.Ra <- sapply(Exp3.Rl,function(x){CR.3[ceiling(x)+1]})
Cases3.Rha <- sapply(Exp3.Rl,function(x){CRh.3[ceiling(x)+1]})

Cases3.Rlb <- sapply(Exp3.R,function(x){CRl.3[ceiling(x)+1]})
Cases3.Rb <- sapply(Exp3.R,function(x){CR.3[ceiling(x)+1]})
Cases3.Rhb <- sapply(Exp3.R,function(x){CRh.3[ceiling(x)+1]})

Cases3.Rlc <- sapply(Exp3.Rh,function(x){CRl.3[ceiling(x)+1]})
Cases3.Rc <- sapply(Exp3.Rh,function(x){CR.3[ceiling(x)+1]})
Cases3.Rhc <- sapply(Exp3.Rh,function(x){CRh.3[ceiling(x)+1]})


Cases.facet <- c(Cases1.Ra, Cases1.Rb, Cases1.Rc, Cases2.Ra, Cases2.Rb, Cases2.Rc, Cases3.Ra,  Cases3.Rb, Cases3.Rc)

CaseRlow <-  c(Cases1.Rla, Cases1.Rlb, Cases1.Rlc, Cases2.Rla, Cases2.Rlb,Cases2.Rlc, Cases3.Rla, Cases3.Rlb, Cases3.Rlc)
  
CaseRhigh <-c(Cases1.Rha, Cases1.Rhb, Cases1.Rhc, Cases2.Rha, Cases2.Rhb, Cases2.Rhc, Cases3.Rha, Cases3.Rhb, Cases3.Rhc)


Freq.facet <- factor(c(rep("28",18),rep("14",18),rep("7",18)),levels=c("28","14","7"))
Incub.f <- factor(c(rep("5.2",6),rep("6.7",6),rep("8.2",6)),levels=c("5.2","6.7","8.2"))
Incub.facet <- c(rep(Incub.f,3))
Test.f<- factor(c("notest", "ind", "npool2", "npool4",  "npool8",  "npool14"),levels=c("notest", "ind", "npool2", "npool4",  "npool8",  "npool14"))
Test.facet <- c(rep(Test.f,9))

data_facet <- data.frame(Test.facet,Cases.facet,CaseRlow,CaseRhigh,Freq.facet,Incub.facet)

label.legend <- c("1"="No screening","2"="Without pooling","3"=expression(n[pool]*"=2"),"4"=expression(n[pool]*"=4"),"5"=expression(n[pool]*"=8"),"6"=expression(n[pool]*"=14"))
label.inc <- c("1"="5.2 days","2"="6.7 days","3"="8.2 days")
label.freq <- c("28"="Low resources","14"="Medium resources","7"="High resources")

 ggplot(data=data_facet,aes(x=Test.facet,y=Cases.facet,ymin=CaseRlow,ymax=CaseRhigh,color=as.factor(Test.facet)))+geom_pointrange()+scale_color_viridis(discrete=TRUE,labels=label.legend)+theme_bw()+
  facet_grid(Freq.facet~Incub.facet,labeller = labeller(Incub.facet=label.inc,Freq.facet=label.freq))+theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10),
        legend.position="bottom",axis.text=element_text(size=10),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title=element_text(size=11))+labs(color="Screening")+xlab("")+ylab("Size of outbreak")#+removeGridX()

#annotate_figure(outbreak,top=text_grob("Mean incubation period",size=10),right=text_grob("Frequency of testing without pooling",size = 10,rot = -90))

 

