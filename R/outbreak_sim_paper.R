#####################################################################
# Determine cumulative cases over time under each testing strategy 
#
#####################################################################



# 1. Renewal equation functions ----------------------------------------------


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


## Renewal equation for syndromic surveillance

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



## Renewal equation for rapid testing


## Renewal equation for syndromic surveillance

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
      GI_j <- tran_curve_rapid((j-1),tau,rel,gamma.sh,gamma.sc,R0)
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



# 2.  Simulation parameters -----------------------------------------------


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

 
# set prevalence equal to 0.02%

sens_HCWp1 <- sens_HCW[sens_HCW$Prevalence==2e-04,] 



# 3.  Cumulative number of cases  -----------------------------------------


## Cumulative number of cases under different R0 and different testing resources


### Scenario 1: testing every 7 days

#### Parameters

R0vect <- c(1.8, 2.5, 3.6)
tau_ind <- 7


#### Define number of cases across time for different R0 values when inc = 7.76, 
#### testing frequency = 28days

# Syndromic
notestarray <- sapply(R0vect,
                      function(x){
                        RE.simul.N.notest(disc,
                                          RCHEpop,
                                          no.inf,
                                          x,
                                          tau_ind,
                                          sh,
                                          sc,
                                          BDrate,
                                          timeframe
                                          )
                        }
                      )

# Unpooled testing 
indarray <- sapply(R0vect,
                   function(x){
                     RE.simul.N(disc, 
                                RCHEpop,
                                no.inf,
                                x,
                                tau_ind,
                                1,
                                sh,
                                sc,
                                BDrate,
                                timeframe
                                )
                     }
                   )

#npool = 2
npool2array <- sapply(R0vect,
                      function(x){
                        RE.simul.N(disc,
                                   RCHEpop,
                                   no.inf,
                                   x,
                                   tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                   sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                   sh,
                                   sc,
                                   BDrate,
                                   timeframe
                                   )
                        }
                      )

#npool = 5
npool5array <- sapply(R0vect,
                      function(x){
                        RE.simul.N(disc,
                                   RCHEpop,
                                   no.inf,
                                   x,
                                   tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                   sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                   sh,
                                   sc,
                                   BDrate,
                                   timeframe
                                   )
                        }
                      )

# npool = 10
npool10array <- sapply(R0vect,
                       function(x){
                         RE.simul.N(disc,
                                    RCHEpop,
                                    no.inf,
                                    x,
                                    tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe
                                    )
                         }
                       )

# combine 
n.array <- c(notestarray,indarray,npool2array,npool5array,npool10array)


## Cumulative cases

# cumulative R low
CRl = RCHEpop - c(notestarray[,1]$S, 
                  indarray[,1]$S, 
                  npool2array[,1]$S,
                  npool5array[,1]$S,
                  npool10array[,1]$S
                  )

# cumulative R mid
CR = RCHEpop - c(notestarray[,2]$S,
                 indarray[,2]$S,
                 npool2array[,2]$S,
                 npool5array[,2]$S,
                 npool10array[,2]$S
                 )

# cumulative R high
CRh = RCHEpop - c(notestarray[,3]$S,
                  indarray[,3]$S,
                  npool2array[,3]$S,
                  npool5array[,3]$S,
                  npool10array[,3]$S
                  )



### Scenario 2: testing every 14 days

tau_ind2 <- 14

# syndromic 
notestarray2 <- sapply(R0vect,
                       function(x){
                         RE.simul.N.notest(disc,
                                           RCHEpop,
                                           no.inf,
                                           x,
                                           tau_ind,
                                           sh,
                                           sc,
                                           BDrate,
                                           timeframe
                                           )
                         }
                       )


# unpooled 
indarray2 <- sapply(R0vect,
                    function(x){
                      RE.simul.N(disc,
                                 RCHEpop,
                                 no.inf,
                                 x,
                                 tau_ind,
                                 1,
                                 sh,
                                 sc,
                                 BDrate,
                                 timeframe
                                 )
                      }
                    )


# npool = 2
npool2array2 <- sapply(R0vect,
                       function(x){
                         RE.simul.N(disc,
                                    RCHEpop,
                                    no.inf,
                                    x,
                                    tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe
                                    )
                         }
                       )

npool5array2 <- sapply(R0vect,
                       function(x){
                         RE.simul.N(disc,
                                    RCHEpop,
                                    no.inf,
                                    x,
                                    tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe
                                    )
                         }
                       )

npool10array2 <- sapply(R0vect,
                        function(x){
                          RE.simul.N(disc,
                                     RCHEpop,
                                     no.inf,
                                     x,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe
                                     )
                          }
                        )


# Cumulative cases

# R low
CRl.2 = RCHEpop - c(notestarray2[,1]$S,
                    indarray2[,1]$S,
                    npool2array2[,1]$S,
                    npool5array2[,1]$S,
                    npool10array2[,1]$S
                    )

# R mid
CR.2 = RCHEpop - c(notestarray2[,2]$S,
                   indarray2[,2]$S,
                   npool2array2[,2]$S,
                   npool5array2[,2]$S,
                   npool10array2[,2]$S
                   )

# R high
CRh.2 = RCHEpop - c(notestarray2[,3]$S,
                    indarray2[,3]$S,
                    npool2array2[,3]$S,
                    npool5array2[,3]$S,
                    npool10array2[,3]$S
                    )



### Scenario 3: testing every 21 days

tau_ind3 <- 21

# syndromic surveillance
notestarray3 <- sapply(R0vect,
                       function(x){
                         RE.simul.N.notest(disc,
                                           RCHEpop,
                                           no.inf,
                                           x,
                                           tau_ind,
                                           sh,
                                           sc,
                                           BDrate,
                                           timeframe
                                           )
                         }
                       )

# unpooled testing 
indarray3 <- sapply(R0vect,
                    function(x){
                      RE.simul.N(disc,
                                 RCHEpop,
                                 no.inf,
                                 x,
                                 tau_ind,
                                 1,
                                 sh,
                                 sc,
                                 BDrate,
                                 timeframe
                                 )
                      }
                    )

# npool = 2
npool2array3 <- sapply(R0vect,
                       function(x){
                         RE.simul.N(disc,
                                    RCHEpop,
                                    no.inf,
                                    x,
                                    tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe
                                    )
                         }
                       )

# npool = 5
npool5array3 <- sapply(R0vect,
                       function(x){
                         RE.simul.N(disc,
                                    RCHEpop,
                                    no.inf,
                                    x,
                                    tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe
                                    )
                         }
                       )


# npool = 10
npool10array3 <- sapply(R0vect,
                        function(x){
                          RE.simul.N(disc,
                                     RCHEpop,
                                     no.inf,
                                     x,
                                     tau_ind3*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe
                                     )
                          }
                        )


## Cumulative cases 

# R low
CRl.3 = RCHEpop - c(notestarray3[,1]$S,
                    indarray3[,1]$S,
                    npool2array3[,1]$S,
                    npool5array3[,1]$S,
                    npool10array3[,1]$S
                    )

# R mid
CR.3 = RCHEpop - c(notestarray3[,2]$S,
                   indarray3[,2]$S,
                   npool2array3[,2]$S,
                   npool5array3[,2]$S,
                   npool10array3[,2]$S
                   )

# R high
CRh.3 = RCHEpop - c(notestarray3[,3]$S,
                    indarray3[,3]$S,
                    npool2array3[,3]$S,
                    npool5array3[,3]$S,
                    npool10array3[,3]$S
                    )




testing <-  c(rep("notest",dim(notest)[1]),
              rep("ind",dim(ind)[1]),
              rep("npool2",dim(npool2)[1]),
              rep("npool5",dim(npool4)[1]),
              rep("npool10",dim(npool8)[1])
              )

time = rep(ind$time,5)



