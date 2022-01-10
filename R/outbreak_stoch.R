#####################################################################
# Determine expected size of the outbreak under each testing strategy 
#
#####################################################################



# 1.  Cumulative number of cases - stochastic  --------------------------------



inf_time <- function(n,result,RCHEpop){
 
mean <- vector()
lower <- vector()
upper <- vector()

pois_mean <- 
  RCHEpop - sapply(1:3,
                   function(x){
                     result[,x]$S
                     }
                   )
#pois_mean <- RCHEpop-sapply(1:3,function(x){result$S})  

  samp <- 
    apply(pois_mean,
          c(1,2),
          rpois,
          n=100
          )
  
  meanRl <- apply(samp[,,1],
                  2,
                  median
                  )
  
  lowerRl <- apply(samp[,,1],
                   2,
                   quantile,
                   probs = 0.025
                   )
    
  upperRl <- apply(samp[,,1],
                   2,
                   quantile,
                   probs = 0.975
                   )

  meanR <- apply(samp[,,2],
                 2,
                 median
                 )
  
  lowerR <- apply(samp[,,2],
                  2,
                  quantile,
                  probs = 0.025
                  )
  
  upperR <- apply(samp[,,2],
                  2,
                  quantile,
                  probs = 0.975
                  )
  
  meanRh <- apply(samp[,,3],
                  2,
                  median
                  )
    
  lowerRh <- apply(samp[,,3],
                   2,
                   quantile,
                   probs = 0.025
                   )
  
  upperRh <- apply(samp[,,3],
                   2,
                   quantile,
                   probs = 0.975
                   )
  
 Rtype <- c(rep("R0=1.8",dim(pois_mean)[1]),
            rep("R0=2.5",dim(pois_mean)[1]),
            rep("R0=3.6",dim(pois_mean)[1])
            )
 
 Rtype <- factor(Rtype,
                 levels = c("R0=1.8",
                            "R0=2.5",
                            "R0=3.6"
                            )
                 )
 
 

 return(
  data.frame(times = rep(1:dim(pois_mean)[1],3),
             c(meanRl,meanR,meanRh),
             c(lowerRl,lowerR,lowerRh),
             c(upperRl,upperR,upperRh),Rtype
             )
  )
 }


# cumulative cases for tau = 7

sched_0 <- inf_time(50000,notestarray,150)
sched_1 <- inf_time(50000,indarray,150)
sched_2 <- inf_time(50000,npool2array,150)
sched_3 <- inf_time(50000,npool5array,150)
sched_4 <- inf_time(50000,npool10array,150)


type <- 
  c(
    rep("synd",dim(sched_0)[1]),
    rep("ind",dim(sched_1)[1]),
    rep("npool2",dim(sched_2)[1]),
    rep("npool5",dim(sched_3)[1]),
    rep("npool10",dim(sched_4)[1])
    )

type <- 
  factor(
    type,
    levels = c("synd","ind","npool2","npool5","npool10")
    )


data_plot1 <- rbind(sched_0,sched_1,sched_2,sched_3,sched_4)
data_plot1 <- cbind(type,data_plot1)



# cumulative cases for tau = 14

sched_5 <- inf_time(50000,notestarray2,150)
sched_6 <- inf_time(50000,indarray2,150)
sched_7 <- inf_time(50000,npool2array2,150)
sched_8 <- inf_time(50000,npool5array2,150)
sched_9 <- inf_time(50000,npool10array2,150)

data_plot2 <- rbind(sched_5, 
                    sched_6, 
                    sched_7, 
                    sched_8, 
                    sched_9
                    )

data_plot2 <- cbind(type,
                    data_plot2
                    )



# cumulative cases for tau = 28

sched_10 <- inf_time(50000,notestarray3,150)
sched_11 <- inf_time(50000,indarray3,150)
sched_12 <- inf_time(50000,npool2array3,150)
sched_13 <- inf_time(50000,npool5array3,150)
sched_14 <- inf_time(50000,npool10array3,150)

data_plot3 <- rbind(sched_10,sched_11,sched_12,sched_13,sched_14)
data_plot3 <- cbind(type,data_plot3)




# 2.  Plot cumulative cases over time ----------------------------------------


## Make facet data

resource <- 
  c(
    rep("High resources",dim(data_plot1)[1]),
    rep("Medium resources",dim(data_plot2)[1]),
    rep("Low resources",dim(data_plot3)[1])
    )

resource <- 
  factor(
    resource,
    levels = c("High resources","Medium resources","Low resources")
    )

data_facetplot <- rbind(data_plot1,data_plot2,data_plot3)
data_facetplot <- cbind(data_facetplot,resource)

names(data_facetplot) <- c("type","times","mean","lower","upper","Rtype","resource")

label.legend <- 
  c("synd" = "Syndromic",
    "ind" = "No pooling",
    "npool2" = expression(n[pool]*"=2"),
    "npool5" = expression(n[pool]*"=5"),
    "npool10" = expression(n[pool]*"=10")
    )

data_facetplot$upper <- 
  ifelse(
    data_facetplot$upper <= 150,
    data_facetplot$upper,
    150
    )

data_facetplot$mean  <- 
  ifelse(
    data_facetplot$mean <= 150,
    data_facetplot$mean,
    150
    )



ggplot(
  data_facetplot, 
  aes(
    x = times,
    fill=type
    )
  ) +
  geom_line(
    data_facetplot,
    mapping = aes(
      y = mean,
      x = times
      ),
    size=0.5
    ) + 
  facet_grid(
    resource~Rtype
    ) +
  geom_ribbon(
    data_facetplot,
    mapping = aes(
      x = times,
      ymin = lower,
      ymax = upper
      ), 
    alpha=0.3
    ) +
  xlab("Time since first infection (days)") + 
  ylab('Cumulative number of cases') + 
  scale_fill_brewer(
    palette = "Dark2",
    labels=label.legend
    ) + 
  theme_bw() + 
  labs(fill = "Screening") + 
  theme(legend.position = "bottom")





# 3.  Expected time of detection T under each testing strategy ------------



det_symp <- function(inc,symp){
  
  time <- (1/symp)*inc[1]+1
  timelow <- (1/symp)*inc[2]+1
  timehigh <- (1/symp)*inc[3]+1
  
  return(c(time,timelow,timehigh))
}


det_symp <- function(inc,symp){
  
  rweibull(1000,1.91,0.11)
  
  time <- (1/symp)*inc[1]+1
  timelow <- (1/symp)*inc[2]+1
  timehigh <- (1/symp)*inc[3]+1
  
  return(c(time,timelow,timehigh))
}

E_synd <- ceiling(det_symp(c(6.7,5.2,8.2),
                           0.65
                           )
                  )



# Low resources

tau_ind <- 21

indlow <- sapply(a,
                 function(x){
                   Tcut_Hlow(x,tau_ind,1)
                   }
                 )

E_ind <- mean(density(indlow)$y
              ) + 1

low_pool2 <- sapply(a,
                    function(x){
                      Tcut_Hlow(x,
                                ceiling(
                                  tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2]
                                  ),
                                se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2]
                      )
                      }
                    )

E_low2 <- mean(density(low_pool2)$y
               ) + 1


low_pool5 <- sapply(a, 
                    function(x){
                      Tcut_Hlow(x,
                                ceiling(
                                  tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5]
                                  ),
                                se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==5]
                      )
                      }
                    )


E_low5 <- mean(density(low_pool5)$y
               ) + 1


low_pool10 <- sapply(a,
                     function(x){
                       Tcut_Hlow(x,
                                 ceiling(
                                   tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10]
                                   ),
                                 se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==10]
                                 )
                       }
                     )


E_low10 <- mean(density(low_pool10)$y
                ) + 1


EtestLa <- c(E_ind,
             E_low2,
             E_low5,
             E_low10
             )

EtestL1 <- c(E_synd[1],
             sapply(EtestLa,
                    function(x){
                      min(E_synd[1], x)
                      }
                    )
             )

EtestL2 <- c(E_synd[2],
             sapply(EtestLa,
                    function(x){
                      min(E_synd[2], x)
                      }
                    )
             )


EtestL3 <- c(E_synd[3],
             sapply(EtestLa,
                    function(x){
                      min(E_synd[3], x)
                      }
                    )
             )


EtestL <- c(EtestL1,
            EtestL2,
            EtestL3
            )



########### Medium resources

tau_ind <- 14

indmed <- sapply(a,
                 function(x){
                   Tcut_Hlow(x,
                             tau_ind,
                             1)
                   }
                 )


E_indM <- mean(density(indmed)$y
               ) + 1

med_pool2 <- sapply(a,
                    function(x){
                      Tcut_Hlow(
                        x,
                        tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
                        se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2])
                      }
                    )


E_med2 <- mean(density(med_pool2)$y
               ) + 1


med_pool5 <- sapply(a,
                    function(x){
                      Tcut_Hlow(x,
                                ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5]),
                                se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==5]
                                )
                      }
                    )

E_med5 <- mean(density(med_pool5)$y
               ) + 1

med_pool10 <- sapply(a,
                     function(x){
                       Tcut_Hlow(x,
                                 ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10]),
                                 se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==10]
                                )
                       }
                     )


E_med10 <- mean(density(med_pool10)$y
                ) + 1


EtestMa <- c(E_indM,
             E_med2,
             E_med5,
             E_med10
             )

EtestM1 <- c(E_synd[1],
             sapply(EtestMa,
                    function(x){
                      min(E_synd[1],x)
                      }
                    )
             )

EtestM2 <- c(E_synd[2],
             sapply(EtestMa,
                    function(x){
                      min(E_synd[2],x)
                      }
                    )
             )

EtestM3 <- c(E_synd[3],
             sapply(EtestMa,
                    function(x){
                      min(E_synd[3],x)
                      }
                    )
             )

EtestM <- c(EtestM1,
            EtestM2,
            EtestM3
            )



########## High resources

tau_ind <- 7


indhigh <- sapply(a,
                  function(x){
                    Tcut_Hlow(x,tau_ind,1)
                    }
                  )

E_indH <- mean(density(indhigh)$y
               ) + 1

high_pool2 <- sapply(a,
                     function(x){
                       Tcut_Hlow(x,
                                 ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2]),
                                 se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2]
                                 )
                       }
                     )
E_high2 <- mean(density(high_pool2)$y
                ) + 1

high_pool5 <- sapply(a,
                     function(x){
                       Tcut_Hlow(x,
                                 ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5]),
                                 se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==5]
                                )
                       }
                     )

E_high5 <- mean(density(high_pool5)$y
                ) + 1

high_pool10 <- sapply(a,
                      function(x){
                        Tcut_Hlow(x,
                                  ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10]),
                                  se_time_H(x)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==10]
                        )
                        }
                      )

E_high10 <- mean(density(high_pool10)$y
                 ) + 1

EtestHa <- c(E_indH,
             E_high2,
             E_high5,
             E_high10
             )

EtestH1 <- c(E_synd[1],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[1],x)
                      }
                    )
             )

EtestH2 <- c(E_synd[2],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[2],x)
                      }
                    )
             )

EtestH3 <- c(E_synd[3],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[3],x)
                      }
                    )
             )

EtestH <- c(EtestH1,
            EtestH2,
            EtestH3
            )


## Expected values and CIs 
meanH <- sapply(EtestH,
                function(x){
                  sched_0$c.meanRl..meanR..meanRh.[sched_0$times==ceiling(x) & 
                                                     sched_0$Rtype=="R0=2.5"]
                  }
                )

lowH <- sapply(EtestH,
               function(x){
                 sched_0$c.lowerRl..lowerR..lowerRh.[sched_0$times==ceiling(x) & 
                                                       sched_0$Rtype=="R0=2.5"]
                 }
               )

uppH <- sapply(EtestH,
               function(x){
                 sched_0$c.upperRl..upperR..upperRh.[sched_0$times==ceiling(x) & 
                                                       sched_0$Rtype=="R0=2.5"]
                 }
               )


meanM <- sapply(EtestM,
                function(x){
                  sched_5$c.meanRl..meanR..meanRh.[sched_5$times==ceiling(x) & 
                                                     sched_5$Rtype=="R0=2.5"]
                  }
                )

lowM <- sapply(EtestM,
               function(x){
                 sched_5$c.lowerRl..lowerR..lowerRh.[sched_5$times==ceiling(x) & 
                                                       sched_5$Rtype=="R0=2.5"]
                 }
               )

uppM <- sapply(EtestM,
               function(x){
                 sched_5$c.upperRl..upperR..upperRh.[sched_5$times==ceiling(x) & 
                                                       sched_5$Rtype=="R0=2.5"]
                 }
               )


meanL <- sapply(EtestL,
                function(x){
                  sched_10$c.meanRl..meanR..meanRh.[sched_10$times==ceiling(x) & 
                                                      sched_10$Rtype=="R0=2.5" ]
                  }
                )

lowL <- sapply(EtestL,
               function(x){
                 sched_10$c.lowerRl..lowerR..lowerRh.[sched_10$times==ceiling(x) & 
                                                        sched_10$Rtype=="R0=2.5"]
                 }
               )

uppL <- sapply(EtestL,
               function(x){
                 sched_10$c.upperRl..upperR..upperRh.[sched_10$times==ceiling(x) & 
                                                        sched_10$Rtype=="R0=2.5"]
                 }
               )


resource <- c(rep("high",length(meanH)),
              rep("med",length(meanM)),
              rep("low",length(meanL))
              )

resource <- factor(resource, 
                   levels = c("high","med","low")
                   )

sched <- c("synd", "ind", "npool2", "npool5", "npool10")
sched <- factor(sched, 
                levels = c("synd",
                           "ind",
                           "npool2",
                           "npool5",
                           "npool10"
                           )
                )

sched <- rep(sched, 9)
incub <- rep(c(rep("medinc", 5), 
               rep("lowinc", 5), 
               rep("highinc", 5)
               ), 
             3)

incub <- factor(incub,
                levels = c("lowinc",
                           "medinc",
                           "highinc"
                           )
                )



# 4.  Plot expected cases  ------------------------------------------------



data_outsize <- data.frame(meanexp = c(meanH,meanM,meanL),
                           low = c(lowH,lowM,lowL),
                           upp = c(uppH,uppM,uppL),
                           resource,
                           sched,
                           incub )

label.legend <- c("synd" = "Syndromic",
                  "ind" = "No pooling",
                  "npool2" = expression(n[pool]*"=2"),
                  "npool5" = expression(n[pool]*"=5"),
                  "npool10" = expression(n[pool]*"=10")
                  )

label.res <- c("low" = "Low resources",
               "med" = "Medium resources",
               "high" = "High resources"
               )

label.inc <- c("lowinc" = "5.2 days",
               "medinc" = "6.7 days",
               "highinc" = "8.2 days"
               )


ggplot(
  data = data_outsize,
  aes(
    x = sched,
    y = meanexp,
    ymin = low,
    ymax = upp, 
    color = as.factor(sched)
    )
  ) + 
  geom_pointrange() +
  scale_color_viridis(
    discrete = TRUE,
    labels = label.legend
    ) + 
  theme_bw()+
  facet_grid(resource~incub, 
             labeller = labeller(resource=label.res,incub=label.inc)
             ) +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), 
        legend.position="bottom",
        axis.text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=11)
        ) + 
  labs(color="") + 
  xlab("") + 
  ylab("Size of outbreak")





# 5.  Sensitivity analysis ------------------------------------------------


## Sensitivity analysis 1 

### Assume Kucirka model of sens

RE.simul.K <- function(N,
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
      GI_j <- tran_curve_K((j-1),tau,rel,gamma.sh,gamma.sc,R0)
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


 ### Assume R0 values and testing frequency
R0vect <- c(1.8,2.5,3.6)
tau_ind <- 14

### Define number of cases across time for different R0 values when inc = 7.76, testing frequency = 7days

notestK <- sapply(R0vect,
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

indarrayK <- sapply(R0vect,
                    function(x){
                      RE.simul.K(disc,
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

npool2K <- sapply(R0vect,
                  function(x){
                    RE.simul.K(disc,
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

npool5K <- sapply(R0vect,
                  function(x){
                    RE.simul.K(disc,
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

npool10K <- sapply(R0vect,
                   function(x){
                     RE.simul.K(disc,
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

sched_0K <- inf_time(50000, notestK, 150)
sched_1K <- inf_time(50000, indarrayK,150)
sched_2K <- inf_time(50000,npool2K, 150)
sched_5K <- inf_time(50000, npool5K, 150)
sched_10K <- inf_time(50000, npool10K, 150)

type <- c(rep("synd",dim(sched_0K)[1]),
          rep("ind",dim(sched_1K)[1]),
          rep("npool2",dim(sched_2K)[1]),
          rep("npool5",dim(sched_5K)[1]),
          rep("npool10",dim(sched_10K)[1])
          )

type <- factor(type,
               levels = c("synd",
                          "ind",
                          "npool2",
                          "npool5",
                          "npool10"
                          )
               )

data_plotK <- rbind(sched_0K,
                    sched_1K,
                    sched_2K,
                    sched_5K,
                    sched_10K
                    )

data_plotK <- cbind(type,
                    data_plotK
                    )

data_sensK <- rbind(data_plot1,
                    data_plotK
                    )

sensfunc <- c(rep("Hellewell", dim(data_plot1)[1]),
              rep("Kucirka", dim(data_plotK)[1])
              )
data_sensK <- cbind(data_sensK, sensfunc)

names(data_sensK) <- c("type",
                       "times",
                       "mean",
                       "lower",
                       "upper",
                       "Rtype",
                       "sensfunc"
                       )

data_sensK$upper <- ifelse(data_sensK$upper<=150,
                           data_sensK$upper,
                           150
                           )

data_sensK$mean  <- ifelse(data_sensK$mean<=150,
                           data_sensK$mean,
                           150
                           )

ggplot(
  data_sensK,
  aes(x = times,
      fill=type
      )
  ) + 
  geom_line(
    data_sensK,
    mapping = aes(y = mean)
    ) +
  geom_ribbon(
    aes(x = times,
        ymin = lower,
        ymax = upper
        ), 
    alpha=0.3
    ) + 
  facet_grid(
    sensfunc~Rtype
    ) +
  xlab("Time since first infection (days)") + 
  ylab('Cumulative number of cases') + 
  scale_fill_brewer(
    palette = "Dark2",
    labels=label.legend
    ) + 
  theme_bw() + 
  labs(fill = "Screening") + 
  theme(legend.position = "bottom")



## Sensitivity analysis 2 

### Vary forward generation time - fig 2

tau_ind <- 7

# Cumulative over time
notestWu <- RE.simul.N.notest(disc,
                              RCHEpop,
                              no.inf,
                              2.5,
                              tau_ind,
                              shWu,
                              scWu,
                              BDrate,
                              timeframe
                              )

indWu <- RE.simul.N(disc,
                    RCHEpop,
                    no.inf,
                    2.5,
                    tau_ind,
                    1,
                    sh,
                    sc,
                    BDrate,
                    timeframe
                    )

npool2Wu <- RE.simul.N(disc,
                       RCHEpop,
                       no.inf,
                       2.5,
                       tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
                       sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                       shWu,
                       scWu,
                       BDrate,
                       timeframe
                       )

npool5Wu <- RE.simul.N(disc,
                       RCHEpop,
                       no.inf,
                       2.5,
                       tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],
                       sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                       shWu,
                       scWu,
                       BDrate,
                       timeframe
                       )

npool10Wu <- RE.simul.N(disc,
                        RCHEpop,
                        no.inf,
                        2.5,
                        tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],
                        sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                        shWu,
                        scWu,
                        BDrate,
                        timeframe
                        )

sched_0Wu <- inf_time(50000,notestWu,150) # Need to change poismean in inf_time
sched_1Wu <- inf_time(50000,indWu,150)
sched_2Wu <- inf_time(50000,npool2Wu,150)
sched_5Wu <- inf_time(50000,npool5Wu,150)
sched_10Wu <- inf_time(50000,npool10Wu,150)

notestLi <- RE.simul.N.notest(disc,
                              RCHEpop,
                              no.inf,
                              2.5,
                              tau_ind,
                              sh_Li,
                              sc_Li,
                              BDrate,
                              timeframe
                              )

indLi <- RE.simul.N(disc,
                    RCHEpop,
                    no.inf,
                    2.5,
                    tau_ind,
                    1,
                    sh_Li,
                    sc_Li,
                    BDrate,
                    timeframe
                    )

npool2Li <- RE.simul.N(disc,
                       RCHEpop,
                       no.inf,
                       2.5,
                       tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
                       sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                       sh_Li,
                       sc_Li,
                       BDrate,
                       timeframe
                       )

npool5Li <- RE.simul.N(disc,
                       RCHEpop,
                       no.inf,
                       2.5,
                       tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],
                       sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                       sh_Li,
                       sc_Li,
                       BDrate,
                       timeframe
                       )

npool10Li <- RE.simul.N(disc,
                        RCHEpop,
                        no.inf,
                        2.5,
                        tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],
                        sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                        sh_Li,
                        sc_Li,
                        BDrate,
                        timeframe
                        )

sched_0Li <- inf_time(50000, notestLi, 150)
sched_1Li <- inf_time(50000, indLi, 150)
sched_2Li <- inf_time(50000, npool2Li, 150)
sched_5Li <- inf_time(50000, npool5Li, 150)
sched_10Li <- inf_time(50000, npool10Li, 150)

EtestHsens <- EtestH[6:10]
  
meanHP <- sapply(EtestHsens,
                 function(x){
                   sched_0$c.meanRl..meanR..meanRh.[sched_0$times==ceiling(x) & 
                                                      sched_0$Rtype=="R0=2.5"]
                   }
                 )

lowHP <- sapply(EtestHsens,
                function(x){
                  sched_0$c.lowerRl..lowerR..lowerRh.[sched_0$times==ceiling(x) & 
                                                        sched_0$Rtype=="R0=2.5"]
                  }
                )

uppHP <- sapply(EtestHsens,
                function(x){
                  sched_0$c.upperRl..upperR..upperRh.[sched_0$times==ceiling(x) & 
                                                        sched_0$Rtype=="R0=2.5"]
                  }
                )

meanHwu <- sapply(EtestHsens,
                  function(x){
                    sched_0Wu$c.meanRl..meanR..meanRh.[sched_0Wu$times==ceiling(x) & 
                                                         sched_0$Rtype=="R0=2.5"]
                    }
                  )

lowHwu <- sapply(EtestHsens,
                 function(x){
                   sched_0Wu$c.lowerRl..lowerR..lowerRh.[sched_0Wu$times==ceiling(x) & 
                                                           sched_0$Rtype=="R0=2.5"]
                   }
                 )

uppHwu <- sapply(EtestHsens,
                 function(x){
                   sched_0Wu$c.upperRl..upperR..upperRh.[sched_0Wu$times==ceiling(x) & 
                                                           sched_0$Rtype=="R0=2.5"]
                   }
                 )

meanHLi <- sapply(EtestHsens,
                  function(x){
                    sched_0Li$c.meanRl..meanR..meanRh.[sched_0Li$times==ceiling(x) &
                                                         sched_0$Rtype=="R0=2.5"]
                    }
                  )

lowHLi <- sapply(EtestHsens,
                 function(x){
                   sched_0Li$c.lowerRl..lowerR..lowerRh.[sched_0Li$times==ceiling(x) & 
                                                           sched_0$Rtype=="R0=2.5"]
                   }
                 )

uppHLi <- sapply(EtestHsens,
                 function(x){
                   sched_0Li$c.upperRl..upperR..upperRh.[sched_0Li$times==ceiling(x) & 
                                                           sched_0$Rtype=="R0=2.5"]
                   }
                 )


sched <- c("synd", "ind", "npool2", "npool5", "npool10")
sched <- factor(sched,
                levels = c("synd",
                           "ind",
                           "npool2",
                           "npool5",
                           "npool10"
                           )
                )

sched <- rep(sched,3)
forwardgen <- c(rep("Park et al.", length(meanHP)),
                rep("Wu et al.",length(meanHwu)),
                rep("Li et al.",length(meanHLi)
                    )
                ) 

data_forward <- data.frame(meanexp = c(meanHP,meanHwu,meanHLi),
                           low = c(lowHP,lowHwu,lowHLi),
                           upp = c(uppHP,uppHwu,uppHLi),
                           sched, 
                           forwardgen)
label.legend <- c("synd" = "No screening",
                  "ind" = "Without pooling",
                  "npool2"=expression(n[pool]*"=2"),
                  "npool5"=expression(n[pool]*"=5"),
                  "npool10"=expression(n[pool]*"=10")
                  )

label.res <- c("low" = "Low resources",
               "med" = "Medium resources",
               "high" = "High resources"
               )

label.inc <- c("lowinc" = "5.2 days",
               "medinc" = "6.7 days",
               "highinc" = "8.2 days"
               )

facet_forgen <- 
  ggplot(
    data = data_forward,
    aes(x = sched,
        y = meanexp,
        ymin = low,
        ymax = upp,
        color = as.factor(sched)
        )
    ) + 
  geom_pointrange() +
  scale_color_viridis(discrete = TRUE,
                      labels = label.legend) + 
  theme_bw() +
  facet_wrap(
    ~forwardgen
    ) + 
  ylim(c(0,20)) +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), 
        legend.position = "bottom",
        axis.text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=11)
        ) + 
  labs(color="") + 
  xlab("") + 
  ylab("Size of outbreak")


forWu <- dgamma(a, shape = shWu, sc = scWu)
forPark <- dgamma(a, shape = sh_Park, sc = sc_Park)
forLi <- dgamma(a, shape = sh_Li, scale = sc_Li)

forgen <- c(forWu,forPark,forLi)
gen <- c(rep("Wu et al. (2020)",length(forWu)),
         rep("Park et al. (2020)",length(forPark)),
         rep("Li et al. (2020)",length(forLi))
         )

timesgen <- rep(a,3)

datagen <- data.frame(timesgen, forgen, gen)

forgen <- 
  ggplot(
    datagen, 
    aes(x = timesgen,
        y = forgen,
        color = gen)
    ) + 
  geom_line() + 
  scale_color_brewer(palette = "Dark2") + 
  xlab("Time since first infection (days)") +
  ylab("Density") + 
  labs(color="Forward \ngeneration \ndistribution") + 
  theme_minimal()


ggarrange(forgen,facet_forgen,ncol=1,labels=c("A","B"))


#### Fig 3

tau_ind <- 7
E_synd <- ceiling(det_symp(c(6.7, 5.2, 8.2),
                           0.85)
                  )

EtestHa <- c(E_indH, 
             E_high2, 
             E_high5, 
             E_high10
             )

EtestH1 <- c(E_synd[1],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[1],x)
                      }
                    )
             )

EtestH2 <- c(E_synd[2],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[2],x)
                      }
                    )
             )

EtestH3 <- c(E_synd[3],
             sapply(EtestHa,
                    function(x){
                      min(E_synd[3],x)
                      }
                    )
             )

EtestHa<- c(EtestH1, EtestH2, EtestH3)

meanHa <- sapply(EtestHa,
                 function(x){
                   sched_0$c.meanRl..meanR..meanRh.[sched_0$times==ceiling(x) & 
                                                      sched_0$Rtype=="R0=2.5"]
                   }
                 )

lowHa <- sapply(EtestHa,
                function(x){
                  sched_0$c.lowerRl..lowerR..lowerRh.[sched_0$times==ceiling(x) & 
                                                        sched_0$Rtype=="R0=2.5"]
                  }
                )

uppHa <- sapply(EtestHa,
                function(x){
                  sched_0$c.upperRl..upperR..upperRh.[sched_0$times==ceiling(x) & 
                                                        sched_0$Rtype=="R0=2.5"]
                  }
                )

sched <- c("synd", "ind", "npool2", "npool5", "npool10")
sched <- factor(sched,
                levels = c("synd","ind","npool2","npool5","npool10")
                )

sched <- rep(sched,3)

propsymp <- c(rep("65%",length(meanH)),
              rep("85%",length(meanHa))
              )
incub <- c(rep("5.2 days",5),
           rep("6.7 days",5),
           rep("8.2 days",5)
           )

rep(incub,2)

data_symp <- data.frame(mean = c(meanH,meanHa),
                        low = c(lowH,lowHa),
                        upp = c(uppH,uppHa),
                        sched,
                        propsymp,
                        incub
                        )

ggplot(
  data = data_symp,
  aes(x = sched,
      y = mean,
      ymin = low,
      ymax = upp,
      color = as.factor(sched))
  ) + 
  geom_pointrange() +
  scale_color_viridis(discrete = TRUE,
                      labels = label.legend) + 
  theme_bw() +
  facet_grid(
    propsymp~incub
    ) +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), 
        legend.position = "bottom",
        axis.text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=11)
        ) +
  labs(color = "") + 
  xlab("") + 
  ylab("Size of outbreak")



## Sensitivity analysis 3 

### Effect of vaccination on testing strategy - Fig 4


RE.simul.vax <- function(N,
                       pop_size, 
                       I.init,
                       R0,
                       tau,
                       rel,
                       gamma.sh,
                       gamma.sc,
                       mu, 
                       horizon,
                       cov,
                       eff){
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
      z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
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


RE.simul.vax.notest <- function(N,
                              pop_size, 
                              I.init,
                              R0,
                              tau,
                              gamma.sh,
                              gamma.sc,
                              mu, 
                              horizon,
                              cov,
                              eff){
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
      z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
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



#########

tau_ind2 <- 14
cov = 0.5
effvect = c(0.3,0.5,0.75)

#### Calculate cases over time for VE=30%, 50%, 75%

R0est = 2.5 # scenario 1

notest_vac <- sapply(effvect,
                     function(x){
                       RE.simul.vax.notest(disc,
                                           RCHEpop,
                                           no.inf,
                                           R0est,
                                           tau_ind,
                                           sh,
                                           sc,
                                           BDrate,
                                           timeframe,
                                           cov,
                                           x
                                           )
                       }
                     )

ind_vac <- sapply(effvect,
                  function(x){
                    RE.simul.vax(disc,
                                 RCHEpop,
                                 no.inf,
                                 R0est,
                                 tau_ind,
                                 1,
                                 sh,
                                 sc,
                                 BDrate,
                                 timeframe,
                                 cov,
                                 x
                                 )
                    }
                  )

npool2_vac <- sapply(effvect,
                     function(x){
                       RE.simul.vax(disc,
                                    RCHEpop,
                                    no.inf,
                                    R0est,
                                    tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe,
                                    cov,
                                    x
                                    )
                       }
                     )

npool5_vac <- sapply(effvect,
                     function(x){
                       RE.simul.vax(disc,
                                    RCHEpop,
                                    no.inf,
                                    R0est,
                                    tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                    sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                    sh,
                                    sc,
                                    BDrate,
                                    timeframe,
                                    cov,
                                    x
                                    )
                       }
                     )

npool10_vac <- sapply(effvect,
                      function(x){
                        RE.simul.vax(disc,
                                     RCHEpop,
                                     no.inf,
                                     R0est,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe,
                                     cov,
                                     x
                                     )
                        }
                      )




CRl_vac = RCHEpop - c(notest_vac[,1]$S,
                      ind_vac[,1]$S,
                      npool2_vac[,1]$S,
                      npool5_vac[,1]$S,
                      npool10_vac[,1]$S
                      )

CR_vac = RCHEpop - c(notest_vac[,2]$S,
                     ind_vac[,2]$S,
                     npool2_vac[,2]$S,
                     npool5_vac[,2]$S,
                     npool10_vac[,2]$S
                     )

CRh_vac = RCHEpop - c(notest_vac[,3]$S,
                      ind_vac[,3]$S,
                      npool2_vac[,3]$S,
                      npool5_vac[,3]$S,
                      npool10_vac[,3]$S
                      )



#####################

R0est = 1.8 # scenario 2

notest_vacl <- sapply(effvect,
                      function(x){
                        RE.simul.vax.notest(disc,
                                            RCHEpop,
                                            no.inf,
                                            R0est,
                                            tau_ind,
                                            sh,
                                            sc,
                                            BDrate,
                                            timeframe,
                                            cov,
                                            x
                                            )
                        }
                      )

ind_vacl <- sapply(effvect,
                   function(x){
                     RE.simul.vax(disc,
                                  RCHEpop,
                                  no.inf,
                                  R0est,
                                  tau_ind,
                                  1,
                                  sh,
                                  sc,
                                  BDrate,
                                  timeframe,
                                  cov,
                                  x
                                  )
                     }
                   )

npool2_vacl <- sapply(effvect,
                      function(x){
                        RE.simul.vax(disc,
                                     RCHEpop,
                                     no.inf,
                                     R0est,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe,
                                     cov,
                                     x
                                     )
                        }
                      )


npool5_vacl <- sapply(effvect,
                      function(x){
                        RE.simul.vax(disc,
                                     RCHEpop,
                                     no.inf,
                                     R0est,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe,
                                     cov,
                                     x
                                     )
                        }
                      )


npool10_vacl <- sapply(effvect,
                       function(x){
                         RE.simul.vax(disc,
                                      RCHEpop,
                                      no.inf,
                                      R0est,
                                      tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                      sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                      sh,
                                      sc,
                                      BDrate,
                                      timeframe,
                                      cov,
                                      x
                                      )
                         }
                       )

CRl_vacl = RCHEpop - c(notest_vacl[,1]$S,
                       ind_vacl[,1]$S,
                       npool2_vacl[,1]$S,
                       npool5_vacl[,1]$S,
                       npool10_vacl[,1]$S
                       )

CR_vacl = RCHEpop - c(notest_vacl[,2]$S,
                      ind_vacl[,2]$S,
                      npool2_vacl[,2]$S,
                      npool5_vacl[,2]$S,
                      npool10_vacl[,2]$S
                      )

CRh_vacl = RCHEpop - c(notest_vacl[,3]$S,
                       ind_vacl[,3]$S,
                       npool2_vacl[,3]$S,
                       npool5_vacl[,3]$S,
                       npool10_vacl[,3]$S
                       )


######################

R0est = 3 # scenario 3

notest_vach <- sapply(effvect,
                      function(x){
                        RE.simul.vax.notest(disc,
                                            RCHEpop,
                                            no.inf,
                                            R0est,
                                            tau_ind,
                                            sh,
                                            sc,
                                            BDrate,
                                            timeframe,
                                            cov,
                                            x
                                            )
                        }
                      )

ind_vach <- sapply(effvect,
                   function(x){
                     RE.simul.vax(disc,
                                  RCHEpop,
                                  no.inf,
                                  R0est,
                                  tau_ind,
                                  1,
                                  sh,
                                  sc,
                                  BDrate,
                                  timeframe,
                                  cov,
                                  x
                                  )
                     }
                   )

npool2_vach <- sapply(effvect,
                      function(x){
                        RE.simul.vax(disc,
                                     RCHEpop,
                                     no.inf,
                                     R0est,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe,
                                     cov,
                                     x
                                     )
                        }
                      )

npool5_vach <- sapply(effvect,
                      function(x){
                        RE.simul.vax(disc,
                                     RCHEpop,
                                     no.inf,
                                     R0est,
                                     tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],
                                     sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                                     sh,
                                     sc,
                                     BDrate,
                                     timeframe,
                                     cov,
                                     x
                                     )
                        }
                      )

npool10_vach <- sapply(effvect,
                       function(x){
                         RE.simul.vax(disc,
                                      RCHEpop,
                                      no.inf,
                                      R0est,
                                      tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                      sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                      sh,
                                      sc,
                                      BDrate,
                                      timeframe,
                                      cov,
                                      x
                                      )
                         }
                       )


CRl_vach = RCHEpop - c(notest_vach[,1]$S,
                       ind_vach[,1]$S,
                       npool2_vach[,1]$S,
                       npool5_vach[,1]$S,
                       npool10_vach[,1]$S
                       )

CR_vach = RCHEpop - c(notest_vach[,2]$S,
                      ind_vach[,2]$S,
                      npool2_vach[,2]$S,
                      npool5_vach[,2]$S,
                      npool10_vach[,2]$S
                      )

CRh_vach = RCHEpop - c(notest_vach[,3]$S,
                       ind_vach[,3]$S,
                       npool2_vach[,3]$S,
                       npool5_vach[,3]$S,
                       npool10_vach[,3]$S
                       )



## Expected time of detection

inc = 6.7

Enotest_vac <- (1/0.64)*inc
Eind_vac <- median(det_symp(tau_ind2,
                            sensind,
                            n,
                            inc)
                   )

Enpool2_vac <- median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==2],
                               sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
                               n,
                               x)
                      )

Enpool5_vac <- median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==5],
                               sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
                               n,
                               x)
                      )

Enpool10_vac <- median(det_symp(tau_ind2*sens_HCWp1$eff[sens_HCWp1$npool==10],
                                sensind*sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
                                n,
                                x)
                       )


Expvac.R <- c(Enotest_vac[1],
              Enotest_vac[1],
              Enpool2_vac[1],
              Enpool5_vac[1],
              Enpool10_vac[1]
              )



### Cases 

## R0 = 2.5
Casesvac.Rla <- sapply(Expvac.R,
                       function(x){
                         CRl_vac[ceiling(x)+1]
                         }
                       )

Casesvac.Ra <- sapply(Expvac.R,
                      function(x){
                        CR_vac[ceiling(x)+1]
                        }
                      )

Casesvac.Rha <- sapply(Expvac.R,
                       function(x){
                         CRh_vac[ceiling(x)+1]
                         }
                       )

# R0 = 1.8
Casesvac.Rlb <- sapply(Expvac.R,
                       function(x){
                         CRl_vacl[ceiling(x)+1]
                         }
                       )

Casesvac.Rb <- sapply(Expvac.R,
                      function(x){
                        CR_vacl[ceiling(x)+1]
                        }
                      )

Casesvac.Rhb <- sapply(Expvac.R,
                       function(x){
                         CRh_vacl[ceiling(x)+1]
                         }
                       )

# R0 = 3.6
Casesvac.Rlc <- sapply(Expvac.R,
                       function(x){
                         CRl_vach[ceiling(x)+1]
                         }
                       )

Casesvac.Rc <- sapply(Expvac.R,
                      function(x){
                        CR_vach[ceiling(x)+1]
                        }
                      )

Casesvac.Rhc <- sapply(Expvac.R,
                       function(x){
                         CRh_vach[ceiling(x)+1]
                         }
                       )


Cases.facet <- c(Casesvac.Rla, Casesvac.Ra, Casesvac.Rha)

CaseRlow <-  c(Casesvac.Rlb, Casesvac.Rb, Casesvac.Rhb)

CaseRhigh <-c(Casesvac.Rlc, Casesvac.Rc, Casesvac.Rhc)

Incub.f <- factor(c(rep("30%", 5),
                    rep("50%", 5),
                    rep("75%", 5)
                    ),
                  levels = c("30%",
                             "50%",
                             "75%")
                  )


Test.f<- factor(c("notest", "ind", "npool2", "npool5",  "npool10"),
                levels = c("notest", "ind", "npool2", "npool5",  "npool10")
                )

Test.facet <- c(rep(Test.f, 3))

data_facet <- data.frame(Test.facet, Cases.facet, Incub.facet, CaseRlow, CaseRhigh)

label.legend <- c("1" = "No screening",
                  "2" = "Without pooling",
                  "3"=expression(n[pool]*"=2"),
                  "4"=expression(n[pool]*"=4"),
                  "5"=expression(n[pool]*"=8"),
                  "6"=expression(n[pool]*"=14")
                  )

label.inc <- c("1" = "VE=30%",
               "2"="VE=50%",
               "3" = "VE=75%")

label.freq <- c("28" = "Low resources",
                "14" = "Medium resources",
                "7" = "High resources"
                )

ggplot(
  data = data_facet, 
  aes(x = Test.facet,
      y = Cases.facet,
      color = as.factor(Test.facet), 
      ymin =  CaseRlow, 
      ymax = CaseRhigh)
  ) +
  geom_pointrange() +
  scale_color_viridis(discrete = TRUE,
                      labels = label.legend) + 
  theme_bw() +
  facet_wrap(~Incub.facet,
             labeller = labeller(Incub.facet = label.inc)
             )+
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        axis.text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size = 11)
        ) + 
  labs(color="") +
  xlab("")+
  ylab("Size of outbreak")+
  ylim(c(0,35))

