######################################################################
# Updated analysis comparing screening in Delta and Omicron outbreaks 
# 
# Date: 2022-01-11
######################################################################


# Functions ---------------------------------------------------------------


#' Renewal Equation
#' with seasonal forcing and
#' vital dynamics 
#' @param N subpartition of time unit. N=1 <=> SIR discrete, N=large <=> SIR continuous.
#' @param R0 Basic reproduction number
#' @param pop_size Population size (=1)
#' @param I.init Initial infectious individuals
#' @param rel relative sensitivity 
#' @param mu Immgration and death rate (balanced)
#' @param prd Period of the seasonal forcing
#' @param amp Amplitude of the seasonal forcing
#' @param horizon Time horizon of the simulation
#' @param cov vaccine coverage
#' @param eff vaccine effectiveness against infection
#' @param testing testing strategy: "syndromic", "molecular", "antigen"
#' 


    RE.simul <- function(N = 1,
                         pop_size = 150, # number of staff and residents 
                         I.init = 1, # one infection initially
                         R0 = 3, # wild type 
                         tau = 7, # weekly testing
                         rel = 1, # unpooled
                         gamma.sh = 1.960923,
                         gamma.sc = 4.334694, 
                         mu = 0, # zero seasonal forcing
                         horizon = 30, # 30 days
                         cov = 0.5, # 50% vax coverage
                         eff = 0.5, # 50% effectiveness vs infection
                         testing = "syndromic",
                         sens_func = "Hellewell",
                         variant = "Delta"){ #assumed pcr sensitivity func
  inc <- vector()
  S <- vector()
  C <- vector()
  inc[1] <- I.init
  S[1] <- pop_size - I.init
  C[1] <- 0
  
  if(variant == "Delta" | variant == "delta"){
    # set generation time distribution
    
  } else 
    if(variant == "Omicron" | variant == "omicron"){
      
    } else { stop("variant not recognised") }
  
  
  for(t in 2:(1 + horizon * N))
  {
    z <- 0
    
    #syndromic
    if(testing == "syndromic"){
   
       for(j in 1:(t-1))
         {
      GI_j <- tran_curve_notest((j-1),tau,gamma.sh,gamma.sc,R0)
      z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
      } 
   
    } else if(testing == "molecular"){
         
      for(j in 1:(t-1))
      {

        if(sens_func == "Hellewell"){
          
        GI_j <- tran_curve_H((j-1),tau,rel,gamma.sh,gamma.sc,R0)
        z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
        
        } else
          if(sens_func == "Kucirka"){
            
            GI_j <- tran_curve_K((j-1),tau,rel,gamma.sh,gamma.sc,R0)
            z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
            
          } else{ stop("sensitivity function not recognised") }
      } 
      
    } else if (testing == "antigen"){ # rapid antigen
        
      for(j in 1:(t-1))
      {
        GI_j <- tran_curve_rapid((j-1), tau, rel, gamma.sh, gamma.sc, R0)
        z <- z + GI_j * inc[t-j] * exp(-mu/N * j) * (cov*(1-eff)+(1-cov))
      } 
       
       } else{ stop("Testing strategy not recognised") }
    
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

    

    
#' Number of cases over time
#' drawing from poisson distribution
#' @param n number of draws
#' @param result mean values over time for different testing strategies
#' @param RCHEpop number of staff and residents in facility
  
    
inf_time <- function(n = 50000,
                     result = NULL, # add results from RE.simul
                     RCHEpop = 150){
      
      mid <- vector()
      lower <- vector()
      upper <- vector()
      
      # use RE.simul results as mean counts
      pois_mean <- RCHEpop - result$S

      
      samp <- purrr::map(.x = pois_mean,
                         .f = rpois,
                          n = 100)
      
      samp_ord <- purrr::map(.x = samp,
                             .f = sort)
      
      
      mid <- purrr::map(.x = samp_ord,
                        .f = median)
      
      lower <- purrr::map(.x = samp_ord,
                          .f = quantile,
                          probs = 0.025)
      
      upper <- purrr::map(.x = samp_ord,
                          .f = quantile,
                          probs = 0.975)
      
      
      
      return(
        tibble(times = 1:length(pois_mean), 
                   mid = unlist(mid), 
                   low = unlist(lower), 
                   upp = unlist(upper))
        )
      
    }
    
    

#' Expected time of detection for each strategy
#' @param inc assumed incubation period
#' @param symp proportion symptomatic


det_symp <- function(inc,symp){
  
  rweibull(1000, 1.91, 0.11)
  
  time <- (1/symp)*inc # why did this have +1 - check! 
  
  return(c(time))
}

    


# Simulation parameters ---------------------------------------------------



## Assume 50 staff in RCHE with 100 beds 
disc <- 1
RCHEpop <- 150
no.inf <- 1
timeframe <- 30
BDrate <- 0


# set prevalence equal to 0.02%

sens_HCWp1 <- sens_HCW[sens_HCW$Prevalence==2e-04,] 






# Compute cumulative cases ------------------------------------------------

#### Parameters

tau_ind = 7

sens_rapid = 1 #0.70
sens_ind = 0.95
tau_rapid = 1

# set according to typical mitigation within healthcare settings
R0_del = 1.9
R0_omi = 1.9
  
# generation gamma params 
## shape
  sh_del = 2.230211 #2.201873
  sh_omi = 2.18759 #0.8889796
  
# scale 
  sc_del = 2.076485 #2.08913
  sc_omi = 0.9681889 #3.712121

test_comp <- c("syndromic", rep("molecular", 4), "antigen")
tau_pool = c(rep(tau_ind, 2), 
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],
             tau_rapid
        )
sens_pool <- c(rep(sens_ind, 2),
               sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
               sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
               sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
               sens_rapid # check parameter isn't used for rapid
               )
               
               
## Run scenarios

# assume same sensitivity over time of tests for delta and omicron? - 
# if generation time is changed this should be changed in line? make prop
# to infectiousness?

### Scenario 1: testing every 7 days

tau_ind = 7

delta_s1 <- purrr::pmap(list(tau = tau_pool,
                             R0 = R0_del, # assumed values for delta
                             gamma.sh = sh_del,  # assumed values for delta
                             gamma.sc = sc_del,  # assumed values for delta
                             rel = sens_pool,
                             testing = test_comp),  
                         RE.simul)

omicron_s1 <- purrr::pmap(list(tau = tau_pool,
                               R0 = R0_omi,  # assumed values for omicron
                               gamma.sh = sh_omi,  # assumed values for omicron
                               gamma.sc = sc_omi,  # assumed values for omicron
                               rel = sens_pool,
                               testing = test_comp,
                               variant = "omicron"),  
                          RE.simul)


### Scenario 2: testing every 14 days

tau_ind = 14
tau_pool = c(rep(tau_ind, 2), 
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2],
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5],
             tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10],
             tau_rapid
)

delta_s2 <- purrr::pmap(list(tau = tau_pool,
                             R0 = R0_del, # assumed values for delta
                             gamma.sh = sh_del,  # assumed values for delta
                             gamma.sc = sc_del,  # assumed values for delta
                             rel = sens_pool,
                             testing = test_comp),  
                        RE.simul)

omicron_s2 <- purrr::pmap(list(tau = tau_pool,
                               R0 = R0_omi,  # assumed values for omicron
                               gamma.sh = sh_omi,  # assumed values for omicron
                               gamma.sc = sc_omi,  # assumed values for omicron
                               rel = sens_pool,
                               testing = test_comp,
                               variant = "omicron"),  
                          RE.simul)



## Stochastic cases

del_s1_res <-
  purrr::map(.x = delta_s1,
             .f = inf_time,
             n = 5000,
             RCHEpop = 150)

del_s2_res <-
  purrr::map(.x = delta_s2,
             .f = inf_time,
             n = 5000,
             RCHEpop = 150)

omi_s1_res <-
  purrr::map(.x = omicron_s1,
             .f = inf_time,
             n = 5000,
             RCHEpop = 150)

omi_s2_res <-
  purrr::map(.x = omicron_s2,
             .f = inf_time,
             n = 5000,
             RCHEpop = 150)





## Plot cumulative cases

# Bind testing strategy data within scenario

res_comb <- 
  c(del_s1_res, del_s2_res, omi_s1_res, omi_s2_res)

results <- 
  res_comb %>%
  bind_rows() %>%
  mutate(
    test = rep(c(rep("notest", timeframe + 1),
                 rep("synd", timeframe + 1),
                 rep("npool2", timeframe + 1),
                 rep("npool5", timeframe + 1),
                 rep("npool10", timeframe + 1),
                 rep("rapid", timeframe + 1)
                  ), 
               4
               ),
    variant = c(rep("delta", (timeframe + 1)*6*2), # number of test strategies = 6
                rep("omicron", (timeframe + 1)*6*2)
                ),
    scenarios = rep(c(rep("tau = 7", (timeframe + 1)*6),
                      rep("tau = 14", (timeframe + 1)*6)
                      ),
                    2
                    )
  )



## Make facet data


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


## Plot

ggplot(
  results, 
  aes(
    x = times,
    fill = test
  )
) +
  geom_line(
    mapping = aes(
      y = mid,
      x = times
    ),
    size=0.5
  ) + 
  facet_grid(
    scenarios~variant
  ) +
  geom_ribbon(
    mapping = aes(
      x = times,
      ymin = low,
      ymax = upp
    ), 
    alpha=0.3
  ) +
  xlab("Time since first infection (days)") + 
  ylab('Cumulative number of cases') + 
  scale_fill_brewer(
    palette = "Dark2",
    #labels=label.legend
  ) + 
  theme_bw() + 
  labs(fill = "Screening") + 
  theme(legend.position = "bottom")





# Expected size of outbreak -----------------------------------------------

## Assume same proportion of symptomatic for omicron and delta - CHECK

inc_del = 5 # check and change
inc_omi = 4 # check and change
  
symp_del = 0.65 # check and change
symp_omi = 0.65 # check and change



# Expected value under syndromic
E_synd <- purrr::pmap(list(inc = c(inc_del, inc_omi),
                               symp = c(symp_del, symp_omi)
                               ),
                     det_symp
                     ) # add ceiling to code


tau_E <- c(tau_ind,
           ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==2]),
           ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==5]),
           ceiling(tau_ind*sens_HCWp1$eff[sens_HCWp1$npool==10]),
           tau_rapid)

sens_E <- c(1,
            se_time_H(a)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==2],
            se_time_H(a)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==5],
            se_time_H(a)*sens_HCWp1$Sensitivity[sens_HCWp1$npool==10],
            se_time_ant(a))


 
surv_del <- 
  purrr::pmap(list(tau = tau_E,
                   rel = sens_E),
              Tcut_Hlow, # edit this to include rapid tests also
              a = 1)


Exp_val <- # need to sort out different lengths
  
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










