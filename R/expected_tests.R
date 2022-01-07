
#################################################################
# Perform the decision analysis to determine the expected number 
# of tests required for different pool sizes 
#
#################################################################



# 1.  Load data -----------------------------------------------------------


# Pooled sensitivity values obtained using the method 1.	Muniesa, A., C, F., Fuertes, 
# H., Halaihel, N. & De Blas, I. Estimation of the relative sensitivity of qPCR analysis 
# using pooled samples. PLoS One e93491 (2014) doi:10.1371/journal.pone.0093491.

  sens_HCW <- read.csv(here("data", "SensPool_HCWscreen.csv")) # read data
  sens_HCW <- sens_HCW[sens_HCW$Prevalence<2,] # restrict to prevalence values of less than 2%
  
  sens_p001 <- sens_HCW[which(sens_HCW$Prevalence==1),] # prevalence = 1%

  

# 2. Determine expected number of tests ------------------------------------------------

  
# Assume specificity 100% for PCR
  sp <- 1
  se = sens_HCW$Sensitivity <- sens_HCW$Sensitivity/100
  sens_low = sens_HCW$SensCIlow <- sens_HCW$SensCIlow/100
  sens_upp = sens_HCW$SensCIupp <- sens_HCW$SensCIupp/100
  npool <- sens_HCW$npool 
  prev = sens_HCW$Prevalence <- sens_HCW$Prevalence/100 # Set prevalence
  N <- 50  # Number of staff for max CA home = 100 beds
    
  
# Run Bernoulli trials       
  set.seed(564)
  prev = unique(prev)
  
  bernest = function(prev,draws){
    
    U <- runif(draws,0,1)
    prob = ifelse(U<=prev,1,0)
    p_est <- sum(prob)/draws
    
    return(p_est)
    }

  p_est <- sapply(prev,function(x){bernest(x,draws=100000)})
  sens_HCW$p_est <- rep(p_est,20)
  
  
# Expected number of tests and 95% CIs
  
  sens_HCW <-  
    sens_HCW %>% 
      mutate(
        Pnp = (1-sp)*(1-p_est)^npool+se*(1-(1-p_est)^npool),
        Z = (N/npool)*((npool+1)*Pnp+(1-Pnp)),
        sig_z = sqrt(N*npool*(((1-p_est)^npool)-(1-p_est)^(2*npool))),
        z_low = Z-1.96*sig_z/sqrt(N*npool),
        z_upp = Z+1.96*sig_z/sqrt(N*npool)
        )
  
  mins <- 
    sens_HCW %>% 
    group_by(Prevalence) %>% 
    filter(Z == min(Z))
  
  
## Plot 
  
 ggtest <- 
   ggplot(
     sens_HCW,
     aes(x = npool,
         y = Z,
         group = as.factor(Prevalence),
         color = as.factor(Prevalence)
         )
     ) +
   geom_line() + 
   geom_ribbon(
     aes(
       ymin = z_low, 
       ymax = z_upp,
       color = as.factor(Prevalence)
       ),
     color=NA,
     alpha=0.1
     ) +
   ylab("Expected number of tests") + 
   geom_point( 
     data = mins,
     aes(x = npool, y = Z),
     col = "black",
     size = 2,
     shape = 4
     ) +
   theme_light() + 
   xlab(expression(n[pool])) + 
   ylim(c(0,N+5)) + 
   geom_hline(
     yintercept = N,
     linetype = "dashed"
     ) +
   scale_color_viridis(
     alpha = 1, 
     begin = 0, 
     end = 1, 
     direction = 1,
     discrete = TRUE, 
     option = "D"
     ) + 
   theme(
     legend.title = element_text(size = 10),
     legend.text = element_text(size = 10),
     axis.text=element_text(size=8),
     axis.title=element_text(size=10)
     ) + 
   labs(color="Prevalence")
  


# 2. Expected number of false negatives  -----------------------------------------

  # Calculate false negatives
 
  sens_HCW <- sens_HCW %>%
     mutate(
       pFN = p_est*(1-Sensitivity)
       ) 

 # Plot
 
  ggplot(
    sens_HCW,
    aes(
      x = npool,
      y = pFN,
      group = as.factor(Prevalence),
      color = as.factor(Prevalence)
      )
    ) +
    geom_line() +
    ylab("Probability of false negative") +
    theme_light() +
    xlab(expression(n[pool])) + 
    ylim(c(0,0.4)) +
    scale_color_viridis( 
      alpha = 1, 
      begin = 0, 
      end = 1, 
      direction = 1,
      discrete = TRUE, 
      option = "D"
      ) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.title=element_text(size=11)
      ) + 
    labs(color = "Prevalence")   



# 3. Efficiency ------------------------------------------------------------

# Calculate efficiency and 95% CIs
  
  sens_HCW <- 
    sens_HCW %>%
    mutate(
      eff = Z/N,
      sig_eff = (1/N)*sqrt(N*npool*(((1-p_est)^npool)-(1-p_est)^(2*npool))),
      eff_low = eff-1.96*sig_eff/sqrt(N*npool),
      eff_upp = eff+1.96*sig_eff/sqrt(N*npool),
      )
  
  mineff <- 
    sens_HCW %>% 
    group_by(Prevalence) %>% 
    filter(eff == min(eff))
  
# Plot
  
 ggeff <- 
   ggplot(
     sens_HCW,
     aes(
       x = npool,
       y = eff,
       group = as.factor(Prevalence),
       color = as.factor(Prevalence)
       )
     ) +
   geom_line() + 
   geom_point(
     data = mineff,
     aes(
       x = npool,
       y = eff
       ),
     col = "black",
     size = 2,
     shape = 4
     ) + 
   ylab("Efficiency") +
   geom_ribbon(
     aes(
       ymin = eff_low, 
       ymax = eff_upp,
       color = as.factor(Prevalence)
       ),
     color = NA,
     alpha = 0.1
     ) +
   theme_light() + 
   xlab(expression(n[pool])) + 
   geom_hline(
     yintercept = 1,
     linetype = "dashed"
     ) +
    scale_color_viridis( 
      alpha = 1, 
      begin = 0, 
      end = 1, 
      direction = 1,
      discrete = TRUE,
      option = "D"
      ) +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text = element_text(size=9),
      axis.title = element_text(size=10)
      ) + 
   labs(color = "Prevalence")

 
 
  # Plot Facet_grid for different pool sizes

  sens_HCW_sub <- sens_HCW[which(sens_HCW$npool %in% c(2,4,6,8,10,12)),]
  sens_HCW_sub <- sens_HCW_sub[which(sens_HCW_sub$Prevalence<0.05),]
  
  multformat <- function(){
    function(x) format(100*x,digits=2)
  }
  
  ggplot(
    sens_HCW_sub, 
    aes(
      x = as.factor(Prevalence),
      y = Z,
      ymin = z_low,
      ymax = z_upp,
      color = as.factor(Prevalence)
      )
    ) +
    geom_pointrange() + 
    facet_wrap(~npool) + 
    scale_color_viridis(
      discrete = TRUE,
      labels=c("1","5","10")
      ) +
    theme_light() +
    theme(
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title=element_text(size=11)
      ) + 
    labs(color="Prevalence") + 
    xlab("") + 
    ylab("Expected number of tests") +
    geom_hline(
      yintercept = N,
      linetype="dotted"
      )
  
  
  
  # Plot facet_grid for expected number of false negatives
  
  sens_HCW_sub <- sens_HCW[which(sens_HCW$npool %in% c(5,10,15,20)),]
  sens_HCW_sub <- sens_HCW_sub[which(sens_HCW_sub$Prevalence<0.2),]
  
  multformat <- function(){
    function(x) format(100*x,digits=2)
  }
  
  
  sens_HCW_sub$npool <- factor(sens_HCW_sub$npool,labels=c("n[pool]==5","n[pool]==10","n[pool]==15","n[pool]==20"))
  labels.sens <- c("0.01"="1","0.005"="0.5","2e-04"="0.02")

   ggpFN <- 
     ggplot(
       sens_HCW_sub, 
       aes(
         x = as.factor(Prevalence),
         y = pFN,
         color = as.factor(Prevalence)
         )
       ) +
    geom_point() + 
    facet_wrap(
      ~npool,
      ncol = 2,
      labeller = label_parsed
      ) + 
     scale_color_viridis(
       discrete = TRUE,
       labels=label.sens
       ) + 
     theme_bw() + 
     scale_y_continuous(
       labels = label_percent(suffix="",accuracy=0.1)
       ) +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      strip.text = element_text(face="bold",size=7),
      axis.text = element_text(size=9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title=element_text(size=9)
      ) + 
     labs(color="Prevalence (%)") + 
     xlab("") + 
     ylab("False negative probability (%)") + 
     removeGridX()
  
 

# 4. Cost savings (%) -----------------------------------------------------


 ## Assume cost of test 500HKD
 
   sens_HCW <- 
     sens_HCW %>%
     mutate(
       costsave = (1-sens_HCW$eff)*100
       )
   
   
   ggcost <- 
     ggplot(
       sens_HCW,
       aes(
         x = npool,
         y = costsave,
         group = as.factor(Prevalence),
         color = as.factor(Prevalence)
         )
       ) +
     geom_line() + 
     ylab("Cost savings (%)") +
     theme_bw() + 
     xlab(expression(n[pool])) + 
     geom_hline( 
       yintercept=100,
       linetype="dashed"
       ) +
     scale_color_viridis(
       alpha = 1, 
       begin = 0, 
       end = 1, 
       direction = 1,
       discrete = TRUE, 
       option = "D"
       ) +
     theme(
       legend.title = element_text(size = 10),
       legend.text = element_text(size = 10),
       axis.text=element_text(size=9),
       axis.title=element_text(size=10)
       ) + 
     labs(color="Prevalence")
   
 
 ggarrange(
   ggsens,
   ggeff,
   ggcost,
   ggpFN,
   heights = c(1,1),
   widths = c(1,1.05),
   common.legend = TRUE,
   legend = "bottom",
   labels = c("A","B","C","D")
   )

   
  