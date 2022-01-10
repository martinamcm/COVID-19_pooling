####################################################
# Estimating Rt under different testing strategies #
#                                                  #
####################################################



# 1.  Forward generation time  --------------------------------------------


## Assign shape/scale to gamma dist for forward generation time f(a)
  
a <- seq(0,30,1) # Time since index infection (days)
  
  # Li et al.
  sh_Li = 4.8685
  sc_Li = 1.821917
  
  # Park et al.
  sh_Park = 1.960923
  sc_Park = 4.334694
    


# Sensitivity functions over time -----------------------------------------

  
  # Sensitivity functions that depend on age of infection 'a'
 
  # Kucirka et al. (2020)
 
  se_time_K <- function(a){
   ifelse(a<=21,
          inv.logit(-29.966+37.713*log(a)-14.452*(log(a))^2+1.721*(log(a))^3),
          inv.logit(6.878-2.436*log(a))
          )
 }
 
 # Hellewell et al. (2020) 

   se_time_H <- function(a){
   ifelse(
     a<=3.272221,
     inv.logit(1.592913+2.168557*(a-3.272221)),
     inv.logit(1.592913-0.2273548*(a-3.272221))
     )
 }



# Survivor functions P(T>tau) ---------------------------------------------


 Tcut_Hlow <- function(a,tau,rel){
   
   if(a<tau){
     se_tau <- rel*se_time_H(a)/tau
     PrTa <- se_tau
   }
   
   else{
     se_tau <- rel*se_time_H(a)/tau
     k <- floor(a-1/tau)
     inputi <- sapply(k,function(x){c(1:x)})
     input <- mapply(function(x,y){
       y-x*tau},x=inputi,y=a  )
     se_a <- sapply(input,function(x){(1-rel*se_time_H(x))})
     se_Ta <- mapply(function(x){prod(x)},x=list(se_a))
     PrTa <- unlist(se_Ta)*se_tau
     
   }
   
   return(PrTa)
 }
 
   
 Tcut_K <- function(a,tau,rel){
   
   if(a<tau){
     se_tau <- rel*se_time_K(a)/tau
     PrTa <- se_tau
   }
   
   else{
     se_tau <- rel*se_time_K(a)/tau
     k <- floor(a-1/tau)
     inputi <- sapply(k,function(x){c(1:x)})
     input <- mapply(function(x,y){
       y-x*tau},x=inputi,y=a  )
     se_a <- sapply(input,function(x){(1-rel*se_time_H(x))})
     se_Ta <- mapply(function(x){prod(x)},x=list(se_a))
     PrTa <- unlist(se_Ta)*se_tau
     
   }
   
   return(PrTa)
 }
 
 
 
 intTa <- function(a,tau,rel){
   
   integ <- 1-unlist(sapply(a,function(x){
            adaptIntegrate(Tcut_Hlow,
                           lower = 0,
                           upper = x,
                           tau=tau,
                           rel=rel
                           )
     }
     )[1,])
   
   return(integ)
 }
 
 
 intTK <- function(a,tau,rel){
   
   integ <- 1-unlist(sapply(a,function(x){
     adaptIntegrate(Tcut_K,
                    lower = 0,
                    upper = x,
                    tau=tau,
                    rel=rel
                    )
     }
     )[1,])
   
   return(integ)
 }
  
 
 


# Transmission curve ------------------------------------------------------

 
tran_curve_K <- function(a,tau,rel,shape,scale,R0){

  lambda = R0*dgamma(a,shape=shape,scale=scale)
  ProbTa = intTK(a,tau,rel)
  
  tran = lambda*ProbTa
  
  return(tran)
    
}


tran_curve_H <- function(a, tau, rel, shape, scale, R0){
  
  lambda = R0*dgamma(a,shape=shape,scale=scale)
  ProbTa = intTa(a,tau,rel)
  
  tran = lambda*ProbTa
  
  return(tran)
  
}


tran_curve_notest <- function(a,tau,shape,scale,R0){
  
  lambda = R0*dgamma(a,shape=shape,scale=scale)
  
  tran = lambda
  
  return(tran)
}


