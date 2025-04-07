
library(dplyr)
library(tidyverse)

######## 
Pre<-read.csv("Heerkou_Pre.csv")

rainsummary<-Pre %>% 
  filter(Month < 10 & Month>3) %>% ###Growth Season April-September
  mutate(Season = "Late", 
         Season = ifelse(Month == 4 | Month == 5 | Month == 6, "Early", Season)) %>% 
  group_by(Year, Season) %>%
  summarize(ppt = sum(Pre)) %>%
  spread(Season, ppt) 

rainsummary<-rainsummary%>%
  mutate(raintype = case_when(
    Early > mean(rainsummary$Early) & Late > mean(rainsummary$Late) ~ 1,##Consistent Wet
    Early > mean(rainsummary$Early) & Late < mean(rainsummary$Late)~ 2,##Post-flood Dry
    Early < mean(rainsummary$Early) & Late > mean(rainsummary$Late)~ 3,##Pre-flood Dry
    Early < mean(rainsummary$Early) & Late < mean(rainsummary$Late)~ 4##Consistent Dry
  ))

rain.pro<-rainsummary %>% 
  group_by(raintype) %>% 
  summarize(count = n()) %>%
  mutate(pro = count / sum(count)) 

post_1000<-readRDS("Output/post_1000.rds")

to_df <- function(mat) {
  as.data.frame(mat)
}

##change the data form
post_list <- lapply(post_1000, function(sublist) {
  lapply(sublist, to_df)
})


sg<-read.csv("s_g_data.csv")

##alpha interaction coefficient
AA.ls<-lapply(1:1000, function(i) ##each sample
  lapply(1:8, function(sp) ##each focal species
    lapply(post_list[[sp]][2:9], function(x) x[i,]))) #each alpha

##lambda intrinsic growth rate
Lambda.ls<-lapply(1:1000, function(i) ##each sample
  lapply(1:8, function(sp) ##each focal species
    lapply(post_list[[sp]][1], function(x) x[i,]))) ##each lambda

s.ls <- lapply(1:4, function(i) sg$s[sg$raintype==i]); s.ls ##seed soil bank survival

g.ls <- lapply(1:4, function(i) sg$g[sg$raintype==i]); g.ls ##germination fraction

pairs <- combn(1:8, 2); pairs


# fucntion
# ------------------------------------------------------------------------------------
# Functions for use in coexistence calcualtions

# Functions for use in invasino growth rate calcualtions, ajusted from Hallet 2019 Ecol Lett.
# ------------------------------------------------------------------------------------

# Determine equilibrium conditions for each species in isolation 
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  N <- s*(1-g)*N0 + N0*g*lambda*exp(-1*(a_intra*g*N0))
  return(N)
}

# invader population growth rate one time step forward 
pop.invade <- function (N0, resident, s, g, g_r, a_inter, lambda) {
  N <- s*(1-g)*N0 + N0*g*lambda*exp(-1*(a_inter*g_r*resident))
  return(N)
}

# resident population growth rate one time step forward 
pop.resident <- function (N0, resident, s, g, g_i, a_intra, a_inter, lambda) {
  N <- s*(1-g)*resident + resident*g*lambda*exp(-1*(a_intra*g*resident+a_inter*g_i*N0))
  return(N)
}

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

##remove infinite value
mean_f <- function(x) {mean(x[is.finite(x)])}

###1000 paramters samples from posterior distribution

time <- length(rainsummary$raintype) ##120year

partitioning_1000<-list()


for (rep in 1:1000) {
  
  A.est<-AA.ls[[rep]]
  lambda.est<-Lambda.ls[[rep]]
  
  partitioning_1000[[rep]]<-list()
  
  partitioning_i<-data.frame()
  partitioning_j<-data.frame()
  
  for (p in 1:ncol(pairs)) {
    
    ####data prep
    # ----------------------------------------------------------------------------------------
      i <- pairs[1,p]; j <- pairs[2,p]
      
      i.dat<-data.frame(raintype=1:4,
                        species=i,
                        a_ii=unlist(A.est[[i]][[i]]),
                        a_ij=unlist(A.est[[i]][[j]]),
                        lambda=unlist(lambda.est[[i]]),
                        s=unlist(lapply(1:4,function(x) s.ls[[x]][i])),
                        g=unlist(lapply(1:4,function(x) g.ls[[x]][i]))) %>%
        mutate(eta=lambda*g/(1-(1-g)*s)) %>% 
        mutate(Nstar_i =log(eta)/(a_ii*g)) ##Equillibrum abundance for i
      
      
      
      j.dat<-data.frame(raintype=1:4,
                        species=j,
                        a_ji=unlist(A.est[[j]][[i]]),
                        a_jj=unlist(A.est[[j]][[j]]),
                        lambda=unlist(lambda.est[[j]]),
                        s=unlist(lapply(1:4,function(x) s.ls[[x]][j])),
                        g=unlist(lapply(1:4,function(x) g.ls[[x]][j]))) %>%
        mutate(eta=lambda*g/(1-(1-g)*s)) %>% 
        mutate(Nstar_j =log(eta)/(a_jj*g))
      
      
      a_ii_weighted <- weighted.mean(i.dat$a_ii,rain.pro$pro)
      a_ij_weighted <- weighted.mean(i.dat$a_ij,rain.pro$pro)
      i_lambda_weighted <- weighted.mean(i.dat$lambda,rain.pro$pro)
      i_g_weighted <- weighted.mean(i.dat$g,rain.pro$pro)
      Nstar_i_weighted<-weighted.mean(i.dat$Nstar_i,rain.pro$pro)
      
      a_jj_weighted <- weighted.mean(j.dat$a_jj,rain.pro$pro)
      a_ji_weighted <- weighted.mean(j.dat$a_ji,rain.pro$pro)
      j_lambda_weighted <- weighted.mean(j.dat$lambda,rain.pro$pro)
      j_g_weighted <- weighted.mean(j.dat$g,rain.pro$pro)
      Nstar_j_weighted<-weighted.mean(j.dat$Nstar_j,rain.pro$pro)
      
      # use the timeseries of environmental conditions for environmental variability
      
      ####Equilibrium abundance
      # for i
      
      N0_i <- Nstar_i_weighted
      
      N_i <- rep(NA, time)
      N_i[1] <- N0_i
      
      for (t in 1:time) {
        par.i <- subset(i.dat, raintype==rainsummary$raintype[t]) ##
        N_i[t+1] <- pop.equilibrium(N0=N_i[t], s=par.i$s, g=par.i$g, 
                                    a_intra=par.i$a_ii, lambda=par.i$lambda)
      }
      
      
      
      
      # for j
      N0_j <- Nstar_j_weighted
      N_j <- rep(NA, time)
      N_j[1] <- N0_j
      
      for (t in 1:time) {
        par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
        N_j[t+1] <- pop.equilibrium(N0=N_j[t], s=par.j$s, g=par.j$g, 
                                    a_intra=par.j$a_jj, lambda=par.j$lambda)
      }
      
      ####Mutual Invasion
      
      # invade i first
      i_invader <- rep (NA, 60) ##
      j_resident <- rep (NA, 60)
      temp <- 1
      ##
      for (t in 60:time) {
        par.i <- subset(i.dat, raintype==rainsummary$raintype[t]) ##
        par.j <- subset(j.dat, raintype==rainsummary$raintype[t]) ##
        i_invader[temp] <- pop.invade(N0=1, resident=ceiling(N_j[t]), s=par.i$s, g=par.i$g, g_r=par.j$g,
                                      a_inter=par.i$a_ij, lambda=par.i$lambda) ##
        
        j_resident[temp] <- pop.resident(N0=1, resident=ceiling(N_j[t]), s=par.j$s, g=par.j$g, g_i=par.i$g,
                              a_intra=par.j$a_jj, a_inter=par.j$a_ji, 
                              lambda=par.j$lambda)/ceiling(N_j[t])
        
        temp  <- temp + 1 
      }
      
      # then have j invade into i
      j_invader <- rep (NA, 60)
      i_resident <- rep (NA, 60)
      temp <- 1
      for (t in 60:time) {
        par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
        par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
        
        j_invader[temp] <- pop.invade(N0=1, resident=ceiling(N_i[t]), s=par.j$s, g=par.j$g, g_r=par.i$g,
                                      a_inter=par.j$a_ji, lambda=par.j$lambda)
        i_resident[temp] <- pop.resident(N0=1, resident=ceiling(N_i[t]), s=par.i$s, g=par.i$g, g_i=par.j$g,
                                         a_intra=par.i$a_ii, a_inter=par.i$a_ij, 
                                         lambda=par.i$lambda)/ceiling(N_i[t])
        temp  <- temp + 1 
      }
      
      ###log transformed##
      i_inv <- log(i_invader) 
      j_inv <- log(j_invader)
      
      i_res <- log(i_resident)
      j_res <- log(j_resident)
    
    # ----------------------------------------------------------------------------------------
    #### Calculating the invasion rate with no variation in both alphas and lambda/g
    # ----------------------------------------------------------------------------------------
    #### Resident equilibriums 
    
    # for i
    N_i_no_var <- rep(NA, time)
    N_i_no_var[1] <- N0_i
    
    for (t in 1:time) {
      N_i_no_var[t+1] <- pop.equilibrium(N0=N_i_no_var[t], s=par.i$s, g=i_g_weighted, 
                                       a_intra=a_ii_weighted, lambda=i_lambda_weighted)
    }
    
    # for j
    N_j_no_var <- rep(NA, time)
    N_j_no_var[1] <- N0_j
    
    for (t in 1:time) {
      N_j_no_var[t+1] <- pop.equilibrium(N0=N_j_no_var[t], s=par.j$s, g=j_g_weighted, 
                                       a_intra=a_jj_weighted, lambda=j_lambda_weighted)
    }
    
    ####Mutual invasion
    
    #  i invade j
    i_invade_no_var <- rep (NA, 60)
    j_resident_no_var <- rep (NA, 60)
    
    temp <- 1
    
    for (t in 60:time) {
      
    i_invade_no_var[temp] <- pop.invade(N0=1, resident=ceiling(N_j_no_var[t]), s=par.i$s, 
                                  g=i_g_weighted, g_r=j_g_weighted, 
                                  a_inter=a_ij_weighted, lambda=i_lambda_weighted)
    
    j_resident_no_var[temp] <- pop.resident(N0=1, resident=ceiling(N_j_no_var[t]), s=par.j$s, 
                                           g=j_g_weighted, g_i=i_g_weighted,
                                           a_intra = a_jj_weighted,
                                           a_inter=a_ji_weighted, lambda=j_lambda_weighted)/ceiling(N_j_no_var[time])
    temp <- temp + 1
    
    }
    

    # j invade i
    j_invade_no_var <- rep (NA, 60)
    i_resident_no_var <- rep (NA, 60)
    
    temp <- 1
    
    for (t in 60:time) {
    j_invade_no_var[temp] <- pop.invade(N0=1, resident=ceiling(N_i_no_var[t]), s=par.j$s, 
                                  g=j_g_weighted, g_r=i_g_weighted, 
                                  a_inter=a_ji_weighted, lambda=j_lambda_weighted)
    
    i_resident_no_var[temp] <- pop.resident(N0=1, resident=ceiling(N_i_no_var[t]), s=par.i$s, 
                                      g=i_g_weighted, g_i=j_g_weighted,
                                      a_intra = a_ii_weighted,
                                      a_inter=a_ij_weighted, lambda=i_lambda_weighted)/ceiling(N_i_no_var[time])
    temp <- temp + 1
    }
    ##log transformed
    
    i_inv_epsilon_0 <- log(i_invade_no_var)
    j_inv_epsilon_0 <- log(j_invade_no_var)
    
    i_res_epsilon_0 <- log(i_resident_no_var)
    j_res_epsilon_0 <- log(j_resident_no_var)
    
    # ----------------------------------------------------------------------------------------
    #### Calculate the invasion rate with variable intrinsic growth rates (lambda) and g , 
    # but with NO variation in alphas
    # ----------------------------------------------------------------------------------------
    
    # find resident equilibrium with variable growth rates
    # for i
    N_i_var_lambda <- rep(NA, time)
    N_i_var_lambda[1] <- N0_i
    
    for (t in 1:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      N_i_var_lambda[t+1] <- pop.equilibrium(N0=N_i_var_lambda[t], s=par.i$s, g=par.i$g, 
                                             a_intra=a_ii_weighted, lambda=par.i$lambda)
    }
    
    
    # for j
    N_j_var_lambda <- rep(NA, time)
    N_j_var_lambda[1] <- N0_j
    
    for (t in 1:time) {
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      N_j_var_lambda[t+1] <- pop.equilibrium(N0=N_j_var_lambda[t], s=par.j$s, g=par.j$g, 
                                             a_intra=a_jj_weighted, lambda=par.j$lambda)
    }
    
    # ----------------------------------------------------------------------------------------
    #### Calculate the invasion rate with variation in lambda/g only
    # ----------------------------------------------------------------------------------------
    # i invade j first
    i_invade_var_lambda <- rep (NA, 60)
    j_resident_var_lambda <- rep (NA, 60)
    temp <- 1
    
    for (t in 60:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      
      #invader
      i_invade_var_lambda[temp] <- pop.invade(N0=1, resident=ceiling(N_j_var_lambda[t]), s=par.i$s, 
                                         g=par.i$g, g_r=par.j$g, 
                                         a_inter=a_ij_weighted, lambda=par.i$lambda)
      
      #resident
      j_resident_var_lambda[temp] <- pop.resident(N0=1, resident=ceiling(N_j_var_lambda[t]), s=par.j$s, 
                                       g=par.j$g, g_i=par.i$g, 
                                       a_intra = a_jj_weighted, a_inter=a_ji_weighted, 
                                       lambda=par.j$lambda)/ceiling(N_j_var_lambda[t])
      
      temp  <- temp + 1 
    }
    
    # then have j invade into i
    j_invade_var_lambda <- rep (NA, 60)
    i_resident_var_lambda <- rep (NA, 60)
    temp <- 1
    
    for (t in 60:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      
      # invader
      j_invade_var_lambda[temp] <- pop.invade(N0=1, resident=ceiling(N_i_var_lambda[t]), s=par.j$s, 
                                         g=par.j$g, g_r=par.i$g, 
                                         a_inter=a_ji_weighted, lambda=par.j$lambda)
      
      #resident
      i_resident_var_lambda[temp] <- pop.resident(N0=1, resident=ceiling(N_i_var_lambda[t]), s=par.i$s, 
                                       g=par.i$g, g_i=par.j$g,
                                       a_intra = a_ii_weighted, a_inter=a_ij_weighted, 
                                       lambda=par.i$lambda)/ceiling(N_i_var_lambda[t])
      

      temp  <- temp + 1 
    }
    
    i_inv_epsilon_lambda <- log(i_invade_var_lambda) - i_inv_epsilon_0
    j_inv_epsilon_lambda <- log(j_invade_var_lambda) - j_inv_epsilon_0
    
    i_res_epsilon_lambda <- log(i_resident_var_lambda) - i_res_epsilon_0
    j_res_epsilon_lambda <- log(j_resident_var_lambda) - j_res_epsilon_0
    
    # ----------------------------------------------------------------------------------------
    #### Calculate the invasion rate with variable alphas, 
    # but with NO variation in intrinsic growth rates and g
    # ----------------------------------------------------------------------------------------
    # find resident equilibrium with variable alphas
    # for i
    N_i_var_alpha <- rep(NA, time)
    N_i_var_alpha[1] <- N0_i
    
    for (t in 1:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      N_i_var_alpha[t+1] <- pop.equilibrium(N0=N_i_var_alpha[t], s=par.i$s, g=i_g_weighted, 
                                            a_intra=par.i$a_ii, lambda=i_lambda_weighted)
    }
    
    
    # for j
    N_j_var_alpha <- rep(NA, time)
    N_j_var_alpha[1] <- N0_j
    
    for (t in 1:time) {
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      N_j_var_alpha[t+1] <- pop.equilibrium(N0=N_j_var_alpha[t], s=par.j$s, g=j_g_weighted, 
                                            a_intra=par.j$a_jj, lambda=j_lambda_weighted)
    }
    
    # Then invade each species into the other at equilibrium with variation in alpha only
    # i invade j first
    i_invade_var_alpha <- rep (NA, 60)
    j_resident_var_alpha <- rep (NA, 60)
    temp <- 1
    
    for (t in 60:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      
      #invader
      i_invade_var_alpha[temp] <- pop.invade(N0=1, resident=ceiling(N_j_var_alpha[t]), s=par.i$s, 
                                        g=i_g_weighted, g_r=j_g_weighted, 
                                        a_inter=par.i$a_ij, lambda=i_lambda_weighted)
      
      #resident
      j_resident_var_alpha[temp] <- pop.resident(N0=1, resident=ceiling(N_j_var_alpha[t]), s=par.j$s, 
                                      g=j_g_weighted, g_i=i_g_weighted, 
                                      a_intra = par.j$a_jj, a_inter=par.j$a_ji, 
                                      lambda=j_lambda_weighted)/ ceiling(N_j_var_alpha[t])
      

      temp  <- temp + 1 
    }
    
    # then have j invade into i
    j_invade_var_alpha <- rep (NA, 60)
    i_resident_var_alpha <- rep (NA, 60)
    temp <- 1
    
    for (t in 60:time) {
      par.i <- subset(i.dat, raintype==rainsummary$raintype[t])
      par.j <- subset(j.dat, raintype==rainsummary$raintype[t])
      
      # invader
      j_invade_var_alpha[temp] <- pop.invade(N0=1, resident=ceiling(N_i_var_alpha[t]), s=par.j$s, 
                                        g=j_g_weighted, g_r=i_g_weighted, 
                                        a_inter=par.j$a_ji, lambda=j_lambda_weighted)
      
      #resident
      i_resident_var_alpha <- pop.resident(N0=1, resident=ceiling(N_i_var_alpha[t]), s=par.i$s, 
                                      g=i_g_weighted, g_i=j_g_weighted, 
                                      a_intra = par.i$a_ii, a_inter=par.i$a_ij, 
                                      lambda=i_lambda_weighted)/ceiling(N_i_var_alpha[t])
      
      temp  <- temp + 1 
    }
    
    
    i_inv_epsilon_alpha <- log(i_invade_var_alpha) - i_inv_epsilon_0
    j_inv_epsilon_alpha <- log(j_invade_var_alpha) - j_inv_epsilon_0
    
    i_res_epsilon_alpha <- log(i_resident_var_alpha) - i_res_epsilon_0
    j_res_epsilon_alpha <- log(j_resident_var_alpha) - j_res_epsilon_0
    # ----------------------------------------------------------------------------------------
    #### Finally calculate the invasion rate with the interaction terms
    # ----------------------------------------------------------------------------------------
    
    ### 
    
    #invaders
    i_inv_epsilon_interaction <- i_inv - 
      (i_inv_epsilon_0 + i_inv_epsilon_alpha + i_inv_epsilon_lambda)
    j_inv_epsilon_interaction <- j_inv - 
      (j_inv_epsilon_0 + j_inv_epsilon_alpha + j_inv_epsilon_lambda)
    
    # residents
    i_res_epsilon_interaction <- i_res - 
      (i_res_epsilon_0 + i_res_epsilon_alpha + i_res_epsilon_lambda)
    j_res_epsilon_interaction <- j_res - 
      (j_res_epsilon_0 + j_res_epsilon_alpha + j_res_epsilon_lambda)
    
    ##
    i_delta_r <- mean_f(i_inv-j_res)
    j_delta_r <- mean_f(j_inv-i_res)
    
    i_delta_0 <- mean_f(i_inv_epsilon_0-j_res_epsilon_0)
    j_delta_0 <- mean_f(j_inv_epsilon_0-i_res_epsilon_0)
    
    i_delta_alpha <- mean_f(i_inv_epsilon_alpha-j_res_epsilon_alpha)
    j_delta_alpha <- mean_f(j_inv_epsilon_alpha-i_res_epsilon_alpha)
    
    i_delta_lambda <- mean_f(i_inv_epsilon_lambda-j_res_epsilon_lambda)
    j_delta_lambda <- mean_f(j_inv_epsilon_lambda-i_res_epsilon_lambda)
    
    i_delta_interaction <- mean_f(i_inv_epsilon_interaction-j_res_epsilon_interaction)
    j_delta_interaction <- mean_f(j_inv_epsilon_interaction-i_res_epsilon_interaction)
    # ----------------------------------------------------------------------------------------    
    #### Compile the data
    # ----------------------------------------------------------------------------------------
    
    # with invader resident comparison
    ir_i_results <- c(i_delta_r, 
                               i_delta_0, 
                               i_delta_alpha, 
                               i_delta_lambda, 
                               i_delta_interaction)
    
    ir_j_results <- c(j_delta_r, 
                               j_delta_0, 
                               j_delta_alpha, 
                               j_delta_lambda, 
                               j_delta_interaction)
    
    partitioning_i<-rbind(partitioning_i,ir_i_results)
    partitioning_j<-rbind(partitioning_j,ir_j_results)
    
  }
  partitioning_1000[[rep]]<-list(partitioning_i,
                        partitioning_j)  
}

saveRDS(partitioning_1000,"Output/partitioning_1000.rds")

library(data.table)

pairs <- combn(1:8, 2)

par1000.i<-  rbindlist(lapply(partitioning_1000, function(x) x[[1]]),use.names = F)

par1000.j <- rbindlist(lapply(partitioning_1000, function(x) x[[2]]),use.names = F)

par1000.i<-data.frame(pair=rep(c(1:28),1000),mark="i",sp=rep(pairs[1,],1000),par1000.i)

colnames(par1000.i)<-c("pair","mark","sp","delta_r","delta_0","delta_alpha","delta_lambda","delta_inter")

par1000.j<-data.frame(pair=rep(c(1:28),1000),mark="j",sp=rep(pairs[2,],1000),par1000.j)

colnames(par1000.j)<-c("pair","mark","sp","delta_r","delta_0","delta_alpha","delta_lambda","delta_inter")

par1000<-rbind(par1000.i,par1000.j)

write.csv(par1000,"Output/par1000.csv")


