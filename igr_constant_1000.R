###Invasion growth rate results for 1000 sample from posterior distribution

library(dplyr)
library(tidyverse)

###Invasion growth rate under Historic rainfall fluctuation
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


# ------------------------------------------------------------------------------------
###Constant Rainfall
# ------------------------------------------------------------------------------------


list.consistent_1000<-list();i.dat<-data.frame();j.dat<-data.frame() 

pairs <- combn(1:8, 2); pairs

for (rep in 1:1000) {
  
  A.est<-AA.ls[[rep]]
  Lambda.est<-Lambda.ls[[rep]]
  
  list.consistent_1000[[rep]]<-list()
  
  for (p in 1:ncol(pairs)) {##for each species-pair
    
    i <- pairs[1,p]; j <- pairs[2,p] ##extract the species from species pair
    
    consistent.out.all<-data.frame()
    
    for(r in 1:4){ ##for each raintype
      
      ##parameters for i
      par.i<-data.frame(raintype=r,
                        species=i,
                        a_ii=unlist(A.est[[i]][[i]][r]),
                        a_ij=unlist(A.est[[i]][[j]][r]),
                        lambda=unlist(Lambda.est[[i]])[r],
                        s=s.ls[[r]][i],
                        g=g.ls[[r]][i]) %>%
        mutate(eta=lambda*g/(1-(1-g)*s)) %>% 
        mutate(Nstar_i =log(eta)/(a_ii*g)) ##Equilibrium abundance for i
      
      ##parameters for j
      par.j<-data.frame(raintype=r,
                        species=j,
                        a_ji=unlist(A.est[[j]][[i]][r]),
                        a_jj=unlist(A.est[[j]][[j]][r]),
                        lambda=unlist(Lambda.est[[j]])[r],
                        s=s.ls[[r]][j],
                        g=g.ls[[r]][j]) %>%
        mutate(eta=lambda*g/(1-(1-g)*s)) %>% 
        mutate(Nstar_j =log(eta)/(a_jj*g))  ##Equilibrium abundance for j
      
      
      #### Equilibrium abundance
      time <- length(rainsummary$raintype) ## 120 year
      
      # for i
      N_i <- rep(NA, time)
      N_i[1] <-  ceiling(par.i$Nstar_i) ##make sure N larger or equal to 1
      for (t in 1:time) {
        N_i[t+1] <- pop.equilibrium(N0=N_i[t], s=par.i$s, g=par.i$g, 
                                    a_intra=par.i$a_ii, lambda=par.i$lambda)
      }
      
      plot(seq(1:(time+1)), N_i, type="l")
      
      # for j
      N_j <- rep(NA, time)
      N_j[1] <- ceiling(par.j$Nstar_j)
      for (t in 1:time) {
        N_j[t+1] <- pop.equilibrium(N0=N_j[t], s=par.j$s, g=par.j$g, 
                                    a_intra=par.j$a_jj, lambda=par.j$lambda)
      }
      
      plot(seq(1:(time+1)), ceiling(N_j), type="l")
      
      ####Mutual Invasion
      
      # i invade into j
      i_invader <- rep (NA, 80) ###select the year 40-120 to calculate the igr
      j_resident <- rep (NA, 80)
      
      temp <- 1
      
      #####for year 40-120
      for (t in 40:time) {
        i_invader[temp] <- pop.invade(N0=1, resident=ceiling(N_j[t]), s=par.i$s, g=par.i$g, g_r=par.j$g,
                                      a_inter=par.i$a_ij, lambda=par.i$lambda) 
        temp  <- temp + 1 
      }
      
      # j invade into i
      j_invader <- rep (NA, 80) ###select the year 40-120 to calculate the igr
      i_resident <- rep (NA, 80)
      
      temp <- 1
      
      for (t in 40:time) {
        
        j_invader[temp] <- pop.invade(N0=1, resident=ceiling(N_i[t]), s=par.j$s, g=par.j$g, g_r=par.i$g,
                                      a_inter=par.j$a_ji, lambda=par.j$lambda)
        
        temp  <- temp + 1 
      }
      
      ###log-transformed##
      i_inv <- log(i_invader) 
      j_inv <- log(j_invader)
      
      
      
      consistent.out <- as.data.frame(cbind( i_inv, j_inv)) %>%
        gather(invader, igr) %>% ##
        mutate(raintype = r) %>%
        group_by(raintype, invader) %>% 
        summarize(mean_igr = mean(igr[is.finite(igr)]), 
                  sd_igr = sd(igr[is.finite(igr)]), 
                  se_igr = calcSE(igr[is.finite(igr)])) %>%
        as.data.frame()%>% mutate(invader=case_when(
          invader=="i_inv"~i,
          invader=="j_inv"~j),
        ) 
      
      ###combine the data for 4 raintype
      consistent.out.all <- rbind(consistent.out.all,consistent.out)
      
    }
    
    list.consistent_1000[[rep]][[p]]<- consistent.out.all
    
  }
}

saveRDS(list.consistent_1000,"Output/list.consistent_1000.rds")


