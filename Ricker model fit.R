###Ricker model fit###

com<-read.csv("Competition.csv")

library(rethinking)
library(dplyr)

com$Tr<-ifelse(com$Raintype=="WW",1,
               ifelse(com$Raintype=="WD",2,
                      ifelse(com$Raintype=="DW",3,4)))

com$Tr<-as.integer(com$Tr)


#############################
##Model fit for each species


######PHAN######
##############
sp <- unique(com$Species)[1];sp

dat.phan <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.PHAN <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.phan,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )


######DAAE######
##############
sp <- unique(com$Species)[2];sp

dat.daae <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.DAAE <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.daae,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######BIPI######
##############
sp <- unique(com$Species)[3];sp

dat.bipi <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.BIPI <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.bipi,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######AGCO######
##############
sp <- unique(com$Species)[4];sp

dat.agco <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.AGCO <-  rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.agco,
  start=list(
    lambda = c(rep(1000,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  constraints=list(lambda="lower=0"),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######POHY######
##############
sp <- unique(com$Species)[5];sp

dat.pohy <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.POHY <-  rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.pohy,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######PRCL######
##############
sp <- unique(com$Species)[6];sp

dat.prcl <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.PRCL <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.prcl,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######OPCO######
##############
sp <- unique(com$Species)[7];sp

dat.opco <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.OPCO <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dexp(1)
  ) ,
  data = dat.opco,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )

######LUHY######
##############
sp <- unique(com$Species)[8];sp

dat.luhy <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.LUHY <- rethinking::ulam( 
  alist(
    num_seeds ~ dnorm( mu , sigma ) ,
    mu <- lambda[Tr]*exp(-(a_AGCO[Tr]*n_AGCO+
                             a_BIPI[Tr]*n_BIPI+
                             a_DAAE[Tr]*n_DAAE+
                             a_LUHY[Tr]*n_LUHY+
                             a_OPCO[Tr]*n_OPCO+
                             a_PHAN[Tr]*n_PHAN+
                             a_POHY[Tr]*n_POHY+
                             a_PRCL[Tr]*n_PRCL
    )),
    lambda[Tr] ~ dnorm( 100 , 10000 ),
    a_AGCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_BIPI[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_DAAE[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_LUHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_OPCO[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PHAN[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_POHY[Tr] ~ dnorm( 0.1 , 10 ) ,
    a_PRCL[Tr] ~ dnorm( 0.1 , 10 ) ,
    sigma ~ dcauchy( 0 , 10 )
  ) ,
  data = dat.luhy,
  start=list(
    lambda = c(rep(100,4)),
    a_AGCO = c(rep(0.1,4)),
    a_BIPI = c(rep(0.1,4)),
    a_DAAE = c(rep(0.1,4)),
    a_LUHY = c(rep(0.1,4)),
    a_OPCO = c(rep(0.1,4)),
    a_PHAN = c(rep(0.1,4)),
    a_POHY = c(rep(0.1,4)),
    a_PRCL = c(rep(0.1,4))
    
  ),
  chains=4 , cores=4,iter=4000 , warmup=1000,log_lik =TRUE )


###Extract the posterior distribution and save the results##

AGCO.par<-extract.samples(Ricker.AGCO,n=1000)
BIPI.par<-extract.samples(Ricker.BIPI,n=1000)
DAAE.par<-extract.samples(Ricker.DAAE,n=1000)
LUHY.par<-extract.samples(Ricker.LUHY,n=1000)
OPCO.par<-extract.samples(Ricker.OPCO,n=1000)
PHAN.par<-extract.samples(Ricker.PHAN,n=1000)
POHY.par<-extract.samples(Ricker.POHY,n=1000)
PRCL.par<-extract.samples(Ricker.PRCL,n=1000)

par<-list(AGCO.par,BIPI.par,DAAE.par,LUHY.par,OPCO.par,PHAN.par,POHY.par,PRCL.par)

saveRDS(par,file = "Output/post_1000.rds")


###create a spreadsheet contain the mean and CI of parameter##

sp_name<-c("AGCO","BIPI","DAAE","LUHY","OPCO","PHAN","POHY","PRCL")

model_name<-ls(pattern="^Ricker\\.")

par_list<-lapply(model_name, function(x){
  data<-get(x)
  alpha<-precis(data,depth = 2,digits=3, prob=0.95)[5:36,1:4]
  colnames(alpha)<-c("alpha","alpha_sd","alpha_lowCI","alpha_highCI")
  lambda<-precis(data,depth = 2,digits=3, prob=0.95)[1:4,1:4]
  colnames(lambda)<-c("lambda","lambda_sd","lambda_lowCI","lambda_highCI")
  combined<-cbind(alpha,lambda)
})

pars_ricker<-do.call(rbind,par_list)

pars_ricker<-cbind(data.frame(focal=rep(sp_name,each=32),competitor=rep(sp_name,each=4),raintype=rep(1:4)),pars_ricker)

rownames(pars_ricker)<-NULL

write.csv(pars_ricker,"Output/pars_ricker.csv")

precis(Ricker.PRCL,depth=2)
