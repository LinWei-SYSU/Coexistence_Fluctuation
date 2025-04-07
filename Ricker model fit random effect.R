### Ricker parameter with Pot 


com<-read.csv("Competition.csv")

library(brms)
library(dplyr)

com$Tr<-ifelse(com$Raintype=="WW",1,
               ifelse(com$Raintype=="WD",2,
                      ifelse(com$Raintype=="DW",3,4)))

com$Tr<-as.integer(com$Tr)

###running your brms model with backend="cmdstanr" if you get error 
###"Sys.setenv(R_MAKEVARS_USER = NULL): 'names' and 'val' are of different lengths"

###PHAN###
sp <- unique(com$Species)[1];sp

dat.phan <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.PHAN.Pot<-brm(data = dat.phan, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.PHAN.Pot)

###DAAE###
sp <- unique(com$Species)[2];sp

dat.daae <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.DAAE.Pot<-brm(data = dat.daae, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.DAAE.Pot)

###BIPI###
sp <- unique(com$Species)[3];sp

dat.bipi <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.BIPI.Pot<-brm(data = dat.bipi, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.BIPI.Pot)

###AGCO###
sp <- unique(com$Species)[4];sp

dat.agco <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.AGCO.Pot<-brm(data = dat.agco, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.AGCO.Pot)

###POHY###
sp <- unique(com$Species)[5];sp

dat.pohy <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.POHY.Pot<-brm(data = dat.pohy, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.POHY.Pot)

###PRCL###
sp <- unique(com$Species)[6];sp

dat.prcl <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.PRCL.Pot<-brm(data = dat.prcl, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL,lb=0),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.PRCL.Pot)

###OPCO###
sp <- unique(com$Species)[7];sp

dat.opco <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.OPCO.Pot<-brm(data = dat.opco, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.OPCO.Pot)

###LUHY###
sp <- unique(com$Species)[8];sp

dat.luhy <- com %>% 
  filter(Species==sp) %>%
  select(c(num_seeds,n_AGCO,n_BIPI,n_OPCO,n_PHAN,n_DAAE,n_POHY,n_PRCL,n_LUHY,Tr,Pot)) %>% 
  mutate(Pot=as.integer(as.factor(as.character(Pot))))

Ricker.LUHY.Pot<-brm(data = dat.luhy, 
                     family = gaussian(),
                     bf(num_seeds ~ lambda * exp(-(aAGCO*n_AGCO+
                                                     aBIPI*n_BIPI+
                                                     aDAAE*n_DAAE+
                                                     aLUHY*n_LUHY+
                                                     aOPCO*n_OPCO+
                                                     aPHAN*n_PHAN+
                                                     aPOHY*n_POHY+
                                                     aPRCL*n_PRCL)),
                        lambda~0+ factor(Tr)+(1|Pot),
                        aAGCO~0 + factor(Tr)+(1|Pot),
                        aBIPI~0 + factor(Tr)+(1|Pot),
                        aDAAE~0 + factor(Tr)+(1|Pot),
                        aLUHY~0 + factor(Tr)+(1|Pot),
                        aOPCO~0 + factor(Tr)+(1|Pot),
                        aPHAN~0 + factor(Tr)+(1|Pot),
                        aPOHY~0 + factor(Tr)+(1|Pot),
                        aPRCL~0 + factor(Tr)+(1|Pot),
                        nl = TRUE
                     ),
                     prior = c(prior(normal(10, 10000),  class = b, nlpar = lambda,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aAGCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aBIPI),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aDAAE),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aLUHY,lb=0),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aOPCO),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPHAN),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPOHY),
                               prior(normal( 0.1 , 10 ), class = b, nlpar = aPRCL),
                               prior(exponential(1), class = sigma)),
                     iter = 4000, warmup = 1000, chains = 4, cores = 4, backend="cmdstanr"
)

summary(Ricker.LUHY.Pot)$fixed[1][1:4,]

agco_pot<-summary(Ricker.AGCO.Pot)

agco_pot$fixed[1][1:4,]
agco_pot$fixed[5:36,1:7]


###create a spreadsheet contain the mean and CI of parameter##

sp_name<-c("AGCO","BIPI","DAAE","LUHY","OPCO","PHAN","POHY","PRCL")

model_name<-ls(pattern="^Ricker\\.")

par_list<-lapply(model_name, function(x){
  data<-get(x)
  alpha<-summary(data)$fixed[5:36,1:7]
  colnames(alpha)<-c("alpha","alpha_sd","alpha_lowCI","alpha_highCI","alpha_Rhat","alpha_Bulk_ESS","alpha_Tail_ESS")
  lambda<-summary(data)$fixed[1:4,1:7]
  colnames(lambda)<-c("lambda","lambda_sd","lambda_lowCI","lambda_highCI","lambda_Rhat","lambda_Bulk_ESS","lambda_Tail_ESS")
  combined<-cbind(alpha,lambda)
})

pars_ricker<-do.call(rbind,par_list)

pars_ricker<-cbind(data.frame(focal=rep(sp_name,each=32),competitor=rep(sp_name,each=4),raintype=rep(1:4)),pars_ricker)

rownames(pars_ricker)<-NULL

write.csv(pars_ricker,"Output/pars_ricker_re.csv")
