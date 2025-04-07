###Parameters visualization

post_1000<-readRDS("Output/post_1000.rds")

##extract Lambda

library(ggplot2)
library(dplyr)
library(ggdist)
library(tidyverse)
library(stringr)


###Lambda

Lambda.ls<-lapply(1:1000, function(i) ##each sample
  lapply(1:8, function(sp) ##each focal species
    lapply(post_1000[[sp]][1], function(x) x[i,]))) ##each lambda

lambda.df<-do.call(rbind, lapply(Lambda.ls, function(y) {
  do.call(rbind, lapply(y, function(x) {
    do.call(rbind, x)
  }))
}))

lambda<-data.frame(rep=rep(1:1000,each=8),species=rep(1:8),lambda.df)

lambda<-lambda %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "lambda")

mapping <- c("1" = "AGCO", "2" = "BIPI", "3" = "DAAE", "4" = "LUHY", "5" = "OPCO", "6" = "PHAN", "7" = "POHY", "8" = "PRCL")

lambda$species <- stringr::str_replace_all(lambda$species, mapping)


ggplot(lambda,aes(x=log(lambda),y=factor(species, 
                                         levels = rev(levels(factor(species)))),
                  color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Intrinsic growth rate", y ="Posterior distribution",subtitle = "(a)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = c(.18, .86),
    legend.background = element_rect(fill = NA),
    legend.title = element_blank()
  ) 

ggsave("Figure/Parameters/Par_lambda_halfeye.pdf",width=6,height = 6)


#######
####alpha##

##AGCO

AGCO_alpha_df<-data.frame()

for (i in 2:9) {

  t<-as.data.frame(post_1000[[1]][i],col.names = names(c("X1","X2","X3","X4")))
  AGCO_alpha_df<-bind_rows(AGCO_alpha_df,t)
                              
}


AGCO_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),AGCO_alpha_df)

AGCO_alpha<-AGCO_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

AGCO_alpha$species <- str_replace_all(AGCO_alpha$species, mapping)


ggplot(AGCO_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                  color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: AGCO)", y ="",subtitle = "(b)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(AGCO)_halfeye.pdf",width=6,height = 6)


##BIPI

BIPI_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[2]][i],col.names = names(c("X1","X2","X3","X4")))
  BIPI_alpha_df<-bind_rows(BIPI_alpha_df,t)
  
}


BIPI_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),BIPI_alpha_df)

BIPI_alpha<-BIPI_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

BIPI_alpha$species <- str_replace_all(BIPI_alpha$species, mapping)


ggplot(BIPI_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: BIPI)", y ="",subtitle = "(c)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(BIPI)_halfeye.pdf",width=6,height = 6)


##DAAE

DAAE_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[3]][i],col.names = names(c("X1","X2","X3","X4")))
  DAAE_alpha_df<-bind_rows(DAAE_alpha_df,t)
  
}


DAAE_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),DAAE_alpha_df)

DAAE_alpha<-DAAE_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

DAAE_alpha$species <- str_replace_all(DAAE_alpha$species, mapping)


ggplot(DAAE_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: DAAE)", y ="Posterior distribution",subtitle = "(d)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(DAAE)_halfeye.pdf",width=6,height = 6)


##LUHY

LUHY_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[4]][i],col.names = names(c("X1","X2","X3","X4")))
  LUHY_alpha_df<-bind_rows(LUHY_alpha_df,t)
  
}


LUHY_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),LUHY_alpha_df)

LUHY_alpha<-LUHY_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

LUHY_alpha$species <- str_replace_all(LUHY_alpha$species, mapping)


ggplot(LUHY_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: LUHY)", y ="",subtitle = "(e)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(LUHY)_halfeye.pdf",width=6,height = 6)


##OPCO

OPCO_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[5]][i],col.names = names(c("X1","X2","X3","X4")))
  OPCO_alpha_df<-bind_rows(OPCO_alpha_df,t)
  
}


OPCO_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),OPCO_alpha_df)

OPCO_alpha<-OPCO_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

OPCO_alpha$species <- str_replace_all(OPCO_alpha$species, mapping)


ggplot(OPCO_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: OPCO)", y ="",subtitle = "(f)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(OPCO)_halfeye.pdf",width=6,height = 6)


##PHAN

PHAN_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[6]][i],col.names = names(c("X1","X2","X3","X4")))
  PHAN_alpha_df<-bind_rows(PHAN_alpha_df,t)
  
}


PHAN_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),PHAN_alpha_df)

PHAN_alpha<-PHAN_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

PHAN_alpha$species <- str_replace_all(PHAN_alpha$species, mapping)


ggplot(PHAN_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: PHAN)", y ="Posterior distribution",subtitle = "(g)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(PHAN)_halfeye.pdf",width=6,height = 6)

##POHY

POHY_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[7]][i],col.names = names(c("X1","X2","X3","X4")))
  POHY_alpha_df<-bind_rows(POHY_alpha_df,t)
  
}


POHY_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),POHY_alpha_df)

POHY_alpha<-POHY_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

POHY_alpha$species <- str_replace_all(POHY_alpha$species, mapping)


ggplot(POHY_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: POHY)", y ="",subtitle = "(h)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(POHY)_halfeye.pdf",width=6,height = 6)

##PRCL

PRCL_alpha_df<-data.frame()

for (i in 2:9) {
  
  t<-as.data.frame(post_1000[[8]][i],col.names = names(c("X1","X2","X3","X4")))
  PRCL_alpha_df<-bind_rows(PRCL_alpha_df,t)
  
}


PRCL_alpha<-data.frame(rep=rep(1:1000),species=rep(1:8,each=1000),PRCL_alpha_df)

PRCL_alpha<-PRCL_alpha %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "raintype", 
               values_to = "alpha")

PRCL_alpha$species <- str_replace_all(PRCL_alpha$species, mapping)


ggplot(PRCL_alpha,aes(x=alpha,y=factor(species, 
                                       levels = rev(levels(factor(species)))),
                      color=raintype,fill = raintype))+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = rel(1)
  ) +
  stat_halfeye(
    slab_alpha = 0.3,
    .width = c(0.5, 0.95),
    point_interval = "median_hdi"
  )+
  labs(x = "Interaction coefficients (i: PRCL)", y ="",subtitle = "(i)")+
  scale_y_discrete(expand = expansion(mult = c(0.03, 0))) +
  scale_color_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                     labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  scale_fill_manual(values=c('X1' ="#040082", 'X2' = "#0080FE",'X3' ="#D2B53A", 'X4' = "#60290C"),
                    labels=c("Consistent Wet", "Post-flood Dry", "Pre-flood Dry", "Consistent Dry"))+
  theme_bw()+
  theme(
    plot.subtitle = element_text(
      hjust = 0,
      vjust = 1,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ggsave("Figure/Parameters/Par_alpha(PRCL)_halfeye.pdf",width=6,height = 6)

