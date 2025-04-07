library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(scatterpie)

###the mutual invasion rate for each species pairs and 1000 replicate
p1000<-read.csv("Output/par1000.csv")

p1000_a<-p1000 %>%
  pivot_longer(cols = starts_with("delta"), names_to = "term", values_to = "value") %>%
  pivot_wider(names_from = mark, values_from = value) %>%
  select(pair,sp, i, j, term)

sp_i<-p1000_a[1:140000,];sp_j<-p1000_a[140001:280000,]


df<-data.frame(rep=rep(1:1000,each=28*5),term=sp_i$term,
               sp_i=sp_i$sp,sp_j=sp_j$sp,
               gr_i=sp_i$i,gr_j=sp_j$j
) %>% 
  mutate(species_pairs=paste0(sp_i,"-",sp_j))

mapping <- c("1" = "AG", "2" = "BI", "3" = "DA", "4" = "LU", "5" = "OP", "6" = "PH", "7" = "PO", "8" = "PR")

df$species_pairs <- stringr::str_replace_all(df$species_pairs, mapping)



df_summary<-df %>% 
  group_by(term,species_pairs) %>% 
  summarise(med_gr_i=median(gr_i,na.rm = T),
            med_gr_j=median(gr_j,na.rm = T),
            CE = mean(gr_i > 0 & gr_j > 0,na.rm = T),
            PR = mean(gr_i < 0 & gr_j < 0,na.rm = T),
            Win_i = mean(gr_i > 0 & gr_j < 0,na.rm = T),
            Win_j = mean(gr_i < 0 & gr_j > 0,na.rm = T)) 

df_summary$term <- factor(df_summary$term, 
                          levels = c("delta_r", "delta_0", "delta_alpha",
                                     "delta_lambda", "delta_inter"))

df_hist<-df_summary %>%filter(term=="delta_r"|term=="delta_0")

#############################

####comparison between invasion growth rate with/without environment fluctuation 

#############################

ggplot(df_hist, aes(x=med_gr_i, y=med_gr_j,col=term,shape=term)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "gray",alpha=0.5) +
  geom_point(size=2) + 
  labs( x = expression(Delta * " of sp1") , y = expression(Delta * " of sp2"))+
  theme_bw() +
  scale_color_manual(values=c('delta_0' ="#818181",
                              'delta_r' = "#B65FCE"),
                     labels=c(expression(bar(Delta)[i]^r),
                              expression(bar(Delta)[i]^0)))+
  scale_shape_manual(values=c('delta_0' =16,
                              'delta_r' = 15),
                     labels=c(expression(bar(Delta)[i]^r),
                              expression(bar(Delta)[i]^0)))+
  geom_hline(yintercept=0, alpha=0.5) + 
  geom_vline(xintercept=0, alpha=0.5) +
  geom_path(data=df_hist,
            aes(x = med_gr_i, y = med_gr_j, group = species_pairs), linetype=5,
            arrow = arrow(length=unit(0.05,"inches"), 
                          ends="last", type = "closed"), size = 0.25,col="black")+
  geom_text_repel(data = df_hist %>% filter(term=="delta_r"), 
                  aes(x=med_gr_i, y=med_gr_j, label = species_pairs), 
                  size = 2,
                  color="black",
                  max.overlaps = Inf,
                  box.padding = 0.3)+
  theme(legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = c(0.86, 0.06),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  guides(color = guide_legend(ncol = 2)) 

ggsave("Figure/Partitioning_visualization/deltar&0 comparison_point_median.pdf",width = 7,height =7)



#############################

##Fluctuation dependent mechanisms

#############################

df_fluc<-df_summary %>%filter(term=="delta_alpha"|term=="delta_lambda"|term=="delta_inter")

df_fluc<- df_fluc%>%
  mutate(term = factor(term,
                       levels = c("delta_alpha", "delta_lambda","delta_inter"),
                       labels = c(expression(italic(bar(Delta)[i]^alpha)), 
                                  expression(italic(bar(Delta)[i]^lambda)),
                                  expression(italic(bar(Delta)[i]^{alpha*lambda})))))
###point graph


ggplot(df_fluc,aes(x=med_gr_i, y=med_gr_j)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "gray",alpha=0.2) +
  geom_hline(yintercept=0, ) + 
  geom_vline(xintercept=0, ) +
  geom_point(aes(color=term))+
  labs( x = expression(Delta * " of sp1") , y = expression(Delta * " of sp2"))+
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 23),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 20))+
  geom_text_repel(aes(x=med_gr_i, y=med_gr_j, label = species_pairs), 
                  size = 3,
                  color="black",
                  max.overlaps = Inf,
                  box.padding = 0.3)+
  facet_wrap(~term,scales="free",labeller = label_parsed)

ggsave("Figure/Partitioning_visualization/fluctuation_scalefree_median.pdf",width = 15,height =6)


##################

#####for 28 species pairs bar chart

##################

##data
df_28<-p1000 %>%
  pivot_longer(cols = starts_with("delta"), names_to = "term", values_to = "value") %>% 
  group_by(pair,sp,mark,term) %>% 
  summarise(
    median=median(value,na.rm = T)
  ) %>% 
  mutate(term=factor(term,levels = c("delta_r", "delta_0", "delta_alpha",
                                     "delta_lambda", "delta_inter"))) %>%
  group_by(pair) %>% 
  mutate(sp_pair=paste(unique(sp),collapse = "-"))  ##get the species pairs

mapping <- c("1" = "AGCO", "2" = "BIPI", "3" = "DAAE", "4" = "LUHY", 
             "5" = "OPCO", "6" = "PHAN", "7" = "POHY", "8" = "PRCL")

df_28$sp <- stringr::str_replace_all(df_28$sp, mapping)
df_28$sp_pair <- stringr::str_replace_all(df_28$sp_pair, mapping)

df_fluc<-df_28 %>%
  filter(term=="delta_alpha"|term=="delta_lambda"|term=="delta_inter")

### all component

ggplot(df_28,aes(x=term,y=median,fill = mark))+
  geom_bar(stat="identity",position="dodge")+
  scale_x_discrete(labels=c(expression(bar(Delta)[i]^r),
                            expression(bar(Delta)[i]^0),
                            expression(bar(Delta)[i]^alpha),
                            expression(bar(Delta)[i]^lambda),
                            expression(bar(Delta)[i]^{alpha*lambda})))+
  ylab("Growth rate when rare")+xlab("")+
  theme_bw()+
  geom_vline(xintercept = 1.5,linetype=5,alpha=0.5)+ 
  facet_wrap(~sp_pair,scales = "free_y",dir="v",nrow = 7)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=25),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title =  element_text(size = 22))

ggsave("Figure/Partitioning_visualization/all_pair_ij_5compo_median.png",width = 18,height =20)
