library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(ggrepel)

setwd("D:/PhD/100_Thesis/110_Chapter_1/annualfluc/Coexistence_Fluctuation")
###the mutual invasion rate for each species pairs and 1000 replicate
p1000<-read.csv("Output/par1000.csv")


p1000_a<-p1000 %>%
  mutate(delta_0a=delta_0+delta_alpha,
         delta_0al=delta_0a+delta_lambda) %>%
  pivot_longer(cols = starts_with("delta"), names_to = "term", values_to = "value") %>%
  pivot_wider(names_from = mark, values_from = value) %>%
  select(pair,sp, i, j, term) %>% 
  filter(term=="delta_0"|term=="delta_0a"|term=="delta_0al"|term=="delta_r")

sp_i<-p1000_a[1:112000,];sp_j<-p1000_a[112001:224000,]


df<-data.frame(rep=rep(1:1000,each=28*4),term=sp_i$term,
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
                          levels = c("delta_0", "delta_0a","delta_0al","delta_r"))

df_mean<-df_summary %>% 
  group_by(term) %>% 
  summarise(mean_CE=mean(CE)) %>% 
  mutate(sp="mean")

###Add speices pairs name

ggplot(df_summary)+
  geom_point(aes(x=term,y=CE,group = species_pairs),col="grey70")+
  geom_line(aes(x=term,y=CE,group = species_pairs),col="grey70")+
  geom_point(data=df_mean,aes(x=term,y=mean_CE),size=2,color="#063900")+
  geom_line(data=df_mean,aes(x=term,y=mean_CE,group = sp),linewidth=1,col="#063900")+
  theme_bw()+
  labs(x="",y=expression(italic("p(Co)")))+
  scale_x_discrete(labels=c(expression(Delta^0),
                            "",
                            "",
                            expression(Delta^r)))+
  #geom_text_repel(data=df_summary %>% filter(term=="delta_r"),aes(x=4, y=CE, label = species_pairs))+
  theme(axis.text = element_text(size=12))+
  annotate("text", x = 1.5, y = -0.02, label = expression(+Delta^alpha), vjust = 1,size = 4)+
  annotate("text", x = 2.5, y = -0.02, label = expression(+Delta^lambda), vjust = 1,size = 4)+
  annotate("text", x = 3.5, y = -0.02, label = expression(+Delta^{alpha*lambda}), vjust = 1,size = 4)

ggsave("Figure/Partitioning_visualization/possibility.pdf",width = 6,height = 6)




