###
setwd("D:/PhD/100_Thesis/110_Chapter_1/annualfluc/Coexistence_Fluctuation")
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(stringr)
library(ggrepel)
library(scatterpie)


###the mutual invasion rate for each species pairs and 1000 replicate
list.consistent_1000<-readRDS("Output/list.consistent_1000.rds")

pairs <- combn(1:8, 2)

##unlist the big list to a data.frame
Consistent_1000 <- do.call(rbind, lapply(list.consistent_1000, function(x) {
  do.call(rbind,lapply(x, function(y) {
    do.call(cbind,y)
  }))
})) 

Data<-Consistent_1000%>%
  as.data.frame() %>% 
  select(raintype ,invader,mean_igr)

Data<-data.frame(rep=rep(1:1000,each=ncol(pairs)*2*4),Data)

###separate odd lines from even lines,odd lines represent the sp_i and even lines represent the sp_j
sp_i <- Data[(1:nrow(Data)) %% 2 != 0,]
sp_j <- Data[(1:nrow(Data)) %% 2 == 0,]

df_4constant<-data.frame(rep=sp_i$rep,raintype=sp_i$raintype,
               sp_i=sp_i$invader,sp_j=sp_j$invader,
               igr_i=sp_i$mean_igr,igr_j=sp_j$mean_igr)

df_4constant%>% 
  mutate(species_pairs=paste0(sp_i,"-",sp_j)) 





###the mutual invasion rate for each species pairs and 1000 replicate
p1000<-read.csv("Output/par1000.csv")

p1000_a<-p1000 %>%
  pivot_longer(cols = starts_with("delta"), names_to = "term", values_to = "value") %>%
  pivot_wider(names_from = mark, values_from = value) %>%
  select(pair,sp, i, j, term)

sp_i<-p1000_a[1:140000,];sp_j<-p1000_a[140001:280000,]


p1000<-data.frame(rep=rep(1:1000,each=28*5),term=sp_i$term,
               sp_i=sp_i$sp,sp_j=sp_j$sp,
               igr_i=sp_i$i,igr_j=sp_j$j
) %>% rename(raintype = term) %>% 
  filter(raintype=="delta_r"|raintype=="delta_0") %>% 
  mutate(raintype = case_when(
    raintype == "delta_r" ~ 6,
    raintype == "delta_0" ~ 5))


df_6<-bind_rows(df_4constant,p1000) 


mapping <- c("1" = "AGCO", "2" = "BIPI", "3" = "DAAE", "4" = "LUHY", 
             "5" = "OPCO", "6" = "PHAN", "7" = "POHY", "8" = "PRCL")  
df_6$sp_i <- str_replace_all(df_6$sp_i, mapping)
df_6$sp_j <- str_replace_all(df_6$sp_j, mapping)


df_6<-df_6 %>% mutate(label=paste0(substr(df_6$sp_i,1,2),"-",substr(df_6$sp_j,1,2)))
                
df_median<-df_6 %>% 
  group_by(raintype,label) %>% 
  summarise(med_igr_i=median(igr_i,na.rm=T),
            med_igr_j=median(igr_j,na.rm=T),
            CE = mean(igr_i > 0 & igr_j > 0,na.rm = T),
            PR = mean(igr_i < 0 & igr_j < 0,na.rm = T),
            Win_i = mean(igr_i > 0 & igr_j < 0,na.rm = T),
            Win_j = mean(igr_i < 0 & igr_j > 0,na.rm = T)) 

df_median$raintype <- factor(df_median$raintype, levels = c(1, 2, 5, 3, 4, 6))

#####pie graph 

lim <- 10 # limit x-axis at 25.
offset <- 1 # the amount to offset the red points from the limit. 

df_cap_median<-df_median %>% 
  mutate(med_igr_i=case_when(
    med_igr_i > lim~lim + offset,
    .default=med_igr_i
  )) %>% 
  mutate(flag=med_igr_i > lim)

pie_median<-ggplot(df_cap_median, aes(x=med_igr_i, y=med_igr_j)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#002B12",alpha=0.2) +
  geom_scatterpie(aes(x=med_igr_i, y=med_igr_j),pie_scale = .8,color=NA,
                  data=df_cap_median, cols=c("CE","PR","Win_i","Win_j"), alpha=.8)+
  geom_hline(yintercept=0, alpha=0.5) + 
  geom_vline(xintercept=0, alpha=0.5) +
  coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10), clip="off") +
  labs( x = "Invasion growth rate of species 1" , y = "Invasion growth rate of species 2")+
  scale_fill_manual(values=c('CE' ="#002B12", 'PR' = "#65676A",'Win_i' ="#0081C1", 'Win_j' = "#FBAF3B"), 
                    name = "Outcome", 
                    labels = c("Coexistence", "Priority Effect", "Species 1 Wins", "Species 2 Wins"))+
  theme_bw() +
  theme(legend.text = element_text(size = 33),
        legend.title = element_text(size = 33),
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 30),
        strip.text = element_text(size = 40))+
  geom_text_repel(data = df_cap_median, 
                  aes(x=med_igr_i, y=med_igr_j, label = label), 
                  size = 8,
                  color="black",
                  max.overlaps = Inf,
                  box.padding = 1)+
  geom_segment(data=df_cap_median %>% filter(flag==T), 
               mapping=aes(x=med_igr_i, y=med_igr_j, 
                           xend=med_igr_i+0.8, yend=med_igr_j),
               linewidth=1.5,
               arrow = arrow(length = unit(0.3, "cm")),
               color="black")+
  facet_wrap(~raintype,labeller = labeller(raintype = c("1" = "Constant Rainfall: Consistent Wet",
                                                        "2" = "Constant Rainfall: Post-flood Dry",
                                                        "3" = "Constant Rainfall: Pre-flood Dry", 
                                                        "4" = "Constant Rainfall: Consistent Dry",
                                                        "5" = "Constant Rainfall: Average",
                                                        "6" = "Variable Rainfall")))

pie_median+theme(plot.margin = margin(0,2,0,1,unit = "cm"),
                 panel.spacing=unit(3,"lines"),
                 legend.position = c(0.07, 0.6))

ggsave("Figure/Coexistence outcome/coexistence-pie_median.pdf",width = 32,height = 22)


df_median%>% 
  group_by(raintype) %>% 
  summarise(CE=mean(CE)) %>% 
  ggplot(aes(x=raintype ,y=CE, fill = raintype))+geom_col()+
  scale_fill_manual(values = c("1" = "lightgrey", "2" = "lightgrey", "3" = "lightgrey", 
                               "4" = "lightgrey", "5" = "lightgrey", "6" = "black")) +
  scale_x_discrete(labels = c("1" = "Consistent Wet", 
                              "2" = "Post-flood Dry", 
                              "3" = "Pre-flood Dry", 
                              "4" = "Consistent Dry", 
                              "5" = "Average", 
                              "6" = "Variable")) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black") +
  ylab("Coexistence probability")+xlab("Rainfall Scenario")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=25),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title =  element_text(size = 22))

ggsave("Figure/Coexistence outcome/coexistence probability comparison.pdf",width=12,height=10)
