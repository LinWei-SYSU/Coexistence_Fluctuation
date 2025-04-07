Pre<-read.csv("Heerkou_Pre.csv")

library(dplyr)
library(tidyverse)

HEK_Y<-Pre %>% 
  group_by(Year) %>% 
  summarise(s=sum(Pre))

avg_Y<-HEK_Y %>% summarise(avg_Y=mean(s))

HEK_M<-Pre %>% 
  group_by(Month) %>% 
  summarise(avg_M=mean(Pre))

###Average Intra-annual Rainfall Pattern

p1<-ggplot(HEK_M, aes(x=Month , y=avg_M/10))  + 
  geom_line()+
  theme_bw() + 
  labs(x="Month", y="Rainfall (mm)",subtitle = "(a)")+
  scale_x_continuous(breaks = seq(1,12,1))+
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#52C9F7",alpha=0.6) +
  annotate("rect", xmin = 6.5, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "#CCF9FF",alpha=0.6) +
  geom_point( size = 3) + 
  annotate("text", x = 5, y = 50, label = "Pre-flood Season", size = 4, fontface = "bold") +
  annotate("text", x = 8, y = 50, label = "Post-flood Season", size = 4, fontface = "bold") +
  theme(text = element_text( size = 7),
        plot.subtitle = element_text(size=12,
                                     hjust = 0,
                                     vjust = 1,
                                     face = "bold"),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

ggsave("Rainfall_month.pdf",path = "Figure/Rain_pattern", width = 12, height = 6.5)


avg_Flood<-HEK_M %>% 
  filter(Month>3&Month<10) %>% 
  summarise(s=sum(avg_M))

Proportion=avg_Flood/avg_Y*100

Range_preflood<-Pre %>% 
  group_by(Year) %>%
  filter(Month>3&Month<7)

Range_proflood<-Pre %>% 
  group_by(Year) %>%
  filter(Month>6&Month<10)

report::report(Range_preflood$Pre)
report::report(Range_proflood$Pre)

###plot

rainsummary<-Pre %>% 
  filter(Month < 10 & Month>3) %>% ###Growth Season April-September
  mutate(Season = "Late", 
         Season = ifelse(Month == 4 | Month == 5 | Month == 6, "Early", Season)) %>% 
  group_by(Year, Season) %>%
  summarize(ppt = sum(Pre)) %>%
  spread(Season, ppt) %>%
  mutate(Total = Early + Late) 

rainsummary<-rainsummary%>%
  mutate(Raintype = case_when(
    Early > mean(rainsummary$Early) & Late > mean(rainsummary$Late) ~ "1",
    Early > mean(rainsummary$Early) & Late < mean(rainsummary$Late)~ "2",
    Early < mean(rainsummary$Early) & Late > mean(rainsummary$Late)~ "3",
    Early < mean(rainsummary$Early) & Late < mean(rainsummary$Late)~ "4"
  ))


###
p2<-ggplot(rainsummary, aes(x=Year, y=Total/10))  + 
  geom_line()+
  geom_point(aes(color = Raintype), size = 3) + 
  theme_bw() + 
  scale_color_manual(values=c('1' ="#040082", '2' = "#CDB8AB",'3' ="#D2B53A", '4' = "#60290C"), name = "Raintype: ", labels = c(" Consistent Wet ", " Post-flood Dry ", " Pre-flood Dry", " Consistent Dry"))+
  labs(x="Year", y="Total rainfall (mm)",subtitle = "(b)")+
  theme(
    text = element_text(size = 7),  
    plot.subtitle = element_text(size = 12,
                                 hjust = 0,
                                 vjust = 1,
                                 face = "bold"),
    axis.title = element_text(size = 13), 
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = c(0.08, 0.85)
  )


ggsave("Raintype_year.pdf",path = "Figure/Rain_pattern", width = 12, height = 6.5)

rainsummary %>% 
  ungroup() %>% 
  count(Raintype) %>% 
  mutate(por=n/120*100)

###combination
library(gridExtra)
g <- grid.arrange(p1,p2) 

ggsave("Rainpattern.pdf",g,path = "Figure/Rain_pattern", width = 10, height = 10)


