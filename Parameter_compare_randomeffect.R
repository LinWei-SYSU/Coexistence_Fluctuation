pars<-read.csv("Output/pars_ricker.csv")
pars_re<-read.csv("Output/pars_ricker_re.csv")

library(dplyr)
library(ggplot2)

library(scales)
show_col(hue_pal()(8))

pars_alpha_1<-pars %>% select(focal,competitor,alpha,alpha_lowCI,alpha_highCI)

pars_alpha_2<-pars_re%>% select(alpha,alpha_lowCI,alpha_highCI)%>%
  rename(
    alpha_re = alpha,
    alpha_lowCI_re = alpha_lowCI,
    alpha_highCI_re = alpha_highCI
  )

compare_alpha<-cbind(pars_alpha_1,pars_alpha_2)

###make figure for each lambda


pars_lambda_1<-pars %>% select(focal,competitor,lambda,lambda_lowCI,lambda_highCI)

pars_lambda_2<-pars_re%>% select(lambda,lambda_lowCI,lambda_highCI)%>%
  rename(
    lambda_re = lambda,
    lambda_lowCI_re = lambda_lowCI,
    lambda_highCI_re = lambda_highCI
  )

compare_lambda<-cbind(pars_lambda_1,pars_lambda_2)


ggplot(compare_lambda,aes(lambda,lambda_re,colour = focal))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for lambda_re
  geom_errorbar(aes(ymin = lambda_lowCI_re, ymax = lambda_highCI_re), width = 0.1,alpha=0.02) +
  # Add horizontal error bars for lambda
  geom_errorbarh(aes(xmin = lambda_lowCI, xmax = lambda_highCI), height = 0.1,alpha=0.02)+
  labs(x = "Intrinsic growth rate", y ="Intrinsic growth rate random effect",subtitle = "(a)")+
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
    legend.position = c(.86, .20),
    legend.background = element_rect(fill = NA),
    legend.title = element_blank()
  ) 
  

ggsave("Figure/Parameter_comparison/compare_lambda.png",width=6,height = 6)


##AGCO

ggplot(compare_alpha %>% filter(focal=="AGCO"),aes(alpha,alpha_re,colour = focal))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2) +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2)+
  labs(x = "Interaction coefficients (i: AGCO)", y ="Interaction coefficients (i: AGCO) ramdom effect",subtitle = "(b)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 
  

ggsave("Figure/Parameter_comparison/compare_alpha_AGCO.png",width=6,height = 6)

##BIPI

ggplot(compare_alpha %>% filter(focal=="BIPI"),aes(alpha,alpha_re))+
  geom_point(colour = '#CD9600')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#CD9600') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#CD9600')+
  labs(x = "Interaction coefficients (i: BIPI)", y ="Interaction coefficients (i: BIPI) ramdom effect",subtitle = "(c)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_BIPI.png",width=6,height = 6)


##DAAE

ggplot(compare_alpha %>% filter(focal=="DAAE"),aes(alpha,alpha_re))+
  geom_point(colour = '#7CAE00')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#7CAE00') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#7CAE00')+
  labs(x = "Interaction coefficients (i: DAAE)", y ="Interaction coefficients (i: DAAE) ramdom effect",subtitle = "(d)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_DAAE.png",width=6,height = 6)


##LUHY

ggplot(compare_alpha %>% filter(focal=="LUHY"),aes(alpha,alpha_re))+
  geom_point(colour = '#00BE67')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#00BE67') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#00BE67')+
  labs(x = "Interaction coefficients (i: LUHY)", y ="Interaction coefficients (i: LUHY) ramdom effect",subtitle = "(f)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_LUHY.png",width=6,height = 6)

##OPCO

ggplot(compare_alpha %>% filter(focal=="OPCO"),aes(alpha,alpha_re))+
  geom_point(colour = '#00BFC4')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#00BFC4') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#00BFC4')+
  labs(x = "Interaction coefficients (i: OPCO)", y ="Interaction coefficients (i: OPCO) ramdom effect",subtitle = "(e)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_OPCO.png",width=6,height = 6)

##PHAN

ggplot(compare_alpha %>% filter(focal=="PHAN"),aes(alpha,alpha_re))+
  geom_point(colour = '#00A9FF')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#00A9FF') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#00A9FF')+
  labs(x = "Interaction coefficients (i: PHAN)", y ="Interaction coefficients (i: PHAN) ramdom effect",subtitle = "(g)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_PHAN.png",width=6,height = 6)


##POHY

ggplot(compare_alpha %>% filter(focal=="POHY"),aes(alpha,alpha_re))+
  geom_point(colour = '#C77CFF')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#C77CFF') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#C77CFF')+
  labs(x = "Interaction coefficients (i: POHY)", y ="Interaction coefficients (i: POHY) ramdom effect",subtitle = "(h)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_POHY.png",width=6,height = 6)


##PRCL

ggplot(compare_alpha %>% filter(focal=="PRCL"),aes(alpha,alpha_re))+
  geom_point(colour = '#FF61CC')+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black")+
  # Add vertical error bars for alpha_re
  geom_errorbar(aes(ymin = alpha_lowCI_re, ymax = alpha_highCI_re), width = 0.1,alpha=0.2,colour = '#FF61CC') +
  # Add horizontal error bars for alpha
  geom_errorbarh(aes(xmin = alpha_lowCI, xmax = alpha_highCI), height = 0.1,alpha=0.2,colour = '#FF61CC')+
  labs(x = "Interaction coefficients (i: PRCL)", y ="Interaction coefficients (i: PRCL) ramdom effect",subtitle = "(i)")+
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
    legend.position = "none",
    legend.background = element_rect(fill = NA),
  ) 


ggsave("Figure/Parameter_comparison/compare_alpha_PRCL.png",width=6,height = 6)


