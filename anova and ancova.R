library(tidyverse)
library(car)
library(sm)
mwwd <- read.csv("data/trial_mww.csv",stringsAsFactors = F)

#test of bmi
qqp(mwwd$bmi_d)
shapiro.test(mwwd$bmi_d)#正态分布检验
bartlett.test(bmi_d~group,data=mwwd)#方差齐性检验
kruskal.test(bmi_d~group, data = mwwd)#KS检验，p=0.19
sm.bmi <- sm.ancova(x=mwwd$gender,y=mwwd$bmi_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.09

#test of body fat
shapiro.test(mwwd$bodyfat_d)#正态分布检验,p<0.05
bartlett.test(bodyfat_d~group,data=mwwd)#方差齐性检验,p<0.05,组间方差不齐用非参数检验
kruskal.test(bodyfat_d~group, data =   mwwd)#KS检验，p=0.42
sm.ancova(x=mwwd$gender,y=mwwd$bodyfat_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.66

#test of liverfat_ave
qqp(mwwd$liverfat_ave_d)
shapiro.test(mwwd$liverfat_ave_d)#正态分布检验,p=0.24
bartlett.test(liverfat_ave_d~group,data=mwwd)#方差齐性检验,p>0.05
fitlf <- aov(liverfat_ave_d~group,data = mwwd)#ANOVA
summary(fitlf)#p=0.46
fitlf_c <- aov(liverfat_ave_d~group+age+gender,data = mwwd)#ancova
summary(fitlf_c)

#test of liverfat_2
qqp(mwwd$liverfat_2_d)
shapiro.test(mwwd$liverfat_2_d)#正态分布检验,p>0.05
bartlett.test(liverfat_2_d~group,data=mwwd)#方差齐性检验,p<0.05,组间方差不齐用非参数检验
fitlf2 <- aov(liverfat_2_d~group,data = mwwd)#ANOVA
summary(fitlf2)#p=0.34
fitlf2_c <- aov(liverfat_2_d~group+age+gender,data = mwwd)#ancova
summary(fitlf2_c)#p=0.35
kruskal.test(liverfat_2_d~group, data = mwwd)#KS检验，p=0.45
sm.ancova(x=mwwd$gender,y=mwwd$liverfat_2_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.56

#test of liverfat-3
qqp(mwwd$liverfat_3_d)
shapiro.test(mwwd$liverfat_3_d)#正态分布检验,p<0.05
bartlett.test(liverfat_3_d~group,data=mwwd)#方差齐性检验,p<0.05,组间方差不齐用非参数检验
fitlf3 <- aov(liverfat_3_d~group,data = mwwd)#ANOVA
summary(fitlf3)#p=0.96
fitlf3_c <- aov(liverfat_3_d~group+age+gender,data = mwwd)#ancova
summary(fitlf3_c)#p=0.96
kruskal.test(liverfat_3_d~group, data = mwwd)#KS检验，p=0.27
sm.ancova(x=mwwd$gender,y=mwwd$liverfat_3_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.555

#test of liverfat-4
qqp(mwwd$liverfat_4_d)
shapiro.test(mwwd$liverfat_4_d)#正态分布检验,p<0.05
bartlett.test(liverfat_4_d~group,data=mwwd)#方差齐性检验,p<0.05,组间方差不齐用非参数检验
fitlf4 <- aov(liverfat_4_d~group,data = mwwd)#ANOVA
summary(fitlf4)#p=0.111
fitlf4_c <- aov(liverfat_4_d~group+age+gender,data = mwwd)#ancova
summary(fitlf4_c)#p=0.112
kruskal.test(liverfat_4_d~group, data = mwwd)#KS检验，p=0.54
sm.ancova(x=mwwd$gender,y=mwwd$liverfat_4_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.44

#test of liverfat-5
qqp(mwwd$liverfat_5_d)
shapiro.test(mwwd$liverfat_5_d)#正态分布检验,p>0.05
bartlett.test(liverfat_5_d~group,data=mwwd)#方差齐性检验,p>0.05
fitlf5 <- aov(liverfat_5_d~group,data = mwwd)#ANOVA
summary(fitlf5)#p=0.32
fitlf5_c <- aov(liverfat_5_d~group+age+gender,data = mwwd)#ancova
summary(fitlf5_c)#p=0.30


#test of liverfat-6
qqp(mwwd$liverfat_6_d)
shapiro.test(mwwd$liverfat_6_d)#正态分布检验,p>0.05
bartlett.test(liverfat_6_d~group,data=mwwd)#方差齐性检验,p>0.05
fitlf6 <- aov(liverfat_6_d~group,data = mwwd)#ANOVA
summary(fitlf6)#p=0.006
fitlf6_c <- aov(liverfat_6_d~group+age+gender,data = mwwd)#ancova
summary(fitlf6_c)#p=0.006

#test of liverfat-7
qqp(mwwd$liverfat_7_d)
shapiro.test(mwwd$liverfat_7_d)#正态分布检验,p>0.05
bartlett.test(liverfat_7_d~group,data=mwwd)#方差齐性检验,p>0.05
fitlf7 <- aov(liverfat_7_d~group,data = mwwd)#ANOVA
summary(fitlf7)#p=0.124
fitlf7_c <- aov(liverfat_7_d~group+age+gender,data = mwwd)#ancova
summary(fitlf7_c)#p=0.128

#test of liverfat-8
qqp(mwwd$liverfat_8_d)
shapiro.test(mwwd$liverfat_8_d)#正态分布检验,p>0.05
bartlett.test(liverfat_8_d~group,data=mwwd)#方差齐性检验,p>0.05
fitlf8 <- aov(liverfat_8_d~group,data = mwwd)#ANOVA
summary(fitlf8)#p=0.322
fitlf8_c <- aov(liverfat_8_d~group+age+gender,data = mwwd)#ancova
summary(fitlf8_c)#p=0.326

#test of score
qqp(mwwd$score_d)
shapiro.test(mwwd$score_d)#正态分布检验,p<0.05
bartlett.test(score_d~group,data=mwwd)#方差齐性检验,p<0.05
kruskal.test(score_d~group, data = mwwd)#KS检验，p=0.83
sm.ancova(x=mwwd$gender,y=mwwd$score_d,group = mwwd$group,model = "equal")#非参数ANCOVA,p=0.86

mwwpe <- read.csv("data/mww_pe.csv",stringsAsFactors = F)
mwwpe81 <- read.csv("data/mww_pe81.csv",stringsAsFactors = F)

ggplot(data=mwwpe81,mapping = aes(x=time,y=bmi,color=as.factor(group)))+
  geom_smooth(method = lm)

ggplot(data=mwwpe81,mapping = aes(x=time,y=whr,color=as.factor(group)))+
  geom_smooth(method = lm)

ggplot(data=mwwpe81,mapping = aes(x=time,y=wc,color=as.factor(group)))+
  geom_smooth(method = lm)

ggplot(data=mwwpe81,mapping = aes(x=time,y=systo,color=as.factor(group)))+
  geom_smooth(method = lm)

ggplot(data=mwwpe81,mapping = aes(x=time,y=diasto,color=as.factor(group)))+
  geom_smooth(method = lm)

ggplot(data=mwwpe81,mapping = aes(x=time,y=wc))+
  geom_smooth(aes(color=as.factor(group)),method = "lm",se=FALSE)+
  geom_errorbar()

#mean of BMI
by.group <- dplyr::group_by(mwwpe81,group,time)
bmi.mean <- summarize(by.group,bmi.mean=mean(bmi,na.rm=T),bmi.sd=sd(bmi,na.rm=T))

bmi.group <- ggplot(data=bmi.mean,mapping = aes(x=time,y=bmi.mean))+
  geom_point(aes(shape=factor(group)),size=10)+
  geom_path(aes(lineend=as.factor(group)),size=1.5)+
  geom_errorbar(aes(ymin=bmi.mean-bmi.sd,ymax=bmi.mean+bmi.sd,width=0.4))+
  scale_x_continuous("时间(周)",breaks = c(0,4,8,12,16),minor_breaks = NULL)+
  scale_y_continuous(quote(BMI (Kg/m^2)))+labs(shape="处理",title="不同处理五次体检BMI变化趋势图")

bmi.style <- bmi.group+theme_light()+
  theme(plot.title = element_text(size=50,hjust=0.5,margin = margin(t=20,b=40)),
        axis.title = element_text(size=40),
        axis.title.x=element_text(margin = margin(t=40)),
        axis.title.y = element_text(margin = margin(r=40)),
        axis.line = element_line(color = "grey50",size = 3),
        axis.text = element_text(size=40),
        legend.text = element_text(size=40),
        legend.title = element_text(size=40))

#mean of WC

wc.mean <- summarize(by.group,wc.mean=mean(wc,na.rm=T),wc.sd=sd(wc,na.rm=T))

wc.group <- ggplot(data=wc.mean,mapping = aes(x=time,y=wc.mean))+
  geom_point(aes(shape=factor(group)),size=10)+
  geom_path(aes(lineend=as.factor(group)),size=1.5)+
  geom_errorbar(aes(ymin=wc.mean-wc.sd,ymax=wc.mean+wc.sd,width=0.4))+
  scale_x_continuous("时间(周)",breaks = c(0,4,8,12,16),minor_breaks = NULL)+
  scale_y_continuous("腰围(cm)")+labs(shape="处理",title="不同处理五次体检腰围变化趋势图")

wc.style <- wc.group+theme_light()+
  theme(plot.title = element_text(size=50,hjust=0.5,margin = margin(t=20,b=40)),
        axis.title = element_text(size=40),
        axis.title.x=element_text(margin = margin(t=40)),
        axis.title.y = element_text(margin = margin(r=40)),
        axis.line = element_line(color = "grey50",size = 3),
        axis.text = element_text(size=40),
        legend.text = element_text(size=40),
        legend.title = element_text(size=40))

#mean of whr

whr.mean <- summarize(by.group,whr.mean=mean(whr,na.rm=T),whr.sd=sd(whr,na.rm=T))

whr.group <- ggplot(data=whr.mean,mapping = aes(x=time,y=whr.mean))+
  geom_point(aes(shape=factor(group)),size=10)+
  geom_path(aes(lineend=as.factor(group)),size=1.5)+
  geom_errorbar(aes(ymin=whr.mean-whr.sd,ymax=whr.mean+whr.sd,width=0.4))+
  scale_x_continuous("时间(周)",breaks = c(0,4,8,12,16),minor_breaks = NULL)+
  scale_y_continuous("腰臀比")+labs(shape="处理",title="不同处理五次体检腰臀比变化趋势图")

whr.style <- whr.group+theme_light()+
  theme(plot.title = element_text(size=50,hjust=0.5,margin = margin(t=20,b=40)),
        axis.title = element_text(size=40),
        axis.title.x=element_text(margin = margin(t=40)),
        axis.title.y = element_text(margin = margin(r=40)),
        axis.line = element_line(color = "grey50",size = 3),
        axis.text = element_text(size=40),
        legend.text = element_text(size=40),
        legend.title = element_text(size=40))

#mean of systo

systo.mean <- summarize(by.group,systo.mean=mean(systo,na.rm=T),systo.sd=sd(systo,na.rm=T))

systo.group <- ggplot(data=systo.mean,mapping = aes(x=time,y=systo.mean))+
  geom_point(aes(shape=factor(group)),size=10)+
  geom_path(aes(lineend=as.factor(group)),size=1.5)+
  geom_errorbar(aes(ymin=systo.mean-systo.sd,ymax=systo.mean+systo.sd,width=0.4))+
  scale_x_continuous("时间(周)",breaks = c(0,4,8,12,16),minor_breaks = NULL)+
  scale_y_continuous("收缩压 (mmHg)")+labs(shape="处理",title="不同处理五次体检收缩压变化趋势图")

systo.style <- systo.group+theme_light()+
  theme(plot.title = element_text(size=50,hjust=0.5,margin = margin(t=20,b=40)),
        axis.title = element_text(size=40),
        axis.title.x=element_text(margin = margin(t=40)),
        axis.title.y = element_text(margin = margin(r=40)),
        axis.line = element_line(color = "grey50",size = 3),
        axis.text = element_text(size=40),
        legend.text = element_text(size=40),
        legend.title = element_text(size=40))

#mean of diasto

diasto.mean <- summarize(by.group,diasto.mean=mean(diasto,na.rm=T),diasto.sd=sd(diasto,na.rm=T))

diasto.group <- ggplot(data=diasto.mean,mapping = aes(x=time,y=diasto.mean))+
  geom_point(aes(shape=factor(group)),size=10)+
  geom_path(aes(lineend=as.factor(group)),size=1.5)+
  geom_errorbar(aes(ymin=diasto.mean-diasto.sd,ymax=diasto.mean+diasto.sd,width=0.4))+
  scale_x_continuous("时间(周)",breaks = c(0,4,8,12,16),minor_breaks = NULL)+
  scale_y_continuous("舒张压 (mmHg)")+labs(shape="处理",title="不同处理五次体检舒张压变化趋势图")

diasto.style <- diasto.group+theme_light()+
  theme(plot.title = element_text(size=50,hjust=0.5,margin = margin(t=20,b=40)),
        axis.title = element_text(size=40),
        axis.title.x=element_text(margin = margin(t=40)),
        axis.title.y = element_text(margin = margin(r=40)),
        axis.line = element_line(color = "grey50",size = 3),
        axis.text = element_text(size=40),
        legend.text = element_text(size=40),
        legend.title = element_text(size=40))
