#####20200513##########
##########决定不做中介模型，做常规模型########
library(tidyverse)
mww_fa.p <- read.csv("mww_FAonly_p.csv",stringsAsFactors = F)#FAonly 仅包含脂肪酸数据
mww_fa_NA.p <- dplyr::filter(mww_fa.p,!is.na(C140)& !is.na(age))#删去没有做脂肪酸检测的观测n=6，以及年龄性别有缺失的观测n=1
mww_fa_NAp2 <- mww_fa_NA.p[,which(colMeans(mww_fa_NA.p!=0)>=0.5)]#删除0值超过50%的脂肪
mww_fa_NAp2$scd16 <- mww_fa_NAp2$C161/mww_fa_NAp2$C160
mww_fa_NAp2$scd18 <- mww_fa_NAp2$C181n9c/mww_fa_NAp2$C180
mww_fa_NAp2$d6d <- mww_fa_NAp2$C183n6/mww_fa_NAp2$C182n6c
mww_fa_NAp2$d5d <- mww_fa_NAp2$C204n6/mww_fa_NAp2$C203n6

#导入认知和炎症数据
mww_ci <- read.csv("mww_cog_imf.csv",stringsAsFactors = F)
mww <- left_join(mww_fa_NAp2,mww_ci,by="id")

#提取脂肪酸数据，做PCA
id <- mww$id
fa_scale <- scale(mww[,4:17],center = F)
fa_scale <- as.data.frame(fa_scale)
fa <- mww[,4:17]
fa_prin <- princomp(fa,scores = T,cor = T)
fascale_prin <- princomp(fa_scale,scores = T,cor = T)

mww_pca <- as.data.frame(fascale_prin$scores)
mww_pca$id <- id
mww_pca <- dplyr::select(mww_pca,id,Comp.1,Comp.2,Comp.3,Comp.4)#用4个pca
mww_pca <- left_join(mww_pca,mww,by="id")
#四个主成分得分取四分位
mww_pca$c1q <- ntile(mww_pca$Comp.1,4)
mww_pca$c2q <- ntile(mww_pca$Comp.2,4)
mww_pca$c3q <- ntile(mww_pca$Comp.3,4)
mww_pca$c4q <- ntile(mww_pca$Comp.4,4)

##加入膳食及其他混杂因子数据
mww_diet <- read.csv("mww_dietary.csv",stringsAsFactors = F)
mww_pca2 <- mww_pca
mww_pca <- left_join(mww_pca2,mww_diet,by="id")
write.csv(mww_pca,file = "mww_pca.csv",row.names = F)
mww_pca <- read.csv(file = "mww_pca.csv",stringsAsFactors = F)

####linear regression
######MMSE
lm_comp1.mmse <- lm(mmse_total~Comp.1+age+gender,data = mww_pca)
summary(lm_comp1.mmse)

lm_comp2.mmse <- lm(mmse_total~Comp.2+age+gender,data = mww_pca)
summary(lm_comp2.mmse)
lm_comp2.mmse2 <- lm(mmse_total~Comp.2+age+gender+bmi+wc,data = mww_pca)
summary(lm_comp2.mmse2)
lm_comp2.mmse3 <- lm(mmse_total~Comp.2+age+gender+bmi+wc+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmse3)
lm_comp2.mmse4 <- lm(mmse_total~Comp.2+age+gender+bmi+wc+smoke,data = mww_pca)
summary(lm_comp2.mmse4)
lm_comp2.mmse5 <- lm(mmse_total~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmse5)

lm_comp1.mmse5 <- lm(mmse_total~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mmse5)

lm_comp3.mmse5 <- lm(mmse_total~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mmse5)

lm_comp4.mmse5 <- lm(mmse_total~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mmse5)

##mmse Orientation
lm_comp1.mmseOri <- lm(MMSEOrientation~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mmseOri)
lm_comp2.mmseOri <- lm(MMSEOrientation~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmseOri)
lm_comp3.mmseOri <- lm(MMSEOrientation~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mmseOri)
lm_comp4.mmseOri <- lm(MMSEOrientation~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mmseOri)

##MMSE Computation
lm_comp1.mmsecomp <- lm(MMSEComputation~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mmsecomp)
lm_comp2.mmsecomp <- lm(MMSEComputation~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmsecomp)
lm_comp3.mmsecomp <- lm(MMSEComputation~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mmsecomp)
lm_comp4.mmsecomp <- lm(MMSEComputation~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mmsecomp)

##MMSE Memory
lm_comp1.mmsemem <- lm(MMSEMemory~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mmsemem)
lm_comp2.mmsemem <- lm(MMSEMemory~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmsemem)
lm_comp3.mmsemem <- lm(MMSEMemory~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mmsemem)
lm_comp4.mmsemem <- lm(MMSEMemory~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mmsemem)

##MMSE Languageskill
lm_comp1.mmselan <- lm(MMSELanguageskill~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mmselan)
lm_comp2.mmselan <- lm(MMSELanguageskill~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mmselan)
lm_comp3.mmselan <- lm(MMSELanguageskill~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mmselan)
lm_comp4.mmselan <- lm(MMSELanguageskill~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mmselan)

##MOCA
lm_comp1.moca <- lm(moca_total~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.moca)
lm_comp2.moca <- lm(moca_total~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.moca)
lm_comp3.moca <- lm(moca_total~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.moca)
lm_comp4.moca <- lm(moca_total~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.moca)

##MOCA naming
lm_comp1.mocaname <- lm(MoCANaming~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocaname)
lm_comp2.mocaname <- lm(MoCANaming~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocaname)
lm_comp3.mocaname <- lm(MoCANaming~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocaname)
lm_comp4.mocaname <- lm(MoCANaming~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocaname)

##MoCA Orientation
lm_comp1.mocao <- lm(MoCAOrientation~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocao)
lm_comp2.mocao <- lm(MoCAOrientation~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocao)
lm_comp3.mocao <- lm(MoCAOrientation~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocao)
lm_comp4.mocao <- lm(MoCAOrientation~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocao)

##MoCA Delayedrecall
lm_comp1.mocade <- lm(MoCADelayedrecall~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocade)
lm_comp2.mocade <- lm(MoCADelayedrecall~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocade)
lm_comp3.mocade <- lm(MoCADelayedrecall~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocade)
lm_comp4.mocade <- lm(MoCADelayedrecall~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocade)

##MoCA Abstractthinking
lm_comp1.mocab <- lm(MoCAAbstractthinking~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocab)
lm_comp2.mocab <- lm(MoCAAbstractthinking~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocab)
lm_comp3.mocab <- lm(MoCAAbstractthinking~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocab)
lm_comp4.mocab <- lm(MoCAAbstractthinking~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocab)

##MoCA Languageskills
lm_comp1.mocal <- lm(MoCALanguageskills~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocal)
lm_comp2.mocal <- lm(MoCALanguageskills~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocal)
lm_comp3.mocal <- lm(MoCALanguageskills~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocal)
lm_comp4.mocal <- lm(MoCALanguageskills~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocal)

##MoCA Visualspatialability
lm_comp1.mocalv <- lm(MoCAVisualspatialability~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocalv)
lm_comp2.mocalv <- lm(MoCAVisualspatialability~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocalv)
lm_comp3.mocalv <- lm(MoCAVisualspatialability~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocalv)
lm_comp4.mocalv <- lm(MoCAVisualspatialability~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocalv)

##MoCA Attention
lm_comp1.mocalva <- lm(MoCAAttention~Comp.1+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp1.mocalva)
lm_comp2.mocalva <- lm(MoCAAttention~Comp.2+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp2.mocalva)
lm_comp3.mocalva <- lm(MoCAAttention~Comp.3+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp3.mocalva)
lm_comp4.mocalva <- lm(MoCAAttention~Comp.4+age+gender+bmi+wc+smoke+energy+crp+NFkB,data = mww_pca)
summary(lm_comp4.mocalva)

hist(mww_pca$mmse_total)
hist(mww_pca$moca_total)

#四分位间mmse均值差异
tapply(mww_pca$mmse_total, mww_pca$c1q, mean,na.rm=T)
tapply(mww_pca$mmse_total, mww_pca$c2q, mean,na.rm=T)
tapply(mww_pca$mmse_total, mww_pca$c3q, mean,na.rm=T)
tapply(mww_pca$mmse_total, mww_pca$c4q, mean,na.rm=T)
#四分位间mmse的p
l.mmse1 <- lm(mmse_total~c1q+gender+age,data = mww_pca)
summary(l.mmse1)
l.mmse2 <- lm(mmse_total~c2q+gender+age,data = mww_pca)
summary(l.mmse2)
l.mmse3 <- lm(mmse_total~c3q+gender+age,data = mww_pca)
summary(l.mmse3)
l.mmse4 <- lm(mmse_total~c4q+gender+age,data = mww_pca)
summary(l.mmse4)

#四分位间moca均值差异
tapply(mww_pca$moca_total, mww_pca$c1q, mean,na.rm=T)
tapply(mww_pca$moca_total, mww_pca$c2q, mean,na.rm=T)
tapply(mww_pca$moca_total, mww_pca$c3q, mean,na.rm=T)
tapply(mww_pca$moca_total, mww_pca$c4q, mean,na.rm=T)
#四分位间moca的p
l.moca1 <- lm(moca_total~c1q+gender+age,data = mww_pca)
summary(l.moca1)
l.mmse2 <- lm(moca_total~c2q+gender+age,data = mww_pca)
summary(l.mmse2)
l.mmse3 <- lm(moca_total~c3q+gender+age,data = mww_pca)
summary(l.mmse3)
l.mmse4 <- lm(moca_total~c4q+gender+age,data = mww_pca)
summary(l.mmse4)

l.lasse <- lm(mmse_total~C182n6c+age+gender,data=mww_pca)
summary(l.lasse)
l.18sse <- lm(mmse_total~C180+age+gender,data=mww_pca)
summary(l.18sse)
l.16sse <- lm(mmse_total~C160+age+gender,data=mww_pca)
summary(l.16sse)

