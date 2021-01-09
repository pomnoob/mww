###粗略看一下跟认知的关系
library(tidyverse)
dbi_score <- read.csv(file = "DBI/DBI膳食质量分数.csv",stringsAsFactors = F)
FA <- read.csv(file = "20200108 数据库(放了血液脂肪酸-重新核对脂肪酸正确版-去掉359兰孝明重复的+计算了认知的小分).csv",stringsAsFactors = F)
FA2 <- FA %>% select(登记号,12:15,年龄,性别,体重指数,腰围,腰臀比,教育,mmse总分,MoCA总分,吸烟史,饮酒史,内毒素,超敏C,能量,
                       MMSEOrientation,MMSEComputation,MMSEMemory,MMSELanguageskill,MoCANaming,MoCAOrientation,
                       MoCADelayedrecall,MoCAAbstractthinking,MoCALanguageskills,MoCAVisualspatialability,MoCAAttention) %>%
  rename(id=登记号,crp=超敏C,lps=内毒素,age=年龄,gender=性别,bmi=体重指数,wc=腰围,wh=腰臀比,
         edu=教育,mmse=mmse总分,moca=MoCA总分,energy=能量,smoke=吸烟史,drink=饮酒史)
                       

mww_dbi <- left_join(FA2,dbi_score,by="id")
fa_p <- read.csv(file = "mww_FAonly_p.csv",stringsAsFactors = F)%>%
  select(-age,-gender)

mww_fa_NA.p <- dplyr::filter(fa_p,!is.na(C140))#删去没有做脂肪酸检测的观测n=6
mww_fa_NAp2 <- mww_fa_NA.p[,which(colMeans(mww_fa_NA.p!=0)>=0.5)]#删除0值超过50%的脂肪酸
mww_dbi <- left_join(mww_dbi,mww_fa_NAp2,by="id")
#可以看看pca的结果，实际上是不显著的，如果要加，则84行需要改一下
mww_pca <- read.csv(file="mww_pca.csv",stringsAsFactors = F) %>% select(id,Comp.1,Comp.2,Comp.3,Comp.4)
mww_dbi <- left_join(mww_dbi,mww_pca,by="id")

#TNFa
library(Hmisc)
describe(mww_dbi$TNFa)#最大值280？
boxplot(mww_dbi$TNFa)
mww_dbi$tnf2 <- mww_dbi$TNFa
mww_dbi$tnf2[mww_dbi$tnf2>280] <- NA
describe(mww_dbi$tnf2)
boxplot(mww_dbi$tnf2)
hist(mww_dbi$tnf2)
hist(sqrt(mww_dbi$tnf2))
mww_dbi$tnf3 <- sqrt(mww_dbi$tnf2)
boxplot(mww_dbi$tnf3)
hist(mww_dbi$tnf3)

#bmi
describe(mww_dbi$bmi)
hist(mww_dbi$bmi)

#IL1β
describe(mww_dbi$IL1β)
hist(mww_dbi$IL1β)
mww_dbi$IL1β2 <- sqrt(mww_dbi$IL1β)
hist(mww_dbi$IL1β2)

#IL10
describe(mww_dbi$IL10)
hist(mww_dbi$IL10)
mww_dbi$IL102 <- sqrt(mww_dbi$IL10)
hist(mww_dbi$IL102)

#Nfkb
describe(mww_dbi$NFkB)
hist(mww_dbi$NFkB)
mww_dbi$NFkB2 <- sqrt(mww_dbi$NFkB)
hist(mww_dbi$NFkB2)

#lps
describe(mww_dbi$lps)
mww_dbi$lps[mww_dbi$lps<0] <- NA
hist(mww_dbi$lps)

#crp
describe(mww_dbi$crp)
hist(mww_dbi$crp)
mww_dbi$crp2 <- mww_dbi$crp
mww_dbi$crp2[mww_dbi$crp2>420] <- NA
hist(mww_dbi$crp2)
mww_dbi$crp3 <- sqrt(mww_dbi$crp2)
hist(mww_dbi$crp3)
describe(mww_dbi$crp3)
#C150
hist(mww_dbi$C150)
describe(mww_dbi$C150)

#C170
describe(mww_dbi$C170)
hist(mww_dbi$C170)

mww_cov <- mww_dbi %>% 
  dplyr::select(age,gender,bmi,energy,smoke,edu)#混杂因子
mww_fa <- mww_dbi%>%
  dplyr::select(c(42:55))
mww_cog <- mww_dbi%>%
  dplyr::select(mmse,moca,c(19:29))
mww_inf <- mww_dbi %>%
  dplyr::select(IL1β2,IL102,NFkB2,tnf3,crp3,lps)
mww_dq <- mww_dbi %>%
  dplyr::select(c(31:41))
mww_ad <- mww_dbi %>%
  select(bmi,wc,wh)
mww_cov2 <- mww_dbi %>% 
  dplyr::select(age,gender,energy,smoke,edu)#混杂因子2,for mww_ad
###################重新计算partial correlation#################################

faColname <- colnames(mww_fa)
cogColname <- colnames(mww_cog)
covColname <- colnames(mww_cov)
dqColname <- colnames(mww_dq)
infColname <- colnames(mww_inf)
adColname <- colnames(mww_ad)
cov2Colname <- colnames(mww_cov2)
###adiposity and cognition
RCol.ac <- list()
PCol.ac <- list()
for (i in 1:length(adColname)) {
  
  y.ac <- adColname[i] 
  crntRcor.ac <- double()
  crntPcor.ac <- double()
  
  for (j in 1:length(cogColname)) {
    
    x.ac <- cogColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.ac <- paste(y.ac,"~")
    xhs.ac <- paste(x.ac,"~")
    
    frma1.ac <- as.formula(paste(yhs.ac,cov_f2))
    mod1.ac <- lm(frma1.ac,data = mww_dbi,na.action = na.exclude)
    r1.ac <- resid(mod1.ac)#residuls for mod1
    
    frma2.ac <- as.formula(paste(xhs.ac,cov_f2))
    mod2.ac <- lm(frma2.ac,data = mww_dbi,na.action = na.exclude)
    r2.ac <- resid(mod2.ac)
    
    rc.ac <- rcorr(r1.ac,r2.ac)
    
    crntRcor.ac[j] <- rc.ac$r[1,2]
    crntPcor.ac[j] <- rc.ac$P[1,2]
    
  }
  
  RCol.ac[[i]] <- crntRcor.ac
  PCol.ac[[i]] <- crntPcor.ac
  
}

RCorMat.ac <- matrix(unlist(RCol.ac), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, adColname))#R值矩阵

PCorMat.ac <- matrix(unlist(PCol.ac), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, adColname))#P值矩阵

###adiposity and inflammation
RCol.ai <- list()
PCol.ai <- list()
for (i in 1:length(adColname)) {
  
  y.ai <- adColname[i] 
  crntRcor.ai <- double()
  crntPcor.ai <- double()
  
  for (j in 1:length(infColname)) {
    
    x.ai <- infColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.ai <- paste(y.ai,"~")
    xhs.ai <- paste(x.ai,"~")
    
    frma1.ai <- as.formula(paste(yhs.ai,cov_f2))
    mod1.ai <- lm(frma1.ai,data = mww_dbi,na.action = na.exclude)
    r1.ai <- resid(mod1.ai)#residuls for mod1
    
    frma2.ai <- as.formula(paste(xhs.ai,cov_f2))
    mod2.ai <- lm(frma2.ai,data = mww_dbi,na.action = na.exclude)
    r2.ai <- resid(mod2.ai)
    
    rc.ai <- rcorr(r1.ai,r2.ai)
    
    crntRcor.ai[j] <- rc.ai$r[1,2]
    crntPcor.ai[j] <- rc.ai$P[1,2]
    
  }
  
  RCol.ai[[i]] <- crntRcor.ai
  PCol.ai[[i]] <- crntPcor.ai
  
}

RCorMat.ai <- matrix(unlist(RCol.ai), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, adColname))#R值矩阵

PCorMat.ai <- matrix(unlist(PCol.ai), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, adColname))#P值矩阵

###adiposity and dbi
RCol.ad <- list()
PCol.ad <- list()
for (i in 1:length(adColname)) {
  
  y.ad <- adColname[i] 
  crntRcor.ad <- double()
  crntPcor.ad <- double()
  
  for (j in 1:length(dqColname)) {
    
    x.ad <- dqColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.ad <- paste(y.ad,"~")
    xhs.ad <- paste(x.ad,"~")
    
    frma1.ad <- as.formula(paste(yhs.ad,cov_f2))
    mod1.ad <- lm(frma1.ad,data = mww_dbi,na.action = na.exclude)
    r1.ad <- resid(mod1.ad)#residuls for mod1
    
    frma2.ad <- as.formula(paste(xhs.ad,cov_f2))
    mod2.ad <- lm(frma2.ad,data = mww_dbi,na.action = na.exclude)
    r2.ad <- resid(mod2.ad)
    
    rc.ad <- rcorr(r1.ad,r2.ad)
    
    crntRcor.ad[j] <- rc.ad$r[1,2]
    crntPcor.ad[j] <- rc.ad$P[1,2]
    
  }
  
  RCol.ad[[i]] <- crntRcor.ad
  PCol.ad[[i]] <- crntPcor.ad
  
}

RCorMat.ad <- matrix(unlist(RCol.ad), 
                     nrow = length(dqColname), 
                     dimnames = list(dqColname, adColname))#R值矩阵

PCorMat.ad <- matrix(unlist(PCol.ad), 
                     nrow = length(dqColname), 
                     dimnames = list(dqColname, adColname))#P值矩阵

###dbi与炎症因子
RCol.di <- list()
PCol.di <- list()
for (i in 1:length(dqColname)) {
  
  y.di <- dqColname[i] 
  crntRcor.di <- double()
  crntPcor.di <- double()
  
  for (j in 1:length(infColname)) {
    
    x.di <- infColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.di <- paste(y.di,"~")
    xhs.di <- paste(x.di,"~")
    
    frma1.di <- as.formula(paste(yhs.di,cov_f))
    mod1.di <- lm(frma1.di,data = mww_dbi,na.action = na.exclude)
    r1.di <- resid(mod1.di)#residuls for mod1
    
    frma2.di <- as.formula(paste(xhs.di,cov_f))
    mod2.di <- lm(frma2.di,data = mww_dbi,na.action = na.exclude)
    r2.di <- resid(mod2.di)
    
    rc.di <- rcorr(r1.di,r2.di)
    
    crntRcor.di[j] <- rc.di$r[1,2]
    crntPcor.di[j] <- rc.di$P[1,2]
    
  }
  
  RCol.di[[i]] <- crntRcor.di
  PCol.di[[i]] <- crntPcor.di
  
}

RCorMat.di <- matrix(unlist(RCol.di), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, dqColname))#R值矩阵

PCorMat.di <- matrix(unlist(PCol.di), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, dqColname))#P值矩阵

###dbi与认知
RCol.dc <- list()#fi means fatty acid and inflammation
PCol.dc <- list()
for (i in 1:length(dqColname)) {
  
  y.dc <- dqColname[i] 
  crntRcor.dc <- double()
  crntPcor.dc <- double()
  
  for (j in 1:length(cogColname)) {
    
    x.dc <- cogColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.dc <- paste(y.dc,"~")
    xhs.dc <- paste(x.dc,"~")
    
    frma1.dc <- as.formula(paste(yhs.dc,cov_f))
    mod1.dc <- lm(frma1.dc,data = mww_dbi,na.action = na.exclude)
    r1.dc <- resid(mod1.dc)#residuls for mod1
    
    frma2.dc <- as.formula(paste(xhs.dc,cov_f))
    mod2.dc <- lm(frma2.dc,data = mww_dbi,na.action = na.exclude)
    r2.dc <- resid(mod2.dc)
    
    rc.dc <- rcorr(r1.dc,r2.dc)
    
    crntRcor.dc[j] <- rc.dc$r[1,2]
    crntPcor.dc[j] <- rc.dc$P[1,2]
    
  }
  
  RCol.dc[[i]] <- crntRcor.dc
  PCol.dc[[i]] <- crntPcor.dc
  
}

RCorMat.dc <- matrix(unlist(RCol.dc), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, dqColname))#R值矩阵

PCorMat.dc <- matrix(unlist(PCol.dc), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, dqColname))#P值矩阵


####脂肪酸与认知
RCol.fc <- list()#fc means fatty acid and cognition
PCol.fc <- list()
for (i in 1:length(faColname)) {
  
  y.fc <- faColname[i] 
  crntRcor.fc <- double()
  crntPcor.fc <- double()
  
  for (j in 1:length(cogColname)) {
    
    x.fc <- cogColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.fc <- paste(y.fc,"~")
    xhs.fc <- paste(x.fc,"~")
    
    frma1.fc <- as.formula(paste(yhs.fc,cov_f))
    mod1.fc <- lm(frma1.fc,data = mww_dbi,na.action = na.exclude)
    r1.fc <- resid(mod1.fc)#residuls for mod1
    
    frma2.fc <- as.formula(paste(xhs.fc,cov_f))
    mod2.fc <- lm(frma2.fc,data = mww_dbi,na.action = na.exclude)
    r2.fc <- resid(mod2.fc)
    
    rc.fc <- rcorr(r1.fc,r2.fc)
    
    crntRcor.fc[j] <- rc.fc$r[1,2]
    crntPcor.fc[j] <- rc.fc$P[1,2]
    
  }
  
  RCol.fc[[i]] <- crntRcor.fc
  PCol.fc[[i]] <- crntPcor.fc
  
}

RCorMat.fc <- matrix(unlist(RCol.fc), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, faColname))#R值矩阵

PCorMat.fc <- matrix(unlist(PCol.fc), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, faColname))#P值矩阵

## Partial correlation between plasma FA and inflammation
RCol.fi <- list()#fi means fatty acid and inflammation
PCol.fi <- list()

for (i in 1:length(faColname)) {
  
  y.fi <- faColname[i] 
  crntRcor.fi <- double()
  crntPcor.fi <- double()
  
  for (j in 1:length(infColname)) {
    
    x.fi <- infColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.fi <- paste(y.fi,"~")
    xhs.fi <- paste(x.fi,"~")
    
    frma1.fi <- as.formula(paste(yhs.fi,cov_f))
    mod1.fi <- lm(frma1.fi,data = mww_dbi,na.action = na.exclude)
    r1.fi <- resid(mod1.fi)#residuls for mod1
    
    frma2.fi <- as.formula(paste(xhs.fi,cov_f))
    mod2.fi <- lm(frma2.fi,data = mww_dbi,na.action = na.exclude)
    r2.fi <- resid(mod2.fi)
    
    rc.fi <- rcorr(r1.fi,r2.fi)
    
    crntRcor.fi[j] <- rc.fi$r[1,2]
    crntPcor.fi[j] <- rc.fi$P[1,2]
    
  }
  
  RCol.fi[[i]] <- crntRcor.fi
  PCol.fi[[i]] <- crntPcor.fi
  
}

RCorMat.fi <- matrix(unlist(RCol.fi), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, faColname))#R值矩阵

PCorMat.fi <- matrix(unlist(PCol.fi), 
                     nrow = length(infColname), 
                     dimnames = list(infColname, faColname))#P值矩阵


## Partial correlation between inflammation and cognition
RCol.ic <- list()#ic means inflammation and cognition
PCol.ic <- list()

for (i in 1:length(infColname)) {
  
  y.ic <- infColname[i] 
  crntRcor.ic <- double()
  crntPcor.ic <- double()
  
  for (j in 1:length(cogColname)) {
    
    x.ic <- cogColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.ic <- paste(y.ic,"~")
    xhs.ic <- paste(x.ic,"~")
    
    frma1.ic <- as.formula(paste(yhs.ic,cov_f))
    mod1.ic <- lm(frma1.ic,data = mww_dbi,na.action = na.exclude)
    r1.ic <- resid(mod1.ic)#residuls for mod1
    
    frma2.ic <- as.formula(paste(xhs.ic,cov_f))
    mod2.ic <- lm(frma2.ic,data = mww_dbi,na.action = na.exclude)
    r2.ic <- resid(mod2.ic)
    
    rc.ic <- rcorr(r1.ic,r2.ic)
    
    crntRcor.ic[j] <- rc.ic$r[1,2]
    crntPcor.ic[j] <- rc.ic$P[1,2]
    
  }
  
  RCol.ic[[i]] <- crntRcor.ic
  PCol.ic[[i]] <- crntPcor.ic
  
}

RCorMat.ic <- matrix(unlist(RCol.ic), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, infColname))#R值矩阵

PCorMat.ic <- matrix(unlist(PCol.ic), 
                     nrow = length(cogColname), 
                     dimnames = list(cogColname, infColname))#P值矩阵
###########################
###########################
p.ic <- as.data.frame(PCorMat.ic) %>% rownames_to_column("cog")
p.fc <- as.data.frame(PCorMat.fc) %>% rownames_to_column("cog")
p.fi <- as.data.frame(PCorMat.fi) %>% rownames_to_column("inf")
r.ic <- as.data.frame(RCorMat.ic) %>% rownames_to_column("cog")
r.fc <- as.data.frame(RCorMat.fc) %>% rownames_to_column("cog")
r.fi <- as.data.frame(RCorMat.fi) %>% rownames_to_column("inf")
write.csv(p.ic,file = "p ic.csv",row.names = F)
write.csv(r.ic,file = "r ic.csv",row.names = F)
write.csv(p.fc,file = "p fc.csv",row.names = F)
write.csv(r.fc,file = "r fc.csv",row.names = F)
write.csv(p.fi,file = "p fi.csv",row.names = F)
write.csv(r.fi,file = "r fi.csv",row.names = F)
#adiposity and cognition
p.ac <- as.data.frame(PCorMat.ac) %>% rownames_to_column("cog")
r.ac <- as.data.frame(RCorMat.ac) %>% rownames_to_column("cog")
write.csv(p.ac,file = "p ac.csv",row.names = F)
write.csv(r.ac,file = "r ac.csv",row.names = F)
#adipositiy and inflammation
p.ai <- as.data.frame(PCorMat.ai) %>% rownames_to_column("cog")
r.ai <- as.data.frame(RCorMat.ai) %>% rownames_to_column("cog")
write.csv(p.ai,file = "p ai.csv",row.names = F)
write.csv(r.ai,file = "r ai.csv",row.names = F)
###dbi与炎症因子
p.di <- as.data.frame(PCorMat.di) %>% rownames_to_column("cog")
r.di <- as.data.frame(RCorMat.di) %>% rownames_to_column("cog")
write.csv(p.di,file = "p di.csv",row.names = F)
write.csv(r.di,file = "r di.csv",row.names = F)
###dbi与cog
p.dc <- as.data.frame(PCorMat.dc) %>% rownames_to_column("cog")
r.dc <- as.data.frame(RCorMat.dc) %>% rownames_to_column("cog")
write.csv(p.dc,file = "p dc.csv",row.names = F)
write.csv(r.dc,file = "r dc.csv",row.names = F)
###dbi与adiposity
p.ad <- as.data.frame(PCorMat.ad) %>% rownames_to_column("cog")
r.ad <- as.data.frame(RCorMat.ad) %>% rownames_to_column("cog")
write.csv(p.ad,file = "p ad.csv",row.names = F)
write.csv(r.ad,file = "r ad.csv",row.names = F)

p.fc.sfa.mmse <- p.fc %>% 
  dplyr::select(cog,C150,C170)%>%
  dplyr::slice(1,3:6)#mmse分数

p.fc.sfa.moca <- p.fc %>% 
  dplyr::select(cog,C140,C150,C160,C170,C180)%>%
  dplyr::slice(2,7:13)#mocA分数

p.fi.sfa <- p.fi %>%
  select(inf,C150,C170)

p.ic.mmse <- p.ic %>% 
  dplyr::slice(1,3:6)#mmse分数

p.ic.moca <- p.ic %>% 
  dplyr::slice(2,7:13)#mmse分数

library(reshape2)
p.fc.sfa.mmse.m <- melt(p.fc.sfa.mmse,value.name = "p")
p.fc.sfa.mmse.m$fdr <- p.adjust(p.fc.sfa.mmse.m$p,method = "fdr")
p.fi.sfa.m <- melt(p.fi.sfa,value.name = "p")
p.fi.sfa.m$fdr <- p.adjust(p.fi.sfa.m$p,method = "fdr")
p.ic.mmse.m <- melt(p.ic.mmse,value.name = "p")
p.ic.mmse.m$fdr <- p.adjust(p.ic.mmse.m$p,method = "fdr")
p.ic.moca.m <- melt(p.ic.moca,value.name = "p")
p.ic.moca.m$fdr <- p.adjust(p.ic.moca.m$p,method = "fdr")

#Mediation

library(lavaan)
model_11 <- ' # direct effect
            mmse ~ c*C170+age+gender
# mediator
NFkB2 ~ a*C170+age+gender
mmse ~ b*NFkB2
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit11 <- sem(model_11, data = mww_dbi)
summary(fit11)

model_12 <- ' # direct effect
            MMSEComputation ~ c*C170+age+gender+bmi+edu+energy+smoke
# mediator
NFkB2 ~ a*C170+age+gender+bmi+edu+energy+smoke
MMSEComputation ~ b*NFkB2
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit12 <- sem(model_12, data = mww_dbi)
summary(fit12)

model_13 <- ' # direct effect
            MoCALanguageskills ~ c*C160+age+gender+bmi+edu+energy+smoke
# mediator
IL102 ~ a*C160+age+gender+bmi+edu+energy+smoke
MoCALanguageskills ~ b*IL102
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit13 <- sem(model_13, data = mww_dbi)
summary(fit13)

model_14 <- ' # direct effect
            MoCALanguageskills ~ c*C180+age+gender+bmi+edu+energy+smoke
# mediator
IL102 ~ a*C180+age+gender+bmi+edu+energy+smoke
MoCALanguageskills ~ b*IL102
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit14 <- sem(model_14, data = mww_dbi)
summary(fit14)

model_15 <- ' # direct effect
            MoCAVisualspatialability ~ c*C160+age+gender+bmi+edu+energy+smoke
# mediator
lps ~ a*C160+age+gender+bmi+edu+energy+smoke
MoCAVisualspatialability ~ b*lps
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit15 <- sem(model_15, data = mww_dbi)
summary(fit15)

###mediation dbi and crp
model_lbc <- ' # direct effect
            crp3 ~ c*lbs+age+gender+edu+energy+smoke
# mediator
bmi ~ a*lbs+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitlbc <- sem(model_lbc, data = mww_dbi)
summary(fitlbc)

model_hbc <- ' # direct effect
            crp3 ~ c*hbs+age+gender+edu+energy+smoke
# mediator
bmi ~ a*hbs+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fithbc <- sem(model_hbc, data = mww_dbi)
summary(fithbc)


model_dbc <- ' # direct effect
            crp3 ~ c*dqd+age+gender+edu+energy+smoke
# mediator
bmi ~ a*dqd+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitdbc <- sem(model_dbc, data = mww_dbi)
summary(fitdbc)

model_d3bc <- ' # direct effect
            crp3 ~ c*dbi3+age+gender+edu+energy+smoke
# mediator
bmi ~ a*dbi3+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitd3bc <- sem(model_d3bc, data = mww_dbi)
summary(fitd3bc)

model_d5bc <- ' # direct effect
            crp3 ~ c*dbi5+age+gender+edu+energy+smoke
# mediator
bmi ~ a*dbi5+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitd5bc <- sem(model_d5bc, data = mww_dbi)
summary(fitd5bc)

model_d6bc <- ' # direct effect
            crp3 ~ c*dbi6+age+gender+edu+energy+smoke
# mediator
bmi ~ a*dbi6+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitd6bc <- sem(model_d6bc, data = mww_dbi)
summary(fitd6bc)

model_d7bc <- ' # direct effect
            crp3 ~ c*dbi7+age+gender+edu+energy+smoke
# mediator
bmi ~ a*dbi7+age+gender+edu+energy+smoke
crp3 ~ b*bmi
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitd7bc <- sem(model_d7bc, data = mww_dbi)
summary(fitd7bc)

hist(mww_dbi$C1702)
mww_dbi$C1702 <- sqrt(mww_dbi$C170)
###从DBI计算的R script文件中计算乳品摄入量及其对应的DBI分数
dbi_dairy <- select(dbi,id,dbi_dairy,dairy)
mww_dbi <- left_join(mww_dbi,dbi_dairy,by="id")
fit_od <- lm(C170~dbi_dairy+age+gender+energy+edu,data=mww_dbi)
summary(fit_od)


#bmi
fit.bmi <- lm(NFkB2~bmi+age+gender,data = mww_dbi)
summary(fit.bmi)
fit.bmi2 <- lm(tnf3~bmi+age+gender,data = mww_dbi)
summary(fit.bmi2)
fit.bmi3 <- lm(crp3~bmi+age+gender,data = mww_dbi)
summary(fit.bmi3)
fit.bmi4 <- lm(IL102~bmi+age+gender,data = mww_dbi)
summary(fit.bmi4)
fit.bmi5 <- lm(IL1β2~bmi+age+gender,data = mww_dbi)
summary(fit.bmi5)


