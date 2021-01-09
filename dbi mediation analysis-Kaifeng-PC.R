###粗略看一下跟认知的关系
library(tidyverse)
dbi_score <- read.csv(file = "DBI/DBI膳食质量分数.csv",stringsAsFactors = F)
FA <- read.csv(file = "20200108 数据库(放了血液脂肪酸-重新核对脂肪酸正确版-去掉359兰孝明重复的+计算了认知的小分).csv",stringsAsFactors = F)
FA2 <- FA %>% select(登记号,12:15,年龄,性别,体重指数,腰围,腰臀比,教育,mmse总分,MoCA总分,吸烟史,饮酒史,内毒素,超敏C,能量,
                       MMSEOrientation,MMSEComputation,MMSEMemory,MMSELanguageskill,MoCANaming,MoCAOrientation,
                       MoCADelayedrecall,MoCAAbstractthinking,MoCALanguageskills,MoCAVisualspatialability,MoCAAttention,) %>%
  rename(id=登记号,crp=超敏C,lps=内毒素,age=年龄,gender=性别,bmi=体重指数,wc=腰围,wh=腰臀比,
         edu=教育,mmse=mmse总分,moca=MoCA总分,energy=能量,smoke=吸烟史,drink=饮酒史)
                       
fa_p <- read.csv(file = "mww_FAonly_p.csv",stringsAsFactors = F)%>%select(-age,-gender)
mww_dbi <- left_join(FA2,dbi_score,by="id")
mww_dbi <- inner_join(FA2,fa_p,by="id")
mww_mplus.dbi <- mww_dbi
mww_mplus.dbi[is.na(mww_mplus.dbi)] <- 999
write.csv(mww_mplus.dbi,file = "mplus/mplus with dbi.csv",row.names = F)

#CRP
boxplot(mww_dbi$crp)
describe(mww_dbi$crp)
mww_dbi$crp[mww_dbi$crp>400] <- NA
hist(mww_dbi$crp)
mww_dbi$crp2 <- sqrt(mww_dbi$crp)
hist(mww_dbi$crp2)

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
describe(mww_dbi$IL102)
hist(mww_dbi$IL102)

#NFkB2
describe(mww_dbi$NFkB)
hist(mww_dbi$NFkB)
mww_dbi$NFkB2 <- sqrt(mww_dbi$NFkB)
hist(mww_dbi$NFkB2)

#lsp
describe(mww_dbi$lps)
hist(mww_dbi$lps)
mww_dbi$lps[mww_dbi$lps<0] <- 0
mww_dbi$lps2 <- mww_dbi$lps
mww_dbi$lps2 <- sqrt(mww_dbi$lps)
hist(mww_dbi$lps)

#C14
describe(mww_dbi$C140)
hist(mww_dbi$C140)
boxplot(mww_dbi$C140)
boxplot(mww_dbi$C160)
boxplot(mww_dbi$C180)

###数据校正后再看一下partial correlation
mww_cov <- mww_dbi %>% 
  dplyr::select(age,gender,bmi,energy,smoke,drink,edu)#混杂因子
mww_fa <- mww_dbi%>%
  dplyr::select(c(30:59))
mww_cog <- mww_dbi%>%
  dplyr::select(mmse,moca,c(19:29))
mww_inf <- mww_dbi%>%
  dplyr::select(IL1β2,IL102,NFkB2,tnf3,crp2,lps2)

RCol.fc <- list()#fc means fatty acid and cognition
PCol.fc <- list()
faColname <- colnames(mww_fa)
cogColname <- colnames(mww_cog)
covColname <- colnames(mww_cov)
infColname <- colnames(mww_inf)

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
#Mediation

library(lavaan)
model_11 <- ' # direct effect
MMSEOrientation ~ d*dqd
# mediator
bmi ~ a*dqd+age+gender
tnf3 ~ b*bmi+age+gender
MMSEOrientation~c*tnf3+age+gender
# indirect effect (a*b*c)
abc := a*b*c
# total effect
total := d + (a*b*c)
'
fit11 <- sem(model_11, data = mww_dbi)
summary(fit11)

library(lavaan)
########C17
modelc17<- ' # direct effect
  mmse ~ c*C170+age+gender
  # mediator
  NFkB2 ~ a*C170+age+gender
  mmse ~ b*NFkB2
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
fitc17 <- sem(modelc17, data = mww_dbi)
summary(fitc17) 
