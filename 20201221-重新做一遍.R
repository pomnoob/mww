library(haven)
library(tidyverse)
#倒入新版数据
mww2 <- read_sav(file = "data/20201201 录入核对数据库 数据清理核对（最终不改版）(1).sav")
# 选择需要分析的变量
mwwSel <- mww2 %>%
  select(登记号,12:15,年龄,性别,BMI,腰围,腰臀比,教育,mmse总分,MoCA总分,吸烟史,饮酒史,内毒素,超敏C,能量,C反应蛋白,
            OGTT,OGTT2,胰岛素,MMSEOrientation,MMSEComputation,MMSEMemory,MMSELanguageskill,MoCANaming,MoCAOrientation,
            MoCADelayedrecall,MoCAAbstractthinking,MoCALanguageskills,MoCAVisualspatialability,MoCAAttention) %>%
  rename(id= 登记号,hscrp=超敏C,crp=C反应蛋白,insulin=胰岛素,lps=内毒素,age=年龄,gender=性别,bmi=BMI,wc=腰围,wh=腰臀比,
         edu=教育,mmse=mmse总分,moca=MoCA总分,energy=能量,smoke=吸烟史,drink=饮酒史)
# 导入新的脂肪酸数据
# MAC 系统下中文读取有些问题，把原文件中的中文名改成了英文字母
fa <- read.csv(file = "data/20201201 录入核对数据库 数据清理核对（最终不改版）-脂肪酸.csv")

# 0值变成 NA
fa[fa==0] <- NA
# 去除NA大于一半的脂肪酸
mww_fa<- fa[,which(colMeans(is.na(fa))<0.5)]#删除0值超过50%的脂肪酸
mww_fa <- mww_fa %>%
  mutate(C140p=C140/TotalFA*100,
         C150p=C150/TotalFA*100,
         C160p=C160/TotalFA*100,
         C161p=C161/TotalFA*100,
         C170p=C170/TotalFA*100,
         C180p=C180/TotalFA*100,
         C181n9cp=C181n9c/TotalFA*100,
         C182n6p=C182n6c/TotalFA*100,
         C183n6p=C183n6/TotalFA*100,
         C183n3p=C183n3/TotalFA*100,
         C203n6p=C203n6/TotalFA*100,
         C204n6p=C204n6/TotalFA*100,
         C205n3p=C205n3/TotalFA*100,
         C226n3p=C226n3/TotalFA*100)

mww_fa <- dplyr::mutate(mww_fa,sfa=C140p+C150p+C160p+C170p+C180p,mufa=C181n9cp,pufa=C182n6p+C183n6p
                             +C183n3p+C203n6p+C204n6p+C205n3p+C226n3p,n3=C183n3p+C205n3p+C226n3p,n6=C182n6p+C183n6p
                             +C203n6p+C204n6p,scd16=C161p/C160p,scd18=C181n9cp/C180p,d6d=C183n6p/C182n6p,d5d=C204n6p/C203n6p)
mww_fa$d5d[is.infinite(mww_fa$d5d)] <- NA
mww_fa$d6d[is.infinite(mww_fa$d6d)] <- NA
mww_fa$scd16[is.infinite(mww_fa$scd16)] <- NA
mww_fa$scd18[is.infinite(mww_fa$scd18)] <- NA
mww_fa <- mww_fa %>%
  rename(id=djh) # djh = 登记号
write.csv(mww_fa,file = "data/脂肪酸比例.csv",row.names=F)

# 导入DBI数据
# MAC 系统里面最好将csv保存为UTF8格式，然后读取的时候添加选项
dbi_score <- read.csv(file = "DBI/DBI膳食质量分数 UTF8.csv",stringsAsFactors = F,
                      fileEncoding = "UTF-8")
# 转化为整数型
mwwSel$id <- as.integer(mwwSel$id)
# 合并数据
mww_dbi <- left_join(mwwSel,dbi_score,by="id")
mww_dbi <- left_join(mww_dbi,mww_fa,by="id")

##筛一遍炎症因子和血糖胰岛素等
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

#hscrp
describe(mww_dbi$hscrp)
hist(mww_dbi$hscrp)
mww_dbi$hscrp2 <- mww_dbi$hscrp
mww_dbi$hscrp2[mww_dbi$hscrp2>420] <- NA
hist(mww_dbi$hscrp2)
mww_dbi$hscrp3 <- sqrt(mww_dbi$hscrp2)
hist(mww_dbi$hscrp3)
describe(mww_dbi$hscrp3)

#crp
describe(mww_dbi$crp)
hist(mww_dbi$crp)
boxplot(mww_dbi$crp)
mww_dbi$crp2 <- mww_dbi$crp
mww_dbi$crp2[mww_dbi$crp>70] <- NA
describe(mww_dbi$crp2)
boxplot(mww_dbi$crp2)
hist(mww_dbi$crp2)
mww_dbi$crp3 <- sqrt(mww_dbi$crp2)
hist(mww_dbi$crp3)

#OGTT
describe(mww_dbi$OGTT)
hist(mww_dbi$OGTT)
boxplot(mww_dbi$OGTT)

#OGTT2
describe(mww_dbi$OGTT2)
hist(mww_dbi$OGTT2)
boxplot(mww_dbi$OGTT2)

#insulin
describe(mww_dbi$insulin)
hist(mww_dbi$insulin)
boxplot(mww_dbi$insulin)
mww_dbi$insulin[mww_dbi$insulin>400] <- NA
mww_dbi$insulin[mww_dbi$insulin<0.1] <- NA
mww_dbi$insulin2 <- mww_dbi$insulin
mww_dbi$insulin2 <- sqrt(mww_dbi$insulin2)
hist(mww_dbi$insulin2)

#计算HOMA-IR和QUICKI
mww_dbi <- dplyr::mutate(mww_dbi,homa_ir=(OGTT+insulin)/22.5,
                         quicki=1/(log(insulin)+log(OGTT*18)))
describe(mww_dbi$homa_ir)
hist(mww_dbi$homa_ir)
is.infinite(mww_dbi$homa_ir)
hist(mww_dbi$quicki)
describe(mww_dbi$quicki)
mww_dbi$homa2 <- log10(mww_dbi$homa_ir)
describe(mww_dbi$homa2)
hist(mww_dbi$homa2)
mww_dbi$quicki2 <- log10(mww_dbi$quicki)
hist(mww_dbi$quicki2)
write.csv(mww_dbi,file = "20201221 分析用样本.csv")
##############
mww_cov <- mww_dbi %>% 
  dplyr::select(age,gender,bmi,energy,smoke,edu,drink)#混杂因子
mww_fa <- mww_dbi%>%
  dplyr::select(C140p,C150p,C160p,C161p,C170p,C180p,C181n9cp,C182n6p,C183n6p,
                C183n3p,C203n6p,C204n6p,C205n3p,C226n3p,sfa,mufa,pufa,n3,n6,
                scd16,scd18,d6d,d5d)

mww_cog <- mww_dbi%>%
  dplyr::select(mmse,moca,MMSEOrientation,MMSEComputation,MMSEMemory,MMSELanguageskill,MoCANaming,MoCAOrientation,
                MoCADelayedrecall,MoCAAbstractthinking,MoCALanguageskills,MoCAVisualspatialability,MoCAAttention)

mww_inf <- mww_dbi %>%
  dplyr::select(IL1β2,IL102,NFkB2,tnf3,hscrp3,crp3,lps,OGTT,OGTT2,insulin,homa2,quicki2)

mww_dq <- mww_dbi %>%
  dplyr::select(dbi1:dbi7,lbs,hbs,dqd)

mww_ad <- mww_dbi %>%
  select(bmi,wc,wh)

mww_cov2 <- mww_dbi %>% 
  dplyr::select(age,gender,energy,smoke,edu,drink)#混杂因子2,for mww_ad

######multiple linear regression########################
faColname <- colnames(mww_fa)
cogColname <- colnames(mww_cog)
covColname <- colnames(mww_cov)
dqColname <- colnames(mww_dq)
infColname <- colnames(mww_inf)
adColname <- colnames(mww_ad)
cov2Colname <- colnames(mww_cov2)

#plasma FA and cognition

RCol.fc <- list()#fc means fatty acid and cognition
PCol.fc <- list()
for (i in 1:length(cogColname)) {
  
  y.fc <- cogColname[i] 
  crntRcor.fc <- double()
  crntPcor.fc <- double()
  
  for (j in 1:length(faColname)) {
    
    x.fc <- faColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.fc <- paste(y.fc,"~")
    xhs.fc <- paste(x.fc,"+")
    
    frma.fc <- as.formula(paste(yhs.fc,xhs.fc,cov_f))
    mod.fc <- lm(frma.fc,data = mww_dbi,na.action = na.exclude)
    coef.fc <- coef(mod.fc)
    names(coef.fc) <- NULL
    p.fc <- anova(mod.fc)
    
    crntRcor.fc[j] <- coef.fc[2]
    crntPcor.fc[j] <- p.fc[1,5]
    
  }
  
  RCol.fc[[i]] <- crntRcor.fc
  PCol.fc[[i]] <- crntPcor.fc
  
}

RCorMat.fc <- matrix(unlist(RCol.fc), 
                     nrow = length(faColname), 
                     dimnames = list(faColname, cogColname))#coef矩阵

PCorMat.fc <- matrix(unlist(PCol.fc), 
                     nrow = length(faColname), 
                     dimnames = list(faColname, cogColname))#P值矩阵

p.fc.mlr <- as.data.frame(PCorMat.fc) %>% rownames_to_column("FA")
e.fc.mlr <- as.data.frame(RCorMat.fc) %>% rownames_to_column("FA")
write.csv(e.fc.mlr,file = "coef of fc using MLR.csv",row.names = F)
write.csv(p.fc.mlr,file = "p value of fc using MLR.csv",row.names = F)

#plasma FA, inflammation and biochemistry

RCol.fib <- list()
PCol.fib <- list()
for (i in 1:length(infColname)) {
  
  y.fib <- infColname[i] 
  crntRcor.fib <- double()
  crntPcor.fib <- double()
  
  for (j in 1:length(faColname)) {
    
    x.fib <- faColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.fib <- paste(y.fib,"~")
    xhs.fib <- paste(x.fib,"+")
    
    frma.fib <- as.formula(paste(yhs.fib,xhs.fib,cov_f))
    mod.fib <- lm(frma.fib,data = mww_dbi,na.action = na.exclude)
    coef.fib <- coef(mod.fib)
    names(coef.fib) <- NULL
    p.fib <- anova(mod.fib)
    
    crntRcor.fib[j] <- coef.fib[2]
    crntPcor.fib[j] <- p.fib[1,5]
    
  }
  
  RCol.fib[[i]] <- crntRcor.fib
  PCol.fib[[i]] <- crntPcor.fib
  
}

RCorMat.fib <- matrix(unlist(RCol.fib), 
                      nrow = length(faColname), 
                      dimnames = list(faColname, infColname))#coef矩阵

PCorMat.fib <- matrix(unlist(PCol.fib), 
                      nrow = length(faColname), 
                      dimnames = list(faColname, infColname))#P值矩阵

p.fib.mlr <- as.data.frame(PCorMat.fib) %>% rownames_to_column("FA")
e.fib.mlr <- as.data.frame(RCorMat.fib) %>% rownames_to_column("FA")
write.csv(e.fib.mlr,file = "coef of fib using MLR.csv",row.names = F)
write.csv(p.fib.mlr,file = "p value of fib using MLR.csv",row.names = F)

#cognition,inflammation and chemistry

RCol.cib <- list()
PCol.cib <- list()
for (i in 1:length(cogColname)) {
  
  y.cib <- cogColname[i] 
  crntRcor.cib <- double()
  crntPcor.cib <- double()
  
  for (j in 1:length(infColname)) {
    
    x.cib <-infColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.cib <- paste(y.cib,"~")
    xhs.cib <- paste(x.cib,"+")
    
    frma.cib <- as.formula(paste(yhs.cib,xhs.cib,cov_f))
    mod.cib <- lm(frma.cib,data = mww_dbi,na.action = na.exclude)
    coef.cib <- coef(mod.cib)
    names(coef.cib) <- NULL
    p.cib <- anova(mod.cib)
    
    crntRcor.cib[j] <- coef.cib[2]
    crntPcor.cib[j] <- p.cib[1,5]
    
  }
  
  RCol.cib[[i]] <- crntRcor.cib
  PCol.cib[[i]] <- crntPcor.cib
  
}

RCorMat.cib <- matrix(unlist(RCol.cib), 
                      nrow = length(infColname), 
                      dimnames = list(infColname, cogColname))#coef矩阵

PCorMat.cib <- matrix(unlist(PCol.cib), 
                      nrow = length(infColname), 
                      dimnames = list(infColname, cogColname))#P值矩阵

p.cib.mlr <- as.data.frame(PCorMat.cib) %>% rownames_to_column("FA")
e.cib.mlr <- as.data.frame(RCorMat.cib) %>% rownames_to_column("FA")
write.csv(e.cib.mlr,file = "coef of cib using MLR.csv",row.names = F)
write.csv(p.cib.mlr,file = "p value of cib using MLR.csv",row.names = F)

#cognition and dbi

RCol.cd <- list()
PCol.cd <- list()
for (i in 1:length(cogColname)) {
  
  y.cd <- cogColname[i] 
  crntRcor.cd <- double()
  crntPcor.cd <- double()
  
  for (j in 1:length(dqColname)) {
    
    x.cd <-dqColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.cd <- paste(y.cd,"~")
    xhs.cd <- paste(x.cd,"+")
    
    frma.cd <- as.formula(paste(yhs.cd,xhs.cd,cov_f))
    mod.cd <- lm(frma.cd,data = mww_dbi,na.action = na.exclude)
    coef.cd <- coef(mod.cd)
    names(coef.cd) <- NULL
    p.cd <- anova(mod.cd)
    
    crntRcor.cd[j] <- coef.cd[2]
    crntPcor.cd[j] <- p.cd[1,5]
    
  }
  
  RCol.cd[[i]] <- crntRcor.cd
  PCol.cd[[i]] <- crntPcor.cd
  
}

RCorMat.cd <- matrix(unlist(RCol.cd), 
                     nrow = length(dqColname), 
                     dimnames = list(dqColname, cogColname))#coef矩阵

PCorMat.cd <- matrix(unlist(PCol.cd), 
                     nrow = length(dqColname), 
                     dimnames = list(dqColname, cogColname))#P值矩阵

p.cd.mlr <- as.data.frame(PCorMat.cd) %>% rownames_to_column("FA")
e.cd.mlr <- as.data.frame(RCorMat.cd) %>% rownames_to_column("FA")
write.csv(e.cd.mlr,file = "coef of cd using MLR.csv",row.names = F)
write.csv(p.cd.mlr,file = "p value of cd using MLR.csv",row.names = F)

#Inflammation, biochemistry and dbi

RCol.ibd <- list()
PCol.ibd <- list()
for (i in 1:length(infColname)) {
  
  y.ibd <- infColname[i] 
  crntRcor.ibd <- double()
  crntPcor.ibd <- double()
  
  for (j in 1:length(dqColname)) {
    
    x.ibd <-dqColname[j] 
    cov_f <- paste(covColname,collapse = "+")
    yhs.ibd <- paste(y.ibd,"~")
    xhs.ibd <- paste(x.ibd,"+")
    
    frma.ibd <- as.formula(paste(yhs.ibd,xhs.ibd,cov_f))
    mod.ibd <- lm(frma.ibd,data = mww_dbi,na.action = na.exclude)
    coef.ibd <- coef(mod.ibd)
    names(coef.ibd) <- NULL
    p.ibd <- anova(mod.ibd)
    
    crntRcor.ibd[j] <- coef.ibd[2]
    crntPcor.ibd[j] <- p.ibd[1,5]
    
  }
  
  RCol.ibd[[i]] <- crntRcor.ibd
  PCol.ibd[[i]] <- crntPcor.ibd
  
}

RCorMat.ibd <- matrix(unlist(RCol.ibd), 
                      nrow = length(dqColname), 
                      dimnames = list(dqColname, infColname))#coef矩阵

PCorMat.ibd <- matrix(unlist(PCol.ibd), 
                      nrow = length(dqColname), 
                      dimnames = list(dqColname, infColname))#P值矩阵

p.ibd.mlr <- as.data.frame(PCorMat.ibd) %>% rownames_to_column("FA")
e.ibd.mlr <- as.data.frame(RCorMat.ibd) %>% rownames_to_column("FA")
write.csv(e.ibd.mlr,file = "coef of ibd using MLR.csv",row.names = F)
write.csv(p.ibd.mlr,file = "p value of ibd using MLR.csv",row.names = F)

#Adiposity and dbi

RCol.add <- list()
PCol.add <- list()
for (i in 1:length(adColname)) {
  
  y.add <- adColname[i] 
  crntRcor.add <- double()
  crntPcor.add <- double()
  
  for (j in 1:length(dqColname)) {
    
    x.add <-dqColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.add <- paste(y.add,"~")
    xhs.add <- paste(x.add,"+")
    
    frma.add <- as.formula(paste(yhs.add,xhs.add,cov_f2))
    mod.add <- lm(frma.add,data = mww_dbi,na.action = na.exclude)
    coef.add <- coef(mod.add)
    names(coef.add) <- NULL
    p.add <- anova(mod.add)
    
    crntRcor.add[j] <- coef.add[2]
    crntPcor.add[j] <- p.add[1,5]
    
  }
  
  RCol.add[[i]] <- crntRcor.add
  PCol.add[[i]] <- crntPcor.add
  
}

RCorMat.add <- matrix(unlist(RCol.add), 
                      nrow = length(dqColname), 
                      dimnames = list(dqColname, adColname))#coef矩阵

PCorMat.add <- matrix(unlist(PCol.add), 
                      nrow = length(dqColname), 
                      dimnames = list(dqColname, adColname))#P值矩阵

p.add.mlr <- as.data.frame(PCorMat.add) %>% rownames_to_column("FA")
e.add.mlr <- as.data.frame(RCorMat.add) %>% rownames_to_column("FA")
write.csv(e.add.mlr,file = "coef of add using MLR.csv",row.names = F)
write.csv(p.add.mlr,file = "p value of add using MLR.csv",row.names = F)

#Adiposity,inflammation and biochemistry marker

RCol.aib <- list()
PCol.aib <- list()
for (i in 1:length(infColname)) {
  
  y.aib <- infColname[i] 
  crntRcor.aib <- double()
  crntPcor.aib <- double()
  
  for (j in 1:length(adColname)) {
    
    x.aib <-adColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.aib <- paste(y.aib,"~")
    xhs.aib <- paste(x.aib,"+")
    
    frma.aib <- as.formula(paste(yhs.aib,xhs.aib,cov_f2))
    mod.aib <- lm(frma.aib,data = mww_dbi,na.action = na.exclude)
    coef.aib <- coef(mod.aib)
    names(coef.aib) <- NULL
    p.aib <- anova(mod.aib)
    
    crntRcor.aib[j] <- coef.aib[2]
    crntPcor.aib[j] <- p.aib[1,5]
    
  }
  
  RCol.aib[[i]] <- crntRcor.aib
  PCol.aib[[i]] <- crntPcor.aib
  
}

RCorMat.aib <- matrix(unlist(RCol.aib), 
                      nrow = length(adColname), 
                      dimnames = list(adColname, infColname))#coef矩阵

PCorMat.aib <- matrix(unlist(PCol.aib), 
                      nrow = length(adColname), 
                      dimnames = list(adColname, infColname))#P值矩阵

p.aib.mlr <- as.data.frame(PCorMat.aib) %>% rownames_to_column("FA")
e.aib.mlr <- as.data.frame(RCorMat.aib) %>% rownames_to_column("FA")
write.csv(e.aib.mlr,file = "coef of aib using MLR.csv",row.names = F)
write.csv(p.aib.mlr,file = "p value of aib using MLR.csv",row.names = F)

#Adiposity and cognition

RCol.ac <- list()
PCol.ac <- list()
for (i in 1:length(cogColname)) {
  
  y.ac <- cogColname[i] 
  crntRcor.ac <- double()
  crntPcor.ac <- double()
  
  for (j in 1:length(adColname)) {
    
    x.ac <-adColname[j] 
    cov_f2 <- paste(cov2Colname,collapse = "+")
    yhs.ac <- paste(y.ac,"~")
    xhs.ac <- paste(x.ac,"+")
    
    frma.ac <- as.formula(paste(yhs.ac,xhs.ac,cov_f2))
    mod.ac <- lm(frma.ac,data = mww_dbi,na.action = na.exclude)
    coef.ac <- coef(mod.ac)
    names(coef.ac) <- NULL
    p.ac <- anova(mod.ac)
    
    crntRcor.ac[j] <- coef.ac[2]
    crntPcor.ac[j] <- p.ac[1,5]
    
  }
  
  RCol.ac[[i]] <- crntRcor.ac
  PCol.ac[[i]] <- crntPcor.ac
  
}

RCorMat.ac <- matrix(unlist(RCol.ac), 
                     nrow = length(adColname), 
                     dimnames = list(adColname, cogColname))#coef矩阵

PCorMat.ac <- matrix(unlist(PCol.ac), 
                     nrow = length(adColname), 
                     dimnames = list(adColname, cogColname))#P值矩阵

p.ac.mlr <- as.data.frame(PCorMat.ac) %>% rownames_to_column("FA")
e.ac.mlr <- as.data.frame(RCorMat.ac) %>% rownames_to_column("FA")
write.csv(e.ac.mlr,file = "coef of ac using MLR.csv",row.names = F)
write.csv(p.ac.mlr,file = "p value of ac using MLR.csv",row.names = F)

#############mediation analysis for Mplus#########
mww_dbimplus <- mww_dbi %>%
  select(age,gender,bmi,smoke,drink,energy,
         mmse,moca,23:33,66:86,sfa,mufa,IL1β2,IL102,NFkB2,tnf3,
         hscrp3,crp3,lps,OGTT,OGTT2,insulin,homa2,quicki2,
         dbi1:dbi7,lbs,hbs,dqd)%>%
  mutate(kcal=energy/1000,
         d6d=d6d*100,
         scd16=scd16*10) # 确保数量级

mww_dbimplus[is.na(mww_dbimplus)] <- 999



write.csv(mww_dbimplus,file = "mplus/mww_dbimplus 20210109.csv",fileEncoding = "UTF-8",
          row.names = F)
# 然后使用MplusAutomation进行中介分析
