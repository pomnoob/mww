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
write.csv(mww2,file = "data/20201201 录入核对数据库 数据清理核对（最终不改版）.csv")
# 导入DBI数据
dbi_score <- read.csv(file = "DBI/DBI膳食质量分数.csv",stringsAsFactors = F)
# 转化为整数型
mwwSel$id <- as.integer(mwwSel$id)
# 合并数据
mww_dbi <- left_join(mwwSel,dbi_score,by="id")
# 倒入脂肪酸数据
fa_p <- read.csv(file = "mww_FAonly_p.csv",stringsAsFactors = F)%>%
  select(-age,-gender)

mww_fa_NA.p <- dplyr::filter(fa_p,!is.na(C140))#删去没有做脂肪酸检测的观测n=6
mww_fa_NAp2 <- mww_fa_NA.p[,which(colMeans(mww_fa_NA.p!=0)>=0.5)]#删除0值超过50%的脂肪酸
mww_fa_NAp2 <- dplyr::mutate(mww_fa_NAp2,sfa=C140+C150+C160+C170+C180,mufa=C181n9c,pufa=C182n6c+C183n6
                             +C183n3+C203n6+C204n6+C205n3+C226n3,n3=C183n3+C205n3+C226n3,n6=C182n6c+C183n6
                             +C203n6+C204n6,scd16=C161/C160,scd18=C181n9c/C180,d6d=C183n6/C182n6c,d5d=C204n6/C203n6)
mww_fa_NAp2$d5d[is.infinite(mww_fa_NAp2$d5d)] <- NA
mww_fa_NAp2$d6d[is.infinite(mww_fa_NAp2$d6d)] <- NA
mww_fa_NAp2$scd16[is.infinite(mww_fa_NAp2$scd16)] <- NA
mww_fa_NAp2$scd18[is.infinite(mww_fa_NAp2$scd18)] <- NA
mww_dbi <- left_join(mww_dbi,mww_fa_NAp2,by="id")

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

summary(mww_dbi$wh)
boxplot(mww_dbi$wh)
mww_dbi$wh[mww_dbi$wh>2] <- NA
mww_dbi$wh[mww_dbi$wh<0.6] <- NA
hist(mww_dbi$wh)
