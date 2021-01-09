######麻老师脂肪酸与认知数据-数据20200108#########


library(tidyverse)

######百分比
mww_fa.p <- read.csv("mww_FAonly_p.csv",stringsAsFactors = F)#FAonly 仅包含脂肪酸数据
mww_fa_NA.p <- dplyr::filter(mww_fa.p,!is.na(C140)& !is.na(age))#删去没有做脂肪酸检测的观测n=6，以及年龄性别有缺失的观测n=1
mww_fa_NAp2 <- mww_fa_NA.p[,which(colMeans(mww_fa_NA.p!=0)>=0.5)]#删除0值超过50%的脂肪酸

#princomp
mww_fa_p.prin <- princomp(mww_fa_NAp2[,c(4:17)],scores = T,cor = T)
screeplot(mww_fa_p.prin,type = "lines")
mww_fa_p.pca$loadings
view(mww_fa_p.pca$scores)

#female
mww_fa_female <- dplyr::filter(mww_fa_NAp2,gender==2)
mww_fa_female.pca <- princomp(mww_fa_female[,c(4:17)],scores = T,cor = T)
screeplot(mww_fa_female.pca,type = "lines")
mww_fa_female.pca$loadings

#male
mww_fa_male <- dplyr::filter(mww_fa_NAp2,gender==1)
mww_fa_male.pca <- princomp(mww_fa_male[,c(4:17)],scores = T,cor = T)
screeplot(mww_fa_male.pca,type = "lines")
mww_fa_male.pca$loadings

#random sampling
set.seed(101)
ran.samp <- sample.int(n=nrow(mww_fa_NAp2),size=floor(0.5*nrow(mww_fa_NAp2)),replace = F)
mww_fa_r1 <- mww_fa_NAp2[ran.samp,]
mww_fa_r2 <- mww_fa_NAp2[-ran.samp,]
  ##r1
  mww_fa_r1.pca <- princomp(mww_fa_r1[,c(4:17)],scores = T,cor = T)
  screeplot(mww_fa_r1.pca,type = "lines")
  mww_fa_r1.pca$loadings
  ##r2
  mww_fa_r2.pca <- princomp(mww_fa_r2[,c(4:17)],scores = T,cor = T)
  screeplot(mww_fa_r2.pca,type = "lines")
  mww_fa_r2.pca$loadings
  
library(mixOmics)
  fa.pca <- pca(mww_fa_NAp2[,c(4:17)],scale = T,ncomp = 4)
  fa.pca$rotation

  fa.pca.f <- pca(mww_fa_female[,c(4:17)],scale = T,ncomp = 4)
  fa.pca.f$rotation

  fa.pca.m <- pca(mww_fa_male[,c(4:17)],scale = T,ncomp = 4)
  fa.pca.m$rotation

  fa.pca.r1 <- pca(mww_fa_r1[,c(4:17)],scale = T,ncomp = 4)
  fa.pca.r1$rotation
  
  fa.pca.r2 <- pca(mww_fa_r2[,c(4:17)],scale = T,ncomp = 4)
  fa.pca.r2$rotation

#对比男女之间脂肪酸差异
  t.test(C140~gender,data=mww_fa_NAp2)
  t.test(C150~gender,data=mww_fa_NAp2)
  t.test(C160~gender,data=mww_fa_NAp2)#不显著
  t.test(C161~gender,data=mww_fa_NAp2)
  t.test(C170~gender,data=mww_fa_NAp2)#不显著
  t.test(C180~gender,data=mww_fa_NAp2)#不显著
  t.test(C181n9c~gender,data=mww_fa_NAp2)#不显著
  t.test(C182n6c~gender,data=mww_fa_NAp2)#不显著
  t.test(C183n6~gender,data=mww_fa_NAp2)
  t.test(C183n3~gender,data=mww_fa_NAp2)#不显著
  t.test(C203n6~gender,data=mww_fa_NAp2)
  t.test(C204n6~gender,data=mww_fa_NAp2)#不显著
  t.test(C205n3~gender,data=mww_fa_NAp2)
  t.test(C226n3~gender,data=mww_fa_NAp2)
  
  ########################################
  ########################################
  ########################################
  ########################################
  ##2020-5-7##############################
  library(tidyverse)
  mww_fa.p <- read.csv("mww_FAonly_p.csv",stringsAsFactors = F)#FAonly 仅包含脂肪酸数据
  mww_fa_NA.p <- dplyr::filter(mww_fa.p,!is.na(C140)& !is.na(age))#删去没有做脂肪酸检测的观测n=6，以及年龄性别有缺失的观测n=1
  mww_fa_NAp2 <- mww_fa_NA.p[,which(colMeans(mww_fa_NA.p!=0)>=0.5)]#删除0值超过50%的脂肪
  mww_fa_NAp2$scd16 <- mww_fa_NAp2$C161/mww_fa_NAp2$C160
  mww_fa_NAp2$scd18 <- mww_fa_NAp2$C181n9c/mww_fa_NAp2$C180
  mww_fa_NAp2$d6d <- mww_fa_NAp2$C183n6/mww_fa_NAp2$C182n6c
  mww_fa_NAp2$d5d <- mww_fa_NAp2$C204n6/mww_fa_NAp2$C203n6
  
  ##
  summary(mww_fa_NAp2$d5d)
  summary(mww_fa_NAp2$d6d)
  summary(mww_fa_NAp2$scd16)
  summary(mww_fa_NAp2$scd18)
  mww_fa_NAp2$d5d[is.infinite(mww_fa_NAp2$d5d)] <- NA
  ##

  
  #导入认知和炎症数据
  mww_ci <- read.csv("mww_cog_imf.csv",stringsAsFactors = F)
  mww <- left_join(mww_fa_NAp2,mww_ci,by="id")
  
  hist(mww$scd16)
  hist(log(mww$scd16))
  mww$logscd16 <- log(mww$scd16)
  mww$sqrtscd18 <- sqrt(mww$scd18)
  mww$logd5d <- log(mww$d5d)
  mww$logd6d <- log(mww$d6d)
  ##d6d存在较多0值
  mww$d6d[mww$d6d==0] <- NA
  mww$logd6d <- log(mww$d6d)
  mww$logIL10 <- log(mww$IL10)
  mww$logmoca <- log(mww$moca_total)
  ##TNF去除0值
  mww$TNFa[mww$TNFa==0] <- NA
  mww$logtnf <- log(mww$TNFa)
  ##IL1b去除0值
  mww$IL1β[mww$IL1β==0] <- NA
  mww$logil1b <- log(mww$IL1β)
  ##
  mww$NFkB[mww$NFkB==0] <- NA
  mww$lognfkb <- log(mww$NFkB)
  
  mww_mplus <- mww
  mww_mplus[is.na(mww_mplus)] <- 999
  write.csv(mww_mplus,file = "mplus/mww_mplus.csv",row.names = F)
  
  ###mediation
  #####Y为mmse_total，M为IL1B，X为desaturase
  
  ########scd16
  model_scd16<- ' # direct effect
             mmse_total ~ c*logscd16
  # mediator
  logil1b ~ a*logscd16
  mmse_total ~ b*logil1b
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit_scd16 <- sem(model_scd16, data = mww)
  summary(fit_scd16) 
  
  ########scd18
  model_scd18<- ' # direct effect
  mmse_total ~ c*sqrtscd18+age+gender
  # mediator
  logil1b ~ a*sqrtscd18+age+gender
  mmse_total ~ b*logil1b
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit_scd18 <- sem(model_scd18, data = mww)
  summary(fit_scd18) 
  
  
  hist(mww$logil1b)
  summary(mww$logil1b)
  summary(mww$IL1β)
  summary(mww$logtnf)
  hist(mww$IL10)
  hist(log(mww$IL10))
  hist(mww$IL1β)
  hist(mww$mmse_total)
  hist(log(mww$mmse_total))
  hist(sqrt(mww$mmse_total))
  hist(mww$moca_total)
  hist(log(mww$moca_total))
  hist(sqrt(mww$moca_total))
  
  
  
  summary(mww$logIL10)
  summary(mww$d6d)
  summary(mww$logd6d)
  summary(mww$logscd16)
  summary(mww$sqrtscd18)
  summary(mww$logd5d)
  summary(mww$logd6d)
  hist(mww$scd18)
  hist(log(mww$scd18))
  hist(sqrt(mww$scd18))
  hist(mww$d5d)
  hist(log(mww$d5d))
  hist(mww$d6d)
  hist(log(mww$d6d))
  
  #rescale FA 
  id <- mww$id
  fa_scale <- scale(mww[,4:17],center = F)
  fa_scale <- as.data.frame(fa_scale)
  fa <- mww[,4:17]
  fa_prin <- princomp(fa,scores = T,cor = T)
  fascale_prin <- princomp(fa_scale,scores = T,cor = T)
  
  mww_noFA <- dplyr::select(mww,-c(4:17))
  mww_scale <- left_join(mww_noFA,fa_scale,by="id")
  mww_mplus <- mww_scale
  mww_mplus[is.na(mww_mplus)] <- 999
  write.csv(mww_mplus,"mplus/mww_mplus.csv",row.names = F)
  
  #pricomp
  mww_s.prin <- princomp(fa_scale,scores = T,cor = T)
  screeplot(mww_s.prin,type = "lines")
  mww_s.prin$loadings
  mww_prin4 <- as.data.frame(mww_s.prin$scores)
  mww_prin4$id <- id
  mww_prin4 <- dplyr::select(mww_prin4,id,Comp.1,Comp.2,Comp.3,Comp.4)
  mww_prin <- left_join(mww_prin4,mww_ci,by="id")
  mww_prin_ <- mww_prin
  mww_prin_[is.na(mww_prin_)] <- 0
  mww_prin_ <- mww_prin_%>%dplyr::select(-c(30:32))
  comp <- mww_prin[2:5]
  inf <- mww_prin[7:11]
  cog <- mww_prin[16:17]
corp <- rcorr(as.matrix(mww_prin_))

model_11 <- ' # direct effect
             MoCANaming ~ c*C180
# mediator
lognfkb ~ a*C180
MoCANaming ~ b*lognfkb
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit11 <- sem(model_11, data = mww)
summary(fit11) 



mww_ <- mww
  mww_[is.na(mww_)] <- 0
  corp <- rcorr(as.matrix(mww_))
  pp <- as.data.frame(corp$P)
  ##cpeptide
  hist(log(mww_prin$cpep))
  mww_prin$logcpep <- log(mww_prin$cpep)
  ##IL10
  hist(log(mww_prin$IL10))
  
  mww_prin$logIL10 <- log(mww_prin$IL10)
  
  hist(mww_prin$logIL10)
  ##mmse_total
  hist(mww_prin$mmse_total)
  hist(mww_prin$moca_total)
  hist(log(mww_prin$moca_total))
  hist(mww_prin$Comp.1)
  summary(mww_prin$IL1β)
  frequency(mww_prin$IL1β)
  
  model_11 <- ' # direct effect
             mmse_total ~ c*Comp.4
  # mediator
  logIL10 ~ a*Comp.4
  mmse_total ~ b*logIL10
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit11 <- sem(model_11, data = mww_prin)
  summary(fit11) 
  
  #mediation
  library(lavaan)
 ##TNF to mmse 
  model_11 <- ' # direct effect
             mmse_total ~ c*Comp.1
           # mediator
             TNFa ~ a*Comp.1
             mmse_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit11 <- sem(model_11, data = mww_prin)
  summary(fit11)
  
  model_21 <- ' # direct effect
             mmse_total ~ c*Comp.2
           # mediator
             TNFa ~ a*Comp.2
             mmse_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit21 <- sem(model_21, data = mww_prin)##显著
  summary(fit21)
  
  hist(mww_prin$logtnf)
  ##TNF去除0值
  mww_prin$TNFa[mww_prin$TNFa==0] <- NA
  mww_prin$logtnf <- log(mww_prin$TNFa)
  mww_prin_ <- mww_prin
  mww_prin$logtnf[is.na(mww_prin$logtnf)] <- 0
  mww_prin$log
  model_211 <- ' # direct effect
             mmse_total ~ c*Comp.2
  # mediator
  logtnf ~ a*Comp.2
  mmse_total ~ b*logtnf
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit211 <- sem(model_211, data = mww_prin)##显著
  summary(fit211)
  
  model_31 <- ' # direct effect
             mmse_total ~ c*Comp.3
           # mediator
             TNFa ~ a*Comp.3
             mmse_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit31 <- sem(model_31, data = mww_prin)
  summary(fit31)
  
  model_41 <- ' # direct effect
             mmse_total ~ c*Comp.4
           # mediator
             TNFa ~ a*Comp.4
             mmse_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit41 <- sem(model_41, data = mww_prin)
  summary(fit41)
  
  #TNF to moca
  model_12 <- ' # direct effect
             moca_total ~ c*Comp.1
           # mediator
             TNFa ~ a*Comp.1
             moca_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit12 <- sem(model_12, data = mww_prin)
  summary(fit12)
  
  model_22 <- ' # direct effect
             moca_total ~ c*Comp.2
           # mediator
             TNFa ~ a*Comp.2
             moca_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit22 <- sem(model_22, data = mww_prin)
  summary(fit22)
  
  model_32 <- ' # direct effect
             moca_total ~ c*Comp.3
           # mediator
             TNFa ~ a*Comp.3
             moca_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit32 <- sem(model_32, data = mww_prin)
  summary(fit32)
  
  model_42 <- ' # direct effect
             moca_total ~ c*Comp.4
           # mediator
             TNFa ~ a*Comp.4
             moca_total ~ b*TNFa
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit42 <- sem(model_42, data = mww_prin)
  summary(fit42)
  
  ##crp to mmse 
  model_13 <- ' # direct effect
             mmse_total ~ c*Comp.1
           # mediator
             crp ~ a*Comp.1
             mmse_total ~ b*crp
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit13 <- sem(model_13, data = mww_prin)
  summary(fit13)
  
  model_23 <- ' # direct effect
             mmse_total ~ c*Comp.2
  # mediator
  crp ~ a*Comp.2
  mmse_total ~ b*crp
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit23 <- sem(model_23, data = mww_prin)#
  summary(fit23)
  
  model_33 <- ' # direct effect
             mmse_total ~ c*Comp.3
  # mediator
  crp ~ a*Comp.3
  mmse_total ~ b*crp
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit33 <- sem(model_33, data = mww_prin)#
  summary(fit33)
  
  model_43 <- ' # direct effect
             mmse_total ~ c*Comp.4
  # mediator
  crp ~ a*Comp.4
  mmse_total ~ b*crp
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit43 <- sem(model_43, data = mww_prin)#
  summary(fit43)
  
  #IL1β to mmse
  model_14 <- ' # direct effect
             mmse_total ~ c*Comp.1
           # mediator
             IL1β ~ a*Comp.1
             mmse_total ~ b*IL1β
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
  fit14 <- sem(model_14, data = mww_prin)
  summary(fit14)
  
  model_24 <- ' # direct effect
             mmse_total ~ c*Comp.2
  # mediator
  IL1β ~ a*Comp.2
  mmse_total ~ b*IL1β
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit24 <- sem(model_24, data = mww_prin)
  summary(fit24)
  
  model_34 <- ' # direct effect
             log(mmse_total) ~ c*Comp.3
  # mediator
  log(IL1β) ~ a*Comp.3
  log(mmse_total) ~ b*log(IL1β)
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit34 <- sem(model_34, data = mww_prin)
  summary(fit34)
  hist(log(mww_prin$IL1β))
  hist(sqrt(mww_prin$mmse_total))
  hist(mww_prin$mmse_total)
  hist(log(mww_prin$moca_total))
  
  
  model_34 <- ' # direct effect
             log(moca_total) ~ c*Comp.3
  # mediator
  log(IL1β) ~ a*Comp.3
  log(moca_total) ~ b*log(IL1β)
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  fit34 <- sem(model_34, data = mww_prin)
  summary(fit34)
  