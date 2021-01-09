library(tidyverse)
drink <- read.csv("drinking.csv",stringsAsFactors = F)
mww_pca <- read.csv(file = "mww_pca.csv",stringsAsFactors = F)
mww_pca <- left_join(mww_pca,drink,by="id")
mww_pca$d5d[is.infinite(mww_pca$d5d)] <- NA
mww_pca$d6d[is.infinite(mww_pca$d6d)] <- NA
mww_pca$scd16[is.infinite(mww_pca$scd16)] <- NA
mww_pca$scd18[is.infinite(mww_pca$scd18)] <- NA

mww_cov <- mww_pca %>% 
  dplyr::select(age,gender,bmi,energy,smoke)#混杂因子
mww_fa <- mww_pca%>%
  dplyr::select(c(8:25))
mww_cog <- mww_pca%>%
  dplyr::select(mmse_total,moca_total,c(39:49))
mww_inf <- mww_pca %>%
  dplyr::select(IL1β,IL10,NFkB,TNFa,crp)

mww_mplus <- mww_pca
mww_mplus[is.na(mww_mplus)] <- 999
write.csv(mww_mplus,file = "mww mplus for moderation.csv")


## Partial correlation between plasma FA and MMSE&MOCA 
library(Hmisc)
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
    mod1.fc <- lm(frma1.fc,data = mww_pca,na.action = na.exclude)
    r1.fc <- resid(mod1.fc)#residuls for mod1
    
    frma2.fc <- as.formula(paste(xhs.fc,cov_f))
    mod2.fc <- lm(frma2.fc,data = mww_pca,na.action = na.exclude)
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
    mod1.fi <- lm(frma1.fi,data = mww_pca,na.action = na.exclude)
    r1.fi <- resid(mod1.fi)#residuls for mod1
    
    frma2.fi <- as.formula(paste(xhs.fi,cov_f))
    mod2.fi <- lm(frma2.fi,data = mww_pca,na.action = na.exclude)
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
    mod1.ic <- lm(frma1.ic,data = mww_pca,na.action = na.exclude)
    r1.ic <- resid(mod1.ic)#residuls for mod1
    
    frma2.ic <- as.formula(paste(xhs.ic,cov_f))
    mod2.ic <- lm(frma2.ic,data = mww_pca,na.action = na.exclude)
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

##Export P value file to csv
write.csv(PCorMat.fc,file ="p value between FA and Cog.csv")
write.csv(PCorMat.fi,file ="p value between FA and Inf.csv")
write.csv(PCorMat.ic,file ="p value between inf and cog.csv")


#Mediation
library(lavaan)
model_11 <- ' # direct effect
             MoCAOrientation ~ c*C180+age+gender+bmi+energy+smoke
# mediator
TNFa ~ a*C180
MoCAOrientation ~ b*TNFa
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'

fit11 <- sem(model_11, data = mww_pca)
summary(fit11) 

##################C16+C18四分位
library(Hmisc
        )
mww_pca <- mww_pca %>% 
  dplyr::mutate(s4q=C160+C180)
mww_pca$q4 <- ntile(mww_pca$s4q,3)

tapply(mww_pca$MoCAOrientation,mww_pca$q4, mean,na.rm=T)
cor(mww_pca$MoCAOrientation,mww_pca$s4q,use = "na.or.complete")
