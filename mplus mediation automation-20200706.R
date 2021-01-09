library(MplusAutomation)

createModels("C:/Users/likai/OneDrive/OneDrive - Business/首医/数据/mww/mww/mplus/auto2/temp.inp")


#批量计算，批量保存结果，批量整理结果
##批量计算
runModels("C:/Users/likai/OneDrive/OneDrive - Business/首医/数据/mww/mww/mplus/auto")
allresults <- readModels("C:/Users/likai/OneDrive/OneDrive - Business/首医/数据/mww/mww/mplus/auto")
total.e <- double()#总效应值
total.p <- double()#总效应p
indirect.e <- double()
indirect.p <- double()
direct.e <- double()
direct.p <- double()
failmodel <- character()
modelnames <- names(allresults)
model.fine <- character()
for (i in 1:length(modelnames)) {
 if (length(allresults[[c(i,11)]])>0) {#长度为零的模型均为中介因子为d6d和scd16
    cc <- allresults[[c(i,11)]]#11为中介模型结果
    model.fine[length(model.fine)+1] <- modelnames[i]
    total.e[length(total.e)+1] <- cc$unstandardized$overall$est[1]
    indirect.e[length(indirect.e)+1]<- cc$unstandardized$overall$est[2]
    direct.e[length(direct.e)+1]<- cc$unstandardized$overall$est[3]
    total.p[length(total.p)+1] <- cc$unstandardized$overall$pval[1]
    indirect.p[length(indirect.p)+1]<- cc$unstandardized$overall$pval[2]
    direct.p[length(direct.p)+1]<- cc$unstandardized$overall$pval[3]
  }
}

mediation.results <- data.frame(model.fine,total.e,total.p,indirect.e,indirect.p,direct.e,direct.p)


#更改了scd16和d6d的数量级，重新批量计算，批量保存结果，批量整理结果
##批量计算
runModels("C:/Users/likai/OneDrive/OneDrive - Business/首医/数据/mww/mww/mplus/auto2")
allresult2 <- readModels("C:/Users/likai/Desktop/auto2")
total.e <- double()#总效应值
total.p <- double()#总效应p
indirect.e <- double()
indirect.p <- double()
direct.e <- double()
direct.p <- double()
failmodel <- character()
modelnames <- names(allresult2)
model.fine <- character()
X <- character()
M <- character()
Y <- character()
for (i in 1:length(modelnames)) {
  if (length(allresult2[[c(i,11)]])>0) {#长度为零的模型均为中介因子为d6d和scd16
    cc <- allresult2[[c(i,11)]]#11为中介模型结果
    X[length(X)+1] <- as.character(cc$unstandardized$specific$pred)
    M[length(M)+1] <- as.character(cc$unstandardized$specific$intervening)
    Y[length(Y)+1] <- as.character(cc$unstandardized$specific$outcome)
    total.e[length(total.e)+1] <- cc$unstandardized$overall$est[1]
    indirect.e[length(indirect.e)+1]<- cc$unstandardized$overall$est[2]
    direct.e[length(direct.e)+1]<- cc$unstandardized$overall$est[3]
    total.p[length(total.p)+1] <- cc$unstandardized$overall$pval[1]
    indirect.p[length(indirect.p)+1]<- cc$unstandardized$overall$pval[2]
    direct.p[length(direct.p)+1]<- cc$unstandardized$overall$pval[3]
  }
  
}


mediation.result2 <- data.frame(X,M,Y,total.e,total.p,indirect.e,indirect.p,direct.e,direct.p)
write.csv(mediation.result2,file = "mediation FA insulin Cog.csv",row.names = F)
save(mediation.result2,file = )


#炎症因子也做一遍
createModels("C:/Users/likai/Desktop/auto3/temp.inp")
runModels("C:/Users/likai/Desktop/auto3")
infresult <- readModels("C:/Users/likai/Desktop/auto3")

total.e <- double()#总效应值
total.p <- double()#总效应p
indirect.e <- double()
indirect.p <- double()
direct.e <- double()
direct.p <- double()
failmodel <- character()
modelnames <- names(infresult)
model.fine <- character()
Xinf <- character()
Minf <- character()
Yinf <- character()

for (i in 1:length(modelnames)) {
  if (length(infresult[[c(i,11)]])>0) {
    cc <- infresult[[c(i,11)]]#11为中介模型结果
    Xinf[length(Xinf)+1] <- as.character(cc$unstandardized$specific$pred)
    Minf[length(Minf)+1] <- as.character(cc$unstandardized$specific$intervening)
    Yinf[length(Yinf)+1] <- as.character(cc$unstandardized$specific$outcome)
    total.e[length(total.e)+1] <- cc$unstandardized$overall$est[1]
    indirect.e[length(indirect.e)+1]<- cc$unstandardized$overall$est[2]
    direct.e[length(direct.e)+1]<- cc$unstandardized$overall$est[3]
    total.p[length(total.p)+1] <- cc$unstandardized$overall$pval[1]
    indirect.p[length(indirect.p)+1]<- cc$unstandardized$overall$pval[2]
    direct.p[length(direct.p)+1]<- cc$unstandardized$overall$pval[3]
  }
}

mediation.result.inf <- data.frame(Xinf,Minf,Yinf,total.e,total.p,indirect.e,indirect.p,direct.e,direct.p)
write.csv(mediation.result.inf,file = "mediation FA inflammation Cog.csv",row.names = F)
