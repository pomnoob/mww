#计算DBI16##

#发现脂肪酸数据中王秀云有两条记录，对比发现可能是重名患者，分别标记为王秀云1、王秀云2
#又发现膳食调查数据中王秀云有两条记录，处理方法同上；黑振成也有两条记录，其中一条热量摄入仅为500多kcal
#因此对此条进行了删除操作
#除外，膳食调查数据中有许多记录的登记号有问题（均为13438），需要用脂肪酸数据中的登记号对其进行校正。
diet <- read.csv(file = "DBI/脂肪酸整理-李亚茹-190917.csv",stringsAsFactor=F)
diet_id <- diet$登记号
for (i in 1:288) {
  name <- diet_id[-i]
  name2 <- diet_id[i]
  if (name2%in%name){
    print(name2)
  }
}
#筛选出登记号有问题的患者
name_id <- diet %>%
  dplyr::select(X,姓名)

#导入脂肪酸数据中正确的登记号
FA <- read.csv(file = "20200108 数据库(放了血液脂肪酸-重新核对脂肪酸正确版-去掉359兰孝明重复的+计算了认知的小分).csv",stringsAsFactors = F)
FA_name <- FA%>%
  dplyr::select(姓名,登记号)
name_id <- inner_join(name_id,FA_name,by="姓名")#大概有10个人脂肪酸文件中没有，排除掉
#导出正确的姓名和对应的登记号
write.csv(name_id,file = "DBI/正确的姓名和登记号.csv",row.names=F)
#仅使用完整的数据，排除掉膳食数据中原有的10个人，因其无其他数据，最后合并成用于分析的膳食数据
diet2 <- dplyr::select(diet,-X,-登记号)
diet_use <- inner_join(name_id,diet2,by="姓名")
diet_use <- dplyr::rename(diet_use,id=登记号)
write.csv(diet_use,file = "DBI/膳食调查数据202005289-正确的姓名+登记号-完整其他数据.csv",row.names=F)

#第一部分 Cereal
#Cereal include rice, wheat , dried legumes(exclude soybean) and tubers. Intake amount means fresh amount. Sweat potato: intake amount divided by 
#3; potato: intake amount divided by 4; yam and yambean: divided by 6; score increased (decreased) 2 with 15g intake increased (decreased) when the 
#energy intake level is 1000 and 1200; score increased (decreased) 2 with 25g intake increased (decreased) when the energy intake level is 
#1400 kcal; 1 with 15g intake for 1600 and 1800 kcal; 1 with 20g for 2000 and 2200 kcal; and 1 with 25g for the energy intake level more than 
#2400 kcal.

#cereal的量等于生的大米、小麦、干豆的量相加；红薯的摄入量除以3，土豆摄入量除以4，芋头摄入量除以6
#实际，本研究所用的FFQ问卷没有区分土豆和红薯，只有一个”薯类“摄入（包含红薯和土豆）；芋头和胡萝卜、白萝卜
#放在一起进行问卷（根茎）。因此计算的时候折中一下：红薯+土豆=薯类/3.5；芋头=根茎/6
library(tidyverse)
diet_use <- read.csv(file = "DBI/膳食调查数据202005289-正确的姓名+登记号-完整其他数据.csv",stringsAsFactor=F)
diet_use[is.na(diet_use)] <- 0##有些受试者膳食有缺失值，全部按照0处理
dbi <- diet_use %>%
  dplyr::mutate(cereal=(大米.两.+面粉.两.+包馅食品.两.+杂粮.两.+油条.两.
                +油饼.两.)*50+方便面.克.+薯类.两.*50/3.5+根茎.两.*50/6,veg=50*(叶菜.两.+
                  包菜.两.+瓜茄.两.+果菜.两.+花菜.两.+菌类.两.),fruit=水果.克.+干果.克.,meat=
                  50*(瘦猪肉.两.+红烧肉.两.+排骨.两.+五花肉.两.+牛肉.两.+牛腩.两.+羊肉.两.+卤煮.两.+
                  炒肝.两.+腰花.两.+大肠.两.+肝脏.两.+肺片.两.+火腿肠.两.+鸡肉.两.+鸡腿.两.+鸡翅.两.),fish=
                  50*(鱼肉.两.+虾.两.),egg=50*鸡蛋.个.,dairy=50*(牛奶.两.+酸奶.克.),soy=豆腐.两.*50/212.5+豆浆.ml.*50/14.6#此处的换算按照2016膳食指南63页的图片而来，豆腐取南北豆腐平均值
                ,oil=食用油.克.,salt=食盐.克.,alcohol=0,sugar=(碳酸.ml.+果汁.ml.)*0.05) %>%
  dplyr::select(id,姓名,cereal,veg,fruit,meat,fish,egg,dairy,soy,salt,oil,alcohol,sugar)

FA <- read.csv(file = "20200108 数据库(放了血液脂肪酸-重新核对脂肪酸正确版-去掉359兰孝明重复的+计算了认知的小分).csv",stringsAsFactors = F)
FA_id <- FA %>% 
  dplyr::select(姓名,登记号,能量,性别) %>%
  dplyr::rename(name=姓名,id=登记号,energy=能量,gender=性别)

dbi <- inner_join(dbi,FA_id,by="id")

##Cereal

dbi$dbi1 <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal==0~-12,
                                   cereal>=75 & cereal<=95~0,
                                   cereal>170~12,
                                   cereal>0 & cereal <75~2*ceiling(cereal/15)-12,#向上取整
                                   cereal>95 & cereal<=170~2*ceiling((cereal-95)/15)))
                    
  } else if(dbi$energy[i]<=1200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<15~-12,
                                 cereal>=90 & cereal<=110~0,
                                 cereal>185~12,
                                 cereal>=15 & cereal <90~2*ceiling((cereal-15)/15)-12,
                                 cereal>110 & cereal<=185~2*ceiling((cereal-110)/15)))
    
  } else if(dbi$energy[i]<=1400 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal==0~-12,
                                 cereal>=125 & cereal<=175~0,
                                 cereal>300~12,
                                 cereal>0 & cereal <125~2*ceiling(cereal/25)-12,
                                 cereal>175 & cereal<=300~2*ceiling((cereal-175)/25)))
    
  } else if(dbi$energy[i]<=1600 & dbi$energy[i]>1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<10~-12,
                                 cereal>=175 & cereal<=225~0,
                                 cereal>390~12,
                                 cereal>=10 & cereal <175~ceiling((cereal-10)/15)-12,
                                 cereal>225 & cereal<=390~ceiling((cereal-225)/15)))
    
  }else if(dbi$energy[i]<=1800 & dbi$energy[i]>1600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<35~-12,
                                 cereal>=200 & cereal<=250~0,
                                 cereal>415~12,
                                 cereal>=35 & cereal <200~ceiling((cereal-35)/15)-12,
                                 cereal>250 & cereal<=415~ceiling((cereal-250)/15)))
  }else if(dbi$energy[i]<=2000 & dbi$energy[i]>1800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<5~-12,
                                 cereal>=225 & cereal<=275~0,
                                 cereal>495~12,
                                 cereal>=5 & cereal <225~ceiling((cereal-5)/20)-12,
                                 cereal>275 & cereal<=495~ceiling((cereal-275)/20)))
  }else if(dbi$energy[i]<=2200 & dbi$energy[i]>2000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<30~-12,
                                 cereal>=250 & cereal<=300~0,
                                 cereal>520~12,
                                 cereal>=30 & cereal <250~ceiling((cereal-30)/20)-12,
                                 cereal>300 & cereal<=520~ceiling((cereal-300)/20)))
  }else if(dbi$energy[i]<=2400 & dbi$energy[i]>2200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal==0~-12,
                                 cereal>=275 & cereal<=325~0,
                                 cereal>600~12,
                                 cereal>=0 & cereal <275~ceiling(cereal/25)-12,
                                 cereal>325 & cereal<=600~ceiling((cereal-325)/25)))
    
  }else if(dbi$energy[i]<=2600 & dbi$energy[i]>2400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<50~-12,
                                 cereal>=325 & cereal<=375~0,
                                 cereal>650~12,
                                 cereal>=50 & cereal <325~ceiling((cereal-50)/25)-12,
                                 cereal>375 & cereal<=650~ceiling((cereal-375)/25)))
    
  }else if(dbi$energy[i]<=2800 & dbi$energy[i]>2600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<75~-12,
                                 cereal>=350 & cereal<=400~0,
                                 cereal>675~12,
                                 cereal>=75 & cereal <350~ceiling((cereal-75)/20)-12,
                                 cereal>400 & cereal<=675~ceiling((cereal-400)/20)))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,]<- dbi[i,] %>% 
      dplyr::mutate(dbi1=case_when(cereal<100~-12,
                                 cereal>=375 & cereal<=425~0,
                                 cereal>700~12,
                                 cereal>=100 & cereal <375~ceiling((cereal-100)/20)-12,
                                 cereal>425 & cereal<=700~ceiling((cereal-425)/20)))
  } 

  }


#计算vegetable和fruit的分数


#Vegetable
dbi$dbi_veg <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=200~0,
                                   veg==0~-6,
                                   veg>=160 & veg<200~-1,
                                   veg>0 & veg <160~-1-ceiling(veg/40)))#向上取整
                                   
    
  } else if(dbi$energy[i]<=1200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=250~0,
                                      veg==0~-6,
                                      veg>=200 & veg<250~-1,
                                      veg>0 & veg <200~-1-ceiling(veg/50)))
    
  } else if(dbi$energy[i]<=1600 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=300~0,
                                      veg==0~-6,
                                      veg>=240 & veg<300~-1,
                                      veg>0 & veg <240~-1-ceiling(veg/60)))
    
  } else if(dbi$energy[i]<=1800 & dbi$energy[i]>1600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=400~0,
                                      veg==0~-6,
                                      veg>=320 & veg<400~-1,
                                      veg>0 & veg <320~-1-ceiling(veg/80)))
    
  }else if(dbi$energy[i]<=2200 & dbi$energy[i]>1800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=450~0,
                                      veg==0~-6,
                                      veg>=360 & veg<450~-1,
                                      veg>0 & veg <360~-1-ceiling(veg/90)))
    
  }else if(dbi$energy[i]<=2800 & dbi$energy[i]>2200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=500~0,
                                      veg==0~-6,
                                      veg>=400 & veg<500~-1,
                                      veg>0 & veg <400~-1-ceiling(veg/100)))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_veg=case_when(veg>=600~0,
                                      veg==0~-6,
                                      veg>=480 & veg<600~-1,
                                      veg>0 & veg <480~-1-ceiling(veg/120)))
    
  }
}

#Fruit
dbi$dbi_fru <- 0
for (i in 1:278) {
  if (dbi$energy[i]<=1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fru=case_when(fruit>=150~0,
                                      fruit==0~-6,
                                      fruit>=120 & fruit<150~-1,
                                      fruit>0 & fruit <120~-1-ceiling(fruit/30)))#向上取整
    
    
  } else if(dbi$energy[i]<=1800 & dbi$energy[i]>1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fru=case_when(fruit>=200~0,
                                      fruit==0~-6,
                                      fruit>=160 & fruit<200~-1,
                                      fruit>0 & fruit <160~-1-ceiling(fruit/40)))
    
  } else if(dbi$energy[i]<=2200 & dbi$energy[i]>1800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fru=case_when(fruit>=300~0,
                                      fruit==0~-6,
                                      fruit>=240 & fruit<300~-1,
                                      fruit>0 & fruit <240~-1-ceiling(fruit/60)))
    
  } else if(dbi$energy[i]<=2800 & dbi$energy[i]>2200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fru=case_when(fruit>=350~0,
                                      fruit==0~-6,
                                      fruit>=280 & fruit<350~-1,
                                      fruit>0 & fruit <280~-1-ceiling(fruit/70)))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fru=case_when(fruit>=400~0,
                                      fruit==0~-6,
                                      fruit>=320 & fruit<400~-1,
                                      fruit>0 & fruit <320~-1-ceiling(fruit/80)))
    
  }
}

#DBI2
dbi$dbi2 <- dbi$dbi_fru+dbi$dbi_veg

#计算肉类、鱼虾、鸡蛋的分数
# meat为畜肉、鸡肉摄入总量，fish为鱼虾摄入总量

dbi$dbi_mp <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-3,
                                     meat>0 & meat<=5~-2,
                                     meat>5 & meat<=10~-1,
                                     meat>10 & meat<=20~0,
                                     meat>20 & meat<=25~1,
                                     meat>25 & meat<=30~2,
                                     meat>30 & meat<=35~3,
                                     meat>35~4))

  } else if(dbi$energy[i]<=1200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-4,
                                     meat>0 & meat<=5~-3,
                                     meat>5 & meat<=10~-2,
                                     meat>10 & meat<=15~-1,
                                     meat>15 & meat<=35~0,
                                     meat>35 & meat<=40~1,
                                     meat>40 & meat<=45~2,
                                     meat>45 & meat<=50~3,
                                     meat>50~4))
    
  } else if(dbi$energy[i]<=1600 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-4,
                                     meat>0 & meat<=10~-3,
                                     meat>10 & meat<=20~-2,
                                     meat>20 & meat<=30~-1,
                                     meat>30 & meat<=50~0,
                                     meat>50 & meat<=60~1,
                                     meat>60 & meat<=70~2,
                                     meat>70 & meat<=80~3,
                                     meat>80~4))
    
  }else if(dbi$energy[i]<=2000 & dbi$energy[i]>1600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-4,
                                     meat>0 & meat<=15~-3,
                                     meat>15 & meat<=30~-2,
                                     meat>30 & meat<=45~-1,
                                     meat>45 & meat<=55~0,
                                     meat>55 & meat<=70~1,
                                     meat>70 & meat<=85~2,
                                     meat>85 & meat<=100~3,
                                     meat>100~4))
    }else if(dbi$energy[i]<=2800 & dbi$energy[i]>2000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-4,
                                     meat>0 & meat<=20~-3,
                                     meat>20 & meat<=40~-2,
                                     meat>40 & meat<=60~-1,
                                     meat>60 & meat<=90~0,
                                     meat>90 & meat<=110~1,
                                     meat>110 & meat<=130~2,
                                     meat>130 & meat<=150~3,
                                     meat>150~4))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_mp=case_when(meat==0~-4,
                                     meat>0 & meat<=25~-3,
                                     meat>25 & meat<=50~-2,
                                     meat>50 & meat<=75~-1,
                                     meat>75 & meat<=125~0,
                                     meat>125 & meat<=150~1,
                                     meat>150 & meat<=175~2,
                                     meat>175 & meat<=200~3,
                                     meat>200~4))
    
  }
}

#Fish & shell
dbi$dbi_fs <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish==0~-4,
                                     fish>0 & fish<5~-3,
                                     fish>=5 & fish<10~-2,
                                     fish>=10 & fish<15~-1,
                                     fish>=15~0))
    
  } else if(dbi$energy[i]<=1200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish<5~-4,
                                     fish>=5 & fish<10~-3,
                                     fish>=10 & fish<15~-2,
                                     fish>=15 & fish<20~-1,
                                     fish>=20~0))
    
  } else if(dbi$energy[i]<=1600 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish<10~-4,
                                     fish>=10 & fish<20~-3,
                                     fish>=20 & fish<30~-2,
                                     fish>=30 & fish<40~-1,
                                     fish>=40~0))
    
  }else if(dbi$energy[i]<=2000 & dbi$energy[i]>1600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish<5~-4,
                                     fish>=5 & fish<20~-3,
                                     fish>=20 & fish<35~-2,
                                     fish>=35 & fish<50~-1,
                                     fish>=50~0))
    
  }else if(dbi$energy[i]<=2600 & dbi$energy[i]>2000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish==0~-4,
                                     fish>0 & fish<25~-3,
                                     fish>=25 & fish<50~-2,
                                     fish>=50 & fish<75~-1,
                                     fish>=75~0))
    
  }else if(dbi$energy[i]<=2800 & dbi$energy[i]>2600){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish<25~-4,
                                     fish>=25 & fish<50~-3,
                                     fish>=50 & fish<75~-2,
                                     fish>=75 & fish<100~-1,
                                     fish>=100~0))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_fs=case_when(fish<50~-4,
                                     fish>=50 & fish<75~-3,
                                     fish>=75 & fish<100~-2,
                                     fish>=100 & fish<125~-1,
                                     fish>=125~0))
  }}

#egg
dbi$dbi_egg <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_egg=case_when(egg==0~-4,
                                      egg>0 & egg<=5~-3,
                                      egg>5 & egg<=10~-2,
                                      egg>10 & egg<=15~-1,
                                      egg>15 & egg<=25~0,
                                      egg>25 & egg<=30~1,
                                      egg>30 & egg<=35~2,
                                      egg>35 & egg<=40~3,
                                      egg>40~4))
    
  } else if(dbi$energy[i]<=1400 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_egg=case_when(egg<5~-4,
                                      egg>=5 & egg<=10~-3,
                                      egg>10 & egg<=15~-2,
                                      egg>15 & egg<=20~-1,
                                      egg>20 & egg<=30~0,
                                      egg>30 & egg<=35~1,
                                      egg>35 & egg<=40~2,
                                      egg>40 & egg<=45~3,
                                      egg>45~4))
    
  } else if(dbi$energy[i]<=1800 & dbi$energy[i]>1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_egg=case_when(egg==0~-4,
                                      egg>0 & egg<=10~-3,
                                      egg>10 & egg<=20~-2,
                                      egg>20 & egg<=30~-1,
                                      egg>30 & egg<=50~0,
                                      egg>50 & egg<=60~1,
                                      egg>60 & egg<=70~2,
                                      egg>70 & egg<=80~3,
                                      egg>80~4))
    
  }else if(dbi$energy[i]>1800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_egg=case_when(egg==0~-4,
                                      egg>0 & egg<=15~-3,
                                      egg>15 & egg<=30~-2,
                                      egg>30 & egg<=45~-1,
                                      egg>45 & egg<=55~0,
                                      egg>55 & egg<=70~1,
                                      egg>70 & egg<=85~2,
                                      egg>85 & egg<=100~3,
                                      egg>100~4))
    
  }}

##DBI4
dbi$dbi4 <- dbi$dbi_mp+dbi$dbi_fs+dbi$dbi_egg

#Dairy and soybean
##Dairy
dbi$dbi_dairy <- 0
for (i in 1:278) {
  if (dbi$energy[i]<=1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_dairy=case_when(dairy>=500~0,
                                        dairy==0~-6,
                                        dairy>0 & dairy<500~-6+ceiling(dairy/100)))
    
  } else if(dbi$energy[i]<=1400 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_dairy=case_when(dairy>=350~0,
                                        dairy==0~-6,
                                        dairy>0 & dairy<350~-6+ceiling(dairy/70)))
    
  } else if(dbi$energy[i]>1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_dairy=case_when(dairy>=300~0,
                                        dairy==0~-6,
                                        dairy>0 & dairy<300~-6+ceiling(dairy/60)))
    
  }}

##Soybean
dbi$dbi_soy <- 0
for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_soy=case_when(soy>=5~0,
                                      soy==0~-6,
                                      soy>0 & soy<5~-6+ceiling(soy)))
    
  } else if(dbi$energy[i]<=2000 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_soy=case_when(soy>=15~0,
                                      soy==0~-6,
                                      soy>0 & soy<15~-6+ceiling(soy/3)))
    
  } else if(dbi$energy[i]>2000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_soy=case_when(soy>=25~0,
                                      soy==0~-6,
                                      soy>0 & soy<25~-6+ceiling(soy/5)))
    
  }}

###dbi3
dbi$dbi3 <- dbi$dbi_dairy+dbi$dbi_soy

#Cooking oil and alcohol
##Cooking oil
dbi$dbi_oil <- 0
for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_oil=case_when(oil<=20~0,
                                      oil>20 & oil<=25~1,
                                      oil>45~6,
                                      oil>25 & oil<=45~1+ceiling((oil-25)/5)))
    
  } else if(dbi$energy[i]<=2200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_oil=case_when(oil<=25~0,
                                      oil>25 & oil<=30~1,
                                      oil>50~6,
                                      oil>30 & oil<=50~1+ceiling((oil-30)/5)))
    
  } else if(dbi$energy[i]<=2800 & dbi$energy[i]>2200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_oil=case_when(oil<=30~0,
                                      oil>30 & oil<=35~1,
                                      oil>55~6,
                                      oil>35 & oil<=55~1+ceiling((oil-35)/5)))
    
  }else if(dbi$energy[i]>2800){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_oil=case_when(oil<=35~0,
                                      oil>35 & oil<=40~1,
                                      oil>60~6,
                                      oil>40 & oil<=60~1+ceiling((oil-40)/5)))
  }}

##Alcohol
###调查问卷里面没有酒类摄入的数据，按照0摄入处理
dbi$dbi_alcohol <- 0
####dbi5
dbi$dbi5 <- dbi$dbi_oil+dbi$dbi_alcohol

#sugar and salt
##sugar:糖类摄入仅通过软饮料计算，因此不可避免地会低估人群中添加糖的摄入量
dbi$dbi_sugar <- 0
dbi <- dplyr::mutate(dbi,dbi_sugar=case_when(sugar<=25~0,
                                             sugar>50~6,
                                             sugar>25 & sugar<=50~ceiling((sugar-25)/5)))

##salt
dbi$dbi_salt <- 0

for (i in 1:278) {
  if (dbi$energy[i]<=1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_salt=case_when(salt<2~0,
                                      salt>=2 & salt<=3~1,
                                      salt>12~6,
                                      salt>3 & salt<=12~1+ceiling((salt-3)/2)))
    
  } else if(dbi$energy[i]<=1200 & dbi$energy[i]>1000){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_salt=case_when(salt<3~0,
                                       salt>=3 & salt<=4~1,
                                       salt>13~6,
                                       salt>4 & salt<=13~1+ceiling((salt-4)/2)))
    
  } else if(dbi$energy[i]<=1400 & dbi$energy[i]>1200){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_salt=case_when(salt<4~0,
                                       salt>=4 & salt<=5~1,
                                       salt>14~6,
                                       salt>5 & salt<=14~1+ceiling((salt-5)/2)))
    
  }else if(dbi$energy[i]>1400){
    dbi[i,] <- dbi[i,] %>% 
      dplyr::mutate(dbi_salt=case_when(salt<6~0,
                                       salt>=6 & salt<=7~1,
                                       salt>16~6,
                                       salt>7 & salt<=16~1+ceiling((salt-7)/2)))
  }}
####dbi6
dbi$dbi6 <- dbi$dbi_sugar+dbi$dbi_salt

#Diet variety
dbi_dv <- diet_use %>%
  dplyr::mutate(cereal=(大米.两.+面粉.两.+包馅食品.两.+杂粮.两.+油条.两.
                          +油饼.两.)*50+方便面.克.+薯类.两.*50/3.5+根茎.两.*50/6,veg=50*(叶菜.两.+
                                                                                包菜.两.+瓜茄.两.+果菜.两.+花菜.两.+菌类.两.),fruit=水果.克.+干果.克.,meat=
                  50*(瘦猪肉.两.+红烧肉.两.+排骨.两.+五花肉.两.+牛肉.两.+牛腩.两.+羊肉.两.+卤煮.两.+
                           炒肝.两.+腰花.两.+大肠.两.+肝脏.两.+肺片.两.+火腿肠.两.+鸡肉.两.+鸡腿.两.+鸡翅.两.),fish=
                  50*(鱼肉.两.+虾.两.),egg=50*鸡蛋.个.,dairy=50*(牛奶.两.+酸奶.克.),soy=豆腐.两.*50/212.5+豆浆.ml.*50/14.6#此处的换算按照2016膳食指南63页的图片而来，豆腐取南北豆腐平均值
                ,oil=食用油.克.,salt=食盐.克.,alcohol=0,sugar=(碳酸.ml.+果汁.ml.)*0.05) %>%
  dplyr::mutate(dv1=case_when(大米.两.*50>=25~0,大米.两.*50<25~-1),
                dv2=case_when((面粉.两.*50+方便面.克.+油饼.两.*50+油条.两.*50)>=25~0,(面粉.两.*50+方便面.克.+油饼.两.*50+油条.两.*50)<25~-1),
                dv3=case_when(杂粮.两.*50>=25~0,杂粮.两.*50<25~-1),
                dv4=case_when(蔬菜.两.*50>=25~0,蔬菜.两.*50<25~-1),
                dv5=case_when(蔬菜.两.*50>=25~0,蔬菜.两.*50<25~-1),
                dv6=case_when(fruit>=25~0,fruit<25~-1),
                dv7=case_when(soy>=5~0,soy<5~-1),
                dv8=case_when(dairy>=25~0,dairy<25~-1),
                dv9=case_when(50*(瘦猪肉.两.+红烧肉.两.+
                                  排骨.两.+五花肉.两.+牛肉.两.+
                                  牛腩.两.+羊肉.两.+卤煮.两.+炒肝.两.+腰花.两.+
                                  大肠.两.+肝脏.两.+肺片.两.+火腿肠.两.)>=25~0,50*(瘦猪肉.两.+红烧肉.两.+
                                                                           排骨.两.+五花肉.两.+牛肉.两.+
                                                                           牛腩.两.+羊肉.两.+卤煮.两.+炒肝.两.+腰花.两.+
                                                                           大肠.两.+肝脏.两.+肺片.两.+火腿肠.两.)<25~-1),
                dv10=case_when(50*(鸡肉.两.+鸡腿.两.+鸡翅.两.)>=25~0,50*(鸡肉.两.+鸡腿.两.+鸡翅.两.)<25~-1),
                dv11=case_when(egg>=25~0,egg<25~-1),
                dv12=case_when(fish>=25~0,fish<25~-1)) %>%
  dplyr::mutate(dbi_dv=rowSums(.[76:87],na.rm = T)) %>% dplyr::select(id,姓名,dbi_dv) %>% dplyr::rename(dbi7=dbi_dv)

dbi_score <- dbi %>%
  dplyr::select(id,dbi1,dbi2,dbi3,dbi4,dbi5,dbi6)
  
dbi_score <- inner_join(dbi_score,dbi_dv,by="id")
##分别计算正端分和负端分
dbi_score <- dbi_score %>%
  dplyr::mutate(dbi1_p=case_when(dbi1>=0~dbi1,dbi1<0~0),
                dbi1_n=case_when(dbi1<=0~-dbi1,dbi1>0~0),
                dbi2_p=case_when(dbi2>=0~dbi2,dbi2<0~0),
                dbi2_n=case_when(dbi2<=0~-dbi2,dbi2>0~0),
                dbi3_p=case_when(dbi3>=0~dbi3,dbi3<0~0),
                dbi3_n=case_when(dbi3<=0~-dbi3,dbi3>0~0),
                dbi4_p=case_when(dbi4>=0~dbi4,dbi4<0~0),
                dbi4_n=case_when(dbi4<=0~-dbi4,dbi4>0~0),
                dbi5_p=case_when(dbi5>=0~dbi5,dbi5<0~0),
                dbi5_n=case_when(dbi5<=0~-dbi5,dbi5>0~0),
                dbi6_p=case_when(dbi6>=0~dbi6,dbi6<0~0),
                dbi6_n=case_when(dbi6<=0~-dbi6,dbi6>0~0),
                dbi7_p=case_when(dbi7>=0~dbi7,dbi7<0~0),
                dbi7_n=case_when(dbi7<=0~-dbi7,dbi7>0~0)) %>%
  dplyr::mutate(lbs=dbi1_n+dbi2_n+dbi3_n+dbi4_n+dbi5_n+dbi6_n+dbi7_n,
                hbs=dbi1_p+dbi2_p+dbi3_p+dbi4_p+dbi5_p+dbi6_p+dbi7_p,
                dqd=lbs+hbs,
                ts=hbs-lbs)
#计算出来发现有5个人lbs为59分，仔细查看后发现这5个人膳食调查数据全部为0，因此从结果中剔除
dbi_score.final <- dbi_score %>%
  filter(lbs!=59) %>%
  select(id,姓名,dbi1:dbi7,lbs,hbs,dqd,ts)

write.csv(dbi_score.final,file = "DBI/DBI膳食质量分数.csv",row.names = F)



