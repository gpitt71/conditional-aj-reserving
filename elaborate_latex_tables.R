rm(list=ls())
# libraries ----

library(data.table)
library(xtable)
library(ggplot2)
library(dplyr)
library(tidyr)

divided.by = 1e+06

#

locn <- "conditional-aj-reserving\\results_csv\\%s"

# Simulations ----

fname_no_trend <- c('simulation_no_features_all_records_4_2023_10_23_11_57.csv',
'simulation_no_features_all_records_5_2023_10_23_11_58.csv',
'simulation_no_features_all_records_6_2023_10_23_12_00.csv',
'simulation_no_features_all_records_7_2023_10_23_15_16.csv')


dt <- data.table()

for(nm in fname_no_trend){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt <- rbind(dt,tmp)
  
  
  
}

dt[['features']] = "x"
dt[['scenario']] = "alpha"

dt[,aj.msep:=sqrt(variance)]


fname_trend <-c('simulation_aytrend4_2024_02_29_13_31.csv',
  'simulation_aytrend5_2024_02_29_13_32.csv',
  'simulation_aytrend6_2024_02_29_13_34.csv',
  'simulation_aytrend7_2024_02_29_13_36.csv')

dt1 <- data.table()

for(nm in fname_trend){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt1 <- rbind(dt1,tmp)}

dt1[['features']] = "v"
dt1[['scenario']] = "beta"


fname_trend_i <- c('simulation_aytrend_interceptmodel_all_records_4_2023_10_24_16_27.csv',
                   'simulation_aytrend_interceptmodel_all_records_5_2023_10_24_17_36.csv',
                   'simulation_aytrend_interceptmodel_all_records_6_2023_10_24_16_31.csv',
                   'simulation_aytrend_interceptmodel_all_records_7_2023_10_24_16_32.csv')
dt2 <- data.table()



for(nm in fname_trend_i){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt2 <- rbind(dt2,tmp)}

dt2[['aj.msep']]=sqrt(dt2[['variance']])
dt2[['features']] = "v2"
dt2[['scenario']] = "beta"

dt.tot <- rbind(dt[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","mcrps")],
                dt1[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","mcrps")],
                dt2[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","mcrps")])

dt.tot<-dt.tot[,.(actual.tot = mean(actual.tot)/divided.by,
          ei_aj=mean(pred.tot/actual.tot-1),
          ei_cl=mean(cl.tot/actual.tot-1),
          aj.msep=mean(aj.msep/pred.tot),
          cl.msep=mean(cl.msep/cl.tot),
          mcrps=mean(mcrps)/divided.by
),
by=.(k,
     scenario,
     features)][order(k),] 

# Table 1 ----

dt.tot[scenario=='alpha',]%>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)

# Table 2 ----

dt11<-dt1[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","mcrps")]
dt11[['relative_mcrps']] <- dt1[["mcrps"]]/dt1[["mcrps"]]
dt21 <- dt2[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","mcrps")]
dt21[['relative_mcrps']] <- dt2[["mcrps"]]/dt1[["mcrps"]]


dt.tot2 <- rbind(dt11[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","relative_mcrps")],
                 dt21[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features","scenario","relative_mcrps")])

dt.tot2[,.(actual.tot = mean(actual.tot)/divided.by,
          ei_aj=mean(pred.tot/actual.tot-1),
          ei_cl=mean(cl.tot/actual.tot-1),
          aj.msep=mean(aj.msep/pred.tot),
          cl.msep=mean(cl.msep/cl.tot),
          mcrps=mean(relative_mcrps)
),
by=.(k,
     scenario,
     features)][order(k),]%>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)


# Real data case studies ----

# Table 4 ----

df_features <- fread(sprintf(locn,'real_data_w_features__2023_11_02_10_42.csv'))
df_features[['features']] <- c("v")
df_features[['relative_crps']] <- df_features[['mean_crps']]/df_features[['mean_crps']] 

df_no_features <-  fread(sprintf(locn,'real_data_no_features__2023_11_02_09_58.csv'))
df_no_features[['features']] <- c("x")
df_no_features[['relative_crps']] <- df_no_features[['mean_crps']]/df_features[['mean_crps']] 


df_tot_rd <- rbind(df_features,
                   df_no_features)[order(k),]


colnames(df_tot_rd)[2] <- 'True_amount'

df_tot_rd[order(k,-features),c("k",
             "features",
             "True_amount",
             "EI_ajr",
             "EI_cl",
             "msep_ajr",
             "msep_cl",
             "relative_crps")] %>%
  reframe(k=k,
    features=features,
    EI_ajr=EI_ajr,
    EI_cl=EI_cl,
    actual.tot = True_amount/divided.by,
    pred.tot=(1+EI_ajr)*True_amount,
    cl.tot=(1+EI_cl)*True_amount,
    aj.msep=msep_ajr/pred.tot,
    cl.msep=msep_cl/cl.tot,
    relative_crps=relative_crps
    #mcrps=mean_crps/divided.by
    
  ) %>%
  select('k','features','actual.tot','EI_ajr','EI_cl','aj.msep', 'cl.msep', 'relative_crps') %>%
  xtable(digits=4)%>%
  print(include.rownames=FALSE)


# Real data application strategy II ----
df_features1 <- fread(sprintf(locn,'real_data_w_features_maxp_6_2023_11_02_11_27.csv'))
df_features1[['features']] <- c("v")
df_features1[['relative_crps']] <- df_features1[['mean_crps']]/df_features1[['mean_crps']]

df_no_features1 <-  fread(sprintf(locn,'real_data_no_features_maxp_6_2023_11_02_09_29.csv'))
df_no_features1[['features']] <- c("x")
df_no_features1[['relative_crps']] <- df_no_features1[['mean_crps']]/df_features1[['mean_crps']]

df_tot_rd <- rbind(df_features1,
                   df_no_features1)[order(k),]

colnames(df_tot_rd)[2] <- 'True_amount'
df_tot_rd[complete.cases(df_tot_rd),][order(k,-features),
          c("k",
             "features",
             "True_amount",
             "EI_ajr",
             "EI_cl",
             "msep_ajr",
             "msep_cl",
             "relative_crps")]  %>%
  reframe(
    k=k,
    EI_ajr=EI_ajr,
    EI_cl=EI_cl,
    features=features,
    actual.tot = True_amount/divided.by,
    pred.tot=(1+EI_ajr)*True_amount,
    cl.tot=(1+EI_cl)*True_amount,
    aj.msep=msep_ajr/pred.tot,
    cl.msep=msep_cl/cl.tot,
    mcrps=relative_crps
    
  ) %>%
  select('k','features','actual.tot','EI_ajr','EI_cl','aj.msep', 'cl.msep', 'mcrps') %>%
  xtable(digits=4)%>%
  print(include.rownames=FALSE) 


# Reserves ----

df_features <- fread(sprintf(locn,'real_data_w_features__2023_11_02_10_42.csv'))
df_features[['features']] <- c("v")


df_no_features <-  fread(sprintf(locn,'real_data_no_features__2023_11_02_09_58.csv'))
df_no_features[['features']] <- c("x")

rbind(df_features,
                   df_no_features)[order(k),]%>%
  group_by(k) %>%
  filter(mean_crps==min(mean_crps)) %>%
  reframe(reserve_aj=reserve_aj/divided.by,
         sd_aj=msep_ajr/divided.by,
         features=features)  %>%
  xtable(digits=4)%>%
  print(include.rownames=FALSE) 
  

# Comparison with other algorithms ----

locn = "conditional-aj-reserving//results_csv//"
dt1 <- fread(paste0(locn,'simulation_hirem_10.csv'))
dt1[['scenario']] <- 'baseline'
dt2 <- fread(paste0(locn,'simulation_hirem_extreme_10.csv'))
dt2[['scenario']] <- 'extreme'
dt3 <- fread(paste0(locn,'simulation_hirem_settlement_10.csv'))
dt3[['scenario']] <- 'settlement'
dt4 <- fread(paste0(locn,'simulation_hirem_claimmix_10.csv'))
dt4[['scenario']] <- 'claimmix'

dt <- rbind(dt1,dt2,dt3,dt4)

dt.ei<-dt[,list(scenario=scenario,
         actual =mean(actual.tot)/divided.by,
         EI_ajr=mean(pred.tot/actual.tot)-1,
         EI_ajr_nf=mean(pred.tot.nf/actual.tot)-1,
         EI_hirem=mean(hirem.tot/actual.tot)-1,
         EI_cl=mean(cl.tot/actual.tot)-1),
   by=.(scenario)] %>%
  pivot_longer(cols= c('EI_ajr','EI_ajr_nf','EI_hirem','EI_cl'),
               names_to = "model",
               values_to="EI") 


dt.ei$model <- recode(dt.ei$model,
                      EI_ajr = "AJ", 
                      EI_ajr_nf = "AJ NF", 
                      EI_hirem="hirem", 
                      EI_cl="CL")


dt.cv <- dt[,list(aj.msep=mean(sqrt(var.aj)/pred.tot),
aj.msep.nf=mean(sqrt(var.aj.nf)/pred.tot.nf),
hirem.msep=mean(sqrt(var.hirem)/hirem.tot),
cl.msep=mean((var.cl)/cl.tot)),
by=.(scenario)] %>%
  pivot_longer(cols= c('aj.msep','aj.msep.nf','hirem.msep','cl.msep'),
               names_to = "model",
               values_to="cv") 

dt.cv$model <- recode(dt.cv$model,
                        aj.msep = "AJ", 
                        aj.msep.nf = "AJ NF", 
                        hirem.msep="hirem", 
                        cl.msep="CL")

dt.crps <- dt[,
list(mcrps.aj=mean(crps.aj/crps.aj),
mcrps.aj.nf=mean(crps.nf/crps.aj),
mcrps.hirem=mean(crps.hirem/crps.aj),
mcrps.cl=NA),
by=.(scenario)] %>%
  pivot_longer(cols= c('mcrps.aj','mcrps.aj.nf','mcrps.hirem','mcrps.cl'),
               names_to = "model",
               values_to="crps") 

dt.crps$model <- recode(dt.crps$model,
                        mcrps.aj = "AJ", 
                        mcrps.aj.nf = "AJ NF", 
                        mcrps.hirem="hirem", 
                        mcrps.cl="CL")

dt <- merge(merge(dt.ei,dt.cv, by = c('scenario','model')),dt.crps, by = c('scenario','model')) 
dt$model <- factor(dt$model, levels = c("AJ",'AJ NF', 'hirem', 'CL'))


setorderv(dt, c('scenario' , 'model')) 

dt[['Type']] <- rep(c('v','x','v','x'),4)

dt %>%
  select('scenario','Type','model','actual','EI','cv','crps') %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE) 

## hirem on our simulator
locn <- "conditional-aj-reserving\\results_csv\\%s"
fname_trend <- c('simulation_aytrend4_2024_02_29_13_31.csv',
                  'simulation_aytrend5_2024_02_29_13_32.csv',
                  'simulation_aytrend6_2024_02_29_13_34.csv',
                  'simulation_aytrend7_2024_02_29_13_36.csv')

dt1 <- data.table()

for(nm in fname_trend){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt1 <- rbind(dt1,tmp)}


dt.ei<-dt1[,list(scenario=k,
                actual =mean(actual.tot)/divided.by,
                EI_ajr=mean(pred.tot/actual.tot)-1,
                EI_hirem=mean(hirem.tot/actual.tot)-1,
                EI_cl=mean(cl.tot/actual.tot)-1),
          by=.(k)] %>%
  pivot_longer(cols= c('EI_ajr','EI_hirem','EI_cl'),
               names_to = "model",
               values_to="EI") 


dt.ei$model <- recode(dt.ei$model,
                      EI_ajr = "AJ", 
                      EI_ajr_nf = "AJ NF", 
                      EI_hirem="hirem", 
                      EI_cl="CL")


dt.cv <- dt1[,list(aj.msep=mean((aj.msep)/pred.tot),
                  # aj.msep.nf=mean(sqrt(var.aj.nf)/pred.tot.nf),
                  hirem.msep=mean((hirem.msep)/hirem.tot),
                  cl.msep=mean((cl.msep)/cl.tot)),
            by=.(k)] %>%
  pivot_longer(cols= c('aj.msep','hirem.msep','cl.msep'),
               names_to = "model",
               values_to="cv") 

dt.cv$model <- recode(dt.cv$model,
                      aj.msep = "AJ", 
                      aj.msep.nf = "AJ NF", 
                      hirem.msep="hirem", 
                      cl.msep="CL")


dt.crps <- dt1[,
              list(mcrps.aj=mean(mcrps/mcrps),
                   # mcrps.aj.nf=mean(crps.nf),
                   mcrps.hirem=mean(crps.hirem/mcrps),
                   mcrps.cl=NA),
              by=.(k)] %>%
  pivot_longer(cols= c('mcrps.aj','mcrps.hirem','mcrps.cl'),
               names_to = "model",
               values_to="crps") 

dt.crps$model <- recode(dt.crps$model,
                        mcrps.aj = "AJ", 
                        mcrps.aj.nf = "AJ NF", 
                        mcrps.hirem="hirem", 
                        mcrps.cl="CL")

dt <- merge(merge(dt.ei,dt.cv, by = c('k','model')),dt.crps, by = c('k','model')) 
dt$model <- factor(dt$model, levels = c("AJ", 'hirem', 'CL'))

# dt <- dt %>% select(-scenario)

dt %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)



