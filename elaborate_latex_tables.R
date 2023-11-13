rm(list=ls())
# libraries ----

library(data.table)
library(xtable)
library(ggplot2)
library(dplyr)
library(tidyr)

divided.by = 1e+06

#

locn <- "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"

# Simulations ----

# Table 1 ----


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


# Table 2 ----


fname_trend <- c('simulation_aytrend4_2023_10_23_12_52.csv',
                 'simulation_aytrend5_2023_10_24_17_35.csv',
                 'simulation_aytrend6_2023_10_23_12_54.csv',
                 'simulation_aytrend7_2023_10_23_12_56.csv')

dt1 <- data.table()

for(nm in fname_trend){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt1 <- rbind(dt1,tmp)}

dt1[['features']] = "v"
dt1[['scenario']] = "beta"


fname_trend_i <- c('simulation_aytrend_interceptmodel_all_records_7_2023_10_24_16_32.csv',
  'simulation_aytrend_interceptmodel_all_records_6_2023_10_24_16_31.csv',
  'simulation_aytrend_interceptmodel_all_records_5_2023_10_24_17_36.csv',
  'simulation_aytrend_interceptmodel_all_records_4_2023_10_24_16_27.csv')
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

dt.tot[,.(actual.tot = mean(actual.tot)/divided.by,
          ei_aj=mean(pred.tot/actual.tot-1),
          ei_cl=mean(cl.tot/actual.tot-1),
          aj.msep=mean(aj.msep/pred.tot),
          cl.msep=mean(cl.msep/cl.tot),
          mcrps=mean(mcrps)/divided.by
),
by=.(k,
     scenario,
     features)][order(k),] %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)


# Real data case studies ----

# Table 4 ----

df_features <- fread(sprintf(locn,'real_data_w_features__2023_11_02_10_42.csv'))
df_features[['features']] <- c("v")


df_no_features <-  fread(sprintf(locn,'real_data_no_features__2023_11_02_09_58.csv'))
df_no_features[['features']] <- c("x")

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
             "mean_crps")] %>%
  reframe(k=k,
    features=features,
    EI_ajr=EI_ajr,
    EI_cl=EI_cl,
    actual.tot = True_amount/divided.by,
    pred.tot=(1+EI_ajr)*True_amount,
    cl.tot=(1+EI_cl)*True_amount,
    aj.msep=msep_ajr/pred.tot,
    cl.msep=msep_cl/cl.tot,
    mcrps=mean_crps/divided.by
    
  ) %>%
  select('k','features','actual.tot','EI_ajr','EI_cl','aj.msep', 'cl.msep', 'mcrps') %>%
  xtable(digits=4)%>%
  print(include.rownames=FALSE)


# Real data application strategy II ----
df_features1 <- fread(sprintf(locn,'real_data_w_features_maxp_6_2023_11_02_11_27.csv'))
df_features1[['features']] <- c("v")

df_no_features1 <-  fread(sprintf(locn,'real_data_no_features_maxp_6_2023_11_02_09_29.csv'))
df_no_features1[['features']] <- c("x")


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
             "mean_crps")]  %>%
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
    mcrps=mean_crps/divided.by
    
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
  





