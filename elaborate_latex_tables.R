rm(list=ls())
# libraries ----

library(data.table)
library(xtable)
library(ggplot2)
library(dplyr)
library(tidyr)

divided.by = 1e+07

#

locn <- "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"

# Simulations ----

# Model selection ----

fname="crpss.csv"

crpss <- fread(sprintf(locn,fname))

crpss[1:9,] %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)

# No features ----

fname_no_trend <- c(
  'simulation_no_features_all_records_4_2023_10_18_17_15.csv',
  # 'simulation_no_features_all_records_7_2023_10_18_17_20.csv',

'simulation_no_features_all_records_5_2023_10_18_17_17.csv',
'simulation_no_features_all_records_6_2023_10_18_17_18.csv')



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

dt[,aj.msep:=sqrt(variance)]


# AY trend ----

fname_trend <- c('simulation_aytrend4_2023_10_23_12_52.csv',
                 'simulation_aytrend5_2023_10_23_12_53.csv',
                 'simulation_aytrend6_2023_10_23_12_54.csv',
                 'simulation_aytrend7_2023_10_23_12_56.csv')


dt1 <- data.table()

for(nm in fname_trend){
  
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  dt1 <- rbind(dt1,tmp)}

dt1[['features']] = "v"

dt.tot <- rbind(dt[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features")],
                dt1[,c("actual.tot","pred.tot","cl.tot","aj.msep","cl.msep","k","features")])


dt.tot[,.(actual.tot = mean(actual.tot)/divided.by,
          ei_aj=mean(pred.tot/actual.tot-1),
          ei_cl=mean(cl.tot/actual.tot-1),
          aj.msep=mean(aj.msep)/divided.by,
          cl.msep=mean(cl.msep)/divided.by
          ),
       by=.(k,
            features)][order(k),] %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)

# Real data case studies ----

# Maximum shapes ---

df_features <- fread(sprintf(locn,'real_data_w_features_2023_10_10_09_57.csv'))
df_features[['features']] <- c("v")

df_no_features <-  fread(sprintf(locn,'real_data_no_features_different_shapes_2023_10_10_19_08.csv'))
df_no_features[['features']] <- c("x")

df_tot_rd <- rbind(df_features,
                   df_no_features)[order(k),]


df_tot_rd[,c("k",
             "features",
             "True amount",
             "EI_ajr",
             "EI_cl",
             "msep_ajr",
             "msep_cl",
             "mean_crps")] %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)


# Real data application strategy II ----


df_features1 <- fread(sprintf(locn,'real_data_w_features_maxp_6_2023_10_10_16_42.csv'))
df_features1[['features']] <- c("v")

df_no_features1 <-  fread(sprintf(locn,'real_data_no_features_maxp_6_2023_10_12_13_09.csv'))
df_no_features1[['features']] <- c("x")


df_tot_rd <- rbind(df_features1,
                   df_no_features1)[order(k),]


df_tot_rd[,c("k",
             "features",
             "True amount",
             "EI_ajr",
             "EI_cl",
             "msep_ajr",
             "msep_cl",
             "mean_crps")] %>%
  xtable(digits=3)%>%
  print(include.rownames=FALSE)










