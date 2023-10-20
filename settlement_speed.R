rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)
source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

df = read.and.pp.data(fname='C:\\Users\\gpitt\\Documents\\Phd\\visiting\\dati\\data_ku\\final_claim_data.csv')

df = df %>% mutate(incPaid=incPaid)

maximum.p <- 8

df =df %>%
  group_by(Claim_number) %>%
  filter(!any(development_period > (maximum.p-1)))

df=df  %>%
  group_by(Claim_number) %>%
  filter(1 %in%delta &(accident_period <= (maximum.p-1-1)))

closed.claims = df%>%
  group_by(Claim_number) %>%
  filter(!any(calendar_period > (maximum.p-1))) %>% as.data.table()

ibnr.claims = df%>%
  group_by(Claim_number) %>%
  filter(all(min(reporting_period)>(maximum.p-1)))%>% as.data.table()

rbns.claims = df%>%
  group_by(Claim_number) %>%
  filter(calendar_period<=(maximum.p-1)) %>%filter(!(1%in%delta))%>%as.data.table()

##
reporting_period_vertical <- rbind(closed.claims,ibnr.claims,rbns.claims)
reporting_period_vertical<-reporting_period_vertical[,reporting_period:=reporting_period-accident_period]
reporting_period_vertical<-reporting_period_vertical[reporting_period<=(maximum.p-1)] #sometimes reporting calendar period is different than developed. Hence we did not filter it before
reporting_period_vertical<-reporting_period_vertical[,.(reporting_counts=.N),by=.(accident_period,reporting_period)]

reporting_period_tr <- clean.lt(as.triangle(reporting_period_vertical,
                                            origin="accident_period",
                                            dev="reporting_period",
                                            value="reporting_counts"))

reporting_period_tr[1,6] <- 0

settled_period_vertical <- rbind(closed.claims,rbns.claims)

settled_period_vertical <- settled_period_vertical %>%
  group_by(Claim_number) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  as.data.table()

settled_period_vertical<-settled_period_vertical[,.(settled=sum(delta)),by=.(accident_period,development_period)]
settled_period_vertical=  rbind(settled_period_vertical,data.table(accident_period=rep(0:6,each=8),
                                           development_period=rep(1:7,8),
                                           settled=0))









