rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

# Construct the model ----

results = matrix(ncol=7)


maximum.p<-4


df = read.and.pp.data(fname='C:\\Users\\gpitt\\Documents\\Phd\\visiting\\dati\\data_ku\\final_claim_data.csv')

df = df %>% mutate(incPaid=incPaid)

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


reporting_period_vertical <- rbind(closed.claims,ibnr.claims,rbns.claims)
reporting_period_vertical<-reporting_period_vertical[,reporting_period:=reporting_period-accident_period]
reporting_period_vertical<-reporting_period_vertical[reporting_period<=(maximum.p-1)] #sometimes reporting calendar period is different than developed. Hence we did not filter it before
reporting_period_vertical<-reporting_period_vertical[,.(reporting_counts=.N),by=.(accident_period,reporting_period)]

reporting_period_tr <- clean.lt(as.triangle(reporting_period_vertical,
                                            origin="accident_period",
                                            dev="reporting_period",
                                            value="reporting_counts"))



closed.claims[,
              events:=encode.cc(development_period,maximum.p),
              by=Claim_number]


closed.claims.fit<-closed.claims[,
                                 .(events=events,
                                   times=incPaid),
                                 by=Claim_number]


rbns.claims.fit <- rbns.claims[,
                               encode.rbns(development_period,
                                           y=incPaid,
                                           ay=last(accident_period),maximum.p=maximum.p),
                               by=Claim_number]


data2fit <- rbind(rbns.claims.fit,
                  closed.claims.fit)

data2fit[,times:=cumsum(times),by=Claim_number]

colnames(data2fit)[colnames(data2fit)=="events"] <- "states"

features.data <- df[,c("Claim_number","Claim_type_key")]
features.data<-setDT(features.data)[,.(Claim_type_key=first(Claim_type_key)),
                                    by=Claim_number]

data2fit<-features.data[data2fit, 
                        on = 'Claim_number']


fit.x <- data2fit[,.(Claim_type_key=unique(Claim_type_key)),by=Claim_number][order(Claim_number),]

data2fit=data2fit[,c("Claim_number", "states", "times","Claim_type_key")]

colnames(data2fit)[4]<-"X"

data2fit[,X:=as.numeric(X)]

data2fit_list <-  split(data2fit,data2fit$Claim_number)

data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1),x[,c('times','states','X')]))})

data2fit <- unname(data2fit)

rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]

out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_resultsX(k,
                                                                         X=as.numeric(fit.x$Claim_type_key),
                                                                         data2fit),
                 by=.(Claim_number)]


