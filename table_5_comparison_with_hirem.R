rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('conditional-aj-reserving\\helper_functions_ajr.R')

# extra for hirem analysis ----

library(hirem)
library(tidyr)

tsh = 1e-08
results <- NULL
maximum.p<-10

# Baseline ----

for(sim in 1:10){
  
  input.data <-hirem::simulate_scenario_baseline(seed=sim,  
                                                 n = 85000,
                                                 max_rep_delay = 1,
                                                 max_set_delay = maximum.p-1) %>%
    as.data.frame()
  

  df <- input.data %>% 
    select(-c(
      "rep.delay",
      "occ.month",
      "rep.date",
      "settlement.date")) %>%
    rename(Claim_number=claim.nr,
           Claim_type_key=type,
           development_period=dev.year,
           accident_period=occ.year,
           reporting_period=rep.year,
           delta=settlement,
           incPaid=size) %>%
    mutate(accident_period=accident_period-min(accident_period),
           calendar_period=accident_period+development_period)

  df$Claim_type_key <- recode(df$Claim_type_key,T1 = "1", T2 = "2", T3="3")
  
  # safety check
  df=df  %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &(first(accident_period) <= (maximum.p-1-1))) %>% #
    as.data.table()
  
  
  closed.claims = df%>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    # filter(calendar_period<=(maximum.p-1)) %>%
    filter(1 %in%delta &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup ()%>% 
    as.data.table()
  dim(closed.claims)
  
  closed.claims<- closed.claims %>%
    group_by(Claim_number) %>%
    mutate(tmp.indicator = diff(c(5000,delta))) %>%
    as.data.table()
  
  closed.claims[payment==0 & tmp.indicator==1, 'incPaid'] <- tsh
  
  closed.claims[payment==0 & tmp.indicator==1, 'payment'] <- 1
  
  closed.claims <- closed.claims %>% filter(payment!=0)
  
  closed.claims <- closed.claims[,-c('tmp.indicator')]
  
  
  
  # safety check -- dim(ibnr.claims)[1] is zero
  ibnr.claims = df%>%
    group_by(Claim_number) %>%
    filter(reporting_period>(maximum.p-1))%>% as.data.table()
  dim(ibnr.claims)
  
  rbns.claims <- df %>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    filter(!(1 %in%delta) &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup() %>%
    filter(payment!=0)%>%
    as.data.table()
  
  # triangles ----
  
  actual_triangle <- incr2cum(as.triangle(df,
                                          origin='accident_period',
                                          dev='development_period',
                                          value='incPaid'))[,1:(maximum.p-1)]
  
  diag_ftr <- sum(ChainLadder::getLatestCumulative(clean.lt(actual_triangle)))
  
  
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
                                             ay=last(accident_period),
                                             maximum.p=maximum.p),
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
  
  data2fit=data2fit[,c("Claim_number", 
                       "states", 
                       "times",
                       "Claim_type_key")]
  
  
  colnames(data2fit)[4]<-"X"
  
  data2fit[,X:=as.numeric(as.character(X))]
  
  data2fit_list <-  split(data2fit,data2fit$Claim_number)
  
  data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1,first(x[,4])),x[,c('times','states','X')]))})
  
  data2fit <- unname(data2fit)
  
  data.list.ref <- lapply(data2fit, function(x){x[[3]]<-unique(x[[3]])
  return(x)
  })
  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                           on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)
  
  
  actual=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    reframe(true.ultimate=sum(incPaid)) %>%
    as.data.table()
  
  rbns.data <-merge(actual,rbns.data,by='Claim_number')
  
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=Claim_type_key ,
                                                                                               models.list), 
                   by=.(Claim_number)]
  
  
  actual.uc=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    summarise(uc=sum(incPaid)) %>% unlist()
  
  rbns.uc <- sum(out$ultimate)
  var.aj <- sum(out$variance)
  closed.uc<-sum(closed.claims$incPaid)
  
  cl.data <- find.t.data(df)
  cl.tot <- cl.calculator(incr2cum(cl.data[,1:(maximum.p-1)]))
  
  rbns.reserve.cl <- sum(cl.tot$ultimate)-diag_ftr
  
  crps.aj <- mean(out$crps_i)
  
  ## hirem 
  
  
  reserving_data <- input.data
  
  upper_triangle <- reserving_data %>% filter(calendar.year <= maximum.p)
  lower_triangle <- reserving_data %>% filter(calendar.year > maximum.p)
  
  model <- hirem(upper_triangle)
  
  model <- hirem(upper_triangle) %>%
    layer_glm('settlement', binomial(link = cloglog)) %>%
    layer_glm('payment', binomial(link = logit)) %>%
    layer_glm('size', Gamma(link = log),
              filter = function(data){data$payment == 1})
  
  model <- fit(model,
               settlement = 'settlement ~  dev.year.fact+factor(type)',
               payment = 'payment ~ dev.year.fact+settlement + factor(type)',
               size = 'size ~ dev.year.fact+settlement + factor(type)')
  
  update <- function(data) {
    data %>%
      dplyr::mutate(dev.year = dev.year + 1,
                    calendar.year = calendar.year + 1)
  }
  
  model <- register_updater(model, update)
  
  
  simul <- simulate(model,
                    nsim = 5,
                    filter = function(data){dplyr::filter(data,
                                                          dev.year <= (maximum.p-1),
                                                          settlement == 0)},
                    data = reserving_data %>% dplyr::filter(calendar.year == (maximum.p-1)))
  
  
  setDT(simul)
  
  rbns_estimate <-simul[,.(rbns = sum(size)) ,by=.(simulation)]
  
  
  rbns_estimate_av <- simul[,.(rbns=sum(size),
                               type=first(type)),
                            by=.(claim.nr,simulation)][,.(rbns=mean(rbns),
                                                          var.hirem=var(rbns),
                                                          type=first(type)),by=.(claim.nr)]
  var.hirem <- sum(rbns_estimate_av$var.hirem)
  
  rbns_estimate_av<- rbns_estimate_av[,-c('var.hirem')]
  
  setDT(reserving_data)
  
  actual.hirem <-  reserving_data[,.(true.ultimate=sum(size)),by='claim.nr']
  
  # find rbns
  commonID<-intersect(lower_triangle$claim.nr, reserving_data$claim.nr)
  
  actual.hirem <- actual.hirem[claim.nr %in% commonID,]
  
  rbns.data.hirem <-merge(actual.hirem,rbns_estimate_av,by='claim.nr')
  
  models.list.hirem <- find.models.list.hirem(hirem.df = rbns_estimate_av)
  
  rbns.data.hirem<-rbns.data.hirem[,c('crps'):=individual_results_hirem(true.ultimate=true.ultimate,
                                                                        X=type,
                                                                        models.list=models.list.hirem),by=.(claim.nr)]
  
  rbns.reserve.hirem <- mean(rbns_estimate$rbns)
  
  crps.hirem <- mean(rbns.data.hirem$crps)
  
  ## unconditional aj
  
  fit1 <- AalenJohansen::aalen_johansen(data.list.ref)
  
  yhat <-  sapply(fit1$p, last)
  
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  x.vals <- fit1$t
  
  yhat[cond]<-1
  # x.vals <- fit$t
  
  cond2 = diff(c(0,yhat))<0
  
  while(sum(cond2)!=0){
    
    yhat <-yhat[!cond2]
    x.vals <-x.vals[!cond2]
    cond2 = diff(c(0,yhat))<0
    
  }
  
  out2 <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_nf(k,
                                                                              true.ultimate=true.ultimate,
                                                                              x.vals=x.vals,
                                                                              yhat=yhat), 
                    by=.(Claim_number)]
  
  
  
  
  rbns.uc.nf <- sum(out2$ultimate)
  var.nf <- sum(out2$variance)
  crps.nf <- mean(out2$crps)
  
  results <- rbind(results,
                   data.frame(pred.tot=(rbns.uc+closed.uc),
                              pred.tot.nf=(rbns.uc.nf+closed.uc),
                              actual.tot=actual.uc,
                              cl.tot=sum(cl.tot$ultimate),
                              hirem.tot=(rbns.reserve.hirem+diag_ftr),
                              var.aj=var.aj,
                              var.aj.nf=var.nf,
                              var.hirem=var.hirem,
                              var.cl=cl.tot$process.error.i,
                              crps.aj=crps.aj,
                              crps.nf=crps.nf,
                              crps.hirem=crps.hirem))
  
 
  
  
  
  
}

finame=paste0("simulation_",
              'hirem_',
              maximum.p,
              ".csv")


fwrite(results, file=sprintf(locn,finame))



# Extreme ----

results <- NULL
maximum.p<-10

for(sim in 1:10){
  
  input.data <-hirem::simulate_scenario_extreme_event(seed=sim,  
                                                      n = 85000,
                                                      max_rep_delay = 1,
                                                      max_set_delay = maximum.p-1) %>%
    as.data.frame()
  

  df <- input.data %>% 
    select(-c(
      "rep.delay",
      "occ.month",
      "rep.date",
      "settlement.date")) %>%
    rename(Claim_number=claim.nr,
           Claim_type_key=type,
           development_period=dev.year,
           accident_period=occ.year,
           reporting_period=rep.year,
           delta=settlement,
           incPaid=size) %>%
    mutate(accident_period=accident_period-min(accident_period),
           calendar_period=accident_period+development_period) 

  df$Claim_type_key <- recode(df$Claim_type_key,T1 = "1", T2 = "2", T3="3")
  
  df %>%
    filter(accident_period==0&development_period == maximum.p-1) %>%
    reframe(sum(incPaid))
  
  # safety check
  df=df  %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &(first(accident_period) <= (maximum.p-1-1))) %>% #
    as.data.table()
  
  
  closed.claims = df%>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    # filter(calendar_period<=(maximum.p-1)) %>%
    filter(1 %in%delta &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup ()%>% 
    as.data.table()
  dim(closed.claims)
  
  closed.claims<- closed.claims %>%
    group_by(Claim_number) %>%
    mutate(tmp.indicator = diff(c(5000,delta))) %>%
    as.data.table()
  
  closed.claims[payment==0 & tmp.indicator==1, 'incPaid'] <- tsh
  
  closed.claims[payment==0 & tmp.indicator==1, 'payment'] <- 1
  
  closed.claims <- closed.claims %>% filter(payment!=0)
  
  
  closed.claims <- closed.claims[,-c('tmp.indicator')]
  
  # safety check -- dim(ibnr.claims)[1] is zero
  ibnr.claims = df%>%
    group_by(Claim_number) %>%
    filter(reporting_period>(maximum.p-1))%>% as.data.table()
  dim(ibnr.claims)
  
  rbns.claims <- df %>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    # filter(calendar_period<=(maximum.p-1)) %>%
    filter(!(1 %in%delta) &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup() %>%
    filter(payment!=0)%>%
    as.data.table()
  
  # triangles ----
  
  actual_triangle <- incr2cum(as.triangle(df,
                                          origin='accident_period',
                                          dev='development_period',
                                          value='incPaid'))[,1:(maximum.p-1)]
  
  diag_ftr <- sum(ChainLadder::getLatestCumulative(clean.lt(actual_triangle)))
  
  
  closed.claims[,
                events:=encode.cc(development_period,maximum.p),
                by=Claim_number]
  
  
  closed.claims.fit<-closed.claims[,
                                   .(events=events,
                                     times=incPaid),
                                   by=Claim_number]
  
  #rbind(rbns.claims1,rbns.claims2)
  rbns.claims.fit <- rbns.claims[,
                                 encode.rbns(development_period,
                                             y=incPaid,
                                             ay=last(accident_period),
                                             maximum.p=maximum.p),
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
  
  data2fit=data2fit[,c("Claim_number", 
                       "states", 
                       "times",
                       "Claim_type_key")]
  
  
  colnames(data2fit)[4]<-"X"
  
  data2fit[,X:=as.numeric(as.character(X))]
  
  data2fit_list <-  split(data2fit,data2fit$Claim_number)
  
  data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1,first(x[,4])),x[,c('times','states','X')]))})
  
  data2fit <- unname(data2fit)
  
  data.list.ref <- lapply(data2fit, function(x){x[[3]]<-unique(x[[3]])
  return(x)
  })
  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                           on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)
  
  
  actual=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    reframe(true.ultimate=sum(incPaid)) %>%
    as.data.table()
  
  rbns.data <-merge(actual,rbns.data,by='Claim_number')
  
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=Claim_type_key ,
                                                                                               models.list), 
                   by=.(Claim_number)]
  
  
  actual.uc=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    summarise(uc=sum(incPaid)) %>% unlist()
  
  rbns.uc <- sum(out$ultimate)
  var.aj <- sum(out$variance)
  closed.uc<-sum(closed.claims$incPaid)
  
  cl.data <- find.t.data(df)
  cl.tot <- cl.calculator(incr2cum(cl.data[,1:(maximum.p-1)]))
  
  rbns.reserve.cl <- sum(cl.tot$ultimate)-diag_ftr
  
  crps.aj <- mean(out$crps_i)
  
  ## hirem 
  
  
  reserving_data <- input.data
  
  upper_triangle <- reserving_data %>% filter(calendar.year <= maximum.p)
  lower_triangle <- reserving_data %>% filter(calendar.year > maximum.p)
  
  model <- hirem(upper_triangle)
  
  model <- hirem(upper_triangle) %>%
    layer_glm('settlement', binomial(link = cloglog)) %>%
    layer_glm('payment', binomial(link = logit)) %>%
    layer_glm('size', Gamma(link = log),
              filter = function(data){data$payment == 1})
  
  model <- fit(model,
               settlement = 'settlement ~  dev.year.fact+factor(type)',
               payment = 'payment ~ dev.year.fact+settlement + factor(type)',
               size = 'size ~ dev.year.fact+settlement + factor(type)')
  
  update <- function(data) {
    data %>%
      dplyr::mutate(dev.year = dev.year + 1,
                    calendar.year = calendar.year + 1)
  }
  
  model <- register_updater(model, update)
  
  
  simul <- simulate(model,
                    nsim = 5,
                    filter = function(data){dplyr::filter(data,
                                                          dev.year <= (maximum.p-1),
                                                          settlement == 0)},
                    data = reserving_data %>% dplyr::filter(calendar.year == (maximum.p-1)))
  
  
  setDT(simul)
  
  rbns_estimate <-simul[,.(rbns = sum(size)) ,by=.(simulation)]
  
  
  rbns_estimate_av <- simul[,.(rbns=sum(size),
                               type=first(type)),
                            by=.(claim.nr,simulation)][,.(rbns=mean(rbns),
                                                          var.hirem=var(rbns),
                                                          type=first(type)),by=.(claim.nr)]
  var.hirem <- sum(rbns_estimate_av$var.hirem)
  
  rbns_estimate_av<- rbns_estimate_av[,-c('var.hirem')]
  
  setDT(reserving_data)
  
  actual.hirem <-  reserving_data[,.(true.ultimate=sum(size)),by='claim.nr']
  
  # find rbns
  commonID<-intersect(lower_triangle$claim.nr, reserving_data$claim.nr)
  
  actual.hirem <- actual.hirem[claim.nr %in% commonID,]
  
  rbns.data.hirem <-merge(actual.hirem,rbns_estimate_av,by='claim.nr')
  
  models.list.hirem <- find.models.list.hirem(hirem.df = rbns_estimate_av)
  
  rbns.data.hirem<-rbns.data.hirem[,c('crps'):=individual_results_hirem(true.ultimate=true.ultimate,
                                                                        X=type,
                                                                        models.list=models.list.hirem),by=.(claim.nr)]
  
  rbns.reserve.hirem <- mean(rbns_estimate$rbns)
  
  crps.hirem <- mean(rbns.data.hirem$crps)
  
  ## unconditional aj
  
  fit1 <- AalenJohansen::aalen_johansen(data.list.ref)
  
  yhat <-  sapply(fit1$p, last)
  
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  x.vals <- fit1$t
  
  yhat[cond]<-1
  
  cond2 = diff(c(0,yhat))<0
  
  while(sum(cond2)!=0){
    
    yhat <-yhat[!cond2]
    x.vals <-x.vals[!cond2]
    cond2 = diff(c(0,yhat))<0
    
  }
  
  out2 <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_nf(k,
                                                                              true.ultimate=true.ultimate,
                                                                              x.vals=x.vals,
                                                                              yhat=yhat), 
                    by=.(Claim_number)]
  
  
  
  
  rbns.uc.nf <- sum(out2$ultimate)
  var.nf <- sum(out2$variance)
  crps.nf <- mean(out2$crps)
  
  results <- rbind(results,
                   data.frame(pred.tot=(rbns.uc+closed.uc),
                              pred.tot.nf=(rbns.uc.nf+closed.uc),
                              actual.tot=actual.uc,
                              cl.tot=sum(cl.tot$ultimate),
                              hirem.tot=(rbns.reserve.hirem+diag_ftr),
                              var.aj=var.aj,
                              var.aj.nf=var.nf,
                              var.hirem=var.hirem,
                              var.cl=cl.tot$process.error.i,
                              crps.aj=crps.aj,
                              crps.nf=crps.nf,
                              crps.hirem=crps.hirem))
  

  
  
  
}

locn = "conditional-aj-reserving\\results_csv\\%s"

finame=paste0("simulation_",
              'hirem_',
              'extreme_',
              maximum.p,
              ".csv")


fwrite(results, file=sprintf(locn,finame))

# Settlement ----

tsh = 1e-08
results <- NULL
maximum.p<-10

for(sim in 1:10){
  
  
  input.data <-hirem::simulate_scenario_change_in_settlement(seed=sim,  
                                                             n = 85000,
                                                             max_rep_delay = 1,
                                                             max_set_delay = maximum.p-1) %>%
    as.data.frame()
  

  df <- input.data %>% 
    select(-c(
      "rep.delay",
      "occ.month",
      "rep.date",
      "settlement.date")) %>%
    rename(Claim_number=claim.nr,
           Claim_type_key=type,
           development_period=dev.year,
           accident_period=occ.year,
           reporting_period=rep.year,
           delta=settlement,
           incPaid=size) %>%
    mutate(accident_period=accident_period-min(accident_period),
           calendar_period=accident_period+development_period) #%>%

  df$Claim_type_key <- recode(df$Claim_type_key,T1 = "1", T2 = "2", T3="3")
  
  df %>%
    filter(accident_period==0&development_period == maximum.p-1) %>%
    reframe(sum(incPaid))
  
  # safety check
  df=df  %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &(first(accident_period) <= (maximum.p-1-1))) %>% #
    as.data.table()
  
  
  closed.claims = df%>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    # filter(calendar_period<=(maximum.p-1)) %>%
    filter(1 %in%delta &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup ()%>% 
    as.data.table()
  dim(closed.claims)
  
  closed.claims<- closed.claims %>%
    group_by(Claim_number) %>%
    mutate(tmp.indicator = diff(c(5000,delta))) %>%
    as.data.table()
  
  closed.claims[payment==0 & tmp.indicator==1, 'incPaid'] <- tsh
  
  closed.claims[payment==0 & tmp.indicator==1, 'payment'] <- 1
  
  closed.claims <- closed.claims %>% filter(payment!=0)
  
  
  closed.claims <- closed.claims[,-c('tmp.indicator')]
  
  # safety check -- dim(ibnr.claims)[1] is zero
  ibnr.claims = df%>%
    group_by(Claim_number) %>%
    filter(reporting_period>(maximum.p-1))%>% as.data.table()
  dim(ibnr.claims)
  
  rbns.claims <- df %>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    filter(!(1 %in%delta) &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup() %>%
    filter(payment!=0)%>%
    as.data.table()

  # triangles ----
  
  actual_triangle <- incr2cum(as.triangle(df,
                                          origin='accident_period',
                                          dev='development_period',
                                          value='incPaid'))[,1:(maximum.p-1)]
  
  diag_ftr <- sum(ChainLadder::getLatestCumulative(clean.lt(actual_triangle)))
  
  
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
                                             ay=last(accident_period),
                                             maximum.p=maximum.p),
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
  
  data2fit=data2fit[,c("Claim_number", 
                       "states", 
                       "times",
                       "Claim_type_key")]
  
  
  colnames(data2fit)[4]<-"X"
  
  data2fit[,X:=as.numeric(as.character(X))]
  
  data2fit_list <-  split(data2fit,data2fit$Claim_number)
  
  data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1,first(x[,4])),x[,c('times','states','X')]))})
  
  data2fit <- unname(data2fit)
  
  data.list.ref <- lapply(data2fit, function(x){x[[3]]<-unique(x[[3]])
  return(x)
  })
  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                           on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)
  
  
  actual=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    reframe(true.ultimate=sum(incPaid)) %>%
    as.data.table()
  
  rbns.data <-merge(actual,rbns.data,by='Claim_number')
  
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=Claim_type_key ,
                                                                                               models.list), 
                   by=.(Claim_number)]
  
  
  actual.uc=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    summarise(uc=sum(incPaid)) %>% unlist()
  
  rbns.uc <- sum(out$ultimate)
  var.aj <- sum(out$variance)
  closed.uc<-sum(closed.claims$incPaid)
  
  cl.data <- find.t.data(df)
  cl.tot <- cl.calculator(incr2cum(cl.data[,1:(maximum.p-1)]))
  
  rbns.reserve.cl <- sum(cl.tot$ultimate)-diag_ftr
  
  crps.aj <- mean(out$crps_i)
  
  ## hirem 
  
  
  reserving_data <- input.data
  
  upper_triangle <- reserving_data %>% filter(calendar.year <= maximum.p)
  lower_triangle <- reserving_data %>% filter(calendar.year > maximum.p)
  
  model <- hirem(upper_triangle)
  
  model <- hirem(upper_triangle) %>%
    layer_glm('settlement', binomial(link = cloglog)) %>%
    layer_glm('payment', binomial(link = logit)) %>%
    layer_glm('size', Gamma(link = log),
              filter = function(data){data$payment == 1})
  
  model <- fit(model,
               settlement = 'settlement ~  dev.year.fact+factor(type)',
               payment = 'payment ~ dev.year.fact+settlement + factor(type)',
               size = 'size ~ dev.year.fact+settlement + factor(type)')
  
  update <- function(data) {
    data %>%
      dplyr::mutate(dev.year = dev.year + 1,
                    calendar.year = calendar.year + 1)
  }
  
  model <- register_updater(model, update)
  
  
  simul <- simulate(model,
                    nsim = 5,
                    filter = function(data){dplyr::filter(data,
                                                          dev.year <= (maximum.p-1),
                                                          settlement == 0)},
                    data = reserving_data %>% dplyr::filter(calendar.year == (maximum.p-1)))
  
  
  setDT(simul)
  
  rbns_estimate <-simul[,.(rbns = sum(size)) ,by=.(simulation)]
  
  rbns_estimate_av <- simul[,.(rbns=sum(size),
                               type=first(type)),
                            by=.(claim.nr,simulation)][,.(rbns=mean(rbns),
                                                          var.hirem=var(rbns),
                                                          type=first(type)),by=.(claim.nr)]
  var.hirem <- sum(rbns_estimate_av$var.hirem)
  
  rbns_estimate_av<- rbns_estimate_av[,-c('var.hirem')]
  
  setDT(reserving_data)
  
  actual.hirem <-  reserving_data[,.(true.ultimate=sum(size)),by='claim.nr']
  
  # find rbns
  commonID<-intersect(lower_triangle$claim.nr, reserving_data$claim.nr)
  
  actual.hirem <- actual.hirem[claim.nr %in% commonID,]
  
  rbns.data.hirem <-merge(actual.hirem,rbns_estimate_av,by='claim.nr')
  
  models.list.hirem <- find.models.list.hirem(hirem.df = rbns_estimate_av)
  
  rbns.data.hirem<-rbns.data.hirem[,c('crps'):=individual_results_hirem(true.ultimate=true.ultimate,
                                                                        X=type,
                                                                        models.list=models.list.hirem),by=.(claim.nr)]
  
  rbns.reserve.hirem <- mean(rbns_estimate$rbns)
  
  crps.hirem <- mean(rbns.data.hirem$crps)
  
  ## unconditional aj
  
  fit1 <- AalenJohansen::aalen_johansen(data.list.ref)
  
  yhat <-  sapply(fit1$p, last)
  
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  x.vals <- fit1$t
  
  yhat[cond]<-1
  # x.vals <- fit$t
  
  cond2 = diff(c(0,yhat))<0
  
  while(sum(cond2)!=0){
    
    yhat <-yhat[!cond2]
    x.vals <-x.vals[!cond2]
    cond2 = diff(c(0,yhat))<0
    
  }
  
  out2 <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_nf(k,
                                                                              true.ultimate=true.ultimate,
                                                                              x.vals=x.vals,
                                                                              yhat=yhat), 
                    by=.(Claim_number)]
  
  
  
  
  rbns.uc.nf <- sum(out2$ultimate)
  var.nf <- sum(out2$variance)
  crps.nf <- mean(out2$crps)
  
  results <- rbind(results,
                   data.frame(pred.tot=(rbns.uc+closed.uc),
                              pred.tot.nf=(rbns.uc.nf+closed.uc),
                              actual.tot=actual.uc,
                              cl.tot=sum(cl.tot$ultimate),
                              hirem.tot=(rbns.reserve.hirem+diag_ftr),
                              var.aj=var.aj,
                              var.aj.nf=var.nf,
                              var.hirem=var.hirem,
                              var.cl=cl.tot$process.error.i,
                              crps.aj=crps.aj,
                              crps.nf=crps.nf,
                              crps.hirem=crps.hirem))
  
}

locn = "conditional-aj-reserving\\results_csv\\%s"

finame=paste0("simulation_",
              'hirem_',
              'settlement_',
              maximum.p,
              ".csv")


fwrite(results, file=sprintf(locn,finame))


# Claim mix ----

tsh = 1e-08
results <- NULL
maximum.p<-10

for(sim in 1:10){
  
  input.data <-hirem::simulate_scenario_claim_mix(seed=sim,  
                                                  n = 85000,
                                                  max_rep_delay = 1,
                                                  max_set_delay = maximum.p-1) %>%
    as.data.frame()
  

  df <- input.data %>% 
    select(-c(
      "rep.delay",
      "occ.month",
      "rep.date",
      "settlement.date")) %>%
    rename(Claim_number=claim.nr,
           Claim_type_key=type,
           development_period=dev.year,
           accident_period=occ.year,
           reporting_period=rep.year,
           delta=settlement,
           incPaid=size) %>%
    mutate(accident_period=accident_period-min(accident_period),
           calendar_period=accident_period+development_period) 
  df$Claim_type_key <- recode(df$Claim_type_key,T1 = "1", T2 = "2", T3="3")
  
  df %>%
    filter(accident_period==0&development_period == maximum.p-1) %>%
    reframe(sum(incPaid))
  
  # safety check
  df=df  %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &(first(accident_period) <= (maximum.p-1-1))) %>% #
    as.data.table()
  
  
  closed.claims = df%>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup ()%>% 
    as.data.table()
  dim(closed.claims)
  
  closed.claims<- closed.claims %>%
    group_by(Claim_number) %>%
    mutate(tmp.indicator = diff(c(5000,delta))) %>%
    as.data.table()
  
  closed.claims[payment==0 & tmp.indicator==1, 'incPaid'] <- tsh
  
  closed.claims[payment==0 & tmp.indicator==1, 'payment'] <- 1
  
  closed.claims <- closed.claims %>% filter(payment!=0)
  
  
  closed.claims <- closed.claims[,-c('tmp.indicator')]
 
  
  # safety check -- dim(ibnr.claims)[1] is zero
  ibnr.claims = df%>%
    group_by(Claim_number) %>%
    filter(reporting_period>(maximum.p-1))%>% as.data.table()
  dim(ibnr.claims)
  
  
  rbns.claims <- df %>%
    filter((calendar_period)<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    filter(!(1 %in%delta) &last(reporting_period)<=(maximum.p-1)) %>%
    ungroup() %>%
    filter(payment!=0)%>%
    as.data.table()
  
  # triangles ----
  
  actual_triangle <- incr2cum(as.triangle(df,
                                          origin='accident_period',
                                          dev='development_period',
                                          value='incPaid'))[,1:(maximum.p-1)]
  
  diag_ftr <- sum(ChainLadder::getLatestCumulative(clean.lt(actual_triangle)))
  
  
  closed.claims[,
                events:=encode.cc(development_period,maximum.p),
                by=Claim_number]
  
  
  closed.claims.fit<-closed.claims[,
                                   .(events=events,
                                     times=incPaid),
                                   by=Claim_number]
  
  #rbind(rbns.claims1,rbns.claims2)
  rbns.claims.fit <- rbns.claims[,
                                 encode.rbns(development_period,
                                             y=incPaid,
                                             ay=last(accident_period),
                                             maximum.p=maximum.p),
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
  
  data2fit=data2fit[,c("Claim_number", 
                       "states", 
                       "times",
                       "Claim_type_key")]
  
  
  colnames(data2fit)[4]<-"X"
  
  data2fit[,X:=as.numeric(as.character(X))]
  
  data2fit_list <-  split(data2fit,data2fit$Claim_number)
  
  data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1,first(x[,4])),x[,c('times','states','X')]))})
  
  data2fit <- unname(data2fit)
  
  data.list.ref <- lapply(data2fit, function(x){x[[3]]<-unique(x[[3]])
  return(x)
  })
  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                           on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)
  
  
  actual=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    group_by(Claim_number) %>%
    reframe(true.ultimate=sum(incPaid)) %>%
    as.data.table()
  
  rbns.data <-merge(actual,rbns.data,by='Claim_number')
  
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=Claim_type_key ,
                                                                                               models.list), 
                   by=.(Claim_number)]
  
  
  actual.uc=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    summarise(uc=sum(incPaid)) %>% unlist()
  
  rbns.uc <- sum(out$ultimate)
  var.aj <- sum(out$variance)
  closed.uc<-sum(closed.claims$incPaid)
  
  cl.data <- find.t.data(df)
  cl.tot <- cl.calculator(incr2cum(cl.data[,1:(maximum.p-1)]))
  
  rbns.reserve.cl <- sum(cl.tot$ultimate)-diag_ftr
  
  crps.aj <- mean(out$crps_i)
  
  ## hirem 
  
  
  reserving_data <- input.data
  
  upper_triangle <- reserving_data %>% filter(calendar.year <= maximum.p)
  lower_triangle <- reserving_data %>% filter(calendar.year > maximum.p)
  
  model <- hirem(upper_triangle)
  
  model <- hirem(upper_triangle) %>%
    layer_glm('settlement', binomial(link = cloglog)) %>%
    layer_glm('payment', binomial(link = logit)) %>%
    layer_glm('size', Gamma(link = log),
              filter = function(data){data$payment == 1})
  
  model <- fit(model,
               settlement = 'settlement ~  dev.year.fact+factor(type)',
               payment = 'payment ~ dev.year.fact+settlement + factor(type)',
               size = 'size ~ dev.year.fact+settlement + factor(type)')
  
  update <- function(data) {
    data %>%
      dplyr::mutate(dev.year = dev.year + 1,
                    calendar.year = calendar.year + 1)
  }
  
  model <- register_updater(model, update)
  
  
  simul <- simulate(model,
                    nsim = 5,
                    filter = function(data){dplyr::filter(data,
                                                          dev.year <= (maximum.p-1),
                                                          settlement == 0)},
                    data = reserving_data %>% dplyr::filter(calendar.year == (maximum.p-1)))
  
  
  setDT(simul)
  
  rbns_estimate <-simul[,.(rbns = sum(size)) ,by=.(simulation)]

  
  rbns_estimate_av <- simul[,.(rbns=sum(size),
                               type=first(type)),
                            by=.(claim.nr,simulation)][,.(rbns=mean(rbns),
                                                          var.hirem=var(rbns),
                                                          type=first(type)),by=.(claim.nr)]
  var.hirem <- sum(rbns_estimate_av$var.hirem)
  
  rbns_estimate_av<- rbns_estimate_av[,-c('var.hirem')]
  
  setDT(reserving_data)
  
  actual.hirem <-  reserving_data[,.(true.ultimate=sum(size)),by='claim.nr']
  
  # find rbns
  commonID<-intersect(lower_triangle$claim.nr, reserving_data$claim.nr)
  
  actual.hirem <- actual.hirem[claim.nr %in% commonID,]
  
  rbns.data.hirem <-merge(actual.hirem,rbns_estimate_av,by='claim.nr')
  
  models.list.hirem <- find.models.list.hirem(hirem.df = rbns_estimate_av)
  
  rbns.data.hirem<-rbns.data.hirem[,c('crps'):=individual_results_hirem(true.ultimate=true.ultimate,
                                                                        X=type,
                                                                        models.list=models.list.hirem),by=.(claim.nr)]
  
  rbns.reserve.hirem <- mean(rbns_estimate$rbns)
  
  crps.hirem <- mean(rbns.data.hirem$crps)
  
  ## unconditional aj
  
  fit1 <- AalenJohansen::aalen_johansen(data.list.ref)
  
  yhat <-  sapply(fit1$p, last)
  
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  x.vals <- fit1$t
  
  yhat[cond]<-1
  # x.vals <- fit$t
  
  cond2 = diff(c(0,yhat))<0
  
  while(sum(cond2)!=0){
    
    yhat <-yhat[!cond2]
    x.vals <-x.vals[!cond2]
    cond2 = diff(c(0,yhat))<0
    
  }
  
  out2 <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_nf(k,
                                                                              true.ultimate=true.ultimate,
                                                                              x.vals=x.vals,
                                                                              yhat=yhat), 
                    by=.(Claim_number)]
  
  
  
  
  rbns.uc.nf <- sum(out2$ultimate)
  var.nf <- sum(out2$variance)
  crps.nf <- mean(out2$crps)
  
  results <- rbind(results,
                   data.frame(pred.tot=(rbns.uc+closed.uc),
                              pred.tot.nf=(rbns.uc.nf+closed.uc),
                              actual.tot=actual.uc,
                              cl.tot=sum(cl.tot$ultimate),
                              hirem.tot=(rbns.reserve.hirem+diag_ftr),
                              var.aj=var.aj,
                              var.aj.nf=var.nf,
                              var.hirem=var.hirem,
                              var.cl=cl.tot$process.error.i,
                              crps.aj=crps.aj,
                              crps.nf=crps.nf,
                              crps.hirem=crps.hirem))
  

}

locn = "conditional-aj-reserving\\results_csv\\%s"

finame=paste0("simulation_",
              'hirem_',
              'claimmix_',
              maximum.p,
              ".csv")


fwrite(results, file=sprintf(locn,finame))





