rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

# Construct the model ----

df = read.and.pp.data(fname='C:\\Users\\gpitt\\Documents\\Phd\\visiting\\dati\\data_ku\\final_claim_data.csv')



for(maximum.p in c(4,5,6,7,8)){

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
  
results = matrix(ncol=7)  

  
  
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
  
data2fit=data2fit[,c("Claim_number", 
                     "states", 
                     "times",
                     "Claim_type_key")]
  
colnames(data2fit)[4]<-"X"
  
data2fit[,X:=as.numeric(as.character(X))]
  
data2fit_list <-  split(data2fit,data2fit$Claim_number)

data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1,x[,4]),x[,c('times','states','X')]))})

data2fit <- unname(data2fit)

data.list.ref <- lapply(data2fit, function(x){x[[3]]<-unique(x[[3]])
                                              return(x)
                                                })
    

# for(k in 4:(maximum.p)){
#   
#   cat(paste0('State space with k '),as.character(k),'\n')
#   
#   if(k!=8){
#   data.list.ref <- lapply(data2fit,reformulate_state_space_X,newk = k)
#   
#   }else{
#     
#     data.list.ref <-data2fit}
  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                          on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)
  
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               X=Claim_type_key,
                                                                                               models.list), by=.(Claim_number)]
  
  
  actual.uc=df %>%
    ungroup() %>%
    filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
    summarise(uc=sum(incPaid)) %>% unlist()
  
  # chain ladder benchmark ----
  cl.data <- find.t.data(df)
  cl.tot <- cl.calculator(incr2cum(cl.data))
  
  # ibnr ----
  
  fj <- NULL
  
  reporting_period_tr[is.na(reporting_period_tr)] <- 0
  reporting_period_tr_cum <- incr2cum(reporting_period_tr)
  
  for(j in 2:(maximum.p-1)){
    
    fj <- c(fj,sum(reporting_period_tr_cum[1:(maximum.p-j),j])/sum(reporting_period_tr_cum[1:(maximum.p-j),j-1]))
  }
  
  
  ultimate <- c()
  diag <- c()
  for(i in 2:(maximum.p-1)){
    
    tmp <- reporting_period_tr_cum[i,maximum.p-i]
    diag<-c(diag,tmp)
    ultimate <- c(ultimate, tmp * last(cumprod(fj[(maximum.p-1-1-i+2):(maximum.p-1-1)])))
    
  }
  
  ibnr_num <- round(sum(ultimate-diag))
  
  # Compute the expected value on the total 
  
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
  
  
  ev.est <- cev(0,x=x.vals,y=yhat)
  
  # estimated uc
  
  ibnr.uc <- ibnr_num*ev.est
  
  # output ----
  
  var.aj <- sum(out$variance)
  crps.aj <- mean(out$crps_i[out$crps_i>0],na.rm=T)
  rbns.uc <- sum(out$ultimate)
  closed.uc<-sum(closed.claims$incPaid)
  
  results = rbind(results,
                  c(k,
                    actual.uc,
                    (closed.uc+ibnr.uc+rbns.uc)/actual.uc-1,
                    sum(cl.tot$ultimate)/actual.uc-1,
                    sqrt(var.aj),
                    cl.tot$msep,
                    crps.aj
                    )
  )
}

results <- results[-1,]
colnames(results) <- c("k",
                       "True amount",
                       "EI_ajr",
                       "EI_cl",
                       "msep_ajr",
                       "msep_cl",
                       "mean_crps")

results %>% xtable::xtable(digits=3) %>%print(include.rownames=FALSE)

fname <- paste0("C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\real_data_w_features_",
                'maxp_',
                maximum.p,
                "_",
                 format(Sys.time(), 
                 "%Y_%m_%d_%H_%M"),
                 ".csv")

fwrite(as.data.table(results),
       fname)


