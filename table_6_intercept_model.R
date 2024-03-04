rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('conditional-aj-reserving\\helper_functions_ajr.R')

# Construct the model ----

results = matrix(ncol=9)

for(maximum.p in c(4,5,6,7)){
  
  df = read.and.pp.data(fname='final_claim_data.csv')
  
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

actual_triangle <- incr2cum(as.triangle(df,
                                        origin='accident_period',
                                        dev='development_period',
                                        value='incPaid'))
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
                                           ay=last(accident_period),maximum.p=maximum.p),
                               by=Claim_number]


data2fit <- rbind(rbns.claims.fit,
                  closed.claims.fit)

data2fit[,times:=cumsum(times),by=Claim_number]

colnames(data2fit)[colnames(data2fit)=="events"] <- "states"

data2fit_list <-  split(data2fit,data2fit$Claim_number)

data2fit <- lapply(data2fit_list, function(x){as.list(rbind(data.frame(times=0.,states=1),x[,c('times','states')]))})

data2fit <- unname(data2fit)

fit1 <- AalenJohansen::aalen_johansen(data2fit)

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

rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]

actual=df %>%
  ungroup() %>%
  filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
  group_by(Claim_number) %>%
  reframe(true.ultimate=sum(incPaid)) %>%
  as.data.table()

rbns.data <-merge(actual,rbns.data,by='Claim_number')

out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_nf(k,
                                                                           true.ultimate=true.ultimate,
                                                                           x.vals=x.vals,
                                                                           yhat=yhat), 
                 by=.(Claim_number)]


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

ev.est <- cev(0,x=x.vals,y=yhat)
ev2.est <- 2*cev2(0,x=x.vals,y=yhat)

# estimated uc

tmp.tr <- clean.lt(reporting_period_tr_cum)
mack.rdel <- cl.calculator(tmp.tr) 

ibnr.uc <- ibnr_num*ev.est
ibnr.var <- ibnr_num*(ev2.est-ev.est^2) + ev.est^2*(mack.rdel$process.error.i^2)

# crps for IBNR

ibnr.claims<-ibnr.claims %>%
  group_by(Claim_number) %>%
  reframe(true.ultimate=floor(sum(incPaid))) %>% 
  as.data.table()

ibnr.claims<-ibnr.claims[,c('crps_i'):=ibnr.crps.computer(true.ultimate=true.ultimate,
                                                          x.vals=x.vals,
                                                          yhat=yhat),by=.(Claim_number)]

rbns.uc <- sum(out$ultimate)
closed.uc<-sum(closed.claims$incPaid)

actual.uc=df %>%
  ungroup() %>%
  filter(accident_period<=(maximum.p-2)&development_period<=(maximum.p-1)) %>%
  summarise(uc=sum(incPaid)) %>% unlist()


cl.data <- find.t.data(df)
cl.tot <- cl.calculator(incr2cum(cl.data))

results = rbind(results,
                c(maximum.p,
                  actual.uc,
                  (closed.uc+ibnr.uc+rbns.uc)/actual.uc-1,
                  sum(cl.tot$ultimate)/actual.uc-1,
                  sqrt(sum(c(out$variance,ibnr.var))),
                  cl.tot$process.error.i,
                  mean(c(out$crps_i, ibnr.claims$crps_i),na.rm=T),
                  (closed.uc+ibnr.uc+rbns.uc)-diag_ftr,
                  sum(cl.tot$ultimate)-diag_ftr
                  
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
                       "mean_crps",
                       "reserve_aj",
                       "reserve_cl")


results %>% xtable::xtable(digits=3) %>%print( include.rownames=FALSE)


fname <- paste0("conditional-aj-reserving\\results_csv\\real_data_no_features_",
                "_",
                format(Sys.time(), 
                       "%Y_%m_%d_%H_%M"),
                ".csv")

fwrite(as.data.table(results),
       fname)


