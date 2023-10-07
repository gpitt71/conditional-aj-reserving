# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')


# Width 7 --- 

input.pbties <-c(-3.97809678024541e-06, 9.51943523703866e-07, 1.61328746890491e-08, 9.85506966190774e-10, 1.3547224236955e-10, 0, 3.00889940264394e-06,
                 0, -6.79208351973386e-06, 1.21175308945881e-06, 3.27015353265806e-08, 5.42628524413257e-09, 0, 5.54220260970434e-06,
                 0, 0, -5.52419512708527e-06, 2.64769885127375e-06, 1.90551855895997e-08, 0, 2.85744109022192e-06,
                 0, 0, 0, -6.40163496831553e-06, 4.21428114058748e-06, 0, 2.18735382772805e-06,
                 0, 0, 0, 0, -2.59150524018556e-06, 1.72767016012371e-06, 8.63835080061855e-07,
                 0, 0, 0, 0, 0, -4.85907232534793e-07, 4.85907232534793e-07,
                 0, 0, 0, 0, 0, 0, 0)

width <- 7

imx <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

jrate <- apply(imx,1,sum)
jmx <-  matrix(rep(jrate,width),byrow = F,ncol=width)

imx[1:(width-1),] <- imx[1:(width-1),]/jmx[1:(width-1),]
apply(imx,1,sum)

lambda <- function(t){
  A <- matrix(c(1*mark_dist(1, t, 0),
                1*mark_dist(2, t, 0),
                1*mark_dist(3, t, 0),
                1*mark_dist(4, t, 0),
                1*mark_dist(5, t, 0),
                1*mark_dist(6, t, 0),
                # 1*mark_dist(7, t, 0),
                # 1*mark_dist(8, t, 0),
                # 1*mark_dist(9, t, 0),
                rep(0, width)),
              nrow = width, ncol = width, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

## w.o. covariates ----
results <- matrix(ncol=5)
trial <- 1
ix2 <- 0


while(trial <20){
  set.seed(trial+ix2)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  
  # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
  # data.setC[,"X"] <- X
  
  fit <- AalenJohansen::aalen_johansen(data.list)
  yhat <-  sapply(fit$p, last)
  # cond <- yhat <=1 & yhat >=0
  # tail(yhat[cond])
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  if(max(yhat[!cond])>0.98){
    # yhat <- yhat[cond]
    
    trial<-trial+1
    yhat[cond]<-1
    x.vals <- fit$t
    # x.vals<- x.vals[cond]
    
    cond2 = diff(c(0,yhat))<0
    
    while(sum(cond2)!=0){
      
      yhat <-yhat[!cond2]
      x.vals <-x.vals[!cond2]
      cond2 = diff(c(0,yhat))<0
      
    }
    
    
    rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
    out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
                        variance=individual_vty_nf(k,x.vals,yhat)),by=.(id)]
    
    actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
    closed <- closed.data[,.(ultimate=last(times)),by=id]
    
    
    actual.tot <- sum(actual$ultimate)
    pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
    cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
    
    crps_data <- full_join(x = out, y = rbns.data, by = "id")[order(id),]
    out_crps <- crps_data[,.(crps_i=crps_computer(k=k,y=ultimate,x=x.vals,cdf_i=yhat)),by=.(id)]
    
    # 1-pred.tot/actual.tot
    # 1-cl.tot$ultimate/actual.tot
    
    if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}
    
    results = rbind(results,c(actual.tot,
                              pred.tot,
                              cl.tot$ultimate,
                              sum(out$variance),
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
  }
  
  
}

results=as.data.table(results[-1,])
colnames(results) <-c('actualtot','pred.tot','cl.tot','variance','crps')
locn = 'C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\%s'
fname = 'results_k7_newkNA.csv'
fwrite(results,sprintf(locn,fname))

mean(results$crps)
#3377406


# newk = 4 ----

results <- matrix(ncol=5)
trial <- 1
ix2 <- 0
newk = 4
while(trial <20){
  
  set.seed(trial+ix2)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  
  # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
  # data.setC[,"X"] <- X
  
  data.list2 <- lapply(data.list,reformulate_state_space,newk = 4)
  
  fit <- AalenJohansen::aalen_johansen(data.list2)
  yhat <-  sapply(fit$p, last)
  # cond <- yhat <=1 & yhat >=0
  # tail(yhat[cond])
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  if(max(yhat[!cond])>0.98){
    # yhat <- yhat[cond]
    
    trial<-trial+1
    yhat[cond]<-1
    x.vals <- fit$t
    # x.vals<- x.vals[cond]
    
    cond2 = diff(c(0,yhat))<0
    
    while(sum(cond2)!=0){
      
      yhat <-yhat[!cond2]
      x.vals <-x.vals[!cond2]
      cond2 = diff(c(0,yhat))<0
      
    }
    
    
    rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
    out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
                        variance=individual_vty_nf(k,x.vals,yhat)),by=.(id)]
    
    actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
    closed <- closed.data[,.(ultimate=last(times)),by=id]
    
    
    actual.tot <- sum(actual$ultimate)
    pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
    cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
    
    crps_data <- full_join(x = out, y = rbns.data, by = "id")[order(id),]
    out_crps <- crps_data[,.(crps_i=crps_computer(k=k,y=ultimate,x=x.vals,cdf_i=yhat)),by=.(id)]
    
    # 1-pred.tot/actual.tot
    # 1-cl.tot$ultimate/actual.tot
    
    # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}
    
    results = rbind(results,c(actual.tot,
                              pred.tot,
                              cl.tot$ultimate,
                              sum(out$variance),
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
  }
  
  
}

results=as.data.table(results[-1,])
colnames(results) <-c('actualtot','pred.tot','cl.tot','variance','crps')
locn = 'C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\%s'
fname = sprintf('results_k7_newk%s.csv',newk)
fwrite(results,sprintf(locn,fname))

mean(results$crps)
#3747450

# newk = 5----
results <- matrix(ncol=5)
trial <- 1
ix2 <- 0

newk = 5

while(trial <20){
  
  set.seed(trial+ix2)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  
  # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
  # data.setC[,"X"] <- X
  
  data.list2 <- lapply(data.list,reformulate_state_space,newk=5)
  
  fit <- AalenJohansen::aalen_johansen(data.list2)
  yhat <-  sapply(fit$p, last)
  # cond <- yhat <=1 & yhat >=0
  # tail(yhat[cond])
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  if(max(yhat[!cond])>0.98){
    # yhat <- yhat[cond]
    
    trial<-trial+1
    yhat[cond]<-1
    x.vals <- fit$t
    # x.vals<- x.vals[cond]
    
    cond2 = diff(c(0,yhat))<0
    
    while(sum(cond2)!=0){
      
      yhat <-yhat[!cond2]
      x.vals <-x.vals[!cond2]
      cond2 = diff(c(0,yhat))<0
      
    }
    
    
    rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
    out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
                        variance=individual_vty_nf(k,x.vals,yhat)),by=.(id)]
    
    actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
    closed <- closed.data[,.(ultimate=last(times)),by=id]
    
    
    actual.tot <- sum(actual$ultimate)
    pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
    cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
    
    crps_data <- full_join(x = out, y = rbns.data, by = "id")[order(id),]
    out_crps <- crps_data[,.(crps_i=crps_computer(k=k,y=ultimate,x=x.vals,cdf_i=yhat)),by=.(id)]
    
    # 1-pred.tot/actual.tot
    # 1-cl.tot$ultimate/actual.tot
    
    # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}
    
    results = rbind(results,c(actual.tot,
                              pred.tot,
                              cl.tot$ultimate,
                              sum(out$variance),
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
  }
  
  
}
results=as.data.table(results[-1,])
colnames(results) <-c('actualtot','pred.tot','cl.tot','variance','crps')
locn = 'C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\%s'
fname = sprintf('results_k7_newk%s.csv',newk)
fwrite(results,sprintf(locn,fname))

mean(results$crps)
#4060944

# newk = 6 ----
results <- matrix(ncol=5)
trial <- 1
ix2 <- 0

newk = 6
while(trial <20){
  
  set.seed(trial+ix2)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  
  # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
  # data.setC[,"X"] <- X
  
  data.list2 <- lapply(data.list,reformulate_state_space,newk = 6)
  
  fit <- AalenJohansen::aalen_johansen(data.list2)
  yhat <-  sapply(fit$p, last)
  # cond <- yhat <=1 & yhat >=0
  # tail(yhat[cond])
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  if(max(yhat[!cond])>0.98){
    # yhat <- yhat[cond]
    
    trial<-trial+1
    yhat[cond]<-1
    x.vals <- fit$t
    # x.vals<- x.vals[cond]
    
    cond2 = diff(c(0,yhat))<0
    
    while(sum(cond2)!=0){
      
      yhat <-yhat[!cond2]
      x.vals <-x.vals[!cond2]
      cond2 = diff(c(0,yhat))<0
      
    }
    
    
    rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
    out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
                        variance=individual_vty_nf(k,x.vals,yhat)),by=.(id)]
    
    actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
    closed <- closed.data[,.(ultimate=last(times)),by=id]
    
    
    actual.tot <- sum(actual$ultimate)
    pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
    cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
    
    crps_data <- full_join(x = out, y = rbns.data, by = "id")[order(id),]
    out_crps <- crps_data[,.(crps_i=crps_computer(k=k,y=ultimate,x=x.vals,cdf_i=yhat)),by=.(id)]
    
    # 1-pred.tot/actual.tot
    # 1-cl.tot$ultimate/actual.tot
    
    # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}
    
    results = rbind(results,c(actual.tot,
                              pred.tot,
                              cl.tot$ultimate,
                              sum(out$variance),
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
  }
  
  
}

results=as.data.table(results[-1,])
colnames(results) <-c('actualtot','pred.tot','cl.tot','variance','crps')
locn = 'C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\%s'
fname = sprintf('results_k7_newk%s.csv',newk)
fwrite(results,sprintf(locn,fname))

mean(results$crps)
# 4390846













