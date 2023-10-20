rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

# Width 4 ----

#yearly
input.pbties <- c(-7.11758383885783e-06, 1.39183828093592e-06, 5.73888289022062e-08, 5.66835672901971e-06,
                  0, -1.16363390439916e-05, 1.9869412608886e-06, 9.64939778310298e-06,
                  0, 0, -4.55533160623821e-06, 4.55533160623821e-06,
                  0, 0 ,0, 0)
width <- 4


lambda <- function(t){
  A <- matrix(c(1*mark_dist(1, t, 0),
                1*mark_dist(2, t, 0),
                1*mark_dist(3, t, 0),
                # 1*mark_dist(4, t, 0),
                # 1*mark_dist(5, t, 0),
                rep(0, width)),
              nrow = width, ncol = width, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

imx <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

jrate <- apply(imx,1,sum)
jmx <-  matrix(rep(jrate,width),byrow = F,ncol=width)

imx[1:(width-1),] <- imx[1:(width-1),]/jmx[1:(width-1),]


# Results for model w.o. covariates ----
results <- matrix(ncol=6)
trial <- 1
ix2 = 0
V=1200

while(trial < 20){
  
  set.seed(trial+ix2)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                    tsh=1e-06,
                                    stdize=1)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter("4"%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!("4"%chin%states)) %>% as.data.table()
  
  # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
  # data.setC[,"X"] <- X
  
  fit <- AalenJohansen::aalen_johansen(data.list)
  # yhat <-  sapply(fit$p, last)
  yhat <-  sapply(fit$p, last)
  
  cond <- (yhat > 1 | yhat <0) | is.na(yhat)
  print(max(yhat[!cond]))
  # if(max(yhat[!cond])>0.0000098){
    # print(max(yhat[!cond]))
    yhat[cond] <- 1
    # tail(yhat)
    x.vals <- fit$t
    # x.vals <-x.vals[cond]
    
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
    1-pred.tot/actual.tot
    1-cl.tot$ultimate/actual.tot
    
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
                              cl.tot$ultimate.se,
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
    
    trial <- trial+1
  # }
  
}

L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

alpha <- c(1,rep(0,width-2))
ev=sum(alpha%*%solve(-L[1:(width-1),1:(width-1)]))

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'variance',
                       'cl.msep',
                       'mcrps')

results$true.ultimate=ev*sum(V*(width-1))

locn = "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'no_features_all_records_',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

fwrite(results, file=sprintf(locn,finame))

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    exact.tot=ev*sum(V*(width-1)),
                    msep= mean(sqrt(variance)),
                    msep.cl=mean(cl.msep),
                    aj.error=1-mean(pred.tot/exact.tot),
                    cl.error=1-mean(cl.tot/exact.tot),
) %>% xtable::xtable()


# Width 5 ----

results <- matrix(ncol=6)

input.pbties <- c(-1.41134747023148e-05, 5.818448310349e-07, 2.64574079945067e-06, 3.21143958772718e-07, 1.05647451130565e-05,
                  0, -1.56578768719176e-05, 9.9084450079893e-06, 3.56570861910867e-07, 5.3928610020174e-06,
                  0, 0, -2.16284547098169e-05, 1.35437786537905e-05, 8.08467605602645e-06,
                  0, 0, 0, -1.3274315082635e-05, 1.3274315082635e-05,
                  0, 0, 0, 0, 0)

width <- 5

imx <- matrix(input.pbties,
              byrow = T,
              ncol=width)

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
                # 1*mark_dist(5, t, 0),
                rep(0, width)),
              nrow = width, ncol = width, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

trial <- 1
ix2 = 0
V=1200
while( trial < 20){
  
  set.seed(trial+ix2)
  # browser()
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                    tsh=1e-06,
                                    stdize=1)
  
  print(trial)
  ix2=ix2+1
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  
  
  fit <- AalenJohansen::aalen_johansen(data.list)
  yhat <-  sapply(fit$p, last)
  
  cond <- (yhat > 1 | yhat <0) | is.na(yhat)
  print(max(yhat[!cond]))
  
  # if(max(yhat[!cond])>0.0000098){
    # print(max(yhat[!cond]))
    yhat[cond] <- 1
    # tail(yhat)
    x.vals <- fit$t
    # x.vals <-x.vals[cond]
    
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
    1-pred.tot/actual.tot
    1-cl.tot$ultimate/actual.tot
    
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
                              cl.tot$ultimate.se,
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
    
    trial <- trial+1
  # }
  
}

L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

alpha <- c(1,rep(0,width-2))
ev=sum(alpha%*%solve(-L[1:(width-1),1:(width-1)]))

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot','pred.tot','cl.tot','variance','cl.msep','mcrps')

results$true.ultimate=ev*sum(V*(width-1))


locn = "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'no_features_all_records_',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

fwrite(results, file=sprintf(locn,finame))

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    exact.tot=ev*sum(V*(width-1)),
                    msep= mean(sqrt(variance)),
                    msep.cl=mean(cl.msep),
                    aj.error=1-mean(pred.tot/exact.tot),
                    cl.error=1-mean(cl.tot/exact.tot),
) %>% xtable::xtable()

# Width 6 -----

# input.pbties <-c(-1.26906163922617e-05, 5.83803794444761e-07, 2.23116423848226e-06, 1.00444539991587e-06, 1.34096859073061e-06, 4.80112100634515e-08, 3.61577265279923e-09, 7.4786073859719e-06,
#                  0, -1.20638161771216e-05, 5.84747900910062e-06, 5.4023737620222e-07, 5.28391796260648e-07, 2.00619513056318e-07, 0, 4.94708848250178e-06,
#                  0, 0, -2.33015116551299e-05, 1.450947061937e-05, 5.3679278175145e-07, 0, 0, 8.25524825400848e-06,
#                  0, 0, 0, -2.01129042710382e-05, 9.74119421745867e-06, 2.3189558073279e-06, 0, 8.05275424625159e-06,
#                  0, 0, 0, 0, -3.23460384127725e-05, 1.99107150573278e-05, 1.94446604962277e-06, 1.04908573058219e-05,
#                  0, 0, 0, 0, 0, -2.72225246947188e-05, 1.90125569296449e-05, 8.20996776507393e-06,
#                  0, 0, 0, 0, 0, 0, -3.23938155023196e-07, 3.23938155023196e-07,
#                  0, 0, 0, 0, 0, 0, 0, 0)


input.pbties <-c(-3.92969143285773e-06, 9.46146090007938e-07, 1.63469251967611e-08, 1.05743452770889e-09, 9.82949116106126e-11, 2.96604268821371e-06,
                 0, -6.26893145831146e-06, 1.25354737215752e-06, 3.00799715378682e-08, 0, 4.98530411461608e-06,
                 0, 0, -4.78554760053172e-06, 2.34236472738395e-06, 0, 2.44318287314777e-06,
                 0, 0, 0, -4.02223209153801e-06, 2.21897636190889e-06, 1.80325572962912e-06,
                 0, 0, 0, 0, -9.71814465069587e-07, 9.71814465069587e-07,
                 0, 0, 0, 0, 0, 0)

width <- 6

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
                # 1*mark_dist(6, t, 0),
                # 1*mark_dist(7, t, 0),
                # 1*mark_dist(8, t, 0),
                rep(0, width)),
              nrow = width, ncol = width, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}


# Results for model w.o. covariates ----
results <- matrix(ncol=6)


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
  
  # if(max(yhat[!cond])>0.0000098){
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
                              cl.tot$ultimate.se,
                              mean(out_crps$crps_i[out_crps$crps_i>0],
                                   na.rm=T)))
  # }
  
  
}


# for(trial in 1:20){
#   set.seed(trial)
#   sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
#
#   input.data <- sim.2.fit
#   data.list <- input.data$sim.2.fit
#
#   rbns.data <- input.data$censored.data
#   rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
#   closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
#   rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
#
#   # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
#   # data.setC[,"X"] <- X
#
#   fit <- AalenJohansen::aalen_johansen(data.list)
#   yhat <-  sapply(fit$p, last)
#   # tail(yhat)
#
#   x.vals <- fit$t
#
#   # cond <- yhat <=1 & yhat >=0 & !is.na(yhat)
#
#   cond <- (yhat >1 | yhat <0) | is.na(yhat)
#   yhat[cond] <- 1
#   # yhat[cond] <- 1
#   # yhat <- yhat[cond]
#   # x.vals<-x.vals[cond]
#
#
#   rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
#   out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
#                       variance=2*cev2(k,x.vals,yhat)-cev(k,x.vals,yhat)^2),by=.(id)]
#
#   actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
#   closed <- closed.data[,.(ultimate=last(times)),by=id]
#
#
#   actual.tot <- sum(actual$ultimate)
#   pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
#   cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
#
#   results = rbind(results,c(actual.tot,
#                             pred.tot,
#                             cl.tot$ultimate,
#                             cl.tot$ultimate.se,
#                             sqrt(sum(out$variance))))
# }

L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

alpha <- c(1,rep(0,width-2))
ev=sum(alpha%*%solve(-L[1:(width-1),1:(width-1)]))

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'variance',
                       'cl.msep',
                       'mcrps')

results$true.ultimate=ev*sum(1200*(width-1))

locn = "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'no_features_all_records_',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

fwrite(results, file=sprintf(locn,finame))

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    exact.tot=ev*sum(V*(width-1)),
                    msep= mean(sqrt(variance)),
                    aj.error=1-mean(pred.tot/exact.tot),
                    cl.error=1-mean(cl.tot/exact.tot),
) %>% xtable::xtable()
1-pred.tot/actual.tot
1-cl.tot$ultimate/actual.tot

# Width 7 ----

# model parameters

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
  
  # if(max(yhat[!cond])>0.0000098){
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
  # }
  
  
}

L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0

alpha <- c(1,rep(0,width-2))
ev=sum(alpha%*%solve(-L[1:(width-1),1:(width-1)]))

results <- data.frame(results[-1,])

colnames(results) <- c('actual.tot','pred.tot','cl.tot','variance','mcrps')

results$true.ultimate=ev*sum(1200*(width-1))

locn = "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'no_features_all_records_',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

fwrite(results, file=sprintf(locn,finame))

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    exact.tot=ev*sum(V*(width-1)),
                    msep= mean(sqrt(variance)),
                    aj.error=mean(pred.tot/exact.tot)-1,
                    cl.error=mean(cl.tot/exact.tot)-1,
) %>% xtable::xtable()


write.table(results, file="C:\\Users\\gpitt\\Documents\\Phd\\visiting\\markov-chains-reserving\\results\\w7nf.csv")

png(filename = "C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\w7nf.png", width=600, height=480, res=72)
understand.the.densities(results)
dev.off()



# Width 8 ----

# input.pbties <-c(-4.02104115693941e-06, 9.41733693345327e-07, 1.45337940215019e-08, 1.0284444599133e-09, 1.28667155196585e-10, 1.41137223345763e-11, 0, 3.06360244423514e-06,
#                  0, -6.79581039063144e-06, 1.38541363530372e-06, 3.02422921453848e-08, 4.66191880121944e-09, 0, 3.70408553883232e-09, 5.37178845884229e-06,
#                  0, 0, -5.97028986443137e-06, 2.59887586380591e-06, 2.77661275734168e-08, 2.3138439644514e-08, 1.19977094453035e-08, 3.30851172396223e-06,
#                  0, 0, 0, -5.23649767140797e-06, 3.19644689000269e-06, 6.99411925618263e-08, 0, 1.97010958884345e-06,
#                  0, 0, 0, 0, -3.16996623129841e-06, 1.64437177740346e-06, 0, 1.52559445389495e-06,
#                  0, 0, 0, 0, 0, -2.10559800765077e-06, 1.13378354258118e-06, 9.71814465069587e-07,
#                  0, 0, 0, 0, 0, 0, -3.23938155023196e-07, 3.23938155023196e-07,
#                  0, 0, 0, 0, 0, 0, 0, 0)

# input.pbties <-c(-1.26906163922617e-05, 5.83803794444761e-07, 2.23116423848226e-06, 1.00444539991587e-06, 1.34096859073061e-06, 4.80112100634515e-08, 3.61577265279923e-09, 7.4786073859719e-06,
#                  0, -1.20638161771216e-05, 5.84747900910062e-06, 5.4023737620222e-07, 5.28391796260648e-07, 2.00619513056318e-07, 0, 4.94708848250178e-06,
#                  0, 0, -2.33015116551299e-05, 1.450947061937e-05, 5.3679278175145e-07, 0, 0, 8.25524825400848e-06,
#                  0, 0, 0, -2.01129042710382e-05, 9.74119421745867e-06, 2.3189558073279e-06, 0, 8.05275424625159e-06,
#                  0, 0, 0, 0, -3.23460384127725e-05, 1.99107150573278e-05, 1.94446604962277e-06, 1.04908573058219e-05,
#                  0, 0, 0, 0, 0, -2.72225246947188e-05, 1.90125569296449e-05, 8.20996776507393e-06,
#                  0, 0, 0, 0, 0, 0, -3.23938155023196e-07, 3.23938155023196e-07,
#                  0, 0, 0, 0, 0, 0, 0, 0)
# 
# width <- 8
# 
# imx <- matrix(input.pbties, byrow = T,ncol=width)
# diag(imx) <- 0
# 
# jrate <- apply(imx,1,sum)
# jmx <-  matrix(rep(jrate,width),byrow = F,ncol=width)
# 
# imx[1:(width-1),] <- imx[1:(width-1),]/jmx[1:(width-1),]
# apply(imx,1,sum)
# 
# lambda <- function(t){
#   A <- matrix(c(1*mark_dist(1, t, 0),
#                 1*mark_dist(2, t, 0),
#                 1*mark_dist(3, t, 0),
#                 1*mark_dist(4, t, 0),
#                 1*mark_dist(5, t, 0),
#                 1*mark_dist(6, t, 0),
#                 1*mark_dist(7, t, 0),
#                 # 1*mark_dist(8, t, 0),
#                 rep(0, width)),
#               nrow = width, ncol = width, byrow = TRUE)
#   diag(A) <- -rowSums(A)
#   A
# }
# 
# 
# # Results for model w.o. covariates ----
# results <- matrix(ncol=5)
# 
# 
# trial <- 1
# ix2 <- 0
# 
# while(trial <20){
#   set.seed(trial+ix2)
#   sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
#   ix2=ix2+1
#   input.data <- sim.2.fit
#   data.list <- input.data$sim.2.fit
#   
#   rbns.data <- input.data$censored.data
#   rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
#   closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
#   rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
#   
#   # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
#   # data.setC[,"X"] <- X
#   
#   fit <- AalenJohansen::aalen_johansen(data.list)
#   yhat <-  sapply(fit$p, last)
#   # cond <- yhat <=1 & yhat >=0
#   # tail(yhat[cond])
#   cond <- (yhat >1 | yhat <0) | is.na(yhat)
#   
#   if(max(yhat[!cond])>0.0000098){
#     # yhat <- yhat[cond]
#     
#     trial<-trial+1
#     yhat[cond]<-1
#     x.vals <- fit$t
#     # x.vals<- x.vals[cond]
#     
#     cond2 = diff(c(0,yhat))<0
#     
#     while(sum(cond2)!=0){
#       
#       yhat <-yhat[!cond2]
#       x.vals <-x.vals[!cond2]
#       cond2 = diff(c(0,yhat))<0
#       
#     }
#     
#     
#     rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
#     out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
#                         variance=individual_vty_nf(k,x.vals,yhat)),by=.(id)]
#     
#     actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
#     closed <- closed.data[,.(ultimate=last(times)),by=id]
#     
#     
#     actual.tot <- sum(actual$ultimate)
#     pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
#     cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
#     
#     crps_data <- full_join(x = out, y = rbns.data, by = "id")[order(id),]
#     out_crps <- crps_data[,.(crps_i=crps_computer(k=k,y=ultimate,x=x.vals,cdf_i=yhat)),by=.(id)]
#     
#     # 1-pred.tot/actual.tot
#     # 1-cl.tot$ultimate/actual.tot
#     
#     if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}
#     
#     results = rbind(results,c(actual.tot,
#                               pred.tot,
#                               cl.tot$ultimate,
#                               sum(out$variance),
#                               mean(out_crps$crps_i[out_crps$crps_i>0],
#                                    na.rm=T)))
#   }
#   
#   
# }


# for(trial in 1:20){
#   set.seed(trial)
#   sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),tsh=1e-06,stdize=1)
#
#   input.data <- sim.2.fit
#   data.list <- input.data$sim.2.fit
#
#   rbns.data <- input.data$censored.data
#   rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
#   closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
#   rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
#
#   # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
#   # data.setC[,"X"] <- X
#
#   fit <- AalenJohansen::aalen_johansen(data.list)
#   yhat <-  sapply(fit$p, last)
#   # tail(yhat)
#
#   x.vals <- fit$t
#
#   # cond <- yhat <=1 & yhat >=0 & !is.na(yhat)
#
#   cond <- (yhat >1 | yhat <0) | is.na(yhat)
#   yhat[cond] <- 1
#   # yhat[cond] <- 1
#   # yhat <- yhat[cond]
#   # x.vals<-x.vals[cond]
#
#
#   rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
#   out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat),
#                       variance=2*cev2(k,x.vals,yhat)-cev(k,x.vals,yhat)^2),by=.(id)]
#
#   actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
#   closed <- closed.data[,.(ultimate=last(times)),by=id]
#
#
#   actual.tot <- sum(actual$ultimate)
#   pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
#   cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
#
#   results = rbind(results,c(actual.tot,
#                             pred.tot,
#                             cl.tot$ultimate,
#                             cl.tot$ultimate.se,
#                             sqrt(sum(out$variance))))
# }
# 
# L <- matrix(input.pbties, byrow = T,ncol=width)
# diag(imx) <- 0
# 
# alpha <- c(1,rep(0,width-2))
# ev=sum(alpha%*%solve(-L[1:(width-1),1:(width-1)]))
# 
# results <- data.frame(results[-1,])
# colnames(results) <- c('actual.tot','pred.tot','cl.tot','variance','mcrps')
# 
# results$true.ultimate=ev*sum(1200*(width-1))
# 
# locn = "C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\results_csv\\%s"
# finame=paste0("simulation_",
#               'no_features_all_records_',
#               width,
#               "_",
#               format(Sys.time(), 
#                      "%Y_%m_%d_%H_%M"),
#               ".csv")
# 
# fwrite(results, file=sprintf(locn,finame))
# 
# results %>% reframe(pred.tot=mean(pred.tot),
#                     simulated.tot=mean(actual.tot),
#                     cl.ultimate=mean(cl.tot),
#                     exact.tot=ev*sum(V*(width-1)),
#                     msep= mean(sqrt(variance)),
#                     aj.error=1-mean(pred.tot/exact.tot),
#                     cl.error=1-mean(cl.tot/exact.tot),
# ) %>% xtable::xtable()
# 1-pred.tot/actual.tot
# 1-cl.tot$ultimate/actual.tot
















