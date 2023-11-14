rm(list=ls())
set.seed(1)
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)

source('conditional-aj-reserving\\helper_functions_ajr.R')

# Width 4 ----
input.pbties <- c(-7.11758383885783e-06, 1.39183828093592e-06, 5.73888289022062e-08, 5.66835672901971e-06,
                  0, -1.16363390439916e-05, 1.9869412608886e-06, 9.64939778310298e-06,
                  0, 0, -4.55533160623821e-06, 4.55533160623821e-06,
                  0, 0 ,0, 0)

results <- matrix(ncol=6)
trial <- 1
ix2 = 0
V=1200
width=4


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


ap.volumes <- rep(V,width-1) -seq(0,100*(width-2),by=100)


{set.seed(1)
  
  X <- sort(1/sample(11:(width-1+10),sum(ap.volumes),replace=T))
  
  }
for(trial in 1:20){
  set.seed(trial)
  
  sim.2.fit <- simulate_mkvian_dataX(ap.volumes=ap.volumes,
                                     tsh=1e-06,
                                     stdize=1,
                                     X=X)
  
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  

  
  
  rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
  rbns.data[["X"]] <- X[rbns.data$id]
  
  models.list <- compute.ajr.models(rbns.data,
                                    X='X',
                                    data.list=data.list)
  
  actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
  
  rbns.data <-merge(actual,rbns.data,by='id')
  colnames(rbns.data)[2] <- 'true.ultimate'
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=X,
                                                                                               models.list), 
                   by=.(id)]
  
  
  
  closed <- closed.data[,.(ultimate=last(times)),by=id]
  
  actual.tot <- sum(actual$ultimate)
  pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
  cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
  
  results = rbind(results, c(actual.tot,
                             pred.tot,
                             cl.tot$ultimate,
                             sqrt(sum(out$variance)),
                             cl.tot$process.se,
                             mean(out$crps_i)))
  
}

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'aj.msep',
                       'cl.msep',
                       'mcrps')

locn = "conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
                'aytrend',
                width,
                "_",
                format(Sys.time(), 
                       "%Y_%m_%d_%H_%M"),
                ".csv")

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    msep= mean(aj.msep),
                    msep.cl=mean(cl.msep),
                    aj.error=mean(pred.tot/simulated.tot)-1,
                    cl.error=mean(cl.tot/simulated.tot)-1,
)

fwrite(results, file=sprintf(locn,finame))

# Width 5 ----




input.pbties <-c(-3.87473643786946e-06, 7.75530168016601e-07, 1.92727077396385e-08, 9.37832380696512e-10, 3.07899572973252e-06,
                 0, -6.0461063538271e-06, 1.40368477593149e-06, 3.24343294024007e-08, 4.60998724849321e-06,
                 0, 0, -4.59531201287676e-06, 2.33954270788095e-06, 2.25576930499581e-06,
                 0, 0, 0, -2.37851748950939e-06, 2.37851748950939e-06,
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

V=1200
results <- matrix(ncol=6)
ap.volumes <- rep(V,width-1) -seq(0,100*(width-2),by=100)


{set.seed(1)
  X <- sort(1/sample(11:(width-1+10),sum(ap.volumes),replace=T))
}
for(trial in 1:20){
  set.seed(trial)
  
  sim.2.fit <- simulate_mkvian_dataX(ap.volumes=ap.volumes,
                                     tsh=1e-06,
                                     stdize=1,
                                     X=X)
  
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  

  
  rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
  rbns.data[["X"]] <- X[rbns.data$id]
  
  
  models.list <- compute.ajr.models(rbns.data,
                                    X='X',
                                    data.list=data.list)
  
  actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
  
  rbns.data <-merge(actual,rbns.data,by='id')
  colnames(rbns.data)[2] <- 'true.ultimate'

  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=X,
                                                                                               models.list), 
                   by=.(id)]
  
  
  
  closed <- closed.data[,.(ultimate=last(times)),by=id]
  
  actual.tot <- sum(actual$ultimate)
  pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
  cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
  
  results = rbind(results, c(actual.tot,
                             pred.tot,
                             cl.tot$ultimate,
                             sqrt(sum(out$variance)),
                             cl.tot$process.se,
                             mean(out$crps_i)))
  
}

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'aj.msep',
                       'cl.msep',
                       'mcrps')

locn = "conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'aytrend',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

results %>% reframe(pred.tot=mean(pred.tot),
                    simulated.tot=mean(actual.tot),
                    cl.ultimate=mean(cl.tot),
                    msep= mean(aj.msep),
                    msep.cl=mean(cl.msep),
                    aj.error=mean(pred.tot/simulated.tot)-1,
                    cl.error=mean(cl.tot/simulated.tot)-1,
)

fwrite(results, file=sprintf(locn,finame))



# Width 6 ----

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

V=1200
results <- matrix(ncol=6)
ap.volumes <- rep(V,width-1)-seq(0,100*(width-2),by=100)


{set.seed(1)
  X <- sort(1/sample(11:(width-1+10),sum(ap.volumes),replace=T))
}
for(trial in 1:20){
  set.seed(trial)
  
  sim.2.fit <- simulate_mkvian_dataX(ap.volumes=ap.volumes,
                                     tsh=1e-06,
                                     stdize=1,
                                     X=X)
  
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()

  
  rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
  rbns.data[["X"]] <- X[rbns.data$id]
  
  models.list <- compute.ajr.models(rbns.data,
                                    X='X',
                                    data.list=data.list)
  
  actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
  
  rbns.data <-merge(actual,rbns.data,by='id')
  colnames(rbns.data)[2] <- 'true.ultimate'

  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=X,
                                                                                               models.list), 
                   by=.(id)]
  
  

  closed <- closed.data[,.(ultimate=last(times)),by=id]
  
  actual.tot <- sum(actual$ultimate)
  pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
  cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
  
  results = rbind(results, c(actual.tot,
                             pred.tot,
                             cl.tot$ultimate,
                             sqrt(sum(out$variance)),
                             cl.tot$process.se,
                             mean(out$crps_i)))
  
}

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'aj.msep',
                       'cl.msep',
                       'mcrps')

locn = "conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'aytrend',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")

fwrite(results, file=sprintf(locn,finame))




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



V=1200
results <- matrix(ncol=6)
ap.volumes <- rep(V,width-1)-seq(0,100*(width-2),by=100)

{set.seed(1)
  X <- sort(1/sample(11:(width-1+10),sum(ap.volumes),replace=T))
}
for(trial in 1:20){
  set.seed(trial)
  
  sim.2.fit <- simulate_mkvian_dataX(ap.volumes=ap.volumes,
                                     tsh=1e-06,
                                     stdize=1,
                                     X=X)
  
  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit
  
  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()
  

  
  
  rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
  rbns.data[["X"]] <- X[rbns.data$id]
  models.list <- compute.ajr.models(rbns.data,
                                    X='X',
                                    data.list=data.list)
  
  actual <- input.data$actual.data[,.(ultimate=last(times)),by=id]
  
  rbns.data <-merge(actual,rbns.data,by='id')
  colnames(rbns.data)[2] <- 'true.ultimate'
  # actual[rbns.data, on = 'id', true.ultimate := ultimate]
  
  out <- rbns.data[,c('ultimate','variance','crps_i'):=individual_results_X_w_pre_saved_models(k,
                                                                                               true.ultimate=true.ultimate,
                                                                                               X=X,
                                                                                               models.list), 
                   by=.(id)]
  
  # 
  closed <- closed.data[,.(ultimate=last(times)),by=id]
  
  actual.tot <- sum(actual$ultimate)
  pred.tot <- sum(out$ultimate) +sum(closed$ultimate)
  cl.tot <- chain_ladder_computer(input.data$actual.data,width=width)
  
  results = rbind(results, c(actual.tot,
                             pred.tot,
                             cl.tot$ultimate,
                             sqrt(sum(out$variance)),
                             cl.tot$process.se,
                             mean(out$crps_i)))
  
}

results <- data.frame(results[-1,])
colnames(results) <- c('actual.tot',
                       'pred.tot',
                       'cl.tot',
                       'aj.msep',
                       'cl.msep',
                       'mcrps')

locn = "conditional-aj-reserving\\results_csv\\%s"
finame=paste0("simulation_",
              'aytrend',
              width,
              "_",
              format(Sys.time(), 
                     "%Y_%m_%d_%H_%M"),
              ".csv")


fwrite(results, file=sprintf(locn,finame))


