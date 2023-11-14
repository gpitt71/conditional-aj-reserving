rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)
library(scales)
source('conditional-aj-reserving\\helper_functions_ajr.R')
unit=1e+06

# width 4 ----

{
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
  results <- matrix(ncol=4)
set.seed(42)
  sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(120,width-1),
                                  tsh=1e-06,
                                  stdize=1)

  input.data <- sim.2.fit
  data.list <- input.data$sim.2.fit

  rbns.data <- input.data$censored.data
  rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
  closed.data <- rbns.data %>% group_by(id) %>% filter("4"%chin%states) %>% as.data.table()
  rbns.data <- rbns.data %>% group_by(id) %>% filter(!("4"%chin%states)) %>% as.data.table()


  fit <- AalenJohansen::aalen_johansen(data.list)
  yhat <-  sapply(fit$p, last)
  x.vals <- fit$t


sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(1200,width-1),
                                  tsh=1e-06,
                                  stdize=1)

input.data <- sim.2.fit
data.list <- input.data$sim.2.fit

rbns.data <- input.data$censored.data
rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
closed.data <- rbns.data %>% group_by(id) %>% filter("4"%chin%states) %>% as.data.table()
rbns.data <- rbns.data %>% group_by(id) %>% filter(!("4"%chin%states)) %>% as.data.table()


fit1 <- AalenJohansen::aalen_johansen(data.list)
yhat1 <-  sapply(fit1$p, last)
x.vals1 <- fit1$t

sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(12000,width-1),
                                  tsh=1e-06,
                                  stdize=1)

input.data <- sim.2.fit
data.list <- input.data$sim.2.fit

rbns.data <- input.data$censored.data
rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
closed.data <- rbns.data %>% group_by(id) %>% filter("4"%chin%states) %>% as.data.table()
rbns.data <- rbns.data %>% group_by(id) %>% filter(!("4"%chin%states)) %>% as.data.table()

fit2 <- AalenJohansen::aalen_johansen(data.list)
yhat2 <-  sapply(fit2$p, last)
x.vals2 <- fit2$t

}


width<-p<-4
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2), 
  y=c(yhat,yhat1,yhat2,ytrue2), 
  lbl=c(rep('Fitted 120',length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        rep('Fitted 12000',length(yhat2)),
        rep('True curve',length(yhat2)))
  )


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2), 
  y=c(yhat,yhat1,ytrue2), 
  lbl=c(rep('Fitted 120',length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        rep('True curve',length(yhat2)))
)

cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x/unit,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")

ggsave(" severity_curve4.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()


# Width 5 ----

{

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

  set.seed(42)

  trial <- 1
  ix2 <- 0
  V <- 120

  while( trial < 2){

    set.seed(trial+ix2)
    # browser()
    sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                      tsh=1e-06,
                                      stdize=1)
    # print(trial)
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

      yhat[cond] <- 1
      x.vals <- fit$t

      rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
      out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat)),by=.(id)]

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

      results = rbind(results,c(actual.tot,pred.tot,cl.tot$ultimate,mean(out_crps$crps_i,
                                                                         na.rm=T)))

      trial <- trial+1

      cond2 = diff(c(0,yhat))<0

      while(sum(cond2)!=0){

        yhat <-yhat[!cond2]
        x.vals <-x.vals[!cond2]
        cond2 = diff(c(0,yhat))<0

      }
    


  }

  trial <- 1
  ix2 <- 0
  V <- 1200

  while( trial < 2){

    set.seed(trial+ix2)

    sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                      tsh=1e-06,
                                      stdize=1)

    ix2=ix2+1
    input.data <- sim.2.fit
    data.list <- input.data$sim.2.fit

    rbns.data <- input.data$censored.data
    rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
    closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
    rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()


    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)

    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)

      yhat1[cond] <- 1

      x.vals1 <- fit1$t

      trial <- trial+1
   

    cond2 = diff(c(0,yhat1))<0

    while(sum(cond2)!=0){

      yhat1 <-yhat1[!cond2]
      x.vals1 <-x.vals1[!cond2]
      cond2 = diff(c(0,yhat1))<0

    }


  }

  trial <- 1
  ix2 <- 0
  V <- 12000

  while( trial < 2){

    set.seed(trial+ix2)
    # browser()
    sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                      tsh=1e-06,
                                      stdize=1)
    # print(trial)
    ix2=ix2+1
    input.data <- sim.2.fit
    data.list <- input.data$sim.2.fit

    rbns.data <- input.data$censored.data
    rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
    closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
    rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()


    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)

    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)

      yhat2[cond] <- 1
      
      x.vals2 <- fit2$t
      

      trial <- trial+1

    cond2 = diff(c(0,yhat2))<0

    while(sum(cond2)!=0){

      yhat2 <-yhat2[!cond2]
      x.vals2 <-x.vals2[!cond2]
      cond2 = diff(c(0,yhat2))<0

    }


  }




}



width<-p<-5
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")

ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),
  y=c(yhat,yhat1,yhat2,ytrue2), 
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))
)


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2),
  y=c(yhat,yhat1,ytrue2), 
  lbl=c(rep("Fitted 120",length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        
        rep('True curve',length(yhat2)))
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x/unit,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")

ggsave(" severity_curve5.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()


# Width 6 ----


{
  
  
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
                  rep(0, width)),
                nrow = width, ncol = width, byrow = TRUE)
    diag(A) <- -rowSums(A)
    A
  }
  
  set.seed(42)
  
  trial <- 1
  ix2 <- 0
  V <- 120
  
  while( trial < 2){
    
    set.seed(trial+ix2)
    
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

      yhat[cond] <- 1
      
      x.vals <- fit$t
      
      
      rbns.data <- rbns.data[,.(k=last(times)),by=.(id)]
      out <- rbns.data[,.(ultimate=k+individual_cost_nf(k,x.vals,yhat)),by=.(id)]
      
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
      
      
      
      
      results = rbind(results,c(actual.tot,pred.tot,cl.tot$ultimate,mean(out_crps$crps_i,
                                                                         na.rm=T)))
      
      trial <- trial+1
      
      cond2 = diff(c(0,yhat))<0
      
      while(sum(cond2)!=0){
        
        yhat <-yhat[!cond2]
        x.vals <-x.vals[!cond2]
        cond2 = diff(c(0,yhat))<0
        
      }
  
    
    
  }
  
  trial <- 1
  ix2 <- 0
  V <- 1200
  
  while( trial < 2){
    
    set.seed(trial+ix2)
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
    
    
    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)
    
    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)
    print(max(yhat1[!cond]))
    if(max(yhat1[!cond])>0.98){
       
      yhat1[cond] <- 1
      
      x.vals1 <- fit1$t
     
      trial <- trial+1
    }
    
    cond2 = diff(c(0,yhat1))<0
    
    while(sum(cond2)!=0){
      
      yhat1 <-yhat1[!cond2]
      x.vals1 <-x.vals1[!cond2]
      cond2 = diff(c(0,yhat1))<0
      
    }
    
    
  }
  
  trial <- 1
  ix2 <- 0
  V <- 12000
  
  while( trial < 2){
    
    set.seed(trial+ix2)
   
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
    
    
    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)
    
    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)
    print(max(yhat2[!cond]))
    
    if(max(yhat2[!cond])>0.94){
      
      yhat2[cond] <- 1
      
      x.vals2 <- fit2$t
      
      
      trial <- trial+1
    }
    
    cond2 = diff(c(0,yhat2))<0
    
    while(sum(cond2)!=0){
      
      yhat2 <-yhat2[!cond2]
      x.vals2 <-x.vals2[!cond2]
      cond2 = diff(c(0,yhat2))<0
      
    }
    
    
  }
  
  
  
  
}


width<-p<-6
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")

ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),
  y=c(yhat,yhat1,yhat2,ytrue2),
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))
)


dt <- data.frame(
  x=c(c(x.vals,last(x.vals)+1),x.vals1,x.vals2),
  y=c(c(yhat,1),yhat1,ytrue2), 
  lbl=c(rep("Fitted 120",length(yhat)+1),
        rep('Fitted 1200',length(yhat1)),
        rep('True curve',length(yhat2)))
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x/unit,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")

ggsave(" severity_curve6.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()



# Width 7 ----


{

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
                  rep(0, width)),
                nrow = width, ncol = width, byrow = TRUE)
    diag(A) <- -rowSums(A)
    A
  }

  set.seed(42)

  trial <- 1
  ix2 <- 0
  V <- 120

  while( trial < 2){

    set.seed(trial+ix2)
    
    sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                      tsh=1e-06,
                                      stdize=1)
   
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

      yhat[cond] <- 1
      x.vals <- fit$t

      trial <- trial+1

      cond2 = diff(c(0,yhat))<0

      while(sum(cond2)!=0){

        yhat <-yhat[!cond2]
        x.vals <-x.vals[!cond2]
        cond2 = diff(c(0,yhat))<0

      }
    


  }

  trial <- 1
  ix2 <- 0
  V <- 1200

  while( trial < 2){

    set.seed(trial+ix2)
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


    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)

    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)

      yhat1[cond] <- 1
      x.vals1 <- fit1$t

      trial <- trial+1

    cond2 = diff(c(0,yhat1))<0

    while(sum(cond2)!=0){

      yhat1 <-yhat1[!cond2]
      x.vals1 <-x.vals1[!cond2]
      cond2 = diff(c(0,yhat1))<0

    }


  }

  trial <- 1
  ix2 <- 0
  V <- 12000

  while( trial < 2){

    set.seed(trial+ix2)
    sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(V,width-1),
                                      tsh=1e-06,
                                      stdize=1)
    ix2=ix2+1
    input.data <- sim.2.fit
    data.list <- input.data$sim.2.fit

    rbns.data <- input.data$censored.data
    rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
    closed.data <- rbns.data %>% group_by(id) %>% filter(as.character(width)%chin%states) %>% as.data.table()
    rbns.data <- rbns.data %>% group_by(id) %>% filter(!(as.character(width)%chin%states)) %>% as.data.table()


    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)

    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)

      yhat2[cond] <- 1
      x.vals2 <- fit2$t
      trial <- trial+1


    cond2 = diff(c(0,yhat2))<0

    while(sum(cond2)!=0){

      yhat2 <-yhat2[!cond2]
      x.vals2 <-x.vals2[!cond2]
      cond2 = diff(c(0,yhat2))<0

    }


  }
  
  }

width<-p<-7
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),
  y=c(yhat,yhat1,yhat2,ytrue2), 
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))
)


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2),
  y=c(yhat,yhat1,ytrue2), 
  lbl=c(rep("Fitted 120",length(yhat)),
        rep('Fitted 1200',length(yhat1)),
       
        rep('True curve',length(yhat2)))
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x/unit,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")

ggsave(" severity_curve7.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()