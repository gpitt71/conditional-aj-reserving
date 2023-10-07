# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)
library(scales)
source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

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

# data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
# data.setC[,"X"] <- X

  fit <- AalenJohansen::aalen_johansen(data.list)
  yhat <-  sapply(fit$p, last)
# tail(yhat)
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

# data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
# data.setC[,"X"] <- X

fit1 <- AalenJohansen::aalen_johansen(data.list)
yhat1 <-  sapply(fit1$p, last)
# tail(yhat)
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

# data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
# data.setC[,"X"] <- X

fit2 <- AalenJohansen::aalen_johansen(data.list)
yhat2 <-  sapply(fit2$p, last)
# tail(yhat)
x.vals2 <- fit2$t


# sim.2.fit <- simulate_mkvian_data(ap.volumes=rep(120000,width-1),
#                                   tsh=1e-06,
#                                   stdize=1)
#
# input.data <- sim.2.fit
# data.list <- input.data$sim.2.fit
#
# rbns.data <- input.data$censored.data
# rbns.data[["states"]] <- as.vector(as.character(rbns.data$states))
# closed.data <- rbns.data %>% group_by(id) %>% filter("4"%chin%states) %>% as.data.table()
# rbns.data <- rbns.data %>% group_by(id) %>% filter(!("4"%chin%states)) %>% as.data.table()
#
# # data.setC <- input.data$censored.data[,.(k=last(times)),by=.(id)]
# # data.setC[,"X"] <- X
#
# fit3 <- AalenJohansen::aalen_johansen(data.list)
# yhat3 <-  sapply(fit$p, last)
# # tail(yhat)
# x.vals3 <- fit$t

}


width<-p<-4
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
# x <- ajr.fit$dy.models[[p-1]]$x
# yhat <- ajr.fit$dy.models[[p-1]]$yhat
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)
# ytrue3 <- True_p(x.vals3)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2), #x.vals3
  y=c(yhat,yhat1,yhat2,ytrue2), #yhat3
  lbl=c(rep('Fitted 120',length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        rep('Fitted 12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
  )


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2), #x.vals3
  y=c(yhat,yhat1,ytrue2), #yhat3
  lbl=c(rep('Fitted 120',length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        # rep('Fitted 12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)

cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          # "Fitted 12000" = "#AAAABC",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")+
  scale_x_continuous(trans = log2_trans(),
                                  breaks = trans_breaks("log2", function(x) 2^x),
                                  labels = trans_format("log2", math_format(2^.x)))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggsave("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\severity_curve4.eps",
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
    if(max(yhat[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat[cond] <- 1
      # tail(yhat)
      x.vals <- fit$t
      # x.vals <-x.vals[cond]

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

      # 1-pred.tot/actual.tot
      # 1-cl.tot$ultimate/actual.tot

      # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}


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


  }

  trial <- 1
  ix2 <- 0
  V <- 1200

  while( trial < 2){

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


    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)

    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)
    print(max(yhat1[!cond]))
    if(max(yhat1[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat1[cond] <- 1
      # tail(yhat)
      x.vals1 <- fit1$t
      # x.vals <-x.vals[cond]

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


    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)

    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)
    print(max(yhat2[!cond]))

    if(max(yhat2[!cond])>0.94){
      # print(max(yhat[!cond]))
      yhat2[cond] <- 1
      # tail(yhat)
      x.vals2 <- fit2$t
      # x.vals <-x.vals[cond]

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



width<-p<-5
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
# x <- ajr.fit$dy.models[[p-1]]$x
# yhat <- ajr.fit$dy.models[[p-1]]$yhat
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),#x.vals3
  y=c(yhat,yhat1,yhat2,ytrue2), #yhat2,yhat3
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2),#x.vals3
  y=c(yhat,yhat1,ytrue2), #yhat2,yhat3
  lbl=c(rep("Fitted 120",length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        # rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          # "Fitted 12000" = "#AAAABC",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggsave("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\severity_curve5.eps",
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
                  # 1*mark_dist(7, t, 0),
                  # 1*mark_dist(8, t, 0),
                  # 1*mark_dist(9, t, 0),
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
    if(max(yhat[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat[cond] <- 1
      # tail(yhat)
      x.vals <- fit$t
      # x.vals <-x.vals[cond]

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

      # 1-pred.tot/actual.tot
      # 1-cl.tot$ultimate/actual.tot

      # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}


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


  }

  trial <- 1
  ix2 <- 0
  V <- 1200

  while( trial < 2){

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


    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)

    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)
    print(max(yhat1[!cond]))
    if(max(yhat1[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat1[cond] <- 1
      # tail(yhat)
      x.vals1 <- fit1$t
      # x.vals <-x.vals[cond]

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


    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)

    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)
    print(max(yhat2[!cond]))

    if(max(yhat2[!cond])>0.94){
      # print(max(yhat[!cond]))
      yhat2[cond] <- 1
      # tail(yhat)
      x.vals2 <- fit2$t
      # x.vals <-x.vals[cond]

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

width<-p<-7
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
# x <- ajr.fit$dy.models[[p-1]]$x
# yhat <- ajr.fit$dy.models[[p-1]]$yhat
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),#x.vals3
  y=c(yhat,yhat1,yhat2,ytrue2), #yhat2,yhat3
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2),#x.vals3
  y=c(yhat,yhat1,ytrue2), #yhat2,yhat3
  lbl=c(rep("Fitted 120",length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        # rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          # "Fitted 12000" = "#AAAABC",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggsave("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\severity_curve7.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()




# Width 8 ----



{


  input.pbties <-c(-4.02104115693941e-06, 9.41733693345327e-07, 1.45337940215019e-08, 1.0284444599133e-09, 1.28667155196585e-10, 1.41137223345763e-11, 0, 3.06360244423514e-06,
                   0, -6.79581039063144e-06, 1.38541363530372e-06, 3.02422921453848e-08, 4.66191880121944e-09, 0, 3.70408553883232e-09, 5.37178845884229e-06,
                   0, 0, -5.97028986443137e-06, 2.59887586380591e-06, 2.77661275734168e-08, 2.3138439644514e-08, 1.19977094453035e-08, 3.30851172396223e-06,
                   0, 0, 0, -5.23649767140797e-06, 3.19644689000269e-06, 6.99411925618263e-08, 0, 1.97010958884345e-06,
                   0, 0, 0, 0, -3.16996623129841e-06, 1.64437177740346e-06, 0, 1.52559445389495e-06,
                   0, 0, 0, 0, 0, -2.10559800765077e-06, 1.13378354258118e-06, 9.71814465069587e-07,
                   0, 0, 0, 0, 0, 0, -3.23938155023196e-07, 3.23938155023196e-07,
                   0, 0, 0, 0, 0, 0, 0, 0)

  width <- 8

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
                  1*mark_dist(7, t, 0),
                  # 1*mark_dist(8, t, 0),
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
    if(max(yhat[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat[cond] <- 1
      # tail(yhat)
      x.vals <- fit$t
      # x.vals <-x.vals[cond]

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

      # 1-pred.tot/actual.tot
      # 1-cl.tot$ultimate/actual.tot

      # if(is.na(mean(out_crps$crps_i))){print(sum(is.na(out_crps$crps_i)))}


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


  }

  trial <- 1
  ix2 <- 0
  V <- 1200

  while( trial < 2){

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


    fit1 <- AalenJohansen::aalen_johansen(data.list)
    yhat1 <-  sapply(fit1$p, last)

    cond <- (yhat1 > 1 | yhat1 <0) | is.na(yhat1)
    print(max(yhat1[!cond]))
    if(max(yhat1[!cond])>0.98){
      # print(max(yhat[!cond]))
      yhat1[cond] <- 1
      # tail(yhat)
      x.vals1 <- fit1$t
      # x.vals <-x.vals[cond]

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


    fit2 <- AalenJohansen::aalen_johansen(data.list)
    yhat2 <-  sapply(fit2$p, last)

    cond <- (yhat2 > 1 | yhat2 <0) | is.na(yhat2)
    print(max(yhat2[!cond]))

    if(max(yhat2[!cond])>0.94){
      # print(max(yhat[!cond]))
      yhat2[cond] <- 1
      # tail(yhat)
      x.vals2 <- fit2$t
      # x.vals <-x.vals[cond]

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


width<-p<-8
L <- matrix(input.pbties, byrow = T,ncol=width)
diag(imx) <- 0
alpha <- c(1,rep(0,width-2))
Lt <- function(t) L
True_p <- Vectorize(function(t) (c(alpha,0)%*%expm::expm(Lt(t)*t))[p],
                    vectorize.args = "t")
# x <- ajr.fit$dy.models[[p-1]]$x
# yhat <- ajr.fit$dy.models[[p-1]]$yhat
ytrue <- True_p(x.vals)
ytrue1 <- True_p(x.vals1)
ytrue2 <- True_p(x.vals2)

dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2,x.vals2),#x.vals3
  y=c(yhat,yhat1,yhat2,ytrue2), #yhat2,yhat3
  lbl=c(rep('fitted120',length(yhat)),
        rep('fitted1200',length(yhat1)),
        rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


dt <- data.frame(
  x=c(x.vals,x.vals1,x.vals2),#x.vals3
  y=c(yhat,yhat1,ytrue2), #yhat2,yhat3
  lbl=c(rep("Fitted 120",length(yhat)),
        rep('Fitted 1200',length(yhat1)),
        # rep('fitted12000',length(yhat2)),
        rep('True curve',length(yhat2)))#yhat3
)


cols <- c("Fitted 120" = "#a71429",
          'Fitted 1200' = "#4169E1",
          # "Fitted 12000" = "#AAAABC",
          "True curve" = "#454555")

dt %>%
  ggplot(aes(x=x,y=y,color=lbl)) +
  geom_line(aes(linetype=lbl), show.legend = FALSE)+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  scale_color_manual(values = cols)+
  xlab("Z")+
  ylab("cdf")+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

ggsave("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\severity_curve8.eps",
       width = 5.5,
       height= 5,
       device=cairo_ps)
dev.off()

