
jump_rate <- function(i, t, u){
  
  jrate[i]
}

mark_dist <- function(i, s, v){
  return(imx[i,])
}

cev <- function(k,x,y){
  "
  It computes E[X|X>k]:conditional on k the expected value of x given the cumulative distrubtion function y.
  "
  # y=y/sum(y)
  
  # if(k<1e-03){return(sum(diff(c(0,x))*(1-y)))}
  if(k<1e-03){return(sum((x-k)*diff(c(0,y))))}
  
  temp.cond <- (x <= k)
  temp.which <- y[temp.cond]
  den <- (1-tail(temp.which,1))
  
  # num <- sum((x[!temp.cond]-k)*diff(c(tail(temp.which,1),
  #                                     y[!temp.cond])))
  
  num <- sum(diff(c(k,x[!temp.cond])) *(1-y[!temp.cond]))
  
  return(num/den)
  
}

cev2 <- function(k,x,y){
  "
  It computes E[X^2|X>k]:conditional on k the expected value of x given the cumulative distrubtion function y.
  "
  # y=y/sum(y)
  
  # if(k<1e-03){return(sum(diff(c(0,x))*(1-y)))}
  if(k<1e-03){return(sum(x*(x-k)*diff(c(0,y))))}
  
  temp.cond <- (x <= k)
  temp.which <- y[temp.cond]
  den <- (1-tail(temp.which,1))
  
  # num <- sum((x[!temp.cond]-k)*diff(c(tail(temp.which,1),
  #                                     y[!temp.cond])))
  
  num <- sum(x[!temp.cond]*diff(c(k,x[!temp.cond])) *(1-y[!temp.cond]))
  
  return(num/den)
  
}


sample.and.censor.ap.strategy <- function(x,
                                          p,
                                          c.val,
                                          tsh=1e-04,
                                          stdize = 1){
  "
  This function takes a state and censors observations. 
  It also adds a small time windown between states until the censoring "
  # browser()
  ap <- x[[3]][1]
  c.val <- p-ap
  
  l <- length(x[[1]])
  times <- x[[1]]/stdize
  states <- x[[2]]
  id <- x[[4]][1]
  
  if(c.val > states[l-1]){ 
    tmp <- x
    tmp[[1]]=times
    tmp[[3]] = rep(ap,l)
    tmp[[4]] = rep(id,l)
    
    names(tmp) <- c("times",
                    "states",
                    "accident_period",
                    "id")
    
    return(tmp)}
  
  exact.cens <- c.val %in% states
  
  if(exact.cens){
    
    cond <- states <= c.val
    v1 <- c(states[cond],last(states[cond]))
    end_time <- times[1:(sum(cond))]
    v2 <- c(end_time,last(end_time)+tsh)
    
    l2 <- length(v2)
    
    v3 = rep(ap,l2)
    v4 = rep(id,l2)
    
    return(list(times=v2,
                states=v1,
                accident_period=v3,
                id=v4
    ))
    
  }else{
    
    
    v1 <- 1:c.val
    cond <- states <= c.val
    if(sum(cond) > 1){
      
      add <- times[1:(sum(cond))]
      
    }else{add <- 0}
    
    # add= switch((sum(cond)==1)+1,times[1:(sum(cond)-1)],NULL) 
    
    v1 <- c(states[cond],last(states[cond]+1):c.val)
    l <- length(v1)+1
    plus.add <- last(add)+cumsum(rep(tsh,l-length(add)))
    
    v2 <- c(add,plus.add)
    
    v3 <- rep(ap,l)
    v4 <- rep(id,l)
    
    return(list(times=v2,
                states=c(v1,c.val),
                accident_period= v3,
                id=v4))
    
  }}

renamer.ap <- function(x,y){
  l <- length(x[[1]])
  x$accident_period <- rep(y,l)
  # past <- names(tmp)
  # names(tmp)<-c("times","states","")
  return(x)
}

renamer.id <- function(x,y){
  l <- length(x[[1]])
  x$id <- rep(y,l)
  # past <- names(tmp)
  # names(tmp)<-c("times","states","")
  return(x)
}

renamer.X <- function(x,y){
  tmp <- append(x,y)
  names(tmp)<-c("times","states","X")
  return(tmp)}


# crps_computer <- function(k,y,x,cdf_i){

  # s_i = 1 - cdf_i
  # # browser()
  # temp.cond <- (x <= k)
  # temp.which <- cdf_i[temp.cond]
  # den <- (1-tail(temp.which,1))
  # # browser()
  # tmp2 <- !temp.cond & (x <= y)
  # tmp3 <- !temp.cond & (x > y)
  # 
  # if(k<1e-03){return(sum(diff(c(0,x))*(s_i)^2))}
  # 
  # num <- sum(diff(c(k,x[tmp2])) *(cdf_i[tmp2])^2) + sum(diff(c(last(x[tmp2]),x[tmp3])) *(s_i[tmp3])^2)
  # 
  # if(is.na(num/den)){
  #   
  #   return(0)
  #   
  #   }else{
  #     return(num/den)
  # }}

crps_computer <- function(k,y,x,cdf_i){
  
  s_i = 1 - cdf_i
  
  x.vals <- c(0,diff(x))

  crps=sum((x.vals*(cdf_i^2))[x<=y])+sum((x.vals*(s_i^2))[x>y])
  
  return(crps)  
  
}

simulate_mkvian_data <- function(rates = jump_rate, 
                                 dists = mark_dist,
                                 ap.volumes,
                                 tsh=1e-04,
                                 stdize=1,
                                 X=NULL){
  
  l <- length(ap.volumes)
  k <- l+1
  
  sim <- list()
  
  tmp.tot <- sum(ap.volumes)
  
  for(i in 1:tmp.tot){
    sim[[i]] <- AalenJohansen:::sim_path(1, 
                                         rates = jump_rate, 
                                         dists = mark_dist,
                                         t=0,
                                         tn=1e+15,
                                         stdize)
  }
  
  ap <- rep(0:(l-1),ap.volumes)
  id <- 1:tmp.tot
  sim <- Map(renamer.ap, sim, ap)
  sim <- Map(renamer.id, sim, id)
  
  actual.data <-rbindlist(sim, fill=TRUE)
  
  # sim <- Map(append, sim, id)
  
  sim.cens <- lapply(sim, 
                     sample.and.censor.ap.strategy,
                     p=k,
                     stdize=stdize,
                     tsh=tsh)
  
  censored.data <- rbindlist(sim.cens, fill=TRUE)
  
  test <- censored.data[,.(any(duplicated(states))),
                        by=.(id)] 
  
  sim.2.fit <- lapply(sim.cens,function(x) x[1:2])
  
  if(!is.null(X)){
    sim.2.fit <- Map(renamer.X, sim.2.fit, X)
  }
  
  return(list(sim.2.fit = sim.2.fit,
              actual.data=actual.data,
              censored.data=censored.data))
  
}

simulate_mkvian_dataX <- function(rates = jump_rate, 
                                  dists = mark_dist,
                                  ap.volumes,
                                  tsh=1e-04,
                                  stdize=1,
                                  X=X){
  
  l <- length(ap.volumes)
  k <- l+1
  
  sim <- list()
  
  tmp.tot <- sum(ap.volumes)
  
  for(i in 1:tmp.tot){
    jrates<- function(i, t, u){return(rates(i, t, u)/X[i])}
    sim[[i]] <- AalenJohansen:::sim_path(1, 
                                         rates = jrates, 
                                         dists = dists,
                                         t=0,
                                         tn=1e+15,
                                         stdize)
  }
  
  ap <- rep(0:(l-1),ap.volumes)
  id <- 1:tmp.tot
  sim <- Map(renamer.ap, sim, ap)
  sim <- Map(renamer.id, sim, id)
  
  actual.data <-rbindlist(sim, fill=TRUE)
  
  # sim <- Map(append, sim, id)
  
  sim.cens <- lapply(sim, 
                     sample.and.censor.ap.strategy,
                     p=k,
                     stdize=stdize,
                     tsh=tsh)
  
  censored.data <- rbindlist(sim.cens, fill=TRUE)
  
  test <- censored.data[,.(any(duplicated(states))),
                        by=.(id)] 
  
  sim.2.fit <- lapply(sim.cens,function(x) x[1:2])
  
  if(!is.null(X)){
    sim.2.fit <- Map(renamer.X, sim.2.fit, X)
  }
  
  return(list(sim.2.fit = sim.2.fit,
              actual.data=actual.data,
              censored.data=censored.data))
  
}


individual_cost_nf <- function(k,x.vals,yhat){
  
  a <- cev(k=k,x=x.vals,y=yhat)
  
  if(is.na(a)){
    return(0)
    }else{
    return(a)
    }
  
  }


individual_vty_nf <- function(k,x.vals,yhat){
  
  a <- cev(k=k,x=x.vals,y=yhat)
  
  if(is.na(a)){
    return(0)
  }else{
    return(2*cev2(k,x.vals,yhat)-cev(k,x.vals,yhat)^2)
  }
  
}



individual_resultsX <- function(k,X,data.list){
  
  fit <- AalenJohansen::aalen_johansen(data.list,x=X)
  yhat <-  sapply(fit$p, last)
  x.vals <- fit$t
  
  cond <- (yhat >1 | yhat <0) | is.na(yhat)
  
  yhat[cond]<-1
  x.vals <- fit$t
  # x.vals<- x.vals[cond]
  
  cond2 = diff(c(0,yhat))<0
  
  while(sum(cond2)!=0){
    
    yhat <-yhat[!cond2]
    x.vals <-x.vals[!cond2]
    cond2 = diff(c(0,yhat))<0
    
  }
  
  ev <- cev(k,x=x.vals,y=yhat)
  if(is.na(ev)){ev <- 0}
  mment2 <- cev2(k,x=x.vals,y=yhat)
  vty <- 2*mment2-ev^2
  if(is.na(mment2)){vty <- 0}
  
  crps_i=crps_computer(k=k,y=k+ev,x=x.vals,cdf_i=yhat)
  
  return(list(ultimate = ev+k,
              variance=vty,
              crps_i=crps_i))
  
  }


compute.ajr.models <- function(dt,
                               true.ultimate,
                               X='Claim_type_key',
                               data.list){
  
  tmp.x <- unique(sort(as.numeric(as.character(dt[[X]]))))
  l <- list()
  
  for(i in tmp.x){
    
    fit <- AalenJohansen::aalen_johansen(data.list,x=i)
    yhat <-  sapply(fit$p, last)
    x.vals <- fit$t
    cond <- (yhat >1 | yhat <0) | is.na(yhat)
    
    yhat[cond]<-1
    x.vals <- fit$t
    # x.vals<- x.vals[cond]
    
    cond2 = diff(c(0,yhat))<0
    
    while(sum(cond2)!=0){
      
      yhat <-yhat[!cond2]
      x.vals <-x.vals[!cond2]
      cond2 = diff(c(0,yhat))<0
      
    }
    
    l[[as.character(i)]] <- list(yhat=yhat,
                                 x.vals=x.vals)
    
  }
  
  return(l)
  
}


individual_results_nf <- function(k,true.ultimate,x.vals,yhat){

  ev <- cev(k,x=x.vals,y=yhat)
  if(is.na(ev)){ev <- 0}
  mment2 <- 2*cev2(k,x=x.vals,y=yhat)
  vty <- mment2-ev^2
  if(is.na(mment2)){vty <- 0}
  
  crps_i=crps_computer(k=k,y=true.ultimate,x=x.vals,cdf_i=yhat)
  
  return(list(ultimate = ev+k,
              variance=vty,
              crps_i=crps_i))
  
}

individual_results_X_w_pre_saved_models <- function(k,true.ultimate,X,models.list){
  xi <- as.character(X)
  
  yhat <- models.list[[xi]]$yhat
  x.vals <- models.list[[xi]]$x.vals
  
  ev <- cev(k,x=x.vals,y=yhat)
  if(is.na(ev)){ev <- 0}
  mment2 <- 2*cev2(k,x=x.vals,y=yhat)
  vty <- mment2-ev^2
  if(is.na(mment2)){vty <- 0}
  
  crps_i=crps_computer(k=k,y=true.ultimate,x=x.vals,cdf_i=yhat)
  
  return(list(ultimate = ev+k,
              variance=vty,
              crps_i=crps_i))
  
}

recode.cl <- function(x){
  if(length(x) == 1){
    # browser()
    return(1)
  }else{
    x[length(x)] <- x[length(x)-1]+1
    return(x)
  }
  
}


cl.calculator<- function(dt){
  
  f <- c()
  sigmas <- c()
  J<-dim(dt)[2]
  
  for(j in 1:(J-1)){
    
    f <- c(f,sum(dt[1:(J-j),j+1])/sum(dt[1:(J-j),j]))
    sigmas <- c(sigmas, sum(dt[1:(J-j),j]*(dt[1:(J-j),j+1]/dt[1:(J-j),j] - last(f))^2)/(J-j-1))
    
    
  }
  
  
  
  dg <- rev(as.numeric(ChainLadder::getLatestCumulative(dt)))
  
  
  tmp.final = (length(dg)-1)
  dt.hat <- dt
  
  for(i in 1:tmp.final){
    
    f.product <- cumprod(f[i:tmp.final])
    fcsts <- dg[i]*f.product
    dt.hat[dim(dt)[1]+1-i,(i+1):dim(dt)[2]] <- fcsts
  }
  
  ultimate <- dt.hat[,dim(dt)[2]]
  
  # variability 
  
  sigmas[length(sigmas)] <- min((sigmas[length(sigmas)-1]^2)/sigmas[length(sigmas)-2],
                                min(sigmas[length(sigmas)-1],sigmas[length(sigmas)-2]))
  
  J <- dim(dt)[2]
  msep.Ri <- c()
  pv.Ri <- c()
  for(i in 2:J){
    
    external.cms <- last(dt.hat[i,])^2
    tmp.sigmas<-sigmas[(J-i+1):(J-1)]
    tmp.fs<-(f[(J-i+1):(J-1)])^2
    
    n1 <- 1/dt.hat[i,(J-i+1):(J-1)]
    n2 <- c()
    for(i2 in (J-i+1):(J-1)){
      n2 <- c(n2,1/sum(dt[1:(J-i2),i2]))
    }
    
    tmp.ri <- external.cms*sum(tmp.sigmas/tmp.fs*(n1+n2))
    
    msep.Ri <- c(msep.Ri,tmp.ri)
    
    pv.Ri <- c(pv.Ri,external.cms*sum(tmp.sigmas/tmp.fs*(n1)))
    
  }
  
  
  
  ## correlations
  
  cors.Ri <- c()
  
  for(i in 2:(J-1)){
    
    last.cms <- last(dt.hat[i,])
    future.cms <- dt.hat[(i+1):J,J]
    
    tmp.sigmas<-sigmas[(J-i+1):(J-1)]
    tmp.fs<-(f[(J-i+1):(J-1)])^2
    
    n2 <- c()
    for(i2 in (J-i+1):(J-1)){
      n2 <- c(n2,1/sum(dt[1:(J-i2),i2]))
    }
    
    tmp.ri2 <- sum(2*tmp.sigmas/tmp.fs*n2)
    
    
    cors.Ri <- c(cors.Ri,last.cms*sum(future.cms)*tmp.ri2)
    
  }
  
  l <- list(
    ultimate = dt.hat[,J],
    process.error.i=sqrt(sum(pv.Ri)),
    msep=sqrt(sum(c(cors.Ri[1:(J-2)],0)+msep.Ri)))
  
  
  
}

chain_ladder_computer <- function(data.set,width){
  
  # cl.data <- data.set[,.(incrementals=diff(c(0,times))),by=.(id,accident_period,states)]
  
  cl.data <- data.set %>%
    group_by(id,accident_period,states) %>%
    filter(states!=1) %>%
    reframe(incrementals=diff(c(0,times))) %>%
    group_by(id) %>% mutate(development_period=recode.cl(states))
  
  
  
  # cl.data<- cl.data[states!=1]
  # cl.data[["states"]]=cl.data[["states"]]-1
  cl.data<-filter(cl.data,(development_period+accident_period)<=(width-1))
  
  
  
  cl.tr <- ChainLadder::as.triangle(cl.data,
                                    origin="accident_period",
                                    dev="development_period",
                                    value="incrementals")
  cl.tr <- ChainLadder::incr2cum(cl.tr)
  
  
  # if(width < 5){
  #   m <- MackChainLadder(cl.tr,est.sigma = 1,mse.method = "Mack")
  #   ultimate.se=0
  # }else{
  #   m <- MackChainLadder(cl.tr,est.sigma = "Mack",mse.method = "Mack")
  #   ultimate.se=summary(m)$Totals['Mack S.E.',]
  # }
  
  cl.output <- cl.calculator(cl.tr)
  
  
  return(list(
    ultimate=sum(cl.output$ultimate),
    ultimate.se=cl.output$msep,
    process.se=cl.output$process.error.i))
  
}

understand.the.densities <- function(data){
  
  ev.teorico <- data$true.ultimate[1]
  
  tmp <- data %>%
    tidyr::pivot_longer(cols=c("pred.tot","cl.tot"),
                        names_to="model",
                        values_to="value")
  
  # Basic density plot with custom color
  ggplot(tmp, aes(x=value, color=model)) + 
    # facet_wrap(~ p)+
    geom_density(size=2)+
    geom_vline(xintercept = ev.teorico, linetype="dotted", 
               color = "black",size=2)+theme_bw()
  
  
  
}


reformulate_state_space <- function(x, newk=4){
  
  tmp<- as.data.table(x)
  setDT(tmp)[states>=newk,states:=newk]
  
  tmp <- tmp[,.(times=sum(times)),by=states]
  
  return(as.list(tmp[,c("times","states")]))
  
  
}

reformulate_state_space_X <- function(x, feature.name="X",newk=4){
  
  tmp<- as.data.table(x)
  setDT(tmp)[states>=newk,states:=newk]
  
  tmp <- tmp[,.(times=sum(times)),by=states]
  
  tmp <- as.list(tmp[,c("times","states")])
  
  tmp[[feature.name]]=x[[feature.name]]
  
  return(tmp)
  
  
}

# Real data ----


read.and.pp.data <- function(fname){
  
  
  df <- read.csv(fname,
                 header = T,
                 colClasses = c(rep('character',5),
                                'factor',
                                'factor',
                                'numeric'))
  
  # df <- df[(df$AM)>=20180101,]
  df <- df[!is.na(df$incPaid),]
  df <- df[df$incPaid>0,]
  df <- df[df$incPaid<1e+06,]
  
  
  df$development_period = floor((as.integer(df$DM)-1)/12)+1
  # df$development_period= df$development_year+1
  baseAM <- as.numeric(substr(df[1,'AM'],1,4))
  acc.st = df$AM
  year.a=as.integer(substr(acc.st ,1,4))-baseAM
  df$accident_period=year.a
  
  rep.st = df$RM
  year.r=as.integer(substr(rep.st ,1,4))-baseAM+1
  df$reporting_period=year.r
  
  df = df %>%
    group_by(Claim_number, Claim_type_key, development_period,accident_period) %>%
    summarise(incPaid=sum(incPaid),
              reporting_period=min(reporting_period)
    ) %>%
    mutate(calendar_period=development_period+accident_period)%>%
    as.data.frame()
  
  df = df %>%
    group_by(Claim_number) %>%
    mutate(delta = c(rep(0,n()-1),1)) %>%
    as.data.frame()
  
  return(df)
  
}

clean.lt <- function(x){
  J=dim(x)[2]
  for(j in 1:J){
    for(i in 1:J){
      
      if((i+j) > (J+1)){x[i,j]=NA}
      
    }}
  return(x)}


find.t.data <- function(x){
  
  setDT(x)
  dt.tr <- x[,.(pays=sum(incPaid)),by=.(development_period,
                    accident_period)]
  
  out <- as.triangle(dt.tr,
              origin = 'accident_period',
              dev='development_period',
              value = 'pays')
  
  return(out)
}

encode.cc <- function(x,maximum.p){
  x=x+1
  x[length(x)] <- maximum.p
  return(x)}

encode.rbns <- function(x,y,ay,maximum.p,tsh=1e-08){
  x=x+1
  lx <- last(x)
  if(lx==(maximum.p-ay)){
    x<-c(x,maximum.p-ay)
  }else{
    x <- c(x, seq(lx+1,maximum.p-ay),maximum.p-ay)
  }
  return(list(events=x,
              times=c(y,rep(tsh,length(x)-length(y)))))
}

input_4_cl_case <- function(maximum.p){
  
  df = read.and.pp.data(fname='C:\\Users\\gpitt\\Documents\\Phd\\visiting\\dati\\data_ku\\final_claim_data.csv')
  
  df = df %>% mutate(incPaid=incPaid)
  
  df =df %>%
    group_by(Claim_number) %>%
    filter(!any(development_period > (maximum.p-1)))
  
  df=df  %>%
    group_by(Claim_number) %>%
    filter(1 %in%delta &(accident_period <= (maximum.p-1-1)))
  
  return(df)
  
}





