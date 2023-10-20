rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ChainLadder)
library(ggplot2)
library(scales)

source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')

# Construct the model ----


for(maximum.p in c(4,5,6,7,8)){
  
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

  
  rbns.data <- rbns.claims.fit[,times:=cumsum(times),by=Claim_number][,.(k=last(times)),by=.(Claim_number)]
  
  
  rbns.data<-features.data[rbns.data, 
                           on = 'Claim_number']
  
  models.list <- compute.ajr.models(rbns.data,
                                    data.list=data.list.ref)

dt <- data.table()

for(nms in names(models.list)){
  
  y <- models.list[[nms]]$yhat
  x<- models.list[[nms]]$x.vals
  
  dt <- rbind(dt,
              data.table(x=x,
                         y=y,
                         claim_type=nms))
  
}

cvi_colours = list(
  cvi_royal = c("#a71429", "#4169E1", "#AAAABC", "#454555",
                  "#FF6A7A"),
  my_favourite_colours = c("#a71429", "#4169E1",    "#AAAABC")
)


cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,
               continuous = grDevices::colorRampPalette(palette)(n),
               discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

cvi_palettes("my_favourite_colours", type = "discrete")


scale_colour_cvi_d = function(name) {
  ggplot2::scale_colour_gradientn(colours = cvi_palettes(name = name,
                                                         type = "discrete"))
}

library(scales)
library(RColorBrewer)
dt %>%
  ggplot(aes(x=x,y=y,color=claim_type)) +
  geom_line()+
  # scale_colour_cvi_d("cvi_royal")+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  # guides(color=guide_legend(title=" "),
         # linetype=guide_legend(title=" "))+
  # scale_color_discrete(breaks=as.character(sort(as.numeric(names(models.list)))))+
  scale_color_brewer(palette="Spectral",
                     breaks=as.character(sort(as.numeric(names(models.list)))))+
  xlab("Z")+
  ylab("cdf")+
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

fname=sprintf("severity_curve_maxp_%s.eps",maximum.p)
ggsave(sprintf("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\%s",fname),
       width = 6.5,
       height= 5,
       device=cairo_ps)
# dev.off()

}



