rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)



source('C:\\Users\\gpitt\\Documents\\GitHub\\conditional-aj-reserving\\helper_functions_ajr.R')


df = read.and.pp.data(fname='C:\\Users\\gpitt\\Documents\\Phd\\visiting\\dati\\data_ku\\final_claim_data.csv')

setDT(df)

dt <- df[accident_period<=5,
   .(final_p=last(development_period)),
   by=.(Claim_number)]

# plot ----

data.frame(cclaim=unname(cumsum(table(dt$final_p)/sum(table(dt$final_p)))),
                    dp=names(table(dt$final_p)))%>%
  mutate(dp=as.numeric(dp)) %>%
  ggplot(aes(x= dp, y= cclaim)) +
  geom_line(linetype="dashed",color="#a71429")+
  geom_point(color="#4169E1")+
  ylim(0.95,1)+
  geom_text(aes(label=round(cclaim,4)),
            vjust = -0.5,color="#454555")+
  theme_bw(base_size=rel(5))+
  theme(plot.title = element_text(size=20),
        legend.text = element_text(size=20))+
  guides(color=guide_legend(title=" "),
         linetype=guide_legend(title=" "))+
  xlab("Development period")+
  ylab(" ")
  
fname=paste0("settlements_",
             format(Sys.time(), 
                    "%Y_%m_%d_%H_%M"),
             ".eps")

ggsave(sprintf("C:\\Users\\gpitt\\Pictures\\markov-chains-reserving\\%s",fname),
       width = 8.5,
       height= 5,
       device=cairo_ps)



