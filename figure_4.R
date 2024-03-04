rm(list=ls())
# libraries ----
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)


unit=1e+06
source('conditional-aj-reserving\\helper_functions_ajr.R')


df = read.and.pp.data(fname='final_claim_data.csv')

setDT(df)


# histogram of ctype

tmp <- df %>%
  group_by(Claim_number) %>%
  summarise(Claim_type_key=unique(Claim_type_key)) 

tmp1 = tmp%>%
  ungroup() %>% 
  group_by(Claim_type_key) %>%
  summarise(counts=n()) 

tmp1 <- tmp1 %>% mutate(counts=counts/dim(tmp)[1])


ggplot(data=tmp1, aes(x=factor(Claim_type_key, 
                               level = as.character(1:19))
                      , y=counts)) +
  geom_bar(stat="identity", fill="#4169E1", color="#454555")+
  theme_bw(base_size=rel(5))+
  ylab("")+
  xlab("claim_type")

dev.off()  



# settlment development period

tmp2 <- df %>%
  group_by(Claim_number) %>%
  filter(delta==1) 

tmp3 = tmp2%>%
  ungroup() %>% 
  group_by(development_period) %>%
  summarise(counts=n()/dim(tmp2)[1]) 


g<- ggplot(data=tmp3, aes(x=development_period
                      , y=counts)) +
  geom_bar(stat="identity", fill="#4169E1", color="#454555")+
  theme_bw(base_size=rel(5))+
  xlim(0.5,8)+
  ylab("")+
  xlab("Settlement DP")
# ggsave(g, file="newsm.eps", device="eps",width = 6,height = 5.5)
# dev.off()  


##

tmp4 <- df %>%
  group_by(Claim_number) %>%
  summarise(csize=sum(incPaid)) 

tmp4 %>%
  ggplot(aes(x=csize))+
  geom_density(fill="#4169E1", color="#454555")+
  theme_bw(base_size=rel(5))+
  scale_x_continuous(trans = log2_trans(),
                                  breaks = trans_breaks("log2", function(x) 2^x),
                                  labels = trans_format("log2", math_format(2^.x)))+
  ylab("")+
  xlab("Y")


dev.off()  

#

tmp5 <- df %>%
  group_by(Claim_number) %>%
  summarise(npays=n()) %>%
  as.data.frame()

tmp6 <- table(tmp5$npays)

tmp6 %>%
  as.data.frame() %>%
  ggplot(aes(x=Var1,y=Freq/dim(tmp5)[1]))+
  geom_bar(fill="#4169E1", color="#454555",stat="identity")+
  theme_bw(base_size=rel(5))+
  ylab("")+
  xlab("Number of payments")


dev.off()  







