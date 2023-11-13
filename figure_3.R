rm(list=ls())
# libraries ----
library(data.table)
library(xtable)
library(ggplot2)
library(dplyr)
library(tidyr)

## 

locn <- "conditional-aj-reserving\\results_csv\\%s"



fname_no_trend <- c('simulation_no_features_all_records_4_2023_10_23_11_57.csv',
                    'simulation_no_features_all_records_5_2023_10_23_11_58.csv',
                    'simulation_no_features_all_records_6_2023_10_23_12_00.csv',
                    'simulation_no_features_all_records_7_2023_10_23_15_16.csv')

cols <- c("CL" = "#AAAABC",
          'AJ' = "#4169E1")


for(nm in fname_no_trend){
  tmp.k <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- tmp.k
  
  tmp %>% 
    reframe(AJ=pred.tot/actual.tot-1,
            CL=cl.tot/actual.tot-1)  %>%
    pivot_longer(cols=c('AJ',
                        'CL'),
                 names_to='name') %>%
    ggplot( aes(x=name, y=value, fill=name))+
    geom_boxplot(show.legend = FALSE) +
    theme_bw(base_size=rel(5))+
    theme(plot.title = element_text(size=20),
          legend.text = element_text(size=20)) +
    guides(color=guide_legend(title="EI"),
           linetype=guide_legend(title=" "))+
    ylim(-0.1,0.5)+
    xlab("")+
    ylab("")+
    # coord_flip()+
    scale_fill_manual(values = cols)
    
  # ggsave(sprintf("boxplot_nofeatures_ei_%s.eps",tmp.k),
  #        width = 5.5,
  #        height= 5,
  #        device=cairo_ps)
  # 
  
}

# ay trend ----

cols <- c("CL" = "#FF6A7A",
          'AJ' = "#a71429")


fname_trend <- c('simulation_aytrend4_2023_10_23_12_52.csv',
                 'simulation_aytrend5_2023_10_24_17_35.csv',
                 'simulation_aytrend6_2023_10_23_12_54.csv',
                 'simulation_aytrend7_2023_10_23_12_56.csv')


for(nm in fname_trend){
  tmp.k <- as.numeric(gsub(".*?([0-9]+).*", "\\1", nm))
  tmp <- fread(sprintf(locn,nm))
  tmp[['k']] <- tmp.k

  tmp %>% 
    reframe(AJ=pred.tot/actual.tot-1,
            CL=cl.tot/actual.tot-1)  %>%
    pivot_longer(cols=c('AJ',
                        'CL'),
                 names_to='name') %>%
    ggplot( aes(x=name, y=value, fill=name))+
    geom_boxplot(show.legend = FALSE) +
    theme_bw(base_size=rel(5))+
    theme(plot.title = element_text(size=20),
          legend.text = element_text(size=20)) +
    guides(color=guide_legend(title="EI"),
           linetype=guide_legend(title=" "))+
    xlab("")+
    ylab("")+
    ylim(-0.1,0.5)+
    # coord_flip()+
    scale_fill_manual(values = cols)
  
  # ggsave(sprintf("boxplot_features_ei_%s.eps",tmp.k),
  #        width = 5.5,
  #        height= 5,
  #        device=cairo_ps)
  
  
}
