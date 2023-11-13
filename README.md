# Individual claims reserving using the Aalen-Johansen estimator

This repository contains the code to replicate the manuscript INDIVIDUAL CLAIMS RESERVING USING THE AALEN-JOHANSEN ESTIMATOR. 

The replication material is organized as described in this file and it can be used for replicating exactly the results of our manuscript.

The script `elaborate_latex_tables.R` contains the code to print the tables in the manuscript from the the results in the `results_csv` folder.

The script `helper_functions_ajr.R` contains functions that we used to perform the calculations.  

The real data that we used in our application will not be shared, but we provide the results that we obtained.

For each section, we describe below the relevant replication material. 

# Section 4. Simulation of RBNS claims

## Section 4.1. Implementation

- Figure 2 can be reproduced with the script `figure_2.R`.

## Section 4.2. Evaluating the models performance

- Figure 3 can be reproduced with the script `figure_3.R`.

## Section 4.3. Empirical analysis

- Table 1 can be reproduced with the script  `table_1.R`.

- Table 2 can be reproduced with the scripts `table_2_intercept_model.R` (results for the intercept model) and `table_2_feature_model.R` (results for the model with feature).

# Section 5. A data application on an insurance portfolio

- We obtained Figure 4 with the script  `figure_4.R`. 

- We obtained Figure 5 with the script  `figure_5.R`. 


## Section 5.1. Model comparison on different datasets

- We obtained Table 4 and Table 5 using the script `table_4_intercept_model.R`(results for the intercept model) and `table_4_feature_model.R` (results for the model with feature).

## Section 5.2. Model comparison on a single dataset

- We obtained Table 6 using the scripts `table_6_intercept_model.R` (results for the intercept model) and `table_6_feature_model.R` (results for the model with feature).

# Content of the results_csv folder

- Results of scenario Alpha: 
  
  - `simulation_no_features_all_records_4_2023_10_23_11_57.csv`
  - `simulation_no_features_all_records_5_2023_10_23_11_58.csv`
  - `simulation_no_features_all_records_6_2023_10_23_12_00.csv`
  - `simulation_no_features_all_records_7_2023_10_23_15_16.csv`
  

- Results of scenario Beta:
  
  - Intercept model
  
    - `simulation_aytrend_interceptmodel_all_records_4_2023_10_24_16_27.csv`
    - `simulation_aytrend_interceptmodel_all_records_5_2023_10_24_17_36.csv`
    - `simulation_aytrend_interceptmodel_all_records_6_2023_10_24_16_31.csv`
    - `simulation_aytrend_interceptmodel_all_records_7_2023_10_24_16_32.csv`
  
  - Model with feature
  
    - `simulation_aytrend4_2023_10_23_12_52.csv`
    - `simulation_aytrend5_2023_10_24_17_35.csv`
    - `simulation_aytrend6_2023_10_23_12_54.csv`
    - `simulation_aytrend7_2023_10_23_12_56.csv`

- Results on the real data (Section 5.1.): 
  
  - `real_data_w_features__2023_11_02_10_42.csv` 
  - `real_data_no_features__2023_11_02_09_58.csv`

- Results on the real data (Section 5.2.): 
  
  - `real_data_w_features_maxp_6_2023_11_02_11_27.csv`
  - `real_data_no_features_maxp_6_2023_11_02_09_29.csv`
  
# Current R Session

Using the command `sessionInfo()`, we report version information about R, the OS and attached or loaded packages.

```
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Italian_Italy.utf8  LC_CTYPE=Italian_Italy.utf8    LC_MONETARY=Italian_Italy.utf8
[4] LC_NUMERIC=C                   LC_TIME=Italian_Italy.utf8    

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] scales_1.2.1       tidyr_1.3.0        xtable_1.8-4       ggplot2_3.4.3      ChainLadder_0.2.18
[6] dplyr_1.1.2        data.table_1.14.8  AalenJohansen_1.0 

loaded via a namespace (and not attached):
 [1] dotCall64_1.0-2    cplm_0.7-11        gtable_0.3.4       spam_2.9-1         biglm_0.9-2.1     
 [6] lattice_0.21-8     quadprog_1.5-8     vctrs_0.6.2        tools_4.3.0        generics_0.1.3    
[11] curl_5.0.0         stats4_4.3.0       parallel_4.3.0     sandwich_3.0-2     tibble_3.2.1      
[16] fansi_1.0.4        xts_0.13.1         pkgconfig_2.0.3    Matrix_1.6-1.1     RColorBrewer_1.1-3
[21] actuar_3.3-2       rootSolve_1.8.2.3  lifecycle_1.0.3    compiler_4.3.0     stringr_1.5.0     
[26] fields_15.2        statmod_1.5.0      munsell_0.5.0      systemfit_1.1-30   carData_3.0-5     
[31] maps_3.4.1         pillar_1.9.0       car_3.1-2          MASS_7.3-58.4      viridis_0.6.3     
[36] StMoMo_0.4.1       abind_1.4-5        nlme_3.1-162       fracdiff_1.5-2     tidyselect_1.2.0  
[41] fanplot_4.0.0      gnm_1.1-5          stringi_1.7.12     reshape2_1.4.4     purrr_1.0.1       
[46] qvcalc_1.0.3       tseries_0.10-54    splines_4.3.0      grid_4.3.0         colorspace_2.1-0  
[51] cli_3.6.1          magrittr_2.0.3     relimp_1.0-5       utf8_1.2.3         withr_2.5.0       
[56] forecast_8.21.1    tweedie_2.3.5      TTR_0.24.3         quantmod_0.4.25    nnet_7.3-18       
[61] gridExtra_2.3      timeDate_4022.108  zoo_1.8-12         expint_0.1-8       urca_1.3-3        
[66] coda_0.19-4        lmtest_0.9-40      viridisLite_0.4.2  rlang_1.1.0        Rcpp_1.0.10       
[71] glue_1.6.2         DBI_1.1.3          rstudioapi_0.14    minqa_1.2.5        R6_2.5.1          
[76] plyr_1.8.8
```