[![DOI](https://zenodo.org/badge/691535587.svg)](https://zenodo.org/doi/10.5281/zenodo.10118895)

# Individual claims reserving using the Aalen-Johansen estimator

This repository contains the code to replicate the manuscript [Individual claims reserving using the Aalen-Johansen estimator](https://arxiv.org/5231629). 

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

## Section 4.4 Comparison with `hirem`

- Table 3 can be reproduced with the script `table_2_feature_model.R`.

- Table 5 can be reproduced with the script `table_5_comparison_with_hirem.R`.

# Section 5. A data application on an insurance portfolio

- We obtained Figure 4 with the script  `figure_4.R`. 

- We obtained Figure 5 with the script  `figure_5.R`. 


## Section 5.1. Model comparison on different datasets

- We obtained Table 6 and Table 7 using the script `table_6_intercept_model.R`(results for the intercept model) and `table_6_feature_model.R` (results for the model with feature).

## Section 5.2. Model comparison on a single dataset

- We obtained Table 8 using the scripts `table_8_intercept_model.R` (results for the intercept model) and `table_8_feature_model.R` (results for the model with feature).

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
  
    - `simulation_aytrend4_2024_02_29_13_31.csv`
    - `simulation_aytrend5_2024_02_29_13_32.csv`
    - `simulation_aytrend6_2024_02_29_13_34.csv`
    - `simulation_aytrend7_2024_02_29_13_36.csv`
    
- Comparison with `hirem` (Section 4.4.):

    - `simulation_hirem_10.csv`
    - `simulation_hirem_extreme_10.csv`
    - `simulation_hirem_settlement_10.csv`
    - `simulation_hirem_claimmix_10.csv`

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
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] hirem_0.1.0        ChainLadder_0.2.18 tidyr_1.3.1        dplyr_1.1.2        ggplot2_3.5.0      xtable_1.8-4      
[7] data.table_1.14.8 

loaded via a namespace (and not attached):
 [1] cplm_0.7-11       gtable_0.3.4      biglm_0.9-2.1     gbm_2.1-06        xfun_0.39         lattice_0.21-8   
 [7] vctrs_0.6.2       tools_4.3.0       generics_0.1.3    stats4_4.3.0      parallel_4.3.0    sandwich_3.0-2   
[13] tibble_3.2.1      fansi_1.0.4       pkgconfig_2.0.3   Matrix_1.6-1.1    actuar_3.3-2      lifecycle_1.0.4  
[19] compiler_4.3.0    stringr_1.5.1     statmod_1.5.0     munsell_0.5.0     systemfit_1.1-30  carData_3.0-5    
[25] htmltools_0.5.5   yaml_2.3.7        pillar_1.9.0      car_3.1-2         MASS_7.3-58.4     abind_1.4-5      
[31] nlme_3.1-162      tidyselect_1.2.0  digest_0.6.31     stringi_1.7.12    reshape2_1.4.4    purrr_1.0.1      
[37] splines_4.3.0     fastmap_1.1.1     grid_4.3.0        colorspace_2.1-0  cli_3.6.1         magrittr_2.0.3   
[43] survival_3.5-5    utf8_1.2.3        withr_3.0.0       scales_1.3.0      bit64_4.0.5       tweedie_2.3.5    
[49] lubridate_1.9.2   timechange_0.2.0  rmarkdown_2.21    bit_4.0.5         AalenJohansen_1.0 zoo_1.8-12       
[55] expint_0.1-8      coda_0.19-4       evaluate_0.21     knitr_1.42        lmtest_0.9-40     rlang_1.1.3      
[61] Rcpp_1.0.10       glue_1.6.2        DBI_1.1.3         rstudioapi_0.15.0 minqa_1.2.5       R6_2.5.1         
[67] plyr_1.8.8    
```