Files in this repository:
- random_forest_feat_selection.R: codes for training the model
- random_forest_feat_selection_test.R: codes for running the trained model on test data
- 3rdChallengeSubmission_LiMaoNguyenTan.tsv: results generated on test data

Dependencies (from SessionInfo()):

R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] xgboost_1.7.8.1      impute_1.78.0        randomForest_4.7-1.2 forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
 [8] readr_2.1.5          tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0      glmnet_4.1-8         Matrix_1.7-1        
[15] lubridate_1.9.3      RSpectra_0.16-2      mogsa_1.38.0         omicade4_1.44.0      ade4_1.7-22         

loaded via a namespace (and not attached):
 [1] DBI_1.2.3                   bitops_1.0-9                GSEABase_1.66.0             rlang_1.1.4                 magrittr_2.0.3             
 [6] matrixStats_1.4.1           compiler_4.4.0              RSQLite_2.3.7               png_0.1-8                   vctrs_0.6.5                
[11] pkgconfig_2.0.3             shape_1.4.6.1               crayon_1.5.3                fastmap_1.2.0               XVector_0.44.0             
[16] svd_0.5.7                   caTools_1.18.3              utf8_1.2.4                  tzdb_0.4.0                  graph_1.82.0               
[21] UCSC.utils_1.0.0            bit_4.5.0                   made4_1.78.0                zlibbioc_1.50.0             cachem_1.1.0               
[26] graphite_1.50.0             GenomeInfoDb_1.40.1         jsonlite_1.8.9              blob_1.2.4                  DelayedArray_0.30.1        
[31] parallel_4.4.0              cluster_2.1.6               R6_2.5.1                    stringi_1.8.4               RColorBrewer_1.1-3         
[36] genefilter_1.86.0           GenomicRanges_1.56.2        Rcpp_1.0.13-1               SummarizedExperiment_1.34.0 iterators_1.0.14           
[41] IRanges_2.38.1              splines_4.4.0               timechange_0.3.0            tidyselect_1.2.1            rstudioapi_0.17.1          
[46] abind_1.4-8                 gplots_3.2.0                codetools_0.2-20            lattice_0.22-6              Biobase_2.64.0             
[51] withr_3.0.2                 KEGGREST_1.44.1             survival_3.7-0              Biostrings_2.72.1           pillar_1.9.0               
[56] MatrixGenerics_1.16.0       KernSmooth_2.23-24          foreach_1.5.2               stats4_4.4.0                generics_0.1.3             
[61] S4Vectors_0.42.1            hms_1.1.3                   munsell_0.5.1               scales_1.3.0                gtools_3.9.5               
[66] xtable_1.8-4                glue_1.8.0                  scatterplot3d_0.3-44        tools_4.4.0                 data.table_1.16.2          
[71] annotate_1.82.0             XML_3.99-0.17               grid_4.4.0                  AnnotationDbi_1.66.0        colorspace_2.1-1           
[76] GenomeInfoDbData_1.2.12     cli_3.6.3                   rappdirs_0.3.3              fansi_1.0.6                 S4Arrays_1.4.1             
[81] corpcor_1.6.10              gtable_0.3.6                BiocGenerics_0.50.0         SparseArray_1.4.8           memoise_2.0.1              
[86] lifecycle_1.0.4             httr_1.4.7                  bit64_4.5.2                 MASS_7.3-61                
