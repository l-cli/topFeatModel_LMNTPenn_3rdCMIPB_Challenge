# use random forest to choose training features for each task
# then use XGboost

library(omicade4)
library(mogsa)
library(RSpectra)
library(lubridate)
library(glmnet)
library(tidyverse)
library(ggplot2)
library(randomForest)
library(impute)  # For KNN imputation
library(xgboost) # For XGBoost model


setwd('~/Desktop/cmipb_2024/data/')
source('https://raw.githubusercontent.com/akonstodata/mcia_mbpca/main/R/MCIA_mbpca_extra.R')

################### 
# define variables
###################
df_source <- readRDS("master_allData_batchCorrected.RDS")

subject_2020 = read_tsv("2020LD_subject.tsv") 
subject_2021 = read_tsv("2021LD_subject.tsv") 
subject_2022 = read_tsv("2022LD_subject.tsv") 
subject_2023 = read_tsv("2023BD_subject.tsv")
subject = bind_rows(subject_2020, subject_2021, subject_2022, subject_2023)

DATASET <- c("2020_dataset", "2021_dataset", "2022_dataset")
TEST_DATASET = "2023_dataset"
TIMEPOINTS <- c(0, 1, 3, 14)

META_COLS <- c("specimen_id", "subject_id", "timepoint", "dataset", 
               "biological_sex", "infancy_vac", "age_at_boost")
ABTITER_COLS <- c("IgG_PT")
RNA_COLS <- c("CCL3")
CELL_COLS <- c("Monocytes")

DEMOGRAPHY_COLS <- c("age_at_boost", "biological_sex", "infancy_vac")

TASK_COLS <- c("Monocytes_D1", "CCL3_D3", "IgG_PT_D14")
TASKS_BASELINES <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
BASE_COLS <- c("Monocytes_D0", "IgG_PT_D0", "CCL3_D0")

LOG_TRANS_COLS <- c("CCL3_D0", "IgG_PT_D0",
                    "CCL3_D3",  "IgG_PT_D14")

################### 
# define variables
###################
# Obtain task df
metaDf <- data.frame(df_source[["subject_specimen"]]) %>% merge(subject)
metaDf["age_at_boost"] <- as.numeric(round(difftime(metaDf$date_of_boost, metaDf$year_of_birth,units="weeks")/52, 2))
metaDf_sel <- metaDf[, META_COLS] %>%
    data.frame()

abtiterDf <- df_source[["plasma_ab_titer"]]$batchCorrected_data %>%
    t() %>%
    data.frame() 

abtiterDf$specimen_id <- as.numeric(rownames(abtiterDf))
abtiterDf_sel <- data.frame(abtiterDf[, c("specimen_id", ABTITER_COLS)])

# Add log2 (11/10/2024)
rnaDf <- df_source[["pbmc_gene_expression"]]$tpm$batchCorrected_data %>%
    t() %>%
    data.frame() %>%
    mutate_all(~ log2(. + 1))

rnaDf$specimen_id <- as.numeric(rownames(rnaDf))
tasks_seq <- c('ENSG00000277632')
for (i in 1:length(tasks_seq)){
    rnaDf_sel <- data.frame(rnaDf %>% rename_at(vars(starts_with(tasks_seq[i])), ~RNA_COLS[i]))
}
rnaDf_sel <- data.frame(rnaDf_sel[, c("specimen_id", RNA_COLS)])


cellDf <- df_source[["pbmc_cell_frequency"]]$batchCorrected_data %>%
    t() %>%
    data.frame()  

cellDf$specimen_id <- as.numeric(rownames(cellDf))
cellDf_sel <- data.frame(cellDf[, c("specimen_id", CELL_COLS)])

list_df <- list(metaDf_sel, cellDf_sel, abtiterDf_sel, rnaDf_sel)
df_merge <- list_df %>% reduce(full_join, by="specimen_id")
df_merge <- df_merge[df_merge$timepoint %in% TIMEPOINTS, ]

df_pivot <- df_merge[, names(df_merge)!="specimen_id"] %>%
    pivot_wider(id_cols=c("subject_id", "dataset", "biological_sex",
                          "infancy_vac", "age_at_boost"),
                names_from = timepoint,
                values_from = all_of(c(CELL_COLS, RNA_COLS, ABTITER_COLS)),
                names_sep = "_D")

#df_pivot <- df_pivot[df_pivot$dataset %in% DATASET, ]
df_pivot <- data.frame(df_pivot %>%
                           mutate(across(everything(),  ~ case_when(.x >=0 ~ .x))))



targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")

fc_cols <- paste(targetY, "FC", sep="_")

ranked_cols <- paste(c("age_at_boost", targetX, targetY, fc_cols), "Rank", sep="_")

df <- df_pivot[, c("subject_id", "dataset", DEMOGRAPHY_COLS, targetX, targetY)]

rankingFunction <- function(x) {
    as.numeric(rank(-x, ties.method = "min", na.last = "keep"))
}

targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")

# Add log2 (11/10/2024)
df[,"Monocytes_D1_FC"] <- log2(df[, "Monocytes_D1"] / df[, "Monocytes_D0"])
df[,"CCL3_D3_FC"] <- log2(df[, "CCL3_D3"] / df[, "CCL3_D0"])
df[,"IgG_PT_D14_FC"] <- log2(df[, "IgG_PT_D14"] / df[, "IgG_PT_D0"])

df[, ranked_cols] <- apply(df[, c("age_at_boost", targetX, targetY, fc_cols)],
                           2, rankingFunction)
df <- data.frame(df)
rownames(df) = df$subject_id

olink = df_source$plasma_cytokine_concentrations_by_olink$batchCorrected_data |>
    t() |>
    data.frame() 

olink$specimen_id <- as.numeric(rownames(olink))

# Function to perform random forest feature selection
select_top_features <- function(train_data, target_column, n_features = 10) {
    rf_model <- randomForest(as.formula(paste(target_column, "~ .")), data = train_data, importance = TRUE)
    feature_importance <- importance(rf_model)
    selected_features <- names(sort(feature_importance[,1], decreasing = TRUE))[1:n_features]
    return(selected_features)
}

make_train_df <- function(task_specific_df, target) {
    # select non-NA rows
    specimen_t0_train <- metaDf %>% filter(timepoint == 0, dataset %in% DATASET) %>% pull(specimen_id)
    specimen_t0_train <- Reduce(intersect, list(specimen_t0_train, task_specific_df$specimen_id))
    print(length(specimen_t0_train))
    task_specific_train = task_specific_df[as.character(specimen_t0_train), ]
    # add metadata
    task_specific_train$age<-metaDf[rownames(task_specific_train),"age_at_boost"]
    task_specific_train$infancy_vac<-as.numeric(as.factor(metaDf[rownames(task_specific_train),'infancy_vac']))
    task_specific_train$biological_sex<-as.numeric(as.factor(metaDf[rownames(task_specific_train),'biological_sex']))
    # join with task_mat
    task_specific_train$subject_id <- metaDf[rownames(task_specific_train), "subject_id"]
    task_mat = df[metaDf$subject_id[metaDf$specimen_id %in% specimen_t0_train], c("subject_id", target)]
    task_specific_train <- merge(task_specific_train, task_mat, by = "subject_id")
    # clean up
    task_specific_train <- task_specific_train %>% select(-specimen_id)
    task_specific_train <- task_specific_train %>% select(-subject_id)
    
    return(task_specific_train)
}

# make_test_df <- function(task_specific_df) {
#     # select non-NA rows
#     all_specimen_t0 = metaDf %>% filter(timepoint == 0, dataset == TEST_DATASET) %>% pull(specimen_id)
#     specimen_t0_test = Reduce(intersect, list(specimen_t0_test, task_specific_df$specimen_id))
#     print(length(specimen_t0_test))
#     task_specific_test = task_specific_df[as.character(specimen_t0_test), ]
#     # add metadata
#     task_specific_test$age <- metaDf$age_at_boost[metaDf$specimen_id %in% specimen_t0_test]
#     task_specific_test$infancy_vac<-as.numeric(as.factor(metaDf$infancy_vac[metaDf$specimen_id %in% specimen_t0_test]))
#     task_specific_test$biological_sex<-as.numeric(as.factor(metaDf$biological_sex[metaDf$specimen_id %in% specimen_t0_test]))
#     task_specific_test <- task_specific_test %>% select(-specimen_id)
#     
#     return(task_specific_test)
# }

make_test_df <- function(task_specific_df) {
    all_specimen_t0 = metaDf %>% filter(timepoint == 0, dataset == TEST_DATASET) %>% pull(specimen_id)
    specimen_t0_test <- Reduce(intersect, list(all_specimen_t0, task_specific_df$specimen_id))
    print(length(specimen_t0_test))
    task_specific_test <- task_specific_df[as.character(specimen_t0_test), ]
    # Identify missing specimen IDs
    missing_specimen_ids <- setdiff(all_specimen_t0, specimen_t0_test)
    # Impute missing data
    if (length(missing_specimen_ids) > 0) {
        print("missing data detected!")
        imputed_data <- sapply(task_specific_test, function(col) {
            if (is.numeric(col)) {
                mean(col, na.rm = TRUE) # Impute numeric columns with their mean
            } else {
                as.character(col[1]) # Impute non-numeric columns with the first value
            }
        })
        imputed_df <- as.data.frame(matrix(
            rep(unlist(imputed_data), length(missing_specimen_ids)), 
            nrow = length(missing_specimen_ids), 
            byrow = TRUE
        ))
        colnames(imputed_df) <- colnames(task_specific_test)
        rownames(imputed_df) <- missing_specimen_ids
        task_specific_test <- rbind(task_specific_test, imputed_df)
    }
    task_specific_test <- task_specific_test[as.character(all_specimen_t0), ]
    task_specific_test$age <- metaDf$age_at_boost[match(all_specimen_t0, metaDf$specimen_id)]
    task_specific_test$infancy_vac <- as.numeric(as.factor(metaDf$infancy_vac[match(all_specimen_t0, metaDf$specimen_id)]))
    task_specific_test$biological_sex <- as.numeric(as.factor(metaDf$biological_sex[match(all_specimen_t0, metaDf$specimen_id)]))
    task_specific_test <- task_specific_test %>% arrange(specimen_id)
    task_specific_test <- task_specific_test %>% select(-specimen_id)
    return(task_specific_test)
}

load("../optim_nfeatures.RData")
load("../optim_features_list.RData")

# IgG_PT task ------------------------------------------------------------------

# train then predict with optimal number of features
targetY <- c("IgG_PT_D14")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)

i = 1
set.seed(1)

test_df <- make_test_df(abtiterDf)
predictions_df <- data.frame(matrix(nrow = length(test_df$IgG_PT), ncol = 6))
rownames(predictions_df) <- rownames(test_df)

for (t in targets) {
    train_df <- make_train_df(abtiterDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeatures = optim_nfeatures[1]
    } else {
        nfeatures = optim_nfeatures[2]
    }
    selected_features <- optim_features_list[[i]]
    dtrain <- train_df[c(selected_features, t)]
    dtest <- test_df[selected_features]
    all_preds<-c()
    
    xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
    xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features]))
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
    xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
    predictions <- predict(xgb_model, xgb_test)
    predictions_df[, i] <- predictions
    colnames(predictions_df)[i] <- t
    i = i + 1
}


# Monocytes Task ---------------------------------------------------------------

# train then predict with optimal # of features
targetY <- c("Monocytes_D1")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
i = 1
set.seed(1)
test_df <- make_test_df(cellDf)

for (t in targets) {
    train_df <- make_train_df(cellDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeatures = optim_nfeatures[3]
    } else {
        nfeatures = optim_nfeatures[4]
    }
    selected_features <- optim_features_list[[i + 2]]
    dtrain <- train_df[c(selected_features, t)]
    dtest <- test_df[selected_features]
    all_preds<-c()
    
    xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
    xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features]))
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
    xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
    predictions <- predict(xgb_model, xgb_test)
    predictions_df[, i + 2] <- predictions
    colnames(predictions_df)[i + 2] <- t
    i = i + 1
}



# CCL3 Task --------------------------------------------------------------------
# train then test with optimal number of features
targetY <- c("CCL3_D3")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)

i = 1
set.seed(1)
test_df <- make_test_df(rnaDf)

for (t in targets) {
    train_df <- make_train_df(rnaDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeatures = optim_nfeatures[5]
    } else {
        nfeatures = optim_nfeatures[6]
    }
    selected_features <- optim_features_list[[i + 4]]
    dtrain <- train_df[c(selected_features, t)]
    dtest <- test_df[selected_features]
    all_preds<-c()
    
    xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
    xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features]))
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
    xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
    predictions <- predict(xgb_model, xgb_test)
    predictions_df[, i + 4] <- predictions
    colnames(predictions_df)[i + 4] <- t
    i = i + 1
}

# get ranking
rownames(predictions_df) <- metaDf$subject_id[metaDf$specimen_id %in% rownames(predictions_df)]
rankings_df <- sapply(predictions_df, rankingFunction)
rownames(rankings_df) <- rownames(predictions_df)
write_tsv(data.frame(rankings_df), file = "../rankings_df.tsv")
