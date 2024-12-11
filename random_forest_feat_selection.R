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

rankingFunction <- function(x) {
    as.numeric(rank(-x, ties.method = "min", na.last = "keep"))
}
targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")

fc_cols <- paste(targetY, "FC", sep="_")

ranked_cols <- paste(c("age_at_boost", targetX, targetY, fc_cols), "Rank", sep="_")

df <- df_pivot[, c("subject_id", "dataset", DEMOGRAPHY_COLS, targetX, targetY)]

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

# helper functions -------------------------------------------------------------
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

make_test_df <- function(task_specific_df, target) {
    # select non-NA rows
    specimen_t0_test = metaDf %>% filter(timepoint == 0, dataset == TEST_DATASET) %>% pull(specimen_id)
    specimen_t0_test = Reduce(intersect, list(specimen_t0_test, task_specific_df$specimen_id))
    print(length(specimen_t0_test))
    task_specific_test = task_specific_df[as.character(specimen_t0_test), ]
    # add metadata
    task_specific_test$age <- metaDf$age_at_boost[metaDf$specimen_id %in% specimen_t0_test]
    task_specific_test$infancy_vac<-as.numeric(as.factor(metaDf$infancy_vac[metaDf$specimen_id %in% specimen_t0_test]))
    task_specific_test$biological_sex<-as.numeric(as.factor(metaDf$biological_sex[metaDf$specimen_id %in% specimen_t0_test]))
    task_specific_test <- task_specific_test %>% select(-specimen_id)
    
    return(task_specific_test)
}


# IgG_PT task ------------------------------------------------------------------

targetY <- c("IgG_PT_D14")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
pred_cor_abtiters <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))

i = 1
set.seed(1)
for (t in targets) {
    train_df <- make_train_df(abtiterDf, t)
    train_df <- na.omit(train_df)
    selected_features <- select_top_features(train_df, t, n_features = 10)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_abtiters[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# try different # of features for abtiter task
n <- 34
nfeatures <- seq(1, n)
pred_cor_abtiter_30 <- data.frame(matrix(nrow=2, ncol=n, dimnames = list(targets, str(nfeatures))))
print(targets)
for (nfeature in nfeatures) {
    i = 1
    set.seed(1)
    for (t in targets) {
        train_df <- make_train_df(abtiterDf, t)
        train_df <- na.omit(train_df)
        selected_features <- select_top_features(train_df, t, n_features = nfeature)
        cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
        train_df <- train_df[c(selected_features, t)]
        all_preds<-c()
        all_true<-c()
        # train all leave-one-out models
        for (j in 1:nrow(train_df)) {
            train <- 1:nrow(train_df)
            train <- train[-c(j)]
            dtrain <- train_df[train,]
            dtest <- train_df[-train,]
            
            xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
            xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
            params <- list(objective = "reg:squarederror", eval_metric = "rmse")
            xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
            predictions <- predict(xgb_model, xgb_test)
            
            true_values <- dtest[[t]]
            all_preds <- c(all_preds, predictions)
            all_true <- c(all_true, true_values)
        }
        pred_cor_abtiter_30[i, nfeature] <- cor(all_preds, all_true)
        i = i + 1
    }
}

df_long <- pred_cor_abtiter_30 %>%
    rownames_to_column("Variable") %>%
    pivot_longer(-Variable, names_to = "Column", values_to = "Value") %>%
    mutate(Column = as.integer(gsub("X", "", Column)))  # Convert Column names to integer

ggplot(df_long, aes(x = Column, y = Value, color = Variable, group = Variable)) +
    geom_line() +
    labs(x = "Column (Integer)", y = "Value", title = "Line Plot of Two Variables Across Columns") +
    theme_minimal()

t1_optim_nfeature <- which.max(abs(pred_cor_abtiter_30[1, ]))
t2_optim_nfeature <- which.max(abs(pred_cor_abtiter_30[2, ]))
optim_nfeatures <- c(t1_optim_nfeature, t2_optim_nfeature)




# Monocytes Task ---------------------------------------------------------------
targetY <- c("Monocytes_D1")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
pred_cor_cell <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))

i = 1
set.seed(1)
for (t in targets) {
    train_df <- make_train_df(cellDf, t)
    train_df <- na.omit(train_df)
    selected_features <- select_top_features(train_df, t, n_features = 15)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_cell[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# try different # of features for cell task
targetY <- c("Monocytes_D1")
fc_cols <- paste(targetY, "FC", sep="_")
n <- 40
nfeatures <- seq(1, n)
pred_cor_cell_30 <- data.frame(matrix(nrow=2, ncol=n, dimnames = list(targets, str(nfeatures))))

for (nfeature in nfeatures) {
    i = 1
    set.seed(1)
    for (t in targets) {
        train_df <- make_train_df(cellDf, t)
        train_df <- na.omit(train_df)
        selected_features <- select_top_features(train_df, t, n_features = nfeature)
        cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
        train_df <- train_df[c(selected_features, t)]
        all_preds<-c()
        all_true<-c()
        # train all leave-one-out models
        for (j in 1:nrow(train_df)) {
            train <- 1:nrow(train_df)
            train <- train[-c(j)]
            dtrain <- train_df[train,]
            dtest <- train_df[-train,]
            
            xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
            xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
            params <- list(objective = "reg:squarederror", eval_metric = "rmse")
            xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
            predictions <- predict(xgb_model, xgb_test)
            
            true_values <- dtest[[t]]
            all_preds <- c(all_preds, predictions)
            all_true <- c(all_true, true_values)
        }
        pred_cor_cell_30[i, nfeature] <- cor(all_preds, all_true)
        i = i + 1
    }
}

df_long <- pred_cor_cell_30 %>%
    rownames_to_column("Variable") %>%
    pivot_longer(-Variable, names_to = "Column", values_to = "Value") %>%
    mutate(Column = as.integer(gsub("X", "", Column)))  # Convert Column names to integer

ggplot(df_long, aes(x = Column, y = Value, color = Variable, group = Variable)) +
    geom_line() +
    labs(x = "Column (Integer)", y = "Value", title = "Line Plot of Two Variables Across Columns") +
    theme_minimal()


t1_optim_nfeature <- which.max(abs(pred_cor_cell_30[1, ]))
t2_optim_nfeature <- which.max(abs(pred_cor_cell_30[2, ]))
optim_nfeatures <- c(optim_nfeatures, t1_optim_nfeature, t2_optim_nfeature)


# CCL3 Task --------------------------------------------------------------------
targetY <- c("CCL3_D3")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
pred_cor_rna <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))

i = 1
set.seed(1)
for (t in targets) {
    train_df <- make_train_df(rnaDf, t)
    train_df <- na.omit(train_df)
    selected_features <- select_top_features(train_df, t, n_features = 20)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_rna[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# check results
rbind(pred_cor_abtiters, pred_cor_cell, pred_cor_rna)


# try different # of features for rna task
n <- 50
nfeatures <- seq(1, n)
pred_cor_rna_30 <- data.frame(matrix(nrow=2, ncol=n, dimnames = list(targets, str(nfeatures))))

for (nfeature in nfeatures) {
    i = 1
    set.seed(1)
    for (t in targets) {
        train_df <- make_train_df(rnaDf, t)
        train_df <- na.omit(train_df)
        selected_features <- select_top_features(train_df, t, n_features = nfeature)
        cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
        train_df <- train_df[c(selected_features, t)]
        all_preds<-c()
        all_true<-c()
        # train all leave-one-out models
        for (j in 1:nrow(train_df)) {
            train <- 1:nrow(train_df)
            train <- train[-c(j)]
            dtrain <- train_df[train,]
            dtest <- train_df[-train,]
            
            xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
            xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
            params <- list(objective = "reg:squarederror", eval_metric = "rmse")
            xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
            predictions <- predict(xgb_model, xgb_test)
            
            true_values <- dtest[[t]]
            all_preds <- c(all_preds, predictions)
            all_true <- c(all_true, true_values)
        }
        pred_cor_rna_30[i, nfeature] <- cor(all_preds, all_true)
        i = i + 1
    }
}

df_long <- pred_cor_rna_30 %>%
    rownames_to_column("Variable") %>%
    pivot_longer(-Variable, names_to = "Column", values_to = "Value") %>%
    mutate(Column = as.integer(gsub("X", "", Column)))  # Convert Column names to integer

ggplot(df_long, aes(x = Column, y = Value, color = Variable, group = Variable)) +
    geom_line() +
    labs(x = "Column (Integer)", y = "Value", title = "Line Plot of Two Variables Across Columns") +
    theme_minimal()

t1_optim_nfeature <- which.max(abs(pred_cor_rna_30[1, ]))
t2_optim_nfeature <- which.max(abs(pred_cor_rna_30[2, ]))
optim_nfeatures <- c(optim_nfeatures, t1_optim_nfeature, t2_optim_nfeature)
save(optim_nfeatures, file = '../optim_nfeatures.RData')


# run with optimal number of features ------------------------------------------

# IgG_PT
load("../optim_nfeatures.RData")
targetY <- c("IgG_PT_D14")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
pred_cor_abtiter <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))
optim_features_list <- list()
i = 1
set.seed(1)
for (t in targets) {
    train_df <- make_train_df(abtiterDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeatures = optim_nfeatures[1]
    } else {
        nfeatures = optim_nfeatures[2]
    }
    selected_features <- select_top_features(train_df, t, n_features = nfeatures)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    optim_features_list <- append(optim_features_list, list(selected_features))
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_abtiter[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# monocyte
targetY <- c("Monocytes_D1")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
i = 1
set.seed(1)
pred_cor_cell <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))

for (t in targets) {
    train_df <- make_train_df(cellDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeature = optim_nfeatures[3]
    } else {
        nfeature = optim_nfeatures[4]
    }
    selected_features <- select_top_features(train_df, t, n_features = nfeature)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    optim_features_list <- append(optim_features_list, list(selected_features))
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_cell[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# CCL3
targetY <- c("CCL3_D3")
fc_cols <- paste(targetY, "FC", sep="_")
targets <- c(targetY, fc_cols)
pred_cor_rna <- data.frame(matrix(nrow=2, ncol=1, dimnames = list(targets, "cor.pred.true")))

i = 1
set.seed(1)
for (t in targets) {
    train_df <- make_train_df(rnaDf, t)
    train_df <- na.omit(train_df)
    if (i == 1) {
        nfeatures = optim_nfeatures[5]
    } else {
        nfeatures = optim_nfeatures[6]
    }
    selected_features <- select_top_features(train_df, t, n_features = nfeatures)
    cat("Selected features for", t, ":", paste(selected_features, collapse = ", "), "\n")
    optim_features_list <- append(optim_features_list, list(selected_features))
    train_df <- train_df[c(selected_features, t)]
    all_preds<-c()
    all_true<-c()
    # train all leave-one-out models
    for (j in 1:nrow(train_df)) {
        train <- 1:nrow(train_df)
        train <- train[-c(j)]
        dtrain <- train_df[train,]
        dtest <- train_df[-train,]
        
        xgb_train <- xgb.DMatrix(data = as.matrix(dtrain[selected_features]), label = dtrain[[t]])
        xgb_test <- xgb.DMatrix(data = as.matrix(dtest[selected_features], label = dtest[[t]]))
        params <- list(objective = "reg:squarederror", eval_metric = "rmse")
        xgb_model <- xgboost(params = params, data = xgb_train, nrounds = 100, verbose = 0)
        predictions <- predict(xgb_model, xgb_test)
        
        true_values <- dtest[[t]]
        all_preds <- c(all_preds, predictions)
        all_true <- c(all_true, true_values)
    }
    pred_cor_rna[i, 1] <- cor(all_preds, all_true)
    i = i + 1
}

# check results
rbind(pred_cor_abtiter, pred_cor_cell, pred_cor_rna)
preds <- rbind(pred_cor_abtiter, pred_cor_cell, pred_cor_rna)
preds <- data.frame(preds)
preds$Task <- rownames(preds)
preds$Model <- "random forest + XGBoost"
preds$Correlation <- preds$cor.pred.true
preds$cor.pred.true <- NULL

# save results
save(optim_features_list, file = "../optim_features_list.RData")
write_tsv(preds, file = "../preds_correlations.tsv")
