library(tidyverse)
library(caret)
library(recipes)
library(pROC)
library(devtools)
library(here)
##
load_all(here("code","rlucas"))  ## -> library(rlucas)
##

##features <- read_csv("../data/training-set.csv")
features <- read_csv(here("data", "training-set.csv"), show_col_types=FALSE)
multinucs <- bins5mb %>% group_by(id) %>%
    summarize(multinucratio = sum(multinucs)/(sum(short+long)))
features <- inner_join(multinucs, features, by="id")

labels <- meta %>% select(id, assigned_group)
labels <- labels %>% filter(assigned_group != 1) %>%
    mutate(type = ifelse(assigned_group==2, "healthy", "cancer")) %>%
    select(-assigned_group)

features <- inner_join(labels, features, by=c("id"="id"))
features <- features  %>% select(-starts_with("cov_"))

features <- features %>% mutate(clinical_smokingstatus=factor(clinical_smokingstatus,
                                                              c("never", "former", "current")),
                                clinical_COPD=as.integer(clinical_COPD))
recipe_seq <- recipe(type ~ ., data=features) %>%
    step_rm(starts_with("clinical_"), multinucratio) %>%
    update_role(id, new_role = "ID") %>%
    step_pca(starts_with("ratio_"), prefix = "ratio_pc_",  threshold=0.90)  %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())


recipe_full_lasso2 <- recipe(type ~ ., data=features) %>%
    update_role(id, new_role = "ID") %>%
    step_rm(clinical_YKL40, clinical_IL6, clinical_CRP, clinical_cfdna_conc,
            clinical_nlratio, multinucratio, clinical_bmi) %>%
    ##step_medianimpute(clinical_packyears) %>%
    step_impute_median(clinical_packyears) %>%
    step_log(clinical_CEA) %>%
    step_dummy(clinical_smokingstatus) %>%
    step_pca(starts_with("ratio_"), prefix = "ratio_pc_",  threshold=0.90) %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())

glmnetGrid <- expand.grid(
    alpha = 1,
    lambda = 10^seq(-5, -1, length.out = 100))
#### Train models
set.seed(1234)
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(features$type, 5, 10),
                     summaryFunction = twoClassSummary)

model_seq <- caret::train(recipe_seq,
                          data = features,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl,
                          metric="ROC")


model_full_lasso2 <- caret::train(recipe_full_lasso2,
                          data = features,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl,
                          metric="ROC")


features <- features %>% mutate(rowIndex = seq_len(nrow(.)))
ids <- inner_join(features %>% select(id, rowIndex), labels, by="id")

pred.full.lasso2 <- model_full_lasso2$pred
pred.full.lasso2 <- pred.full.lasso2 %>%
    group_by(rowIndex) %>%
    summarize(score.full.lasso2 = mean(cancer))

pred.seq <- model_seq$pred
pred.seq <- pred.seq %>%
    group_by(rowIndex) %>%
    summarize(score.seq = mean(cancer))

preds <- inner_join(pred.seq, pred.full.lasso2, by="rowIndex")

preds <- inner_join(ids, preds, by="rowIndex")

preds <- preds %>% select(-rowIndex)
##save(preds, file="prediction_lucas_c1.rda")

##saveRDS(model_seq, "../models_c1/model_seq_glm.rds")
saveRDS(model_seq, here("code", "models_c1/model_seq_glm.rds"))
saveRDS(model_full_lasso2, here("code", "models_c1/model_full_lasso2.rds"))

tmp <- readRDS(here("code", "models_c1", "model_seq_glm.rds"))
pr <- prep(tmp$recipe)
traindat <- model_seq$trainingData %>%
    mutate(id=factor(id),
           type=factor(type))
model_input <- bake(pr, new_data = traindat)
loadings <- data.frame(tidy(pr, number = 2))
