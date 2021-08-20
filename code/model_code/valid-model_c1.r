library(devtools)
library(here)
load_all(here("code","rlucas"))  ## -> library(rlucas)
library(tidyverse)
library(caret)

features <- read_csv(here("data","testing-set.csv"))

model <- readRDS("../models_c1/model_seq_glm.rds")

features2 <- features  %>% select(-starts_with("cov_"), -id)
modelpreds <- predict(model, newdata=features2, type="prob")
modelpreds <- modelpreds %>% mutate(id = features$id)
write_csv(modelpreds, "../models_c1/validation_preds.csv")

joinedmeta <- left_join(valid_meta, modelpreds, by="id")
