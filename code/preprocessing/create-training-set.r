library(data.table)
library(devtools)
load_all("../code/rlucas")
library(readxl)
setDTthreads(1)
library(tidyverse)
library(dplyr)

#### Baseline covariates
setDT(meta)
# meta[,IL6:=ifelse(is.na(IL6), median(IL6, na.rm=TRUE), IL6)]
# meta[,YKL40:=ifelse(is.na(YKL40), median(YKL40, na.rm=TRUE), YKL40)]
meta[,smokestatus:=fcase(grepl("never smoker", Smoking_status), "never",
                          grepl("over", Smoking_status), "former",
                          !grepl("never|over", Smoking_status), "current")]

setnames(meta, "neutrophile_wbc_ratio", "nlratio")
features.clinical <- meta[,c("id","Patient", "nlratio", "CRP", "cfdna_conc",
                            "age", "IL6", "YKL40",  "CEA", "bmi", "Packyears",
                            "smokestatus", "COPD")]
setnames(features.clinical, c("nlratio", "CRP", "cfdna_conc", "age", "IL6",
                              "YKL40", "CEA", "bmi", "Packyears",
                              "smokestatus", "COPD"),
         c("clinical_nlratio", "clinical_CRP", "clinical_cfdna_conc",
           "clinical_age", "clinical_IL6", "clinical_YKL40", "clinical_CEA",
           "clinical_bmi", "clinical_packyears", "clinical_smokingstatus",
           "clinical_COPD"))


### Combine data

features.nonclinical <- setDT(read_csv("allfeatures.csv"))
setkey(features.clinical, id)
setkey(features.nonclinical, id)

features.full <- features.clinical[features.nonclinical,nomatch=NULL]
features.full$Patient <- NULL
features.full <- as_tibble(features.full)
#### Save training data

###
#old_trainingdata <- read_csv("lucas_wflow_training-set_to_reproduce.csv")
#old_trainingdata <- old_trainingdata %>% filter(id != "PGDX26935P")

#dplyr::all_equal(old_trainingdata,features.full,convert=TRUE)


write_csv(features.full,"Hiseq_55_training-set.csv")

