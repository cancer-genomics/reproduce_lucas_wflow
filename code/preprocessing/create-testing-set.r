library(data.table)
library(devtools)
load_all("../code/rlucas")
library(readxl)
setDTthreads(1)
library(tidyverse)
library(dplyr)

features.full <- read_csv("allfeatures_valid.csv")
refids <- valid_meta[ref_panel==FALSE][,id]
features.full <- setDT(features.full)[id %in% refids]

#### Save training data

###
#old_trainingdata <- read_csv("lucas_wflow_testing-set_to_reproduce.csv")
#old_trainingdata <- old_trainingdata %>% filter(id != "PGDX26935P")

#dplyr::all_equal(old_trainingdata,features.full,convert=TRUE)


write_csv(features.full,"Hiseq_55_testing-set.csv")

