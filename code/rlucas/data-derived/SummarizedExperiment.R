library(tidyverse)
library(readxl)
library(magrittr)
library(stringr)
library(purrr)
library(kableExtra)
library(SummarizedExperiment)
library(here)
library(devtools)
devtools::load_all("../../rlucas")
load(here("data", "prediction_lucas.rda"))
load(here("data", "metadata.rda"))
data(bins1kb, package="svfilters.hg19")
load(here("data", "lucas_5mb.rda"))

##
## Clinical and demographic data
##
meta <- meta %>%
    set_colnames(tolower(colnames(.))) %>%
    set_colnames(str_replace_all(colnames(.), "\\.", "_"))
colnames(meta)[1:2] <- c("lab_id", "pgdx_id")
coldat <- meta %>%
    select(-lab_id) 

labid_to_pgdxid <- select(meta, lab_id, pgdx_id)

path <- here("inst", "extdata")
tfile <- file.path(path, "Updated_treatment_LUCAS.xlsx")
treatment <- read_excel(tfile, sheet=1) %>%
    set_colnames(tolower(colnames(.))) %>%
    set_colnames(str_replace_all(colnames(.), " ", "_"))
colnames(treatment)[c(1, 6)] <- c("lab_id", "first_line")
trt_levels <- c("Palliative therapy", "CR with curative intent",
                "Surgery", "Surgery+CR")
treatment2 <- treatment %>%
    "["(, !colnames(treatment) %in% colnames(coldat)) %>%
    left_join(labid_to_pgdxid,
              by="lab_id") %>%
    select(pgdx_id, first_line, survival_status) %>%
    mutate(first_line=factor(first_line, levels=trt_levels))
coldat2 <- left_join(coldat, treatment2, by="pgdx_id") %>%
    mutate(include_coxph=pgdx_id %in% treatment2$pgdx_id)

##
## Age group
coldat3 <- coldat2 %>%
    mutate(histology=factor(histology,
                            levels=c("SCLC",
                                     "Adenocarcinoma",
                                     "Squamous",
                                     "Metastasis from another primary cancer"))) %>%
    mutate(age2=round(age, 0),
           age2=ifelse(age2 < 55, "<55", age2),
           age2=ifelse(age2 >= 55 & age2 <= 74, "55-74",
                       age2),
           age2=ifelse(age2 >= 75, ">=75", age2),
           age2=factor(age2, levels=c("<55", "55-74", ">=75")),
           age_group=age2) %>%
    select(-age2) %>%
    mutate(stage=factor(stage, levels=c("I", "II", "III", "IV")))
if(any(is.na(coldat3$age_group))) stop('nas')

## gender
st <- read_excel(file.path(path, "Supplementary tables.xlsx"), sheet=1,
                 skip=1) %>%
    set_colnames(tolower(colnames(.))) %>%
    set_colnames(str_replace_all(colnames(.), " ", "_")) %>%
    select(patient, sex) %>%
    set_colnames(c("lab_id", "sex")) %>%
    left_join(labid_to_pgdxid, by="lab_id") %>%
    filter(!is.na(pgdx_id)) %>%
    select(pgdx_id, sex)
autoimmune <- read_excel(file.path(path, "COPD_autoimmunity_ data.xlsx"),
                      sheet=1) %>%
    set_colnames(c("lab_id", "copd", "autoimmune_dis")) %>%
    left_join(labid_to_pgdxid, by="lab_id") %>%
    filter(!is.na(pgdx_id)) %>%
    select(pgdx_id, autoimmune_dis)

coldat3 <- left_join(coldat3, st, by=c("pgdx_id", "sex")) %>%
    left_join(autoimmune, by="pgdx_id")

predictions <- preds %>%
    set_colnames(str_replace_all(colnames(.), "\\.", "_")) %>%
    select(-type) %>%
    mutate(delfi_group=ifelse(score_seq >= 0.5, "DELFI >= 0.5", "DELFI < 0.5"),
           delfi_group=factor(delfi_group, levels=c("DELFI < 0.5", "DELFI >= 0.5")))    
colnames(predictions)[1] <- "pgdx_id"
coldat4 <- coldat3 %>%
    left_join(predictions, by="pgdx_id")
stopifnot(identical(coldat4$pgdx_id, labid_to_pgdxid$pgdx_id))

coldat5 <- as(coldat4, "DataFrame")
rownames(coldat5) <- labid_to_pgdxid$lab_id

getAssayData <- function(colname, bins5mb){
    vars <- c("bin", "id", colname)
    bins5mb %>%
        select(all_of(vars)) %>%
        pivot_wider(names_from=id, values_from=starts_with(colname)) %>%
        select(-bin) %>%
        as.matrix()    
}
assaydat <- as.list(c("short.cor", "long.cor", "multinucs.cor",
                      "frag.gc")) %>%
    map(getAssayData, bins5mb) %>%
    setNames(c("short", "long", "multinuc", "frag_gc"))
assaycols <- tibble(pgdx_id=colnames(assaydat[["short"]])) %>%
    left_join(labid_to_pgdxid, by="pgdx_id")  %>%
    filter(!is.na(lab_id))
assaydat2 <- assaydat %>%
    map(function(x, id, id2){
        x <- x[, id]
        colnames(x) <- id2
        return(x)
    }, id=assaycols$pgdx_id, id2=assaycols$lab_id) %>%
    as("SimpleList")

bindata <- filter(bins5mb, !duplicated(bin)) %>%
    select(chr, start, end, bin, gc, map, arm)
gr <- select(bindata, chr, start, end) %>%
    set_colnames(c("seqnames", "start", "end")) %>%
    as_tibble() %>%
    as("GRanges")
mcols(gr) <- as(select(bindata, -c(chr, start, end)), "DataFrame")
## just use this package to get hg19 seqinfo
seqinfo(gr) <- seqinfo(bins1kb)[levels(seqnames(gr))]
names(gr) <- 1:473

se <- SummarizedExperiment(assays=assaydat2,
                           rowRanges=gr,
                           colData=coldat5)

## Add iteration 2 predictions
descr <- tibble("model"=c("score_seq",
                 "score_full_lasso1",
                 "score_full_lasso2",
                 "score_clinical",
                 "score_full1",
                 "score_full2"),
       seq_features=c("PCs of S/L ratios and z-scores",
                      "PCs of S/L ratios and z-scores",
                      "PCs of S/L ratios and z-scores",
                      "none",
                      "PCs of S/L ratios and z-scores",
                      "PCs of S/L ratios and z-scores"),
       nonseq_features=c("none",
                         "bmi, ykl40, il6, packyears, sm status, cea, copd(?)",
                         "cea, packyears, sm status, copd(?)",
                         "bmi, ykl40, il6, packyears, sm status, cea, copd(?)",
                         "bmi, ykl40, il6, packyears, sm status, cea, copd(?)",
                         "cea, packyears, sm status, copd(?)"),
       ML=rep("Logistic regression", 6),
       penalty=c(rep("LASSO", 4), c("none", "none")))
kbl(descr, "html") %>%
    kable_styling()
##
## Add metadata from Supplementary Table S1
##
s1 <- read_excel(here("inst", "extdata", "Supplementary tables.xlsx"),
                 sheet=1, skip=1) %>%
    set_colnames(tolower(colnames(.))) %>%
    set_colnames(str_replace_all(colnames(.), " ", "_")) %>%
    select(patient, histological_diagnosis, biopsy_type) %>%
    set_colnames(c("lab_id", "hist_dx", "biopsy_type"))

##
## Add score from sequence model excluding previous cancers
##
load(here("data", "prediction_lucas_v2.rda"))
dat <- preds %>%
    as_tibble() %>%
    select(id, score.seq) %>%
    set_colnames(c("pgdx_id", "score_no_priorcancer"))



##
## Model trained only on patients 50-80, packyears >20
##
## add score for high risk smokers (score_hr_smokers)
dat2 <- read_excel(here("inst", "extdata", "50_80\ yo_20py.xlsx"),
                   sheet=1) %>%
    set_colnames(tolower(colnames(.))) %>%
    set_colnames(str_replace_all(colnames(.), " ", "_")) %>%
    select(id, score.seq) %>%
    set_colnames(c("pgdx_id", "score_hr_smokers"))
dat3 <- left_join(dat, dat2, by="pgdx_id")


coldat <- colData(se) %>%
    as_tibble() %>%
    mutate(lab_id=colnames(se)) %>%
    left_join(dat3, by="pgdx_id") %>%
    left_join(s1, by="lab_id") %>%
    as(., "DataFrame")
rownames(coldat) <- colnames(se)
colData(se) <- coldat
save(se, file=here("data", "se.rda"), compression_level=9)

### New Stuff Adding for Cohort Descriptions
coldat$Cohort3=!(is.na(coldat$score_hr_smokers))
coldat$Cohort2=!(is.na(coldat$score_no_priorcancer))
coldat$Cohort1=!(is.na(coldat$score_seq))
###

colData(se) <- coldat
save(se, file=here("data", "se.rda"), compression_level=9)
