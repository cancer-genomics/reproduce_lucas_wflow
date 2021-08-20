
args <- commandArgs(trailingOnly = TRUE)
zpath <- args[1]
covpath <- args[2]
outfile <- args[3]

#collect all the bins
library(data.table)
suppressMessages(library(devtools))
library(dplyr)
library(tidyverse)

zscores<- read_csv(zpath)
cov <- read_csv(covpath)

setDT(zscores)
setDT(cov)
setkey(cov, id)
setkey(zscores, id)
full_features <- zscores[cov]


write_csv(full_features,outfile)
