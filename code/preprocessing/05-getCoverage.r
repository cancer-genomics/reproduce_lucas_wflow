
args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
datadir <- args[2]

#collect all the bins
library(data.table)
suppressMessages(library(devtools))
library(dplyr)
library(tidyverse)

print("loaded!")
filenames<-list.files(datadir,pattern=".csv",full.names=TRUE)

x <- rbindlist(lapply(filenames, fread))

x_f <- x %>% filter(chr != "chrX")


setkey(x_f, id)
x_f[,cov:=short+long]
x_f[,ratio.cor:=short.cor/long.cor]
x_f[,ratio.scaled:=scale(ratio.cor), by=id]
x_f[,ratiovar:=factor(paste0("ratio_", bin), paste0("ratio_", 1:.N)), by=id]

features.ratios <- dcast(x_f, id ~ ratiovar, value.var="ratio.scaled")
###

## Create coverage ##
x_f[,cov.cor:=short.cor+long.cor]
x_f[,cov.scaled:=scale(cov.cor), by=id]
x_f[,covvar:=factor(paste0("cov_", bin), paste0("cov_", 1:.N)), by=id]

features.covs <- dcast(x_f, id ~ covvar, value.var="cov.scaled")
setkey(features.ratios, id)
setkey(features.covs, id)
features.full <- features.covs[features.ratios]








#x_f <- x_f %>% group_by(id) %>%
 # mutate(ratio.cor = short.cor/ long.cor,
  #       ratio.scaled = scale(ratio.cor),
   #      ratiovar = factor(paste0("ratio_", bin), paste0("ratio_",1:473)),
    #     cov.cor = short.cor+long.cor,
     #    cov.scaled = scale(cov.cor),
      #   covvar = factor(paste0("cov_", bin), paste0("cov_",1:473)),
       #  short.scaled = scale(short.cor))

#x_features <- x_f %>% ungroup() %>%  select(id, bin, ratio.scaled, ratiovar,short.scaled, cov.scaled, covvar)

#x_features <- x_features %>%
 #pivot_wider(id, names_from=c(bin),
              #values_from=c(ratio.scaled,short.scaled, cov.scaled)) %>% rowwise() 

write_csv(features.full,outfile)
