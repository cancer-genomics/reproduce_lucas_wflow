
args <- commandArgs(trailingOnly = TRUE)
seqtype <- args[1]
outfile <- args[2]
datadir <- args[3]

#collect all the bins
library(data.table)
suppressMessages(library(devtools))
library(dplyr)
library(tidyverse)

print("loaded!")
filenames<-list.files(datadir,patter=".csv",full.names=TRUE)
#filenames<-list.files(datadir,pattern=".rds",full.names=TRUE)

x <- rbindlist(lapply(filenames, fread))
#x <- rbindlist(lapply(filenames, readRDS))


x<-setDT(x %>% filter(chr != "chrX"))
#x<-setDT(x %>% filter(chromosome != "chrX"))

#setnames(x, "sample", "id")

setkey(x, id)
x[,cov:=short+long]

#binding the reference info to the samples
if (seqtype=="nova") {
	load_all("/dcs04/scharpf/data/annapragada/Novaseq_stuff/target_distributions/PlasmaToolsNovaseq.hg19")

	} else if (seqtype == "hi") {
	load_all("/dcs04/scharpf/data/annapragada/Novaseq_stuff/target_distributions/PlasmaToolsHiseq.hg19")
	}
library(PlasmaTools.lucas)  
#suppressMessages(load_all("/dcl01/scharpf/data/pipeline-hub/pipeline-code/pipeline-devel/PlasmaTools"))

refbins <- ref54_bins
refbins <- data.table(refbins)
refids <- refbins$id
setkey(refbins, id)
refbins[,cov:=short+long]

binsforzscores <- rbind(refbins, x)
binsforzscores<-countAndNormalize(binsforzscores, measure="cov")
armmeansdt_all <- PlasmaTools:::getArmMeans(binsforzscores)
armmeansdt_ref <- armmeansdt_all %>% filter(id %in% refids)
armmeansdt_lucas <- armmeansdt_all %>% filter(!id %in% refids)
stats <- armmeansdt_ref %>% group_by(arm) %>% dplyr::summarize(Mean=mean(armmean)) 
stats1 <- armmeansdt_ref %>% group_by(arm)  %>% dplyr::summarize(STD=sd(armmean))
stats <- inner_join(stats,stats1)

zscores <- rep(NA, length(armmeansdt_lucas$armmean))
for (i in 1:length(armmeansdt_lucas$armmean)){
  m <- armmeansdt_lucas[i]$armmean
  a <- armmeansdt_lucas[i]$arm
  mu <- (stats %>% filter(arm==a))$Mean
  sigma <- (stats %>% filter(arm==a))$STD
  zscores[i] <- (m - mu)/sigma
}
armmeansdt_lucas[,zscore := zscores]

#reframe the data-table shape
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
armmeansdt_lucas[,arm:=factor(arm, armlevels)]
armmeansdt_lucas[,armvar:=factor(paste0("zscore_", arm),  paste0("zscore_", armlevels))]
features.zscores_lucas <- dcast(armmeansdt_lucas, id  ~ armvar, value.var="zscore")

write_csv(features.zscores_lucas,outfile)
