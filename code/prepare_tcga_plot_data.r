library(tidyverse)
library(devtools)
library(data.table)
library(GenomicRanges)
library(here)
load_all(here("code","rlucas"))

bins <- unique(GRanges(bins5mb %>% select(chr, start, end, arm, bin)))


#--------------------------------------------------------#
get_arm_ranges <- function(assembly){

	mySession <- rtracklayer::browserSession()
	GenomeInfoDb::genome(mySession) <- assembly

	gaps <- rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track="cytoBand"))
	gaps$arm <- substr(gaps$name, 1, 1)
	gaps$id <- apply(gaps[,c('chrom', 'arm')], 1, function(x) gsub('chr', '', paste0(x[1], x[2])))

	arm.ranges <- plyr::ddply(gaps, plyr::.(id), plyr::summarize, chrom = unique(chrom), start = min(chromStart), end = max(chromEnd))
	clean.arm.ranges <- subset(arm.ranges,! id %in% c('13p', '14p', '15p', '21p', '22p', 'Xp', 'Xq', 'Yp', 'Yq'))
	clean.arm.ranges <- clean.arm.ranges[,c('chrom', 'start', 'end', 'id')]
	colnames(clean.arm.ranges) <- c('seqnames', 'start', 'end', 'arm')
	return(sort(GRanges(clean.arm.ranges)))
}


find_covered_intervals <- function(intervals, data, min.coverage){
	o <- findOverlaps(intervals, data)
	ow <- width(pintersect(intervals[queryHits(o)], data[subjectHits(o)]))
	base <- intervals[queryHits(o)]
	base$covered.fraction <- ow / width(intervals[queryHits(o)])
	total.feature.coverage <- unlist(lapply(split(base, f = base$bin), FUN = function(x) sum(x$covered.fraction)))
	pos <- intervals[intervals$bin %in% names(which(total.feature.coverage >= min.coverage))]
	return(pos)
}

find_overlapping_genes <- function(ranges, genes, within = 0){
	if (within == 1){
		hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'within')$gene_name))
	}
	if (within == 0){
		hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'any')$gene_name))
	}
	return(paste(hits, collapse = ','))
}


vectorize <- function(bins, intervals){
  pos.ids <- find_covered_intervals(bins, intervals, 0.9)$bin
  out <- rep(0, length(bins))
  out[bins$bin %in% pos.ids] <- 1
  names(out) <- bins$bin
  return(out)
}
#--------------------------------------------------------#
# determine the cn profile of LUAD tumors
# threshold to call level changes is based on Davoli et al., Science, supp table S1
date()
luad.thresh <- 0.175

luad <- fread(here("data","TCGA_Lung","LUAD.tsv"))
luad$gain <- ifelse(luad$Segment_Mean > luad.thresh, 1, 0)
luad$loss <- ifelse(luad$Segment_Mean < (-1) * luad.thresh, 1, 0)

luad <- split(luad, f = luad$Sample)
luad <- lapply(luad, function(x) {y = x
                                  y$Sample = NULL
                                  y <- subset(y, Chromosome %in% seq(1,22))
                                  y$Chromosome <- sapply(y$Chromosome, function(x) paste0('chr',x))
                                  return(y)})
luad <- lapply(luad, GRanges)

luad.losses <- lapply(luad, function(x) vectorize(bins, subset(x, loss == 1)))
luad.losses <- do.call(rbind, luad.losses)
aliquot <- sapply(rownames(luad.losses), function(x) strsplit(as.character(x), split = '-')[[1]][4])
luad.loss <- apply(luad.losses[aliquot %in% c('01A', '01B', '02A'),], 2, mean)

# 518 luad
luad.gains <- lapply(luad, function(x) vectorize(bins, subset(x, gain == 1)))
luad.gains <- do.call(rbind, luad.gains)
aliquot <- sapply(rownames(luad.gains), function(x) strsplit(as.character(x), split = '-')[[1]][4])
luad.gain <- apply(luad.gains[aliquot %in% c('01A', '01B', '02A'),], 2, mean)
#--------------------------------------------------------#
date()
lusc.thresh <- 0.225

lusc <- fread(here("data","TCGA_Lung","lusc.tsv"))
lusc$gain <- ifelse(lusc$Segment_Mean > lusc.thresh, 1, 0)
lusc$loss <- ifelse(lusc$Segment_Mean < (-1) * lusc.thresh, 1, 0)

lusc <- split(lusc, f = lusc$Sample)
lusc <- lapply(lusc, function(x) {y = x
                                  y$Sample = NULL
                                  y <- subset(y, Chromosome %in% seq(1,22))
                                  y$Chromosome <- sapply(y$Chromosome, function(x) paste0('chr',x))
                                  return(y)})
lusc <- lapply(lusc, GRanges)

# 501 LUSC

lusc.losses <- lapply(lusc, function(x) vectorize(bins, subset(x, loss == 1)))
lusc.losses <- do.call(rbind, lusc.losses)
aliquot <- sapply(rownames(lusc.losses), function(x) strsplit(as.character(x), split = '-')[[1]][4])
lusc.loss <- apply(lusc.losses[aliquot %in% c('01A', '01B'),], 2, mean)

lusc.gains <- lapply(lusc, function(x) vectorize(bins, subset(x, gain == 1)))
lusc.gains <- do.call(rbind, lusc.gains)
aliquot <- sapply(rownames(lusc.gains), function(x) strsplit(as.character(x), split = '-')[[1]][4])
lusc.gain <- apply(lusc.gains[aliquot %in% c('01A', '01B'),], 2, mean)

bins$`LUAD Loss` = luad.loss
bins$`LUAD Gain` = luad.gain
bins$`LUSC Loss` = lusc.loss
bins$`LUSC Gain` = lusc.gain

bins <- data.frame(bins)
bins$pos = apply(bins[,c('start', 'end')], 1, mean)
bins[,c('start', 'end', 'width', 'strand')] = NULL
b = melt(bins, id.vars = c('seqnames', 'pos', 'arm' ,'bin'))
b$arm <- factor(b$arm, levels = unique(bins$arm))
b$disease <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][1])
b$change <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][2])

b[which(b$change == 'Loss'), 'value'] = b[which(b$change == 'Loss'), 'value'] * (-1)
#--------------------------------------------------------#

fig.data <- b
saveRDS(fig.data, file = here('data', 'TCGA_Lung', 'fig2c_p2_data.rds'))
