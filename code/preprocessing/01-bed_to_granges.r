library(getopt)

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
bedfile <- options.args[1]
outdir <- options.args[2]
names <- gsub(".bed", "", basename(bedfile))

print(names)
print(bedfile)

suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(devtools))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))

## Read beed file
bed <- fread(bedfile)
setnames(bed, c("chr", "start", "end", "mapq"))

## bed files are 0 indexed at start position and 1 indexed at end.
## For consistency in R, add 1 to start position but leave end alone.
## Filter on mapq to retain only high quality alignments.
bed <- bed[,start:=start+1][mapq >= 30]
frags <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)

refdir <- "/dcl01/scharpf/data/data-warehouse/references/genome"
cytosines <- readRDS(file.path(refdir, "cytosine_ref.rds"))

## filter outliers. Very large fragments are likely alignment artifacts.
w.all <- width(frags)
# q.all <- quantile(w.all, c(0.999))
frags <- frags[which(w.all < 1000)]

## Get fragment level GC. Use lookup for where all Cs and Gs are in hg19.
## cytosines contains C locations on both strands, so can
## use ignore.strand=TRUE to get G locations also
frags$gc_count <- countOverlaps(frags, cytosines, ignore.strand=TRUE)
rm(cytosines)
gc()

frags$w <- width(frags)
frags$gc <- frags$gc_count/frags$w

saveRDS(frags, file.path(outdir, paste0(names, ".rds")) )
q('no')
