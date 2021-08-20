suppressMessages(library(getopt))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
data.table::setDTthreads(1)

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

fragfile <- options.args[1]
outdir <- options.args[2]


sample <- gsub(".rds", "", basename(fragfile))

DT <- as.data.table(readRDS(fragfile))
DT[,id:=sample]
DT <- DT[seqnames != "chrM" & width >= 100 & width <= 220]
DT[, gc := round(gc, 2)]

DT.gc <- DT[,.(n=.N), by=.(id, gc, seqnames)]
DT.gc <- DT.gc[gc >= .20 & gc <= .80]
DT.gc <- DT.gc[order(seqnames, gc)]

filename <- file.path(outdir, paste0(sample, "_gc.csv"))
fwrite(DT.gc, file=filename)
q('no')
