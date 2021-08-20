### TODO:
#

## Read in data.table of fragments, return data.table of density
frag.density <- function(fragments, x="w", ...,  groups=NULL) {
    d.list <- function(x, ...) {
        d <- density(x, ...)
        list(x=d$x, y=d$y)
    }

    ## Add grouping arguments (gc etc) to pass into keyby
    dens <- fragments[,c(d.list(get(x), ...)), keyby=groups]
    dens[,ynorm:=y/max(y)][]
}

frag.stats <- function(fragments, cutoff=150, groups=NULL) {
    Mode <- function(x){
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    fragments[,.(mode=Mode(w), median=as.integer(median(w)),
                 mean=mean(w), size_iqr=iqr(w),
                 nfrags=.N, mononucs=sum(w>=100 & w<=250),
                 multinucs = sum(w>=250), ultrashort = sum(w<100), 
                 slratio = sum(w<=cutoff)/sum(w > cutoff & w <= 250),
                 meangc=mean(gc), mediangc=median(gc),
                 gc_iqr=iqr(gc), chrMrep=-log(sum(chr=="chrM")/.N)),
            keyby=groups][]
}

gcCorrectTarget <- function(fragments, ref, bychr=TRUE){
    fragments[, gc := round(gc, 2)]
    if(bychr) {
        DT.gc <- fragments[,.(n=.N), by=.(gc, chr)]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc, chr)]
    } else {
        DT.gc <- fragments[,.(n=.N), by=gc]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc)]
    }
# setkey(mediandt, gc, seqnames)

    if(bychr) {
        setkey(DT.gc, gc, chr)
        setkey(ref, gc, chr)
    } else {
        setkey(DT.gc, gc)
        setkey(ref, gc)
    }
#     DT.gc <- DT.gc[ref][order(chr, gc)]
    DT.gc <- DT.gc[ref]
    DT.gc[,w:=target/n]
    if(bychr) { 
        fragments[DT.gc, on= .(chr, gc), weight := i.w]
    }
    else fragments[DT.gc, on= .(gc), weight := i.w]
    fragments <- fragments[!is.na(weight)]
    fragments[]
}


## Be sure to setname(fragments, c("seqnames", "gc"), c("chr", "fraggc"))
binFrags <- function(fragments, bins, cutoff=150,
                     chromosomes=paste0("chr",c(1:22, "X"))) {
    setkey(bins, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)
    bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                         short = sum(w >= 100 & w <= cutoff ),
                         long = sum(w > cutoff & w <= 250),
                         short.cor = sum(weight[w >= 100 & w <= cutoff]),
                         long.cor = sum(weight[w > cutoff & w <= 250]),
                         ultrashort = sum(w < 100),
                         ultrashort.cor = sum(weight[w < 100]),
                         multinucs = sum(w > 250),
                         multinucs.cor = sum(weight[w > 250]),
                         mediansize = as.integer(median(w)),
                         frag.gc = mean(fraggc)),
            by=.(chr, start, end)]

    setkey(bins2, chr, start, end)
    bins2 <- bins2[bins]
    bins2 <- bins2[is.na(i.gc), which(grepl("short|long|multi", colnames(bins2))):=0]
    bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
    bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
    bins2[, bin:=1:.N]
    setcolorder(bins2, c("chr", "start", "end", "bin"))
    bins2[]
}
