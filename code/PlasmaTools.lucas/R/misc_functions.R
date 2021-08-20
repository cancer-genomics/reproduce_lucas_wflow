#### General functions for operations on the genome
#### Functions in here likely can be sped up with compiled code

#input a bed file, get out a wig file
bedtowig <- function(pathToBed){
    
    input <- pathToBed
    out <- gsub(".bed", ".wig", input)
    
    library(data.table)
    library(rtracklayer)
    tilesdt <- fread("../bins_500kb.csv")
    
    fragbed <- fread(input)
    setnames(fragbed, paste0("V",1:4), c("chr", "start", "end", "mapq"))
    fragbed <- fragbed[mapq >= 20][,`:=`(start=start+1, end=end+1)]
    
    setkey(tilesdt, chr, start, end)
    
    bins <- foverlaps(fragbed, tilesdt, type="within", nomatch=0)
    bins2 <- bins[,.(score=.N), by=c("chr", "start", "end")]
    setkey(bins2, chr, start, end)
    bins2 <- bins2[tilesdt]
    bins2[is.na(score),score:=0]
    bins2[,`:=`(gc=NULL, map=NULL)]
    
    gr <- makeGRangesFromDataFrame(bins2, keep.extra.columns=TRUE)
    
    gr$score <- as.numeric(gr$score)*2
    
    export.wig(gr, out, genome="hg19")
    q('no')
    
    print("bed to wig file conversion is done")
}


### function for getting running average
filterv <- function(x, filter, ...){
 res <- stats::filter(x, filter,...)
 ##now fill out the NAs
 M <- length(filter)
 N <- (M- 1)/2
 L <- length(x)
 for(i in 1:N){
   w <- filter[(N-i+2):M]
   y <- x[1:(M-N+i-1)]
   res[i] <- sum(w*y)/sum(w)
   w <- rev(w)
   ii <- (L-(i-1))
   y <- x[(ii-N):L]
   res[ii] <- sum(w*y)/sum(w)
 }
 return(res)
}

# Kolmogorov-Zurbenko filter (kative moving average)
# x: numeric vector
# w: window width (must be an odd integer)
# k: number of iterations
# optimal: should the optimal number of iterations be used?
# tolerance: to control optimal number of iterations
# Credit to: Jean-Philippe Fortin

kz <- function(x, w=3, k=1, na.rm=TRUE, optimal=FALSE, tolerance=0.05, verbose=TRUE){

    check.integer <- function(N){
        !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
    }
    if (!check.integer(w) | w<=0 | (w %% 2)==0){
        stop("w must be a strictly positive odd integer")
    }
    stopifnot(is.numeric(x))

    w <- (w-1)/2
    .movingAverage <- function(x, w=1, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)

        window.mean <- function(x, j, w, na.rm=na.rm){
            if (w>=1){
                return(mean(x[(j-(w+1)):(j+w)], na.rm=na.rm))
            } else {
                return(x[j])
            }    
        }

        for (i in (w+1):(n-w)){
            y[i] <- window.mean(x,i,w, na.rm)
        }
        for (i in 1:w){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-w+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }

    .iterativeMovingAverage <- function(x, w, k=1, na.rm=TRUE){
        for (i in 1:k){
            x <- .movingAverage(x, w=w, na.rm=na.rm)
        } 
        x
    }


    .optimalMovingAverage <- function(x, w, na.rm=TRUE){
        continue <- TRUE
        k <- 0
        xs <- list()
        xs[[1]] <- x
        while (continue){
            k <- k+1
            xs[[k+1]] <- .movingAverage(xs[[k]], w=w, na.rm = na.rm)
            diff <- sum((xs[[k+1]]-xs[[k]])^2)/sum(xs[[k]]^2)
            if (diff < tolerance){
                continue  <- FALSE
            }
        }
        if (verbose){
            cat(paste0("Optimal number of iterations: ", k-1, "\n"))
        }
        return(xs[[k]])
    }



    # Iterative part:
    if (w >=1){
        
        if (!optimal){
            x <- .iterativeMovingAverage(x, w=w, na.rm=na.rm, k=k)
        } else {
            x <- .optimalMovingAverage(x, w=w, na.rm)
        }
        
    }
    x
}


mappability <- function(map, bins){
  ## summarize mappability by bin
  score <- rep(NA, length(bins))
  o <- findOverlaps(bins, map)
  if(length(o)==0) return(score)
  j <- subjectHits(o)
  i <- queryHits(o)
  subjectIndex <- split(j, i)
  w <- width(map)
  s <- map$score
  ## weight the mappability score by the width of the mappability interval
  avg <- sapply(subjectIndex, function(i, score, weight){
    (sum(score[i]*weight[i],na.rm=TRUE))/sum(weight[i], na.rm=TRUE)
  }, score=s, weight=w)
  score[as.integer(names(avg))] <- avg
  score
}

### Replace GCcontent with biostrings?
bin.genome <- function(binsize, build="hg19", chrX=TRUE) {
    library(BSgenome.Hsapiens.UCSC.hg19)

    ## +1 only if granges, not if bed -- should handle with more care.
    .countFilteredWidth <- function(start, end, x1, x2) {
        x11 <- pmax(start, x1)
        x22 <- pmin(end, x2)
        w <- x22 - x11 + 1
        return(w)
    }

    tiles <- tileGenome(keepStandardChromosomes(seqinfo(Hsapiens)),
                        tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
    tiles$gc <- GCcontent(Hsapiens, tiles)[,1]
    tilesdt <- as.data.table(tiles)
    tilesdt <- tilesdt[width == binsize]
    setnames(tilesdt, "seqnames", "chr")

    armsdt <- getArms(build, chrX)

    setkey(armsdt, chr, start, end)
    tilesdt <- foverlaps(tilesdt, armsdt, type="within", nomatch=NULL)
    tilesdt[,`:=`(start=NULL, end=NULL, width=NULL)]
    setnames(tilesdt, c("i.start", "i.end"), c("start", "end"))

    map <- rtracklayer::import.bw("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig")
    mapdt <- as.data.table(map)
    setnames(mapdt, "seqnames", "chr")

    setkey(tilesdt, chr, start, end)
    mapbins <- foverlaps(mapdt, tilesdt, type="any", nomatch=0)
    mapbins[,width2:=.countFilteredWidth(start, end, i.start, i.end)]

    mapbins2 <- mapbins[,.(map=weighted.mean(score, width2)),
                           by=c("chr", "start", "end")]
    setkey(mapbins2, chr, start, end)
    mapbins2 <- mapbins2[tilesdt]
    mapbins2 <- mapbins2[is.na(map), map:=0]
    mapbins2[,strand:=NULL]

    ## Count filtered bases (maybe should go in binFragments?)
    filtersdt <- as.data.table(filters.hg19)
    setnames(filtersdt, "seqnames", "chr")
    fov <- foverlaps(filtersdt, mapbins2, type="any", nomatch=NULL)
    fov[,width2:=.countFilteredWidth(start, end, i.start, i.end)]
    fov <- fov[,.(filtered.bases=sum(width2)), by=.(chr, start, end)]
    setkey(fov, chr, start, end)
    mapbins2 <- fov[mapbins2]
    setcolorder(mapbins2, c("chr", "start", "end", "arm", "gc", "map"))

    mapbins2[is.na(filtered.bases), filtered.bases:=0][]

}

getArms <- function(build="hg19", chrX=FALSE) {
    library(BSgenome.Hsapiens.UCSC.hg19)
    chromosomes <- GRanges(paste0("chr", c(1:22, "X")),
                           IRanges(1, seqlengths(Hsapiens)[1:23]),
                           seqinfo=seqinfo(Hsapiens))

    tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]
    # tcmeres is 0-indexed? TODO: Fix. Trim for now.

    arms <- GenomicRanges::setdiff(chromosomes, trim(tcmeres))
    arms <- arms[-c(25,27,29,41,43)]

    armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                   "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                   "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                   "19p", "19q","20p","20q","21q","22q", "Xp", "Xq")

    if(!chrX) arms <- arms[1:39]
    arms$arm <- armlevels
    armsdt <- as.data.table(arms)
    setnames(armsdt, "seqnames", "chr")
    armsdt[,`:=`(width=NULL, strand=NULL)]
    armsdt[,arm:=factor(arm, armlevels)][]
}
