#' Counts number of reads overlapping bins and performs GC-correction
#'
#' Uses LOESS to return GC-adjusted read-depth for desired bins
#'
#' @param reads a \code{GRanges} object, obtained from \function{filterReads}
#' @param bins a \code{GRanges} object, obtained from \function{makeBins}
#'
#' @return a \code{GRanges} object, with columns representing the raw counts
#'  for each bin (counts), the log2 weighted count (counts.mult),
#'  the predicted log2 coverage based on LOESS (loess.pred), and the
#'  normalized coverage (adjusted). The adjusted column is used for further
#'  z-score and PA-score analysis
#' 
#' @export
countAndNormalize <- function(bins, measure) {
    ### normalize for filtered bases
    ### weight is to get count per expected width of basepairs
    ## Maybe 
    bins[,weight:=(end - start + 1)/(end - start + 1 - filtered.bases)]
    bins[,counts.mult:=log2((get(measure) * weight + 1))]
    bins[,loess.pred:=gcCorrectLoess(counts.mult, gc),by=id]
    bins[,adjusted:=counts.mult - loess.pred]
    bins[]
}

getArmMeans <- function(bins) {
        bins[,adjusted.cent:=adjusted-median(adjusted, na.rm=TRUE), by=id]
        arms <- bins[,.(armmean=mean(adjusted.cent, na.rm=TRUE)),
                     by=.(id, arm)] 
        arms <- arms[,armmean:=armmean - mean(armmean, na.rm = TRUE)]
        arms[]
}

#' Aneuploidy z-scores
#'
#' Returns the z-scores based on normalized genome coverage.
#'
#' @param bincounts.list a list, with each object being a \code{GRanges}
#'  object obtained from \function{countAndNormalize}
#' @param normal.indices a integer vector with the indices specifying
#'  which elements in bincounts.list are to be used as controls
#' @param loo logical, indicating whether the z-scores for the normal samples
#'  should be calculated leave-one-out, as described in the paper. Default is
#'  FALSE
#'
#' @return a list, with each member being a numeric vector with the arm-level
#'  z-scores for that sample
#'
#' @export
getZscores <- function(bins, normal.ids, loo = FALSE, measure="cov") {
    bins2 <- copy(bins)
    countAndNormalize(bins2, measure=measure)
    armmeansdt <- getArmMeans(bins2)

    armmeansdt[,normal.indices:=createNormalIndex(id, normal.ids)]
#     armmeans[,normal.means:=mean(armmean[normal.indices]),by=.(arm)]
    armmeansdt[,zscore:=.loo.zscore(armmean, normal.indices), by=arm]
    armmeansdt[]
}

## either figure out how to vectorize or rewrite in cpp
.loo.zscore <-  function(x, indices) {
    zscores <- rep(NA, length(x))
    for(i in seq_along(indices)) {
        if(indices[i]) {
            keep <- indices==1
            keep[i] <- FALSE
        } else {
            keep <- indices==1
        }
        mu <- mean(x[keep])
        sigma <- sd(x[keep])
        zscores[i] <- (x[i] - mu)/sigma
    }
    zscores
}
#' Aneuploidy PA-scores
#'
#' Returns the PA-scores based on normalized genome coverage.
#'
#' @param bincounts.list a list, with each object being a \code{GRanges}
#'  object obtained from \function{countAndNormalize}
#' @param normal.indices a integer vector with the indices specifying
#'  which elements in bincounts.list are to be used as controls
#' @param loo a logical scalar, indicating whether the PA-scores for the normal samples
#'  should be calculated leave-one-out, as described in the paper. Default is
#'  FALSE
#'
#' @return a vector with the PA-scores for each sample
#'
#' @export
getPAscores <- function(zscores, loo = FALSE) {

    topzs <- zscores[, .(zscore=zscore[order(abs(zscore), decreasing=TRUE)[1:5]],
                         arm=arm[order(abs(zscore), decreasing=TRUE)[1:5]]),
                     by=.(id, normal.indices)]

    pvals <-  topzs[,.(logP=sum(-log(2*pt(abs(zscore), 3,  lower=FALSE)))),
                    by=.(id, normal.indices)]


    pvals[,pascore:=abs(.loo.zscore(logP, normal.indices))]
    pvals[]
}


gcCorrectLoess <- function(counts.mult, gc) {
    lower <- 0
    upper <- quantile(counts.mult, .99, na.rm=TRUE)
    trend.points <- counts.mult > lower & counts.mult < upper
    trend.counts <- counts.mult[trend.points]
    trend.gc <- gc[trend.points]
    num.sample.points <- min(length(trend.counts), 10E3L)
    samp <- sample(1:length(trend.counts), num.sample.points)
    #    pad counts and GC
    include <- c(which(gc==min(gc)), which(gc==max(gc)))
    trend.counts <- c(counts.mult[include], trend.counts[samp])
    trend.gc <- c(gc[include], trend.gc[samp])
    initial.trend <- loess(trend.counts ~ trend.gc)
    i <- seq(min(gc, na.rm=TRUE), max(gc, na.rm=TRUE), by = 0.001)
    final.trend <- loess(predict(initial.trend, i) ~ i)
    counts.pred <- predict(final.trend, gc)
    return(counts.pred)
}

createNormalIndex <- function(ids, normal.ids){
    ifelse(ids %in% normal.ids, 1L, 0L)
}
