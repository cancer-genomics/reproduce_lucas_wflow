#' @export
countFragmentOverlaps <- function(fragments, tiles,
                                  path="../pcr.data/inst/extdata/fragments_gc"){
  chr19indices <- which(as.character(seqnames(tiles)) == "chr19")
  counts <- matrix(NA, length(tiles), ncol(bamviews))
  count.list <- vector("list", ncol(bamviews))
  for(j in seq_len(ncol(bamviews))){
    cat(".")
    rdsname <- paste0(colnames(bamviews)[j], ".rds")
    fragments <- readRDS(file.path(path, rdsname))
    counts <- findOverlaps(tiles, fragments) %>%
      as_tibble %>%
      mutate(gc=fragments$gc[subjectHits],
             is_chr19=queryHits %in% chr19indices) %>%
      group_by(queryHits) %>%
      summarize(gc=mean(gc),
                count=n(),
                is_chr19=all(is_chr19))
    count.list[[j]] <- counts
  }
  counts <- do.call(rbind, count.list)
  counts
}

shortCounts <- function(se){
  assays(se)[["short"]]
}

longCounts <- function(se){
  assays(se)[["long"]]
}

#' @export
ratios <- function(se){
  shortCounts(se)/longCounts(se)
}
