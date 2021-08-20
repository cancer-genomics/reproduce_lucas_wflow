#' @include functions.R
NULL

setGeneric("short", function(x) standardGeneric("short"))
setGeneric("long", function(x) standardGeneric("long"))
setGeneric("ratios", function(x) standardGeneric("ratios"))
setGeneric("total", function(x) standardGeneric("total"))

setMethod("short", "SummarizedExperiment", function(x) assays(x)[["short"]] )
setMethod("long", "SummarizedExperiment", function(x) assays(x)[["long"]] )
setMethod("ratios", "SummarizedExperiment", function(x){
    R <- short(x)/long(x)
    R
})
setMethod("total", "SummarizedExperiment", function(x){
    T <- short(x) + long(x)
    T
})
