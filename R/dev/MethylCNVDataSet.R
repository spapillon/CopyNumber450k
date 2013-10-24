setClass("MethylCNVDataSet", representation(summary = "list"), contains = "eSet")

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
}) 
