setClass("MethylCNVDataSet", representation(summary = "list"), contains = "eSet")

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
})

# TODO: Add validation for summary

# TODO: Add validation for (phenoData) sexes

# TODO: Add validation for (phenoData) groups

setGeneric("predictSampleSexes", function(object, threshold = -3) standardGeneric("predictSampleSexes"))

setMethod("predictSampleSexes", signature("MethylCNVDataSet"), function(object, threshold) {
    cnQuantiles <- object@summary$cnQuantiles
    diffs <- log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ])
    predicted_sexes <- ifelse(diffs <= threshold, "Female", "Male")
}) 
