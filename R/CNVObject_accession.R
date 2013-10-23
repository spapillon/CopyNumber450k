#
# Accession methods
#

setMethod("intensityMatrix",
          signature("CNVObject"),
          function(object) object@intensity_matrix)

setMethod("RGSetSummary",
          signature("CNVObject"),
          function(object) object@RGSetSummary)

setMethod("probesAnnotation",
          signature("CNVObject"),
          function(object) object@probe_annotation)

setMethod("segments",
          signature("CNVObject"),
          function(object) object@segments)

setMethod("filters",
          signature("CNVObject"),
          function(object) object@filters)

setMethod("sampleGroups",
          signature("CNVObject"),
          function(object) object@sample_groups)

setMethod("sampleSexes",
          signature("CNVObject"),
          function(object) object@sample_sexes)

setMethod("sampleNames",
          signature("CNVObject"),
          function(object) object@sample_names)

setMethod("usedProbes",
          signature("CNVObject"),
          function(object) object@used_probes)

setMethod("isNormalized",
          signature("CNVObject"),
          function(object) object@is_normalized)
      
setMethod("sampleChipRows",
          signature("CNVObject"),
          function(object) object@sample_chip_rows)
      
setMethod("sampleChipColumns",
          signature("CNVObject"),
          function(object) object@sample_chip_columns)
      
setMethod("sampleChipIDs",
          signature("CNVObject"),
          function(object) object@sample_chip_ids)
