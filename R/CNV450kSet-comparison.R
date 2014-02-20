################################################################################ 

setMethod("findCNV", signature("CNV450kSet"), function(object, gene_names, type) {
    if (!is.character(gene_names)) {
        stop("Argument gene_names must be a character vector.")
    } else if (!type %in% c("gain", "loss", "both")) {
        stop("Argument type must be {gain, loss, both}.")
    }
    
    segments_list <- getSegments(object)
    
    if (length(segments_list) == 0) {
        stop("Object has not been segmentized yet.")
    }
    
    if (type == "both") {
        op <- function(a) a[, "isSignificant"]
    } else if (type == "gain") {
        op <- function(a) (a[, "isSignificant"] & as.numeric(a[, "seg.mean"]) > 0)
    } else if (type == "loss") {
        op <- function(a) (a[, "isSignificant"] & as.numeric(a[, "seg.mean"]) < 0)
    }
    x <- sapply(1:length(segments_list), function(i) {
        used_CNVs <- op(segments_list[[i]])
        genes <- unique(unlist(strsplit(x = segments_list[[i]][used_CNVs, "genes"], 
            ";")))
        gene_names %in% genes
    })

    ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
    rownames(x) <- names(segments_list)
    colnames(x) <- gene_names
    x
})

################################################################################ 

setMethod("intersectCNV", signature("CNV450kSet"), function(object, sample_indices, 
    type) {
    segments_list <- getSegments(object)
    
    if (length(segments_list) == 0) {
        stop("Object has not been segmentized yet.")
    }
    
    if (missing(sample_indices)) {
        sample_indices <- 1:length(segments_list)
    }

    if (!is.numeric(sample_indices)) {
        stop("Argument sample_indices must be a numeric vector.")
    } else if (length(setdiff(sample_indices, 1:length(segments_list))) > 0) {
        stop("Argument sample_indices elements must exist in object.")
    } else if (!type %in% c("gain", "loss", "both")) {
        stop("Argument type must be {gain, loss, both}.")
    }
    
    segments_list <- segments_list[sample_indices]
    
    if (type == "both") {
        op <- function(a) a[, "isSignificant"]
    } else if (type == "gain") {
        op <- function(a) (a[, "isSignificant"] & as.numeric(a[, "seg.mean"]) > 0)
    } else if (type == "loss") {
        op <- function(a) (a[, "isSignificant"] & as.numeric(a[, "seg.mean"]) < 0)
    }
    
    x <- unlist(sapply(1:length(segments_list), function(i) {
        used_CNVs <- op(segments_list[[i]])
        unique(unlist(strsplit(x = segments_list[[i]][used_CNVs, "genes"], ";")))
    }))
    
    sort(table(x), decreasing = TRUE)
})

################################################################################ 

# Internal method for CNVObject::subgroupDifference
subgroupDifferenceCNVByType <- function(group1_CNVs, group2_CNVs, group1_size, group2_size) {
    target_CNVs <- unique(c(names(group1_CNVs), names(group2_CNVs)))
    x <- t(sapply(target_CNVs, function(cnv) {
        if (cnv %in% names(group1_CNVs)) {
            group1_hit <- group1_CNVs[cnv]
        } else {
            group1_hit <- 0
        }
        group1_fail <- group1_size - group1_hit
        
        if (cnv %in% names(group2_CNVs)) {
            group2_hit <- group2_CNVs[cnv]
        } else {
            group2_hit <- 0
        }
        group2_fail <- group2_size - group2_hit
        
        pvalue <- fisher.test(matrix(c(group1_hit, group2_hit, group1_fail, group2_fail), 
            2))$p.value
        c(group1_hit, group2_hit, pvalue)
    }))
    
    x <- x[order(x[, 3]), ]
}

################################################################################ 

setMethod("subgroupDifference", signature("CNV450kSet"), function(object, group1_indices, 
    group2_indices) {
    segments_list <- getSegments(object)
    if (length(segments_list) == 0) {
        stop("Object has not been segmentized yet.")
    }
    sample_count <- length(segments_list)
          
    if (!is.numeric(group1_indices)) {
        stop("Argument group1_indices must be a numeric vector.")
    } else if (length(setdiff(group1_indices, 1:sample_count)) > 0) {
        stop("Argument group1_indices' elements must exist in object.")
    }
    
    if (!is.numeric(group2_indices)) {
        stop("Argument group2_indices must be a numeric vector.")
    } else if (length(setdiff(group2_indices, 1:sample_count)) > 0) {
        stop("Argument group2_indices' elements must exist in object.")
    }
    
    if (length(intersect(group1_indices, group2_indices)) > 0) {
        warning("Some samples are present in both subgroups.")
    }
    
    group1_size <- length(group1_indices)
    group2_size <- length(group2_indices)
    
    # Gains
    group1_CNVs <- intersectCNV(object, group1_indices, type = "gain")
    group2_CNVs <- intersectCNV(object, group2_indices, type = "gain")
    gains <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, group2_size)
    colnames(gains) <- c("Group 1", "Group 2", "pvalue")
    
    # Losses
    group1_CNVs <- intersectCNV(object, group1_indices, type = "loss")
    group2_CNVs <- intersectCNV(object, group2_indices, type = "loss")
    losses <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, 
        group2_size)
    colnames(losses) <- c("Group 1", "Group 2", "pvalue")
    
    list(gains = gains, losses = losses)
})

################################################################################  
