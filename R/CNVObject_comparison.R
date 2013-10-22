#
# Comparison Methods
#

setMethod("findCNV",
          signature("CNVObject"),
          function(object,
                   CNVs,
                   type) {
                
	if (!is.character(CNVs)) {
		stop("Expected argument CNVs to be a character vector.")
  }
  
	segments_list <- segments(object)
	filters_list  <- filters(object)

  # TODO: Highly inefficient loop because static variable "type" is
  # evaluated every iteration
  # Suggestion
  # define three functions and pass them by parameter
  # if (type == "both") op <- function (i) filters_list[[i]]
  # else if (type == "gain") op <- ...
  # else if (type == "loss") op <- ...
  # else stop("Expected argument type to be {both, gain, loss}.")
  # x <- sapply(..., function(i) { used_CNVs <- op(i); ... }
	x <- sapply(1:length(segments_list), function(i) {
        if (type == "both")
					used_CNVs <- filters_list[[i]]
				else if(type == "gain")
					used_CNVs <- segments_list[[i]][,'logratio'] > 0 & filters_list[[i]]
				else if(type == "loss")
					used_CNVs <- segments_list[[i]][,'logratio'] < 0 & filters_list[[i]]
						
				genes <- unique(unlist(strsplit(x=segments_list[[i]][used_CNVs,'genes'], ";")))
				return(CNVs %in% genes)
			})
	ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
	rownames(x) <- names(segments_list)
	colnames(x) <- CNVs
	x
})



setMethod("intersectCNV",
          signature("CNVObject"),
          function(object,
                   sample_indices,
                   type) {

	sample_count <- length(segments(object))
  
	if (missing(sample_indices)) {
		sample_indices <- 1:sample_count
  }
  
	if (!is.numeric(sample_indices) ||
      length(setdiff(sample_indices, 1:sample_count)) > 0) {
		stop("Expected argument sample_indices to be a integer vector with values corresponding to existing samples indices.")
  }
    
	segments_list <- segments(object)[sample_indices]
	filters_list  <- filters(object)[sample_indices]
  
  # TODO: See above function (findCNV)
	x <- unlist(sapply(1:length(segments_list), function(i) {
						if(type == "both")
							used_CNVs <- filters_list[[i]]
						else if(type == "gain")
							used_CNVs <- segments_list[[i]][,'logratio'] > 0 & filters_list[[i]]
						else if(type == "loss")
							used_CNVs <- segments_list[[i]][,'logratio'] < 0 & filters_list[[i]]
						unique(unlist(strsplit(x=segments_list[[i]][used_CNVs,'genes'], ";")))
					}))
  
	sort(table(x), decreasing = TRUE)
})



setMethod("subgroupDifference",
          signature("CNVObject"),
          function(object,
                   group1_indices,
                   group2_indices) {
                 
	sample_count <- length(segments(object))
  
	if (!is.numeric(group1_indices) ||
      length(setdiff(group1_indices, 1:sample_count)) > 0 ) {
		stop("Expected group1_indices argument to be a integer vector with values corresponding to existing samples indices.")
  }
  
	if(!is.numeric(group2_indices) ||
      length(setdiff(group1_indices, 1:sample_count)) > 0 ) {
    stop("Expected group2_indices argument to be a integer vector with values corresponding to existing samples indices.")
  }

	if (length(intersect(group1_indices, group2_indices)) > 0) {
		warning("Some samples are present in both subgroups.")
  }
	
	group1_size <- length(group1_indices)
	group2_size <- length(group2_indices)
			
	# Gains
	group1_CNVs <- intersectCNV(object, group1_indices, type = "gain")
	group2_CNVs <- intersectCNV(object, group2_indices, type = "gain")
	gains <- subgroupDifferenceCNVByType(group1_CNVs,
                                       group2_CNVs,
                                       group1_size,
                                       group2_size) 
	colnames(gains) <- c("Group 1", "Group 2", "pvalue")
			
	# Losses
	group1_CNVs <- intersectCNV(object, group1_indices, type = "loss")
	group2_CNVs <- intersectCNV(object, group2_indices, type = "loss")
	losses <- subgroupDifferenceCNVByType(group1_CNVs,
                                        group2_CNVs,
                                        group1_size,
                                        group2_size) 
	colnames(losses) <- c("Group 1", "Group 2", "pvalue")
			
	list("gains"  = gains,
       "losses" = losses)
})
