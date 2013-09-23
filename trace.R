# TODO: Add comment
# 
# Author: spapillo
###############################################################################


# load everything (RGSet)
library(minfi)
source('~/git/CopyNumber450k/R/generics.R')
source('~/git/CopyNumber450k/R/extractFromRGSet450k.R')
source('~/git/CopyNumber450k/R/CNVObject.R')
source('~/git/CopyNumber450k/R/formatSegments.R')
source('~/git/CopyNumber450k/R/FunNormCN450k.R')


path <- '~/Documents/iChange/data_ETMR'
RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))

CNVobj <- CNVObject(RGset)

CNVobj <- filterSNPProbes(CNVobj)

CNVobj <- normalize(CNVobj)

CNVobj <- buildSegments(CNVobj)

CNVobj <- createFilters(CNVobj)



common_gains <- lapply(common, function(name) {
			f_sample <- fishing[[name]]
			f_calls <- as.numeric(f_sample[,'p.adj.holm']) <= 0.01 &  as.numeric(f_sample[,'seg.mean']) > 0
			f_genes <- unique(unlist(strsplit(x=as.character(f_sample[f_calls, 'genes.list']), ";")))
			
			m_sample <- formatted.CNVobject[[name]]
			m_calls <- filter[[name]] & m_sample[,'logratio'] > 0
			m_genes <- unique(unlist(strsplit(x=as.character(m_sample[m_calls, 'genes']), ";")))
			
			return(list(exome=f_genes, methylation=m_genes))
			
		})

common_losses <- lapply(common, function(name) {
			f_sample <- fishing[[name]]
			f_calls <- as.numeric(f_sample[,'p.adj.holm']) <= 0.01 &  as.numeric(f_sample[,'seg.mean']) < 0
			f_genes <- unique(unlist(strsplit(x=as.character(f_sample[f_calls, 'genes.list']), ";")))
			
			m_sample <- formatted.CNVobject[[name]]
			m_calls <- filter[[name]] & m_sample[,'logratio'] < 0
			m_genes <- unique(unlist(strsplit(x=as.character(m_sample[m_calls, 'genes']), ";")))
			
			return(list(exome=f_genes, methylation=m_genes))
			
		})

lapply(1:length(common), function(i) {
			v <- VennFromSets(common_gains[[i]])
			pdf(paste(common[i], "_gains.pdf", sep=""))
			plot(v, doWeights=T)

			dev.off()
			pdf(paste(common[i], "_losses.pdf", sep=""))
			v <- VennFromSets(common_losses[[i]])
			plot(v, doWeights=T)
			dev.off()
		})

lapply(1:length(formatted.CNVobject), function(i) {
			
		})