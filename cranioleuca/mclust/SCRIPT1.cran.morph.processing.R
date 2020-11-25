##################################################
#	Description: Species delimitation using Normal Mixture Models for scenarios of simulated extinction in Cranioleuca antisiensis
#	Script 1 - Bayesian PCA imputation of missing values in morphological data  
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################


#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/mclust/')

#load packages
library(pcaMethods) #installation	https://www.bioconductor.org/packages/release/bioc/html/pcaMethods.html
library(plyr)

#load raw morphological data and modify certain columns
data = read.delim('cran.morph.raw.txt',stringsAsFactors=F)
data[data == 'damaged'] = NA
data[data == 'molt'] = NA
cols.to.numeric = c('mass','b.l','b.w','b.d','w.l','t.max','t.min','t.w','ts.l','h.l')
data[,cols.to.numeric] = apply(data[,cols.to.numeric],2,as.numeric)

#remove rows with juveniles
data = data[!(data$age %in% c('j')), ]

#			bPCA imputation of missing values
# 			Brown et al. 2012 Syst. Biol. 61(6):941–954
tmp = data
#retain samples that were measured (i.e. rows that are not all NA)
measured.samples = which(!apply(is.na(tmp[,-c(1:5)]),1,all))

ID = tmp[measured.samples,1]
mass = tmp[measured.samples,5]
pca.data = tmp[measured.samples,-c(1:5)]

## Perform Bayesian PCA with x components, where x is 1 less than the number of variables
pc = pca(pca.data, method="bpca", nPcs=ncol(pca.data)-1)
imputed = completeObs(pc)
data.morph = data.frame(ID,mass,imputed) 

write.table(data.morph,'cran.morph.txt',row.names=F,col.names=T,quote=F,sep='\t')


