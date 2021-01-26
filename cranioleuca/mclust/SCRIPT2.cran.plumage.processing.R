##################################################
#	Description: Species delimitation using Normal Mixture Models for scenarios of simulated extinction in Cranioleuca antisiensis
#	Script 2 - PCA of plumage color coordinates into PCA for each plumage patch   
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/mclust/')

#Load Packages
library(pavo)
library(plyr)

#VERSION 1

#raw color coordinates for each plumage patch
data = read.delim('cran.plumage.raw.txt',stringsAsFactors=F)

#PCA of raw color coordinates for each patch separately
	patch.names.new = c('crow','back','rump','tail','bell','brea','thro','wing','cheek')
	patch.data = data.frame(data$ID) ; colnames(patch.data) = 'ID'

i = 1
for(i in 1:length(patch.names.new)){
	patch.name = patch.names.new[i]

	tmp = data[,grep(paste0('ID|',patch.name),colnames(data))]
	pca = prcomp(tmp[,-1],scale=T,center=T)
	col.data = data.frame(ID=tmp$ID, PC1=pca$x[,1])
	colnames(col.data) = c('ID',paste0('PC1.',patch.name))
	
	patch.data = merge(patch.data,col.data,by='ID',all.x=T)
	
}

write.table(patch.data,'cran.plumage.txt',row.names=F,col.names=T,quote=F,sep='\t')


# #VERSION 2
# #import spec data
# myspecs = get(load('cran.spectra.Rdata'))
# patch.names = c('crow','back','rump','tail','bell','brea','thro','wing','cheek')

# patch.data = list()

# i = 1
# for(i in 1:length(patch.names)){
	# patch = myspecs[,c(1,grep(patch.names[i],colnames(myspecs)))]
	# patch = procspec(patch, opt='smooth', span = 0.2)
	# sum = summary(patch)
	# redness = apply(patch[,-1],2,function(x) lm(x~patch$wl)$coefficients[2] )
	# col.data = data.frame(ID=gsub(paste0('.',patch.names[i]),'',rownames(sum)),brightness=sum$B2,redness=redness)
	# colnames(col.data) = c('ID',paste0(c('brightness','redness'),'.',patch.names[i]))
		
	# patch.data[[patch.names[i]]] = col.data

# }

# patch.data = join_all(patch.data,by='ID')

# write.table(patch.data,'cran.plumage.2.txt',row.names=F,col.names=T,quote=F,sep='\t')




