##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	- code adapted from Cadena et al. 2018
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/')

#load packages
library(mclust)

##################################################
# 1.2) Read data

#read the morphological data and examine the resulting data frame
data = read.table("Geospiza.data.csv", header=T, sep=",")

d = data
#d = data[data$New_Taxonomy %in% c('G. fortis','G. fuliginosa','G. magnirostris'), ]


#empty lists for mclust model and results 
mclust.data = list()
results = c()


#extinction scenarios
g = 1
extinction.ranges = c()
for(g in 1:5){

	gap = g	
	floor = floor(min(d$Bdepth))
	ceiling = ceiling(max(d$Bdepth))
	min = floor:(ceiling-gap)
	max = (floor+gap):ceiling
	
	extinction.ranges = c(extinction.ranges,paste0(min,'-',max))
}
extinction.ranges = c('none',extinction.ranges)


#create empty results table
colnames = c('gap','extinction.range','min','max','max.BIC.G')
results = data.frame(matrix(nrow=length(extinction.ranges),ncol=length(colnames)))
colnames(results) = colnames
results$extinction.range = extinction.ranges
results$min = as.numeric(sapply(strsplit(extinction.ranges,'-'),'[',1))
results$max = as.numeric(sapply(strsplit(extinction.ranges,'-'),'[',2))
results$gap = results$max - results$min



#Full Dataset: no extinction
extant = d
extant.ln = log(extant[,c(9:14)])
colnames(extant.ln) = c("LnWing", "LnTail", "LnBlength", "LnBdepth", "LnBwidth", "LnTarsus")
#PCA using the covariance matrix
extant.ln.pca <- prcomp(extant.ln, center = T, scale = F)
#change default as needed
mclust.options(hcUse="VARS")
#Mclust analysis
#PC1-4 were used following the variable selection procedure of Cadena et al. 2018
Mcluster.extant.ln.pca.subset = Mclust(extant.ln.pca$x[,c('PC1','PC2','PC3','PC4')], G=1:10)
#extract BIC values for the best model conditional on the number of groups
BIC.Best.Model.Per.G = apply(Mcluster.extant.ln.pca.subset$BIC, 1, max, na.rm=T)
max.BIC.G = which.max(BIC.Best.Model.Per.G)    
max.BIC = max(BIC.Best.Model.Per.G)

mclust.data[['no.extinction']] = list(mclust = Mcluster.extant.ln.pca.subset, extant = d)
save(mclust.data,file='mclust.data.rda')

results[results$extinction.range %in% 'none','max.BIC.G'] = max.BIC.G
save(results,file='results.rda')

	
i = 2
for(i in 2:nrow(results)){
	
	cat(' scenario ',i,'/',nrow(results),' \n',sep='')
	
	extant = d[d$Bdepth <= results[i,'min'] | d$Bdepth >= results[i,'max'], ]	

	extant.ln = log(extant[,c(9:14)])
	colnames(extant.ln) = c("LnWing", "LnTail", "LnBlength", "LnBdepth", "LnBwidth", "LnTarsus")
	#PCA using the covariance matrix
	extant.ln.pca = prcomp(extant.ln, center = T, scale = F)
	#change default as needed
	mclust.options(hcUse="VARS")
	#Mclust analysis
	#PC1-4 were used following the variable selection procedure of Cadena et al. 2018
	Mcluster.extant.ln.pca.subset = Mclust(extant.ln.pca$x[,c('PC1','PC2','PC3','PC4')], G=1:10)
	#extract BIC values for the best model conditional on the number of groups
	BIC.Best.Model.Per.G = apply(Mcluster.extant.ln.pca.subset$BIC, 1, max, na.rm=T)
	max.BIC.G = which.max(BIC.Best.Model.Per.G)    
	max.BIC = max(BIC.Best.Model.Per.G)
	results[i,'max.BIC.G'] = max.BIC.G
	
	name = results[i,'extinction.range']
	mclust.data[[name]]$mclust = Mcluster.extant.ln.pca.subset
	mclust.data[[name]]$extant = extant[,'Museum.Number']

	save(mclust.data,file='mclust.data.rda')

	save(results,file='results.rda')

}





# run NMM across the entire dataset and for individual islands
# remove individuals based on bill depth or bill size (PC1) for a few different scenarios

# run analysis on all species and then only for G. fuliginosa, fortis, and magnirostris where the relationship between bill size is most tightly correlated with the adaptive landscape, for other species, specialized foraging niches may obscure relationship 



# # #For each island


# islands = c('ALL',unique(data$Island))
# res = c()

# f = islands[2]
# i = 1
# for(f in islands)	
	
	# if(f %in% 'ALL'){ d = data }else{ d = data[data$Island %in% f, ] }
	
	# floor = floor(min(d$Bdepth))
	# ceiling = ceiling(max(d$Bdepth))
	# min = floor:(ceiling-1)
	# max = (floor+1):ceiling
	# colnames = c('Island','extinction','max.G','min','max','max.BIC.G')
	# foo = data.frame(matrix(nrow=length(min)+1,ncol=length(colnames)))
	# colnames(foo) = colnames
	# foo$min = c(NA,min)
	# foo$max = c(NA,max)
	# foo$Island = f
	# foo$extinction = c('none',paste0(min,'-',max))
	# i = 2
	# for(i in 1:nrow(foo)){
		# cat(f,'(',which(islands %in% f),'/',19,')',' scenario',i,'/',nrow(foo),sep='')
		

		# if(foo$extinction[i] %in% 'none'){
			# tmp = d
			# G = 1:15
		# }else{
			# tmp = d[	d$Bdepth <= foo[i,'min'] | d$Bdepth >= foo[i,'max'], ]	
			# max.G = (length(unique(d$New_Taxonomy))*2)
			# if(max.G < 3) max.G = 3 
			# G = 1:3	
		# }

		# foo$max.G[i] = max(G)
		
		# table(data[,c('New_Taxonomy','Island')])
		
		# tmp.ln = log(tmp[,c(9:14)])
		# colnames(tmp.ln) = c("LnWing", "LnTail", "LnBlength", "LnBdepth", "LnBwidth", "LnTarsus")
		# #PCA using the covariance matrix
		# tmp.ln.pca <- prcomp(tmp.ln, center = T, scale = F)
		# #change default as needed
		# mclust.options(hcUse="VARS")
		# #Mclust analysis
		# Mcluster.tmp.ln.pca.subset = Mclust(tmp.ln.pca$x[,c('PC1','PC2','PC3','PC4')], G=G)
		# #extract BIC values for the best model conditional on the number of groups
		# BIC.Best.Model.Per.G = apply(Mcluster.tmp.ln.pca.subset$BIC, 1, max, na.rm=T)
		# max.BIC.G = which.max(BIC.Best.Model.Per.G)    
		# max.BIC = max(BIC.Best.Model.Per.G)
		# foo[i,'max.BIC.G'] = max.BIC.G
	# }
	# res = rbind(res,foo)
# }

    