##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	Script 1 - mclust modelling
#	- code adapted from Cadena et al. 2018
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/')


#load packages
library(mclust)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


##################################################
# 1.1) Read data

#read the morphological data and examine the resulting data frame
	data = read.table("Geospiza.data.csv", header=T, sep=",")
	data$ID = paste0(data$Institution,data$Museum.Number)

trait = 'Bdepth'	#trait selected against
d = data[data$New.Taxonomy %in% c('G. fortis','G. fuliginosa','G. magnirostris'), ]
mclust.data.file = 'mclust.Bdepth.SML.groundfinches.rda'
results.file = 'results.Bdepth.SML.groundfinches.rda'
G=1:10

#analyses for all Geospiza finches (LEGACY)
# trait = 'Bdepth'	#trait selected against
# d = data
# mclust.data.file = 'mclust.Bdepth.rda'
# results.file = 'results.Bdepth.rda'
# G=1:14



##################################################
# 1.2) Create Extinction Scenarios 

#extinction scenarios
	extinction.ranges = c()
	floor = floor(min(d[,trait]))
	ceiling = ceiling(max(d[,trait]))
	min.gap = 0.25
	max.gap = 5

	for(i in seq(floor,ceiling-min.gap,min.gap)){
		min = i
		max = seq(min+min.gap,min+max.gap,min.gap)
		max = max[max <= ceiling]
		#format to 2 decimal places
		min = trimws(format(round(min, 2), nsmall = 2))
		max = trimws(format(round(max, 2), nsmall = 2))

		extinction.ranges = c(extinction.ranges,paste0(min,'-',max))
	}	
	extinction.ranges = c('no.extinction',extinction.ranges)


#extinction.ranges = c('no.extinction','7.00-8.25')

##################################################
# 1.3) Run Extinction Scenarios

#empty lists for mclust model and results 
	mclust.data = list()
i = 1
#run extinction scenarios in parrellel
#parrellelization necessary for laptops
#runs in about 5 minutes
mclust.data = foreach(i=1:length(extinction.ranges),.packages = 'mclust') %dopar% {
		
		#cat(' scenario ',i,'/',length(extinction.ranges),' \n',sep='')
		
		if(extinction.ranges[i] %in% 'no.extinction'){
			#Ordering by Blength helps ensure that corresponding pre- and post-extinction morphospecies have same name 
			#pre-extinction morphospecies 2 corresponds to post-extinction morphospecies 2 assuming no splitting or lumping  
			extant = d[order(d$Blength), ]	
			#extant = d
			extant.ln = log(extant[,c(9:14)])
			
			colnames(extant.ln) = c("LnWing", "LnTail", "LnBlength", "LnBdepth", "LnBwidth", "LnTarsus")
			#PCA using the covariance matrix
			extant.ln.pca <- prcomp(extant.ln, center = T, scale = F)
			mclust.options(hcUse="VARS")
			Mcluster.extant.ln.pca = Mclust(extant.ln.pca$x[,c('PC1','PC2','PC3','PC4')], G=G)
			l = list()
			l$extinction.range = extinction.ranges[i]
			l$mclust = Mcluster.extant.ln.pca
			l$extant = extant[,'ID']
			l
		
		}else{
					
			min = as.numeric(strsplit(extinction.ranges[i],'-')[[1]][1])
			max = as.numeric(strsplit(extinction.ranges[i],'-')[[1]][2])
			
			extant = d[d[,trait] <= min | d[,trait] >= max, ]	
			#Ordering by Blength helps ensure that corresponding pre- and post-extinction morphospecies have same name 
			#pre-extinction morphospecies 2 corresponds to post-extinction morphospecies 2 assuming no splitting or lumping  
			extant = extant[order(extant$Blength), ]

			extant.ln = log(extant[,c(9:14)])
			colnames(extant.ln) = c("LnWing", "LnTail", "LnBlength", "LnBdepth", "LnBwidth", "LnTarsus")
			#PCA using the covariance matrix
			extant.ln.pca = prcomp(extant.ln, center = T, scale = F)
			#change default as needed
			#Mclust analysis
			#PC1-4 were used following the variable selection procedure of Cadena et al. 2018
			mclust.options(hcUse="VARS")
			Mcluster.extant.ln.pca = Mclust(extant.ln.pca$x[,c('PC1','PC2','PC3','PC4')], G=G)
			#extract BIC values for the best model conditional on the number of groups
			BIC.Best.Model.Per.G = apply(Mcluster.extant.ln.pca$BIC,1, function(x) max(x, na.rm=T))
			max.BIC.G = which.max(BIC.Best.Model.Per.G)    
			max.BIC = max(BIC.Best.Model.Per.G)
			
			l = list()
			l$extinction.range = extinction.ranges[i]
			l$mclust = Mcluster.extant.ln.pca
			l$extant = extant[,'ID']
			
			l
		}	
	}

names(mclust.data) = unlist(lapply(mclust.data,function(x) x$extinction.range))

save(mclust.data,file=mclust.data.file)
		

