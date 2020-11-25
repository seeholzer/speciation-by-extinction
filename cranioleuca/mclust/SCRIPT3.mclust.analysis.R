##################################################
#	Description: Species delimitation using Normal Mixture Models for scenarios of simulated extinction in Cranioleuca antisiensis
#	Script 3 - mclust analysis of phenotypic data with simulated extinction 
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#mclust analysis of phenotypic data 
#Seeholzer & Brumfield 20XX - Speciation-by-Extinction
#Author: Glenn F. Seeholzer
#Last Update: 21 June 2020


#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/mclust/')

#load packages
library(plyr)
library(mclust)
library(RColorBrewer)
options(scipen=999)

################################################################
####	load and transform data
################################################################
#Meta
meta = read.delim('cran.meta.txt',stringsAsFactors=F)
#Plumage
plumage = read.delim('cran.plumage.txt')
#Morphology
morph = read.delim('cran.morph.txt',stringsAsFactors=F)
#Merge meta and phenotypic data
data = join_all(list(meta,morph,plumage),by='ID',type='full')
#Remove non-adults, non-antisiensis
data = data[!(data$age %in% c('j')) & !(data$pop %in% 'curtata'), ]
#Remove populations not in original Seeholzer and Brumfield 2018 dataset
data = data[!(data$pop %in% c('Zarate','SanDamien','Pariacoto','Tayabamba')), ]
#Trait sets
morph.traits = c('b.l','b.w','b.d','w.l','t.max','t.min','t.w','ts.l','h.l')
plumage.traits = c('PC1.crow','PC1.back','PC1.rump','PC1.tail','PC1.bell','PC1.brea','PC1.thro','PC1.wing','PC1.cheek')
all.traits = c(morph.traits,plumage.traits)
traits = list(morphology=morph.traits,plumage=plumage.traits,combined=all.traits)


#Simulated Extinction Scenarios
popdata = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]


#Simulated extinction scenario x trait combinations
x = expand.grid(colnames(popdata)[-grep('pop|long|lat',colnames(popdata))],names(traits))
scenario.trait = apply(x,1,function(x) paste(x,collapse='.'))


################################################################
####	mclust analysis
################################################################

colnames = c('scenario.name','trait.set','bestK_BIC','BIC.G1','BIC.G2')
sum = data.frame(matrix(nrow=length(scenario.trait),ncol=length(colnames)))
colnames(sum) = colnames
sum[,c('scenario.name','trait.set')] = t(sapply(strsplit(scenario.trait,'\\.'),'['))

results.list = list()
i = 2
for(i in 1:length(scenario.trait)){
		scenario.name = strsplit(scenario.trait[i],'\\.')[[1]][1]
		trait.set = strsplit(scenario.trait[i],'\\.')[[1]][2]
		cat(i,'of',length(scenario.trait),' ',trait.set,' ', scenario.name,'\n')		
		scenario = popdata[,c('pop',scenario.name)]
		ext = scenario[scenario[,2] %in% 'extinct','pop'] 
		d = data[!(data$pop %in% ext),]

		complete.indices = which(complete.cases(d[, traits[[trait.set]] ]))
		subdata = d[complete.indices, c('ID', traits[[trait.set]])]
		rownames(subdata) = subdata[,1] #set rownames to ID
		#Transformation
		#mclust default is to transform variables using "scaled singular value decomposition (SVD) transformation"
		subdata = subdata[,-1]								
		
		#Alternative: scale variables and use directly without transformation
		#subdata = scale(subdata[,-1]) #scale data and remove ID column
		#mclust.options(hcUse="VARS") 
		
		#constrain models to just those comparing 1 vs. 2 clusters
		G = 2
		mclust = Mclust(subdata, G=1:G)
		sum[i,'bestK_BIC'] = summary(mclust)$G
		sum[i,c('BIC.G1','BIC.G2')] = apply(mclust$BIC, 1, max, na.rm=T)
		
		#classification summary
		class = mclust$classification; names(class) = dimnames(mclust$data)[[1]]
		foo = data[data$ID %in% names(class), c('ID','pop')]
		foo$class = mapvalues(foo$ID,names(class),class)
				
		p = popdata[,c('pop','long','lat', scenario.name)]
	
		Ntotal	= aggregate(class ~ pop,data=foo,function(x) length(x));colnames(Ntotal)[2] = 'Ntotal'
		p		= merge(p,Ntotal,all=T) 
		for(g in 1:max(class)){
			x = aggregate(class ~ pop,data=foo,function(x) length(which(x==g)))
			colnames(x)[2] = paste0('N',g)
			p = merge(p,x,all=T) 
		}
		
		results.list[[scenario.trait[i]]]$scenario.name = scenario.name
		results.list[[scenario.trait[i]]]$trait.set = trait.set
		results.list[[scenario.trait[i]]]$ext = ext
		results.list[[scenario.trait[i]]]$subdata = subdata
		results.list[[scenario.trait[i]]]$mclust = mclust
		results.list[[scenario.trait[i]]]$p = p

}#scenario loop

#positive values indicates support for two species
#negative values indicates support for one species species 
sum$deltaBIC = sum$BIC.G2 - sum$BIC.G1

#summary table of results
write.table(sum,'results.mclust.txt',sep='\t',col.names=T,row.names=F,quote=F)

#raw data used in plotting Figure 3
save(results.list,file='results.list.rda')

