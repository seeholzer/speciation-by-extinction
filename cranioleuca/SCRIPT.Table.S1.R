##################################################
#	Description: Table S1
#
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github')

#load packages
library(plyr)

library(RColorBrewer)
library(colorspace)
library(stringi)

#load extinction scenarios used in mclust
#Extinction Scenarios
popdata = read.delim('cranioleuca/cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]

#load summary tables for bfd
bfd = read.delim('cranioleuca/results.bfd.txt',stringsAsFactors=F)

#load summary tables for mclust
mclust = read.delim('cranioleuca/results.mclust.txt',stringsAsFactors=F)
colnames(mclust)[1] = 'scenario'
#only use combined dataset
mclust = mclust[mclust$trait.set %in% 'combined', ]


#get the number of extinction events for each scenario 
scenarios = c(paste0('FULL_v',1:8),paste0('EXT_v',1:28))
tmp = popdata[,scenarios]
nExtinctions = apply(tmp,2,function(x){
	foo = data.frame(group=popdata$group,scenario=x)
	foo = unique(foo)
	nExtinctions = length(which(foo$scenario %in% 'extinct'))
	return(nExtinctions)
})
nExtinctions = data.frame(scenario=names(nExtinctions),nExtinctions=nExtinctions)


#get the midpoint of the extinction event 
# or
#point of the transition from antisiensis to baroni for the full dataset (nExtinctions=0)

x = tmp[,1]
extMidpoint = apply(tmp,2,function(x){
	foo = data.frame(group=popdata$group,scenario=x)
	foo = unique(foo)
	extMidpoint = mean(foo[foo$scenario %in% 'extinct','group'])
	if(is.na(extMidpoint)){	#transition point from ant to bar
		ant = max(foo[foo$scenario %in% 'ant','group'])
		bar = min(foo[foo$scenario %in% 'bar','group'])
		extMidpoint = mean(c(ant,bar))
	}
	return(extMidpoint)
})
extMidpoint = data.frame(scenario=names(extMidpoint), extMidpoint =extMidpoint)
extMidpoint$extMidpoint[is.nan(extMidpoint$extMidpoint)] = NA


#merge nExtinctionEvents and extMidpoint with bfd results
bfd = join_all(list(bfd,nExtinctions,extMidpoint),by='scenario',match='all')

#merge nExtinctionEvents and extMidpoint with mclust results
mclust = join_all(list(mclust,nExtinctions,extMidpoint),by='scenario',match='all')
mclust[mclust$scenario == 'FULL','nExtinctions'] = 0 

#make a key for the latitude of each extMidpoint
foo = aggregate(lat ~ group, data=popdata, mean)
colnames(foo)[1] = 'extMidpoint'
coo = data.frame(extMidpoint=foo$extMidpoint[-length(foo$extMidpoint)] + .5,lat= (diff(foo$lat)/2) + foo$lat[-length(foo$lat)])
lat = rbind(foo,coo)
key  = lat[order(lat[,1]),]


######################################################################
######################################################################
######################################################################
#	Table S1
######################################################################
######################################################################
######################################################################



#combine bfd and mclust
foo = merge(bfd[,c('scenario','nExtinctions','extMidpoint','MLE_1sp','MLE_2sp','BF')], mclust[,c('scenario','BIC.G1','BIC.G2','deltaBIC')],by='scenario',all=T)


full = foo[grep('FULL',foo$scenario), ]
ext = foo[grep('EXT',foo$scenario), ]

head(ext)

#reorder ext
ext = ext[order(as.numeric(gsub('EXT_v','', ext$scenario))),]

#combine Pre and Post-Extinction scenarios
table = rbind(full,ext)

table$extMidpointLat = mapvalues(table$extMidpoint,key[,'extMidpoint'],key[,'lat'])

table$PostPre = NA
table$PostPre[grep('FULL',table$scenario)] = 'Pre'
table$PostPre[grep('EXT',table$scenario)] = 'Post'

reformat = c('extMidpointLat','BF','deltaBIC')
table[,reformat] = apply(table[,reformat],2,function(x) round(x,2))

reformat = c('MLE_1sp','MLE_2sp','BIC.G1','BIC.G2')
table[,reformat] = apply(table[,reformat],2,function(x) round(x))


table$scenario = gsub('FULL_v|EXT_v','',table$scenario)
table$scenario = gsub('FULL','1',table$scenario)

table$dataset = mapvalues(table$PostPre,c('Pre','Post'),c('G','G+P'))
table$dataset[is.na(table$nExtinctions)] = 'P'

table = table[,c('dataset','PostPre','scenario','nExtinctions','extMidpoint','extMidpointLat','MLE_1sp','MLE_2sp','BF','BIC.G1','BIC.G2','deltaBIC')]

write.table(table,'cranioleuca/Table.S1.txt',quote=F,sep='\t',col.names=T,row.names=F)
